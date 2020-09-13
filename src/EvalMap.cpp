/* Copyright (C) 2012-2020 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
#include <helib/EvalMap.h>
#include <helib/apiAttributes.h>

// needed to get NTL's TraceMap functions...needed for ThinEvalMap
#include <NTL/lzz_pXFactoring.h>
#include <NTL/GF2XFactoring.h>

namespace helib {

// Forward declarations
static BlockMatMul1D* buildStep1Matrix(const EncryptedArray& ea,
                                       std::shared_ptr<CubeSignature> sig,
                                       const NTL::Vec<long>& reps,
                                       long dim,
                                       long cofactor,
                                       bool invert,
                                       bool normal_basis);
static MatMul1D* buildStep2Matrix(const EncryptedArray& ea,
                                  std::shared_ptr<CubeSignature> sig,
                                  const NTL::Vec<long>& reps,
                                  long dim,
                                  long cofactor,
                                  bool invert);
static void init_representatives(NTL::Vec<long>& representatives,
                                 long dim,
                                 const NTL::Vec<long>& mvec,
                                 const PAlgebra& zMStar);

// Constructor: initializing tables for the evaluation-map transformations

EvalMap::EvalMap(const EncryptedArray& _ea,
                 bool minimal,
                 const NTL::Vec<long>& mvec,
                 bool _invert,
                 bool build_cache,
                 bool normal_basis) :
    ea(_ea), invert(_invert)
{
  const PAlgebra& zMStar = ea.getPAlgebra();

  long p = zMStar.getP();
  long d = zMStar.getOrdP();

  // FIXME: we should check that ea was initialized with
  // G == factors[0], but this is a slight pain to check
  // currently

  // NOTE: this code is derived from a more general setting, and
  // could certainly be greatly simplified

  nfactors = mvec.length();
  assertTrue(nfactors > 0, "Invalid argument: mvec must not be empty");

  for (long i = 0; i < nfactors; i++) {
    for (long j = i + 1; j < nfactors; j++) {
      assertEq(NTL::GCD(mvec[i], mvec[j]),
               1l,
               "Invalid argument: mvec elements must be pairwise co-prime");
    }
  }

  long m = computeProd(mvec);
  assertEq(m,
           (long)zMStar.getM(),
           "Invalid argument: Product of mvec elements does not match "
           "ea.zMStar.getM()");

  NTL::Vec<long> phivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    phivec[i] = phi_N(mvec[i]);
  long phim = computeProd(phivec);

  NTL::Vec<long> dprodvec(NTL::INIT_SIZE, nfactors + 1);
  dprodvec[nfactors] = 1;

  for (long i = nfactors - 1; i >= 0; i--)
    dprodvec[i] =
        dprodvec[i + 1] *
        multOrd(NTL::PowerMod(p % mvec[i], dprodvec[i + 1], mvec[i]), mvec[i]);

  NTL::Vec<long> dvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    dvec[i] = dprodvec[i] / dprodvec[i + 1];

  long nslots = phim / d;
  assertEq(d, dprodvec[0], "dprodvec must start with d");
  assertEq(nslots,
           (long)zMStar.getNSlots(),
           "Slot count mismatch between ea and phi(m)/d");

  long inertPrefix = 0;
  for (long i = 0; i < nfactors && dvec[i] == 1; i++) {
    inertPrefix++;
  }

  if (inertPrefix != nfactors - 1)
    throw LogicError("EvalMap: case not handled: bad inertPrefix");

  NTL::Vec<NTL::Vec<long>> local_reps(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    init_representatives(local_reps[i], i, mvec, zMStar);

  NTL::Vec<long> crtvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    crtvec[i] = (m / mvec[i]) * NTL::InvMod((m / mvec[i]) % mvec[i], mvec[i]);

  NTL::Vec<long> redphivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    redphivec[i] = phivec[i] / dvec[i];

  CubeSignature redphisig(redphivec);

  NTL::Vec<std::shared_ptr<CubeSignature>> sig_sequence;
  sig_sequence.SetLength(nfactors + 1);
  sig_sequence[nfactors] = std::make_shared<CubeSignature>(phivec);

  NTL::Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors - 1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = std::make_shared<CubeSignature>(reduced_phivec);
  }

  long dim = nfactors - 1;
  std::unique_ptr<BlockMatMul1D> mat1_data;
  mat1_data.reset(buildStep1Matrix(ea,
                                   sig_sequence[dim],
                                   local_reps[dim],
                                   dim,
                                   m / mvec[dim],
                                   invert,
                                   normal_basis));
  mat1.reset(new BlockMatMul1DExec(*mat1_data, minimal));

  matvec.SetLength(nfactors - 1);
  for (dim = nfactors - 2; dim >= 0; --dim) {
    std::unique_ptr<MatMul1D> mat_data;

    mat_data.reset(buildStep2Matrix(ea,
                                    sig_sequence[dim],
                                    local_reps[dim],
                                    dim,
                                    m / mvec[dim],
                                    invert));
    matvec[dim].reset(new MatMul1DExec(*mat_data, minimal));
  }

  if (build_cache)
    upgrade();
}

void EvalMap::upgrade()
{
  mat1->upgrade();
  for (long i = 0; i < matvec.length(); i++)
    matvec[i]->upgrade();
}

// Applying the evaluation (or its inverse) map to a ciphertext
void EvalMap::apply(Ctxt& ctxt) const
{
  if (!invert) { // forward direction
    mat1->mul(ctxt);

    for (long i = matvec.length() - 1; i >= 0; i--)
      matvec[i]->mul(ctxt);
  } else { // inverse transformation
    for (long i = 0; i < matvec.length(); i++)
      matvec[i]->mul(ctxt);

    mat1->mul(ctxt);
  }
}

static void init_representatives(NTL::Vec<long>& representatives,
                                 long dim,
                                 const NTL::Vec<long>& mvec,
                                 const PAlgebra& zMStar)
{
  assertInRange(dim,
                0l,
                mvec.length(),
                "Invalid argument: dim must be between 0 and mvec.length()");

  // special case
  if (dim >= LONG(zMStar.numOfGens())) {
    representatives.SetLength(1);
    representatives[0] = 1;
    return;
  }

  long m = mvec[dim];
  long D = zMStar.OrderOf(dim);
  long g = NTL::InvMod(zMStar.ZmStarGen(dim) % m, m);

  representatives.SetLength(D);
  for (long i = 0; i < D; i++)
    representatives[i] = NTL::PowerMod(g, i, m);
}

// The callback interface for the matrix-multiplication routines.

//! \cond FALSE (make doxygen ignore these classes)
template <typename type>
class Step2Matrix : public MatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<RX> A;

public:
  // constructor
  Step2Matrix(const EncryptedArray& _ea,
              std::shared_ptr<CubeSignature> _sig,
              const NTL::Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert = false) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig->getDim(dim) must equal reps.length()");

    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++)
      points[j] = RX(reps[j] * cofactor, 1) % G;

    A.SetDims(sz, sz);
    for (long j = 0; j < sz; j++)
      A[0][j] = 1;

    for (long i = 1; i < sz; i++)
      for (long j = 0; j < sz; j++)
        A[i][j] = (A[i - 1][j] * points[j]) % G;

    if (invert) {
      REBak ebak;
      ebak.save();
      ea.restoreContextForG();

      mat_RE A1, A2;
      conv(A1, A);

      long p = _ea.getAlMod().getZMStar().getP();
      long r = _ea.getAlMod().getR();

      ppInvert(A2, A1, p, r);
      conv(A, A2);
    }
  }

  bool get(RX& out, long i, long j, UNUSED long k) const override
  {
    out = A[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static MatMul1D* buildStep2Matrix(const EncryptedArray& ea,
                                  std::shared_ptr<CubeSignature> sig,
                                  const NTL::Vec<long>& reps,
                                  long dim,
                                  long cofactor,
                                  bool invert)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new Step2Matrix<PA_GF2>(ea, sig, reps, dim, cofactor, invert);

  case PA_zz_p_tag:
    return new Step2Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor, invert);

  default:
    return 0;
  }
}

template <typename type>
class Step1Matrix : public BlockMatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<mat_R> A;

public:
  // constructor
  Step1Matrix(const EncryptedArray& _ea,
              std::shared_ptr<CubeSignature> _sig,
              const NTL::Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert,
              bool normal_basis) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();
    long d = deg(G);

    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig->getDim(dim) must equal reps.length()");
    assertEq(dim,
             sig->getNumDims() - 1,
             "Invalid argument: dim must be one less than sig->getNumDims()");
    assertEq(sig->getSize(), ea.size(), "sig and ea do not have matching size");

    // so sz == phi(m_last)/d, where d = deg(G) = order of p mod m

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++)
      points[j] = RX(reps[j] * cofactor, 1) % G;

    NTL::Mat<RX> AA(NTL::INIT_SIZE, sz * d, sz);
    for (long j = 0; j < sz; j++)
      AA[0][j] = 1;

    for (long i = 1; i < sz * d; i++)
      for (long j = 0; j < sz; j++)
        AA[i][j] = (AA[i - 1][j] * points[j]) % G;

    A.SetDims(sz, sz);
    for (long i = 0; i < sz; i++)
      for (long j = 0; j < sz; j++) {
        A[i][j].SetDims(d, d);
        for (long k = 0; k < d; k++)
          VectorCopy(A[i][j][k], AA[i * d + k][j], d);
      }

    if (invert) {
      mat_R A1, A2;
      A1.SetDims(sz * d, sz * d);
      for (long i = 0; i < sz * d; i++)
        for (long j = 0; j < sz * d; j++)
          A1[i][j] = A[i / d][j / d][i % d][j % d];

      long p = _ea.getAlMod().getZMStar().getP();
      long r = _ea.getAlMod().getR();

      ppInvert(A2, A1, p, r);

      for (long i = 0; i < sz * d; i++)
        for (long j = 0; j < sz * d; j++)
          A[i / d][j / d][i % d][j % d] = A2[i][j];

      if (normal_basis) {
        const NTL::Mat<R>& CB = ea.getNormalBasisMatrix();

        // multiply each entry of A on the right by CB
        for (long i = 0; i < sz; i++)
          for (long j = 0; j < sz; j++)
            A[i][j] = A[i][j] * CB;
      } // if (normal_basis)
    }   // if (invert)
  }     // constructor

  bool get(mat_R& out, long i, long j, UNUSED long k) const override
  {
    out = A[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static BlockMatMul1D* buildStep1Matrix(const EncryptedArray& ea,
                                       std::shared_ptr<CubeSignature> sig,
                                       const NTL::Vec<long>& reps,
                                       long dim,
                                       long cofactor,
                                       bool invert,
                                       bool normal_basis)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new Step1Matrix<PA_GF2>(ea,
                                   sig,
                                   reps,
                                   dim,
                                   cofactor,
                                   invert,
                                   normal_basis);

  case PA_zz_p_tag:
    return new Step1Matrix<PA_zz_p>(ea,
                                    sig,
                                    reps,
                                    dim,
                                    cofactor,
                                    invert,
                                    normal_basis);

  default:
    return 0;
  }
}
//! \endcond

//=============== ThinEvalMap stuff

// needed to make generic programming work

void RelaxedInv(NTL::Mat<NTL::zz_p>& x, const NTL::Mat<NTL::zz_p>& a)
{
  relaxed_inv(x, a);
}

void RelaxedInv(NTL::Mat<NTL::GF2>& x, const NTL::Mat<NTL::GF2>& a)
{
  inv(x, a);
}

void TraceMap(NTL::GF2X& w,
              const NTL::GF2X& a,
              long d,
              const NTL::GF2XModulus& F,
              const NTL::GF2X& b)

{
  if (d < 0)
    throw InvalidArgument("TraceMap: d is negative");

  NTL::GF2X y, z, t;

  z = b;
  y = a;
  clear(w);

  while (d) {
    if (d == 1) {
      if (IsZero(w))
        w = y;
      else {
        CompMod(w, w, z, F);
        add(w, w, y);
      }
    } else if ((d & 1) == 0) {
      Comp2Mod(z, t, z, y, z, F);
      add(y, t, y);
    } else if (IsZero(w)) {
      w = y;
      Comp2Mod(z, t, z, y, z, F);
      add(y, t, y);
    } else {
      Comp3Mod(z, t, w, z, y, w, z, F);
      add(w, w, y);
      add(y, t, y);
    }

    d = d >> 1;
  }
}

// Forward declarations
static MatMul1D* buildThinStep1Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor);
static MatMul1D* buildThinStep2Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor,
                                      bool invert,
                                      bool inflate = false);
static void init_representatives(NTL::Vec<long>& representatives,
                                 long dim,
                                 const NTL::Vec<long>& mvec,
                                 const PAlgebra& zMStar);

// Constructor: initializing tables for the evaluation-map transformations

ThinEvalMap::ThinEvalMap(const EncryptedArray& _ea,
                         bool minimal,
                         const NTL::Vec<long>& mvec,
                         bool _invert,
                         bool build_cache) :
    ea(_ea), invert(_invert)
{
  const PAlgebra& zMStar = ea.getPAlgebra();

  long p = zMStar.getP();
  long d = zMStar.getOrdP();
  long sz = zMStar.numOfGens();

  // FIXME: we should check that ea was initialized with
  // G == factors[0], but this is a slight pain to check
  // currently

  // NOTE: this code is derived from a more general setting, and
  // could certainly be greatly simplified

  nfactors = mvec.length();
  assertTrue(nfactors > 0, "Invalid argument: mvec must have positive length");

  for (long i = 0; i < nfactors; i++) {
    for (long j = i + 1; j < nfactors; j++) {
      assertEq(NTL::GCD(mvec[i], mvec[j]),
               1l,
               "Invalid argument: mvec must have pairwise-disjoint entries");
    }
  }

  long m = computeProd(mvec);
  assertEq(m,
           (long)zMStar.getM(),
           "Invalid argument: mvec's product does not match ea's m");

  NTL::Vec<long> phivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    phivec[i] = phi_N(mvec[i]);
  long phim = computeProd(phivec);

  NTL::Vec<long> dprodvec(NTL::INIT_SIZE, nfactors + 1);
  dprodvec[nfactors] = 1;

  for (long i = nfactors - 1; i >= 0; i--)
    dprodvec[i] =
        dprodvec[i + 1] *
        multOrd(NTL::PowerMod(p % mvec[i], dprodvec[i + 1], mvec[i]), mvec[i]);

  NTL::Vec<long> dvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    dvec[i] = dprodvec[i] / dprodvec[i + 1];

  long nslots = phim / d;
  assertEq(d, dprodvec[0], "d must match the first entry of dprodvec");
  assertEq(nslots,
           (long)zMStar.getNSlots(),
           "Invalid argument: mismatch of number of slots");

  long inertPrefix = 0;
  for (long i = 0; i < nfactors && dvec[i] == 1; i++) {
    inertPrefix++;
  }

  if (inertPrefix != nfactors - 1)
    throw LogicError("ThinEvalMap: case not handled: bad inertPrefix");

  NTL::Vec<NTL::Vec<long>> local_reps(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    init_representatives(local_reps[i], i, mvec, zMStar);

  NTL::Vec<long> crtvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    crtvec[i] = (m / mvec[i]) * NTL::InvMod((m / mvec[i]) % mvec[i], mvec[i]);

  NTL::Vec<long> redphivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    redphivec[i] = phivec[i] / dvec[i];

  CubeSignature redphisig(redphivec);

  NTL::Vec<std::shared_ptr<CubeSignature>> sig_sequence;
  sig_sequence.SetLength(nfactors + 1);
  sig_sequence[nfactors] = std::make_shared<CubeSignature>(phivec);

  NTL::Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors - 1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = std::make_shared<CubeSignature>(reduced_phivec);
  }

  matvec.SetLength(nfactors);

  if (invert) {
    long dim = nfactors - 1;
    std::unique_ptr<MatMul1D> mat1_data;
    mat1_data.reset(buildThinStep1Matrix(ea,
                                         sig_sequence[dim],
                                         local_reps[dim],
                                         dim,
                                         m / mvec[dim]));
    matvec[dim].reset(new MatMul1DExec(*mat1_data, minimal));
  } else if (sz == nfactors) {
    long dim = nfactors - 1;
    std::unique_ptr<MatMul1D> mat1_data;
    mat1_data.reset(buildThinStep2Matrix(ea,
                                         sig_sequence[dim],
                                         local_reps[dim],
                                         dim,
                                         m / mvec[dim],
                                         invert,
                                         /*inflate=*/true));
    matvec[dim].reset(new MatMul1DExec(*mat1_data, minimal));
  }

  for (long dim = nfactors - 2; dim >= 0; --dim) {
    std::unique_ptr<MatMul1D> mat_data;

    mat_data.reset(buildThinStep2Matrix(ea,
                                        sig_sequence[dim],
                                        local_reps[dim],
                                        dim,
                                        m / mvec[dim],
                                        invert));
    matvec[dim].reset(new MatMul1DExec(*mat_data, minimal));
  }

  if (build_cache)
    upgrade();
}

void ThinEvalMap::upgrade()
{
  for (long i = 0; i < matvec.length(); i++)
    if (matvec[i])
      matvec[i]->upgrade();
}

// Applying the evaluation (or its inverse) map to a ciphertext
void ThinEvalMap::apply(Ctxt& ctxt) const
{
  if (!invert) { // forward direction
    for (long i = matvec.length() - 1; i >= 0; i--)
      if (matvec[i])
        matvec[i]->mul(ctxt);
  } else { // inverse transformation
    for (long i = 0; i < matvec.length(); i++)
      matvec[i]->mul(ctxt);
    traceMap(ctxt);
  }
}

// The callback interface for the matrix-multiplication routines.

//! \cond FALSE (make doxygen ignore these classes)
template <typename type>
class ThinStep2Matrix : public MatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<RX> A;

public:
  // constructor
  ThinStep2Matrix(const EncryptedArray& _ea,
                  std::shared_ptr<CubeSignature> _sig,
                  const NTL::Vec<long>& reps,
                  long _dim,
                  long cofactor,
                  bool invert,
                  bool inflate) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig and reps have inconsistent "
             "dimension");

    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();
    long d = deg(G);

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++) {
      points[j] = RX(reps[j] * cofactor, 1) % G;
      if (inflate)
        points[j] = NTL::PowerMod(points[j], d, G);
    }

    A.SetDims(sz, sz);
    for (long j = 0; j < sz; j++)
      A[0][j] = 1;

    for (long i = 1; i < sz; i++)
      for (long j = 0; j < sz; j++)
        A[i][j] = (A[i - 1][j] * points[j]) % G;

    if (invert) {
      REBak ebak;
      ebak.save();
      ea.restoreContextForG();

      mat_RE A1, A2;
      conv(A1, A);

      long p = _ea.getAlMod().getZMStar().getP();
      long r = _ea.getAlMod().getR();

      ppInvert(A2, A1, p, r);
      conv(A, A2);
    }
  }

  bool get(RX& out, long i, long j, UNUSED long k) const override
  {
    out = A[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static MatMul1D* buildThinStep2Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor,
                                      bool invert,
                                      bool inflate)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new ThinStep2Matrix<PA_GF2>(ea,
                                       sig,
                                       reps,
                                       dim,
                                       cofactor,
                                       invert,
                                       inflate);

  case PA_zz_p_tag:
    return new ThinStep2Matrix<PA_zz_p>(ea,
                                        sig,
                                        reps,
                                        dim,
                                        cofactor,
                                        invert,
                                        inflate);

  default:
    return 0;
  }
}

template <typename type>
class ThinStep1Matrix : public MatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<RX> A_deflated;

public:
  // constructor
  ThinStep1Matrix(const EncryptedArray& _ea,
                  std::shared_ptr<CubeSignature> _sig,
                  const NTL::Vec<long>& reps,
                  long _dim,
                  long cofactor) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RXModulus G(ea.getG());
    long d = deg(G);

    long p = _ea.getAlMod().getZMStar().getP();
    long r = _ea.getAlMod().getR();

    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig and reps have inconsistent "
             "dimension");
    assertEq(dim,
             sig->getNumDims() - 1,
             "Invalid argument: dim must be one less than "
             "sig->getNumDims()");
    assertEq(sig->getSize(), ea.size(), "sig and ea do not have matching size");

    // so sz == phi(m_last)/d, where d = deg(G) = order of p mod m

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++)
      points[j] = RX(reps[j] * cofactor, 1) % G;

    NTL::Mat<RX> AA(NTL::INIT_SIZE, sz * d, sz);
    for (long j = 0; j < sz; j++)
      AA[0][j] = 1;

    for (long i = 1; i < sz * d; i++)
      for (long j = 0; j < sz; j++)
        AA[i][j] = (AA[i - 1][j] * points[j]) % G;

    NTL::Mat<mat_R> A;
    A.SetDims(sz, sz);
    for (long i = 0; i < sz; i++)
      for (long j = 0; j < sz; j++) {
        A[i][j].SetDims(d, d);
        for (long k = 0; k < d; k++)
          VectorCopy(A[i][j][k], AA[i * d + k][j], d);
      }

    // if (invert) {
    // NOTE: this version is only used for the inverse matrix
    mat_R A1, A2;
    A1.SetDims(sz * d, sz * d);
    for (long i = 0; i < sz * d; i++)
      for (long j = 0; j < sz * d; j++)
        A1[i][j] = A[i / d][j / d][i % d][j % d];

    ppInvert(A2, A1, p, r);

    for (long i = 0; i < sz * d; i++)
      for (long j = 0; j < sz * d; j++)
        A[i / d][j / d][i % d][j % d] = A2[i][j];
    // }

    A_deflated.SetDims(sz, sz);
    vec_R v, w;
    v.SetLength(d);
    w.SetLength(d);

    RX h; // set h = X^p mod G
    PowerXMod(h, p, G);

    NTL::Vec<R> trace_vec;
    trace_vec.SetLength(2 * d - 1);
    // set trace_vec[i] = trace(X^i mod G)
    for (long i = 0; i < 2 * d - 1; i++) {
      RX trace_val;
      TraceMap(trace_val, (RX(i, 1) % G), d, G, h);
      assertTrue(deg(trace_val) <= 0, "trace_val is positive");
      trace_vec[i] = ConstTerm(trace_val);
    }

    NTL::Mat<R> trace_mat;
    trace_mat.SetDims(d, d);
    // set trace_mat[i][j] = trace(X^{i+j} mod G)
    for (long i = 0; i < d; i++)
      for (long j = 0; j < d; j++)
        trace_mat[i][j] = trace_vec[i + j];

    NTL::Mat<R> trace_mat_inv;
    RelaxedInv(trace_mat_inv, trace_mat);

    for (long i = 0; i < sz; i++)
      for (long j = 0; j < sz; j++) {
        for (long i1 = 0; i1 < d; i1++)
          v[i1] = A[i][j][i1][0];
        mul(w, v, trace_mat_inv);
        conv(A_deflated[i][j], w);
      }
  } // constructor

  bool get(RX& out, long i, long j, UNUSED long k) const override
  {
    out = A_deflated[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static MatMul1D* buildThinStep1Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new ThinStep1Matrix<PA_GF2>(ea, sig, reps, dim, cofactor);

  case PA_zz_p_tag:
    return new ThinStep1Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor);

  default:
    return 0;
  }
}
//! \endcond

} // namespace helib
