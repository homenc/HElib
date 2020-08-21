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

#include <helib/PAlgebra.h>
#include <helib/hypercube.h>
#include <helib/timing.h>
#include <helib/range.h>

#include <NTL/ZZXFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>
#include <NTL/BasicThreadPool.h>

#include <algorithm> // defines count(...), min(...)
#include <cmath>
#include <mutex> // std::mutex, std::unique_lock

namespace helib {

// polynomials are sorted lexicographically, with the
// constant term being the "most significant"

template <typename RX>
bool poly_comp(const RX& a, const RX& b);

bool less_than(NTL::GF2 a, NTL::GF2 b) { return rep(a) < rep(b); }
bool less_than(NTL::zz_p a, NTL::zz_p b) { return rep(a) < rep(b); }

bool less_than(const NTL::GF2X& a, const NTL::GF2X& b)
{
  return poly_comp(a, b);
}
bool less_than(const NTL::zz_pX& a, const NTL::zz_pX& b)
{
  return poly_comp(a, b);
}

bool less_than(const NTL::GF2E& a, const NTL::GF2E& b)
{
  return less_than(rep(a), rep(b));
}
bool less_than(const NTL::zz_pE& a, const NTL::zz_pE& b)
{
  return less_than(rep(a), rep(b));
}

bool less_than(const NTL::GF2EX& a, const NTL::GF2EX& b)
{
  return poly_comp(a, b);
}
bool less_than(const NTL::zz_pEX& a, const NTL::zz_pEX& b)
{
  return poly_comp(a, b);
}

template <typename RX>
bool poly_comp(const RX& a, const RX& b)
{
  long na = deg(a) + 1;
  long nb = deg(b) + 1;

  long i = 0;
  while (i < na && i < nb && coeff(a, i) == coeff(b, i))
    i++;

  if (i < na && i < nb)
    return less_than(coeff(a, i), coeff(b, i));
  else
    return na < nb;
}

bool PAlgebra::operator==(const PAlgebra& other) const
{
  if (m != other.m)
    return false;
  if (p != other.p)
    return false;

  return true;
}

long PAlgebra::exponentiate(const std::vector<long>& exps,
                            bool onlySameOrd) const
{
  if (isDryRun())
    return 1;
  long t = 1;
  long n = std::min(exps.size(), gens.size());
  for (long i = 0; i < n; i++) {
    if (onlySameOrd && !SameOrd(i))
      continue;
    long g = NTL::PowerMod(gens[i], exps[i], m);
    t = NTL::MulMod(t, g, m);
  }
  return t;
}

void PAlgebra::printout(std::ostream& out) const
{
  out << "m = " << m << ", p = " << p;
  if (isDryRun()) {
    out << " (dry run)" << std::endl;
    return;
  }

  out << ", phi(m) = " << phiM << std::endl;
  out << "  ord(p) = " << ordP << std::endl;
  out << "  normBnd = " << normBnd << std::endl;
  out << "  polyNormBnd = " << polyNormBnd << std::endl;

  std::vector<long> facs;
  factorize(facs, m);
  out << "  factors = " << facs << std::endl;

  for (std::size_t i = 0; i < gens.size(); i++)
    if (gens[i]) {
      // FIXME: is it really possible that gens[i] can be 0?
      // There is very likely some code here and there that
      // would break if that happens.

      out << "  generator " << gens[i] << " has order (";
      if (FrobPerturb(i) == 0)
        out << "=";
      else if (FrobPerturb(i) > 0)
        out << "!";
      else
        out << "!!";
      out << "= Z_m^*) of ";
      out << OrderOf(i) << std::endl;
    }

  if (cube.getSize() < 40) {
    out << "  T = [ ";
    for (auto const& t : T)
      out << t << " ";
    out << "]" << std::endl;
  }
}

void PAlgebra::printAll(std::ostream& out) const
{
  printout(out);
  if (cube.getSize() < 40) {
    out << "  Tidx = [ ";
    for (const auto& x : Tidx)
      out << x << " ";
    out << "]\n";
    out << "  zmsIdx = [ ";
    for (const auto& x : zmsIdx)
      out << x << " ";
    out << "]\n";
    out << "  zmsRep = [ ";
    for (const auto& x : zmsRep)
      out << x << " ";
    out << "]\n";
  }
}

static double cotan(double x) { return 1 / tan(x); }

half_FFT::half_FFT(long m) : fft(m / 2)
{
  typedef std::complex<double> cmplx_t;
  typedef long double ldbl;

  pow.resize(m / 2);
  for (long i : range(m / 2)) {
    // pow[i] = 2^{2*pi*I*(i/m)}
    ldbl angle = -((2.0L * PI) * (ldbl(i) / ldbl(m)));
    pow[i] = cmplx_t(std::cos(angle), std::sin(angle));
  }
}

quarter_FFT::quarter_FFT(long m) : fft(m / 4)
{
  typedef std::complex<double> cmplx_t;
  typedef long double ldbl;

  pow1.resize(m / 4);
  pow2.resize(m / 4);
  for (long i : range(m / 2)) {
    // pow[i] = 2^{2*pi*I*(i/m)}
    ldbl angle = -((2.0L * PI) * (ldbl(i) / ldbl(m)));
    if (i % 2)
      pow1[i >> 1] = cmplx_t(std::cos(angle), std::sin(angle));
    else
      pow2[i >> 1] = cmplx_t(std::cos(angle), std::sin(angle));
  }
}

static inline std::complex<double> MUL(std::complex<double> a,
                                       std::complex<double> b)
{
  double x = a.real(), y = a.imag(), u = b.real(), v = b.imag();
  return std::complex<double>(x * u - y * v, x * v + y * u);
}

static inline double ABS(std::complex<double> a)
{
  double x = a.real(), y = a.imag();
  return std::sqrt(x * x + y * y);
}

double calcPolyNormBnd(long m)
{
  assertTrue(m >= 1, "m >= 1");

  typedef std::complex<double> cmplx_t;
  typedef long double ldbl;

  // first, remove 2's
  while (m % 2 == 0)
    m /= 2;

  if (m == 1) {
    return 1;
  }

  std::vector<long> fac;
  factorize(fac, m);

  long radm = 1;
  for (long p : fac)
    radm *= p;

  if (fac.size() == 1) {
    long u = fac[0];
    return 2.0L * cotan(PI / (2.0L * u)) / u;
  }

  m = radm;

  long n = phi_N(m);

  NTL::ZZX PhiPoly = Cyclotomic(m);

  std::vector<double> a(n);
  for (long i : range(n))
    conv(a[i], PhiPoly[i]);
  // a does not include the leading coefficient 1
  // NOTE: according to the Arnold and Monogan paper
  // (Table 6) the least m such that the coefficients of Phi_m
  // do not fit in 53-bits is m=43,730,115.

  std::vector<cmplx_t> roots(m);
  std::vector<cmplx_t> x(n);

  for (long i : range(m)) {
    ldbl re = std::cos(2.0L * PI * (ldbl(i) / ldbl(m)));
    ldbl im = std::sin(2.0L * PI * (ldbl(i) / ldbl(m)));
    roots[i] = cmplx_t(re, im);
  }

  std::vector<long> res_tab(n);

  long row_num = 0;
  for (long i : range(1, m)) {
    if (NTL::GCD(i, m) != 1)
      continue;
    x[row_num] = roots[i];
    res_tab[row_num] = i;
    row_num++;
  }

  std::vector<double> dist_tab_vec(2 * m - 1);
  std::vector<int> dist_exp_tab_vec(2 * m);

  double* dist_tab = &dist_tab_vec[m - 1];
  int* dist_exp_tab = &dist_exp_tab_vec[m - 1];

  const double sqrt2_inv = 1.0 / std::sqrt(ldbl(2));
  constexpr long FREXP_ITER = 1600;

  for (long i : range(1, m)) {
    dist_tab[i] = std::frexp(double(2.0L * std::sin(PI * (ldbl(i) / ldbl(m)))),
                             &dist_exp_tab[i]);

    if (dist_tab[i] < sqrt2_inv) {
      dist_tab[i] *= 2.0;
      dist_exp_tab[i]--;
    }

    dist_tab[-i] = dist_tab[i];
    dist_exp_tab[-i] = dist_exp_tab[i];
  }

  dist_tab[0] = 1;
  dist_exp_tab[0] = 0;

  std::vector<double> global_norm_col(n);
  for (long i : range(n))
    global_norm_col[i] = 0;
  std::mutex global_norm_col_mutex;

  NTL_EXEC_RANGE(n, first, last)

  std::vector<double> norm_col(n);
  for (long i : range(n))
    norm_col[i] = 0;

  long j = first;

  for (; j <= last - 2; j += 2) {
    // process columns j and j+1 of inverse matrix
    // NOTE: processing columns two at a time gives an almost 2x speedup

    long res_j = res_tab[j];
    long res_j_1 = res_tab[j + 1];

    double prod = 1;
    double prod_1 = 1;
    long e_total = 0;
    long e_total_1 = 0;
    int e;
    int e_1;

    {
      long i = 0;
      while (i <= n - FREXP_ITER) {
        for (long k = 0; k < FREXP_ITER; k++) {
          long res_i = res_tab[i + k];
          prod *= dist_tab[res_i - res_j];
          prod_1 *= dist_tab[res_i - res_j_1];
          e_total += dist_exp_tab[res_i - res_j];
          e_total_1 += dist_exp_tab[res_i - res_j_1];
        }
        prod = std::frexp(prod, &e);
        prod_1 = std::frexp(prod_1, &e_1);
        e_total += e;
        e_total_1 += e_1;

        i += FREXP_ITER;
      }
      while (i < n) {
        long res_i = res_tab[i];
        prod *= dist_tab[res_i - res_j];
        prod_1 *= dist_tab[res_i - res_j_1];
        e_total += dist_exp_tab[res_i - res_j];
        e_total_1 += dist_exp_tab[res_i - res_j_1];
        i++;
      }
    }

    prod = std::ldexp(prod, e_total);
    prod_1 = std::ldexp(prod_1, e_total_1);

    double inv_prod = 1.0 / prod;
    double inv_prod_1 = 1.0 / prod_1;

    cmplx_t xj = x[j];
    cmplx_t xj_1 = x[j + 1];
    cmplx_t q = 1;
    cmplx_t q_1 = 1;

    norm_col[0] += (inv_prod + inv_prod_1);

    for (long i : range(1, n)) {
      q = MUL(q, xj) + a[n - i];
      q_1 = MUL(q_1, xj_1) + a[n - i];
      norm_col[i] += (ABS(q) * inv_prod + ABS(q_1) * inv_prod_1);
    }
  }

  if (j == last - 1) {
    // process column j of inverse matrix

    long res_j = res_tab[j];

    double prod = 1;
    long e_total = 0;
    int e;

    {
      long i = 0;
      while (i <= n - FREXP_ITER) {
        for (long k = 0; k < FREXP_ITER; k++) {
          long res_i = res_tab[i + k];
          prod *= dist_tab[res_i - res_j];
          e_total += dist_exp_tab[res_i - res_j];
        }
        prod = std::frexp(prod, &e);
        e_total += e;
        i += FREXP_ITER;
      }
      while (i < n) {
        long res_i = res_tab[i];
        prod *= dist_tab[res_i - res_j];
        e_total += dist_exp_tab[res_i - res_j];
        i++;
      }
    }

    prod = std::ldexp(prod, e_total);

    double inv_prod = 1.0 / prod;

    cmplx_t xj = x[j];
    cmplx_t q = 1;
    norm_col[0] += inv_prod;
    for (long i : range(1, n)) {
      q = MUL(q, xj) + a[n - i];
      norm_col[i] += ABS(q) * inv_prod;
    }
  }

  std::lock_guard<std::mutex> guard(global_norm_col_mutex);

  for (long i : range(n))
    global_norm_col[i] += norm_col[i];

  NTL_EXEC_INDEX_END

  double max_norm = 0;
  for (long i : range(n)) {
    if (max_norm < global_norm_col[i])
      max_norm = global_norm_col[i];
  }

  return max_norm;
}

PAlgebra::PAlgebra(long mm,
                   long pp,
                   const std::vector<long>& _gens,
                   const std::vector<long>& _ords) :
    m(mm), p(pp), cM(1.0) // default value for the ring constant
{
  assertInRange<InvalidArgument>(mm,
                                 2l,
                                 NTL_SP_BOUND,
                                 "mm is not in [2, NTL_SP_BOUND)");
  if (pp == -1) // pp==-1 signals using the complex field for plaintext
    pp = m - 1;
  else {
    assertTrue<InvalidArgument>((bool)NTL::ProbPrime(pp),
                                "Modulus pp is not prime (nor -1)");
    assertNeq<InvalidArgument>(mm % pp, 0l, "Modulus pp divides mm");
  }

  long k = NTL::NextPowerOfTwo(mm);
  if (static_cast<unsigned long>(mm) == (1UL << k)) // m is a power of two
    pow2 = k;
  else if (p != -1) // is not power of two, set to zero (even if m is even!)
    pow2 = 0;
  else // CKKS requires m to be a power of two.  Throw if not.
    throw InvalidArgument("CKKS scheme only supports m as a power of two.");

  // For dry-run, use a tiny m value for the PAlgebra tables
  if (isDryRun())
    m = (p == 3) ? 4 : 3;

  // Compute the generators for (Z/mZ)^* (defined in NumbTh.cpp)

  std::vector<long> tmpOrds;
  if (_gens.size() > 0 && _gens.size() == _ords.size() && !isDryRun()) {
    // externally supplied generator,orders
    tmpOrds = _ords;
    this->gens = _gens;
    this->ordP = multOrd(pp, mm);
  } else
    // treat externally supplied generators (if any) as candidates
    this->ordP = findGenerators(this->gens, tmpOrds, mm, pp, _gens);

  // Record for each generator gi whether it has the same order in
  // ZM* as in Zm* /(p,g1,...,g_{i-1})

  resize(native, lsize(tmpOrds));
  resize(frob_perturb, lsize(tmpOrds));
  std::vector<long> p_subgp(mm);
  for (long i : range(mm))
    p_subgp[i] = -1;
  long pmodm = pp % mm;
  p_subgp[1] = 0;
  for (long i = 1, p2i = pmodm; p2i != 1; i++, p2i = NTL::MulMod(p2i, pmodm, m))
    p_subgp[p2i] = i;
  for (long j : range(tmpOrds.size())) {
    tmpOrds[j] = std::abs(tmpOrds[j]);
    // for backward compatibility, a user supplied
    // ords value could be negative, but we ignore that here.
    // For testing and debugging, we may want to not ignore this...

    long i = NTL::PowerMod(this->gens[j], tmpOrds[j], m);

    native[j] = (i == 1);
    frob_perturb[j] = p_subgp[i];
  }

  cube.initSignature(tmpOrds); // set hypercube with these dimensions

  phiM = ordP * getNSlots();

  NTL::Vec<NTL::Pair<long, long>> factors;
  factorize(factors, mm);
  nfactors = factors.length();

  radm = 1;
  for (long i : range(nfactors))
    radm *= factors[i].a;

  normBnd = 1.0;
  for (long i : range(nfactors)) {
    long u = factors[i].a;
    normBnd *= 2.0L * cotan(PI / (2.0L * u)) / u;
  }

  polyNormBnd = calcPolyNormBnd(mm);

  // Allocate space for the various arrays
  resize(T, getNSlots());
  Tidx.assign(mm, -1);   // allocate m slots, initialize them to -1
  zmsIdx.assign(mm, -1); // allocate m slots, initialize them to -1
  resize(zmsRep, phiM);
  long i, idx;
  for (i = idx = 0; i < mm; i++) {
    if (NTL::GCD(i, mm) == 1) {
      zmsIdx[i] = idx++;
      zmsRep[zmsIdx[i]] = i;
    }
  }

  // Now fill the Tidx translation table. We identify an element t \in T
  // with its representation t = \prod_{i=0}^n gi^{ei} mod m (where the
  // gi's are the generators in gens[]) , represent t by the vector of
  // exponents *in reverse order* (en,...,e1,e0), and order these vectors
  // in lexicographic order.

  // FIXME: is the comment above about reverse order true?
  // It doesn't seem like it to me, VJS.
  // The comment about reverse order is correct, SH.

  // buffer is initialized to all-zero, which represents 1=\prod_i gi^0
  std::vector<long> buffer(gens.size()); // temporary holds exponents
  i = idx = 0;
  long ctr = 0;
  do {
    ctr++;
    long t = exponentiate(buffer);

    // sanity check for user-supplied gens
    assertEq(NTL::GCD(t, mm), 1l, "Bad user-supplied generator");
    assertEq(Tidx[t], -1l, "Slot at index t has already been assigned");

    T[i] = t;      // The i'th element in T it t
    Tidx[t] = i++; // the index of t in T is i

    // increment buffer by one (in lexicographic order)
  } while (nextExpVector(buffer)); // until we cover all the group

  // sanity check for user-supplied gens
  assertEq(ctr, getNSlots(), "Bad user-supplied generator set");

  PhimX = Cyclotomic(mm); // compute and store Phi_m(X)
  //  pp_factorize(mFactors,mm); // prime-power factorization from NumbTh.cpp

  if (mm % 2 == 0)
    half_fftInfo = std::make_shared<half_FFT>(mm);
  else
    fftInfo = std::make_shared<PGFFT>(mm);

  // fftInfo = std::make_shared<PGFFT>(mm); // Need this for some
  // debugging/timing

  if (mm % 4 == 0)
    quarter_fftInfo = std::make_shared<quarter_FFT>(mm);
}

bool comparePAlgebra(const PAlgebra& palg,
                     unsigned long m,
                     unsigned long p,
                     UNUSED unsigned long r,
                     const std::vector<long>& gens,
                     const std::vector<long>& ords)
{
  if (static_cast<unsigned long>(palg.getM()) != m ||
      static_cast<unsigned long>(palg.getP()) != p ||
      static_cast<std::size_t>(palg.numOfGens()) != gens.size() ||
      static_cast<std::size_t>(palg.numOfGens()) != ords.size())
    return false;

  for (long i = 0; i < (long)gens.size(); i++) {
    if (long(palg.ZmStarGen(i)) != gens[i])
      return false;

    if ((palg.SameOrd(i) && palg.OrderOf(i) != ords[i]) ||
        (!palg.SameOrd(i) && palg.OrderOf(i) != -ords[i]))
      return false;
  }
  return true;
}

long PAlgebra::frobeniusPow(long j) const
{
  return NTL::PowerMod(mcMod(p, m), j, m);
  // Don't forget to reduce p mod m!!
}

long PAlgebra::genToPow(long i, long j) const
{
  long sz = gens.size();

  if (i == sz) {
    assertTrue(j == 0, "PAlgebra::genToPow: i == sz but j != 0");
    return 1;
  }

  assertTrue(i >= -1 && i < LONG(gens.size()), "PAlgebra::genToPow: bad dim");

  long res;
  if (i == -1)
    res = frobeniusPow(j);
  else
    res = NTL::PowerMod(gens[i], j, m);

  return res;
}

/***********************************************************************

  PAlgebraMod stuff....

************************************************************************/

PAlgebraModBase* buildPAlgebraMod(const PAlgebra& zMStar, long r)
{
  long p = zMStar.getP();

  if (p == -1) // complex plaintext space
    return new PAlgebraModCx(zMStar, r);

  assertTrue<InvalidArgument>(p >= 2,
                              "Modulus p is less than 2 (nor -1 for CKKS)");
  assertTrue<InvalidArgument>(r > 0, "Hensel lifting r is less than 1");
  if (p == 2 && r == 1)
    return new PAlgebraModDerived<PA_GF2>(zMStar, r);
  else
    return new PAlgebraModDerived<PA_zz_p>(zMStar, r);
}

template <typename T>
void PAlgebraLift(const NTL::ZZX& phimx,
                  const T& lfactors,
                  T& factors,
                  T& crtc,
                  long r);

// Missing NTL functionality

void EDF(NTL::vec_zz_pX& v, const NTL::zz_pX& f, long d)
{
  EDF(v, f, PowerXMod(NTL::zz_p::modulus(), f), d);
}

NTL::zz_pEX FrobeniusMap(const NTL::zz_pEXModulus& F)
{
  return PowerXMod(NTL::zz_pE::cardinality(), F);
}

template <typename type>
PAlgebraModDerived<type>::PAlgebraModDerived(const PAlgebra& _zMStar, long _r) :
    zMStar(_zMStar), r(_r)

{
  long p = zMStar.getP();
  long m = zMStar.getM();

  // For dry-run, use a tiny m value for the PAlgebra tables
  if (isDryRun())
    m = (p == 3) ? 4 : 3;

  assertTrue<InvalidArgument>(r > 0l, "Hensel lifting r is less than 1");

  NTL::ZZ BigPPowR = NTL::power_ZZ(p, r);
  assertTrue((bool)BigPPowR.SinglePrecision(),
             "BigPPowR is not SinglePrecision");
  pPowR = to_long(BigPPowR);

  long nSlots = zMStar.getNSlots();

  RBak bak;
  bak.save();
  SetModulus(p);

  // Compute the factors Ft of Phi_m(X) mod p, for all t \in T

  RX phimxmod;

  conv(phimxmod, zMStar.getPhimX()); // Phi_m(X) mod p

  vec_RX localFactors;

  EDF(localFactors, phimxmod, zMStar.getOrdP()); // equal-degree factorization

  RX* first = &localFactors[0];
  RX* last = first + lsize(localFactors);
  RX* smallest =
      std::min_element(first,
                       last,
                       static_cast<bool (*)(const RX&, const RX&)>(less_than));
  swap(*first, *smallest);

  // We make the lexicographically smallest factor have index 0.
  // The remaining factors are ordered according to their representatives.

  RXModulus F1(localFactors[0]);
  for (long i = 1; i < nSlots; i++) {
    long t = zMStar.ith_rep(i);      // Ft is minimal poly of x^{1/t} mod F1
    long tInv = NTL::InvMod(t, m);   // tInv = t^{-1} mod m
    RX X2tInv = PowerXMod(tInv, F1); // X2tInv = X^{1/t} mod F1
    NTL::IrredPolyMod(localFactors[i], X2tInv, F1);
    // IrredPolyMod(X,P,Q) returns in X the minimal polynomial of P mod Q
  }
  /* Debugging sanity-check #1: we should have Ft= GCD(F1(X^t),Phi_m(X))
  for (i=1; i<nSlots; i++) {
    long t = T[i];
    RX X2t = PowerXMod(t,phimxmod);  // X2t = X^t mod Phi_m(X)
    RX Ft = GCD(CompMod(F1,X2t,phimxmod),phimxmod);
    if (Ft != localFactors[i]) {
      cout << "Ft != F1(X^t) mod Phi_m(X), t=" << t << endl;
      exit(0);
    }
  }*******************************************************************/

  if (r == 1) {
    build(PhimXMod, phimxmod);
    factors = localFactors;
    pPowRContext.save();

    // Compute the CRT coefficients for the Ft's
    resize(crtCoeffs, nSlots);
    for (long i = 0; i < nSlots; i++) {
      RX te = phimxmod / factors[i];        // \prod_{j\ne i} Fj
      te %= factors[i];                     // \prod_{j\ne i} Fj mod Fi
      InvMod(crtCoeffs[i], te, factors[i]); // \prod_{j\ne i} Fj^{-1} mod Fi
    }
  } else {
    PAlgebraLift(zMStar.getPhimX(), localFactors, factors, crtCoeffs, r);
    RX phimxmod1;
    conv(phimxmod1, zMStar.getPhimX());
    build(PhimXMod, phimxmod1);
    pPowRContext.save();
  }

  // set factorsOverZZ
  resize(factorsOverZZ, nSlots);
  for (long i = 0; i < nSlots; i++)
    conv(factorsOverZZ[i], factors[i]);

  genCrtTable();
  genMaskTable();
}

// Assumes current zz_p modulus is p^r
// computes S = F^{-1} mod G via Hensel lifting
void InvModpr(NTL::zz_pX& S,
              const NTL::zz_pX& F,
              const NTL::zz_pX& G,
              long p,
              long r)
{
  NTL::ZZX ff, gg, ss, tt;

  ff = to_ZZX(F);
  gg = to_ZZX(G);

  NTL::zz_pBak bak;
  bak.save();
  NTL::zz_p::init(p);

  NTL::zz_pX f, g, s, t;
  f = to_zz_pX(ff);
  g = to_zz_pX(gg);
  s = InvMod(f, g);
  t = (1 - s * f) / g;
  assertTrue(static_cast<bool>(s * f + t * g == 1l),
             "Arithmetic error during Hensel lifting");
  ss = to_ZZX(s);
  tt = to_ZZX(t);

  NTL::ZZ pk = NTL::to_ZZ(1);

  for (long k = 1; k < r; k++) {
    // lift from p^k to p^{k+1}
    pk = pk * p;

    assertTrue((bool)divide(ss * ff + tt * gg - 1, pk),
               "Arithmetic error during Hensel lifting");

    NTL::zz_pX d = to_zz_pX((1 - (ss * ff + tt * gg)) / pk);
    NTL::zz_pX s1, t1;
    s1 = (s * d) % g;
    t1 = (d - s1 * f) / g;
    ss = ss + pk * to_ZZX(s1);
    tt = tt + pk * to_ZZX(t1);
  }

  bak.restore();

  S = to_zz_pX(ss);

  assertTrue(static_cast<bool>((S * F) % G == 1),
             "Hensel lifting failed to find solutions");
}

// FIXME: Consider changing this function to something non-templated.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template <typename T>
void PAlgebraLift(const NTL::ZZX& phimx,
                  const T& lfactors,
                  T& factors,
                  T& crtc,
                  long r)
{
  throw LogicError("Uninstantiated version of PAlgebraLift");
}
#pragma GCC diagnostic pop

// This specialized version of PAlgebraLift does the hensel
// lifting needed to finish off the initialization.
// It assumes the zz_p modulus is initialized to p
// when called, and leaves it set to p^r

template <>
void PAlgebraLift(const NTL::ZZX& phimx,
                  const NTL::vec_zz_pX& lfactors,
                  NTL::vec_zz_pX& factors,
                  NTL::vec_zz_pX& crtc,
                  long r)
{
  long p = NTL::zz_p::modulus();
  long nSlots = lsize(lfactors);

  NTL::vec_ZZX vzz; // need to go via ZZX

  // lift the factors of Phi_m(X) from mod-2 to mod-2^r
  if (lsize(lfactors) > 1)
    MultiLift(vzz, lfactors, phimx, r); // defined in NTL::ZZXFactoring
  else {
    resize(vzz, 1);
    vzz[0] = phimx;
  }

  // Compute the zz_pContext object for mod p^r arithmetic
  NTL::zz_p::init(NTL::power_long(p, r));

  NTL::zz_pX phimxmod = to_zz_pX(phimx);
  resize(factors, nSlots);
  for (long i = 0; i < nSlots; i++) // Convert from ZZX to zz_pX
    conv(factors[i], vzz[i]);

  // Finally compute the CRT coefficients for the factors
  resize(crtc, nSlots);
  for (long i = 0; i < nSlots; i++) {
    NTL::zz_pX& fct = factors[i];
    NTL::zz_pX te = phimxmod / fct;   // \prod_{j\ne i} Fj
    te %= fct;                        // \prod_{j\ne i} Fj mod Fi
    InvModpr(crtc[i], te, fct, p, r); // \prod_{j\ne i} Fj^{-1} mod Fi
  }
}

// Returns a vector crt[] such that crt[i] = p mod Ft (with t = T[i])
template <typename type>
void PAlgebraModDerived<type>::CRT_decompose(std::vector<RX>& crt,
                                             const RX& H) const
{
  long nSlots = zMStar.getNSlots();

  if (isDryRun()) {
    crt.clear();
    return;
  }
  resize(crt, nSlots);
  for (long i = 0; i < nSlots; i++)
    rem(crt[i], H, factors[i]); // crt[i] = H % factors[i]
}

template <typename type>
void PAlgebraModDerived<type>::embedInAllSlots(
    RX& H,
    const RX& alpha,
    const MappingData<type>& mappingData) const
{
  if (isDryRun()) {
    H = RX::zero();
    return;
  }
  HELIB_TIMER_START;
  long nSlots = zMStar.getNSlots();

  std::vector<RX> crt(nSlots); // allocate space for CRT components

  // The i'th CRT component is (H mod F_t) = alpha(maps[i]) mod F_t,
  // where with t=T[i].

  if (IsX(mappingData.G) || deg(alpha) <= 0) {
    // special case...no need for CompMod, which is
    // is not optimized for this case

    for (long i = 0; i < nSlots; i++) // crt[i] = alpha(maps[i]) mod Ft
      crt[i] = ConstTerm(alpha);
  } else {
    // general case...

    // FIXME: should update this to use matrix_maps, but this routine
    // isn't actually used anywhere

    for (long i = 0; i < nSlots; i++) // crt[i] = alpha(maps[i]) mod Ft
      CompMod(crt[i], alpha, mappingData.maps[i], factors[i]);
  }

  CRT_reconstruct(H, crt); // interpolate to get H
  HELIB_TIMER_STOP;
}

template <typename type>
void PAlgebraModDerived<type>::embedInSlots(
    RX& H,
    const std::vector<RX>& alphas,
    const MappingData<type>& mappingData) const
{
  if (isDryRun()) {
    H = RX::zero();
    return;
  }
  HELIB_TIMER_START;

  long nSlots = zMStar.getNSlots();
  // assert(lsize(alphas) == nSlots);
  assertEq(
      lsize(alphas),
      nSlots,
      "Cannot embed in slots: alphas size is different than number of slots");

  long d = mappingData.degG;
  for (long i = 0; i < nSlots; i++)
    assertTrue(deg(alphas[i]) < d,
               "Bad alpha element at index i: its degree is greater or "
               "equal than mappingData.degG");

  std::vector<RX> crt(nSlots); // allocate space for CRT components

  // The i'th CRT component is (H mod F_t) = alphas[i](maps[i]) mod F_t,
  // where with t=T[i].

  if (IsX(mappingData.G)) {
    // special case...no need for CompMod, which is
    // is not optimized for this case

    for (long i = 0; i < nSlots; i++) // crt[i] = alpha(maps[i]) mod Ft
      crt[i] = ConstTerm(alphas[i]);
  } else {
    // general case...still try to avoid CompMod when possible,
    // which is the common case for encoding masks

    HELIB_NTIMER_START(CompMod);

#if 0
    for (long i: range(nSlots)) {
      if (deg(alphas[i]) <= 0)
        crt[i] = alphas[i];
      else
        CompMod(crt[i], alphas[i], mappingData.maps[i], factors[i]);
    }
#else
    vec_R in, out;

    for (long i : range(nSlots)) {
      if (deg(alphas[i]) <= 0)
        crt[i] = alphas[i];
      else {
        VectorCopy(in, alphas[i], d);
        mul(out, in, mappingData.matrix_maps[i]);
        conv(crt[i], out);
      }
    }
#endif
  }

  CRT_reconstruct(H, crt); // interpolate to get p

  HELIB_TIMER_STOP;
}

template <typename type>
void PAlgebraModDerived<type>::CRT_reconstruct(RX& H,
                                               std::vector<RX>& crt) const
{
  if (isDryRun()) {
    H = RX::zero();
    return;
  }
  HELIB_TIMER_START;
  long nslots = zMStar.getNSlots();

  const std::vector<RX>& ctab = crtTable;

  clear(H);
  RX tmp1, tmp2;

  bool easy = true;
  for (long i = 0; i < nslots; i++)
    if (!IsZero(crt[i]) && !IsOne(crt[i])) {
      easy = false;
      break;
    }

  if (easy) {
    for (long i = 0; i < nslots; i++)
      if (!IsZero(crt[i]))
        H += ctab[i];
  } else {
    std::vector<RX> crt1;
    resize(crt1, nslots);
    for (long i = 0; i < nslots; i++)
      MulMod(crt1[i], crt[i], crtCoeffs[i], factors[i]);

    evalTree(H, crtTree, crt1, 0, nslots);
  }
  HELIB_TIMER_STOP;
}

template <typename type>
void PAlgebraModDerived<type>::mapToFt(RX& w,
                                       const RX& G,
                                       long t,
                                       const RX* rF1) const
{
  if (isDryRun()) {
    w = RX::zero();
    return;
  }
  long i = zMStar.indexOfRep(t);
  if (i < 0) {
    clear(w);
    return;
  }

  if (rF1 == nullptr) { // Compute the representation "from scratch"
    // special case
    if (G == factors[i]) {
      SetX(w);
      return;
    }

    // special case
    if (deg(G) == 1) {
      w = -ConstTerm(G);
      return;
    }

    // the general case: currently only works when r == 1
    assertEq(r, 1l, "Bad Hensel lifting value in general case: r is not 1");

    REBak bak;
    bak.save();
    RE::init(factors[i]); // work with the extension field GF_p[X]/Ft(X)
    REX Ga;
    conv(Ga, G); // G as a polynomial over the extension field

    vec_RE roots;
    FindRoots(roots, Ga); // Find roots of G in this field
    RE* first = &roots[0];
    RE* last = first + lsize(roots);
    RE* smallest = std::min_element(
        first,
        last,
        static_cast<bool (*)(const RE&, const RE&)>(less_than));
    // make a canonical choice
    w = rep(*smallest);
    return;
  }
  // if rF1 is set, then use it instead, setting w = rF1(X^t) mod Ft(X)
  RXModulus Ft(factors[i]);
  //  long tInv = InvMod(t,m);
  RX X2t = PowerXMod(t, Ft);  // X2t = X^t mod Ft
  w = CompMod(*rF1, X2t, Ft); // w = F1(X2t) mod Ft

  /* Debugging sanity-check: G(w)=0 in the extension field (Z/2Z)[X]/Ft(X)
  RE::init(factors[i]);
  REX Ga;
  conv(Ga, G); // G as a polynomial over the extension field
  RE ra;
  conv(ra, w);         // w is an element in the extension field
  eval(ra,Ga,ra);  // ra = Ga(ra)
  if (!IsZero(ra)) {// check that Ga(w)=0 in this extension field
    cout << "rF1(X^t) mod Ft(X) != root of G mod Ft, t=" << t << endl;
    exit(0);
  }*******************************************************************/
}

template <typename type>
void PAlgebraModDerived<type>::mapToSlots(MappingData<type>& mappingData,
                                          const RX& G) const
{
  assertTrue<InvalidArgument>(
      deg(G) > 0,
      "Polynomial G is constant (has degree less than one)");
  assertEq(zMStar.getOrdP() % deg(G),
           0l,
           "Degree of polynomial G does not divide zMStar.getOrdP()");
  assertTrue<InvalidArgument>(static_cast<bool>(LeadCoeff(G) == 1l),
                              "Polynomial G is not monic");
  mappingData.G = G;
  mappingData.degG = deg(mappingData.G);
  long d = deg(G);
  long ordp = zMStar.getOrdP();

  long nSlots = zMStar.getNSlots();
  long m = zMStar.getM();

  resize(mappingData.maps, nSlots);

  mapToF1(mappingData.maps[0], mappingData.G); // mapping from base-G to base-F1
  for (long i = 1; i < nSlots; i++)
    mapToFt(mappingData.maps[i],
            mappingData.G,
            zMStar.ith_rep(i),
            &(mappingData.maps[0]));

  // create matrices to streamline CompMod operations
  resize(mappingData.matrix_maps, nSlots);
  for (long i : range(nSlots)) {
    mat_R& mat = mappingData.matrix_maps[i];
    mat.SetDims(d, ordp);
    RX pow;
    pow = 1;
    for (long j : range(d)) {
      VectorCopy(mat[j], pow, ordp);
      if (j < d - 1)
        MulMod(pow, pow, mappingData.maps[i], factors[i]);
    }
  }

  REBak bak;
  bak.save();
  RE::init(mappingData.G);
  mappingData.contextForG.save();

  if (deg(mappingData.G) == 1)
    return;

  resize(mappingData.rmaps, nSlots);

  if (G == factors[0]) {
    // an important special case

    for (long i = 0; i < nSlots; i++) {
      long t = zMStar.ith_rep(i);
      long tInv = NTL::InvMod(t, m);

      RX ct_rep;
      PowerXMod(ct_rep, tInv, G);

      RE ct;
      conv(ct, ct_rep);

      REX Qi;
      SetCoeff(Qi, 1, 1);
      SetCoeff(Qi, 0, -ct);

      mappingData.rmaps[i] = Qi;
    }
  } else {
    // the general case: currently only works when r == 1

    assertEq(r, 1l, "Bad Hensel lifting value in general case: r is not 1");

    vec_REX FRts;
    for (long i = 0; i < nSlots; i++) {
      // We need to lift Fi from R[Y] to (R[X]/G(X))[Y]
      REX Qi;
      long t, tInv = 0;

      if (i == 0) {
        conv(Qi, factors[i]);
        FRts = EDF(Qi, FrobeniusMap(Qi), deg(Qi) / deg(G));
        // factor Fi over GF(p)[X]/G(X)
      } else {
        t = zMStar.ith_rep(i);
        tInv = NTL::InvMod(t, m);
      }

      // need to choose the right factor, the one that gives us back X
      long j;
      for (j = 0; j < lsize(FRts); j++) {
        // lift maps[i] to (R[X]/G(X))[Y] and reduce mod j'th factor of Fi

        REX FRtsj;
        if (i == 0)
          FRtsj = FRts[j];
        else {
          REX X2tInv = PowerXMod(tInv, FRts[j]);
          IrredPolyMod(FRtsj, X2tInv, FRts[j]);
        }

        // FRtsj is the jth factor of factors[i] over the extension field.
        // For j > 0, we save some time by computing it from the jth factor
        // of factors[0] via a minimal polynomial computation.

        REX GRti;
        conv(GRti, mappingData.maps[i]);
        GRti %= FRtsj;

        if (IsX(rep(ConstTerm(GRti)))) { // is GRti == X?
          Qi = FRtsj;                    // If so, we found the right factor
          break;
        } // If this does not happen then move to the next factor of Fi
      }

      assertTrue(j < lsize(FRts),
                 "Cannot find the right factor Qi. Loop did not "
                 "terminate before visiting all elements");
      mappingData.rmaps[i] = Qi;
    }
  }
}

template <typename type>
void PAlgebraModDerived<type>::decodePlaintext(
    std::vector<RX>& alphas,
    const RX& ptxt,
    const MappingData<type>& mappingData) const
{
  long nSlots = zMStar.getNSlots();
  if (isDryRun()) {
    alphas.assign(nSlots, RX::zero());
    return;
  }

  // First decompose p into CRT components
  std::vector<RX> CRTcomps(nSlots); // allocate space for CRT component
  CRT_decompose(CRTcomps, ptxt);    // CRTcomps[i] = p mod facors[i]

  if (mappingData.degG == 1) {
    alphas = CRTcomps;
    return;
  }

  resize(alphas, nSlots);

  REBak bak;
  bak.save();
  mappingData.contextForG.restore();

  for (long i = 0; i < nSlots; i++) {
    REX te;
    conv(te, CRTcomps[i]); // lift i'th CRT component to mod G(X)
    te %= mappingData
              .rmaps[i]; // reduce CRTcomps[i](Y) mod Qi(Y), over (Z_2[X]/G(X))

    // the free term (no Y component) should be our answer (as a poly(X))
    alphas[i] = rep(ConstTerm(te));
  }
}

template <typename type>
void PAlgebraModDerived<type>::buildLinPolyCoeffs(
    std::vector<RX>& C,
    const std::vector<RX>& L,
    const MappingData<type>& mappingData) const
{
  REBak bak;
  bak.save();
  mappingData.contextForG.restore();

  long d = RE::degree();
  long p = zMStar.getP();

  assertEq(lsize(L), d, "Vector L size is different than RE::degree()");

  vec_RE LL;
  resize(LL, d);

  for (long i = 0; i < d; i++)
    conv(LL[i], L[i]);

  vec_RE CC;
  ::helib::buildLinPolyCoeffs(CC, LL, p, r);

  resize(C, d);
  for (long i = 0; i < d; i++)
    C[i] = rep(CC[i]);
}

// code for generating mask tables
// the tables are generated "on demand"

template <typename type>
void PAlgebraModDerived<type>::genMaskTable()
{
  // This is only called by the constructor, which has already
  // set the zz_p context and the crtTable
  resize(maskTable, zMStar.numOfGens());
  for (long i = 0; i < (long)zMStar.numOfGens(); i++) {
    long ord = zMStar.OrderOf(i);
    resize(maskTable[i], ord + 1);
    maskTable[i][ord] = 0;
    for (long j = ord - 1; j >= 1; j--) {
      // initialize mask that is 1 whenever the ith coordinate is at least j
      // Note: maskTable[i][0] = constant 1, maskTable[i][ord] = constant 0
      maskTable[i][j] = maskTable[i][j + 1];
      for (long k = 0; k < (long)zMStar.getNSlots(); k++) {
        if (zMStar.coordinate(i, k) == j) {
          add(maskTable[i][j], maskTable[i][j], crtTable[k]);
        }
      }
    }
    maskTable[i][0] = 1;
  }
}

// code for generating crt tables
// the tables are generated "on demand"

template <typename type>
void PAlgebraModDerived<type>::genCrtTable()
{
  // This is only called by the constructor, which has already
  // set the zz_p context

  long nslots = zMStar.getNSlots();
  resize(crtTable, nslots);
  for (long i = 0; i < nslots; i++) {
    RX allBut_i = PhimXMod / factors[i]; // = \prod_{j \ne i }Fj
    allBut_i *= crtCoeffs[i]; // = 1 mod Fi and = 0 mod Fj for j \ne i
    crtTable[i] = allBut_i;
  }

  buildTree(crtTree, 0, nslots);
}

template <typename type>
void PAlgebraModDerived<type>::buildTree(std::shared_ptr<TNode<RX>>& res,
                                         long offset,
                                         long extent) const
{
  if (extent == 1)
    res = buildTNode<RX>(nullTNode<RX>(), nullTNode<RX>(), factors[offset]);
  else {
    long half = extent / 2;
    std::shared_ptr<TNode<RX>> left, right;
    buildTree(left, offset, half);
    buildTree(right, offset + half, extent - half);
    RX data = left->data * right->data;
    res = buildTNode<RX>(left, right, data);
  }
}

template <typename type>
void PAlgebraModDerived<type>::evalTree(RX& res,
                                        std::shared_ptr<TNode<RX>> tree,
                                        const std::vector<RX>& crt1,
                                        long offset,
                                        long extent) const
{
  if (extent == 1)
    res = crt1[offset];
  else {
    long half = extent / 2;
    RX lres, rres;
    evalTree(lres, tree->left, crt1, offset, half);
    evalTree(rres, tree->right, crt1, offset + half, extent - half);
    RX tmp1, tmp2;
    mul(tmp1, lres, tree->right->data);
    mul(tmp2, rres, tree->left->data);
    add(tmp1, tmp1, tmp2);
    res = tmp1;
  }
}

// Explicit instantiation

template class PAlgebraModDerived<PA_GF2>;
template class PAlgebraModDerived<PA_zz_p>;

} // namespace helib
