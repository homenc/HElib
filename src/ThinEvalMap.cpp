#include <NTL/BasicThreadPool.h>

//=================== ThinEvalMap.h ================//

#include "EncryptedArray.h"
#include "matmul.h"

//! @class ThinEvalMap
//! @brief Class that provides the functionality for the 
//! linear transforms used in boostrapping.
//! The constructor is invoked with three arguments:
//!   - an EncryptedArray object ea
//!   - an integer vector mvec
//!   - a boolean flag invert
//! The mvec vector specifies the factorization of m
//! to use in the "powerful basis" decomposition.
//!
//! If the invert flag is false, the forward transformation is
//! used.  This transformation views the slots as being packed
//! with powerful-basis coefficients and performs a multi-point
//! polynomial evaluation.  This is the second transformation
//! used in bootstrapping.
//!
//! If invert flag is true, the inverse transformation is used.
//! In addition, the current implementation folds into the
//! inverse transformation a transformation that moves
//! the coefficients in each slot into a normal-basis representation,
//! which helps with the unpacking procedure.
//! 
//! The constructor precomputes certain values, and the linear
//! transformation itself is effected using the apply method.
//!
//! Note that the factorization in mvec must correspond to the
//! generators used in PAlgebra.  The best way to ensure this is
//! to directly use the output of the program in  params.cpp: that
//! program computes values for mvec (to be used here), and gens
//! and ords (to be used in initialize the FHEcontext).

class ThinEvalMap {
private:
  const EncryptedArray& ea;
  bool invert;   // apply transformation in inverser order?
  long nfactors; // how many factors of m
  NTL::Vec<std::unique_ptr<MatMulExecBase>>  matvec; // regular matrices

public:
  ThinEvalMap(const EncryptedArray& _ea, 
          bool minimal, 
          const Vec<long>& mvec, 
          bool _invert,
          bool build_cache);

  void upgrade();
  void apply(Ctxt& ctxt) const;
};

//=================== ThinEvalMap.cpp ================//

//FIXME
//#include "thinEvalMap.h"

// to get TraceMap
#include <NTL/lzz_pXFactoring.h>
#include <NTL/GF2XFactoring.h>


// needed to make generic programming work

void
RelaxedInv(Mat<zz_p>& x, const Mat<zz_p>& a)
{
   relaxed_inv(x, a);
}

void
RelaxedInv(Mat<GF2>& x, const Mat<GF2>& a)
{
   inv(x, a);
}


void TraceMap(GF2X& w, const GF2X& a, long d, const GF2XModulus& F, 
              const GF2X& b)

{
   if (d < 0) LogicError("TraceMap: bad args");

   GF2X y, z, t;

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
      }
      else if ((d & 1) == 0) {
         Comp2Mod(z, t, z, y, z, F);
         add(y, t, y);
      }
      else if (IsZero(w)) {
         w = y;
         Comp2Mod(z, t, z, y, z, F);
         add(y, t, y);
      }
      else {
         Comp3Mod(z, t, w, z, y, w, z, F);
         add(w, w, y);
         add(y, t, y);
      }

      d = d >> 1;
   }
}



// Forward declerations
static MatMul1D*
buildThinStep1Matrix(const EncryptedArray& ea, shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps, long dim, long cofactor);
static MatMul1D*
buildThinStep2Matrix(const EncryptedArray& ea, shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps, long dim, long cofactor,
                 bool invert, bool inflate=false);
static void
init_representatives(Vec<long>& representatives, long dim, 
                     const Vec<long>& mvec, const PAlgebra& zMStar);


// Constructor: initializing tables for the evaluation-map transformations

ThinEvalMap::ThinEvalMap(const EncryptedArray& _ea, 
                 bool minimal,
                 const Vec<long>& mvec, 
                 bool _invert,
                 bool build_cache)

  : ea(_ea), invert(_invert)
{
  const FHEcontext& context = ea.getContext();
  const PAlgebra& zMStar = ea.getPAlgebra();
  
  long p = zMStar.getP();
  long d = zMStar.getOrdP();

  // FIXME: we should check that ea was initilized with 
  // G == factors[0], but this is a slight pain to check
  // currently

  // NOTE: this code is derived from a more general setting, and
  // could certainly be greatly simplified

  nfactors = mvec.length();
  assert(nfactors > 0);

  for (long i = 0; i < nfactors; i++)
    for (long j = i+1; j < nfactors; j++)
      assert(GCD(mvec[i], mvec[j]) == 1);

  long m = computeProd(mvec);
  assert(m == long(zMStar.getM()));

  Vec<long> phivec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)  phivec[i] = phi_N(mvec[i]);
  long phim = computeProd(phivec);

  Vec<long> dprodvec(INIT_SIZE, nfactors+1);
  dprodvec[nfactors] = 1;
  
  for (long i = nfactors-1; i >= 0; i--)
    dprodvec[i] = dprodvec[i+1] *
      multOrd(PowerMod(p % mvec[i], dprodvec[i+1], mvec[i]), mvec[i]);

  Vec<long> dvec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    dvec[i] = dprodvec[i] / dprodvec[i+1];

  long nslots = phim/d;
  assert(d == dprodvec[0]);
  assert(nslots == long(zMStar.getNSlots()));

  long inertPrefix = 0;
  for (long i = 0; i < nfactors && dvec[i] == 1; i++) {
    inertPrefix++;
  }

  if (inertPrefix != nfactors-1)
    Error("ThinEvalMap: case not handled: bad inertPrefix");

  Vec< Vec<long> > local_reps(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    init_representatives(local_reps[i], i, mvec, zMStar);

  Vec<long> crtvec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++) 
    crtvec[i] = (m/mvec[i]) * InvMod((m/mvec[i]) % mvec[i], mvec[i]);

  Vec<long> redphivec(INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    redphivec[i] = phivec[i]/dvec[i];

  CubeSignature redphisig(redphivec);

  Vec< shared_ptr<CubeSignature> > sig_sequence;
  sig_sequence.SetLength(nfactors+1);
  sig_sequence[nfactors] = shared_ptr<CubeSignature>(new CubeSignature(phivec));

  Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors-1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = 
      shared_ptr<CubeSignature>(new CubeSignature(reduced_phivec));
  }

  matvec.SetLength(nfactors);

  if (invert) {
     long dim = nfactors - 1;
     unique_ptr<MatMul1D> mat1_data;
     mat1_data.reset(buildThinStep1Matrix(ea, sig_sequence[dim],
		     local_reps[dim], dim, m/mvec[dim]));
     matvec[dim].reset(new MatMul1DExec(*mat1_data, minimal));
  }
  else {
     long dim = nfactors - 1;
     unique_ptr<MatMul1D> mat1_data;
     mat1_data.reset(buildThinStep2Matrix(ea, sig_sequence[dim],
		     local_reps[dim], dim, m/mvec[dim], invert, /*inflate=*/true));
     matvec[dim].reset(new MatMul1DExec(*mat1_data, minimal));
  }

  for (long dim=nfactors-2; dim>=0; --dim) {
    unique_ptr<MatMul1D> mat_data;

    mat_data.reset(buildThinStep2Matrix(ea, sig_sequence[dim], local_reps[dim],
				       dim, m/mvec[dim], invert));
    matvec[dim].reset(new MatMul1DExec(*mat_data, minimal));
  }

  if (build_cache) upgrade();
}

void ThinEvalMap::upgrade()
{
  for (long i = 0; i < matvec.length(); i++)
    matvec[i]->upgrade();
}

// Applying the evaluation (or its inverse) map to a ciphertext
void ThinEvalMap::apply(Ctxt& ctxt) const
{
  if (!invert) { // forward direction
    for (long i = matvec.length()-1; i >= 0; i--)
      matvec[i]->mul(ctxt);
  }
  else {         // inverse transformation
    for (long i = 0; i < matvec.length(); i++)
      matvec[i]->mul(ctxt);
    traceMap(ctxt);
  }
}


static void
init_representatives(Vec<long>& representatives, long dim, 
                     const Vec<long>& mvec, const PAlgebra& zMStar)
{
  assert(dim >= 0 && dim < mvec.length());

  // special case
  if (dim >= LONG(zMStar.numOfGens())) {
    representatives.SetLength(1);
    representatives[0] = 1;
    return;
  }
  
  long m = mvec[dim];
  long D = zMStar.OrderOf(dim);
  long g = InvMod(zMStar.ZmStarGen(dim) % m, m);

  representatives.SetLength(D);
  for (long i = 0; i < D; i++)
    representatives[i] = PowerMod(g, i, m);
}

// The callback interface for the matrix-multiplication routines.

//! \cond FALSE (make doxygen ignore these classes)
template<class type> class ThinStep2Matrix : public MatMul1D_derived<type> 
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  shared_ptr<CubeSignature> sig;
  long dim;
  Mat<RX> A;

public:
  // constructor
  ThinStep2Matrix(const EncryptedArray& _ea,
              shared_ptr<CubeSignature> _sig, const Vec<long>& reps,
              long _dim, long cofactor, bool invert, bool inflate)
    : base_ea(_ea), sig(_sig), dim(_dim)
  {
    long sz = sig->getDim(dim);
    assert(sz == reps.length());

    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();
    long d = deg(G);

    Vec<RX> points(INIT_SIZE, sz);
    for (long j = 0; j < sz; j++) {
      points[j] = RX(reps[j]*cofactor, 1) % G;
      if (inflate) points[j] = PowerMod(points[j], d, G);
    }

    A.SetDims(sz, sz);
    for (long j = 0; j < sz; j++)
      A[0][j] = 1;

    for (long i = 1; i < sz; i++)
      for (long j = 0; j < sz; j++)
	A[i][j] = (A[i-1][j] * points[j]) % G;

    if (invert) {
      REBak ebak; ebak.save(); ea.restoreContextForG();

      mat_RE A1, A2;
      conv(A1, A);

      long p = _ea.getAlMod().getZMStar().getP();
      long r = _ea.getAlMod().getR();

      ppInvert(A2, A1, p, r);
      conv(A, A2);
    }
  }

  bool get(RX& out, long i, long j, long k) const override {
    out = A[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static MatMul1D*
buildThinStep2Matrix(const EncryptedArray& ea, shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps, long dim, long cofactor,
                 bool invert, bool inflate)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: 
    return new ThinStep2Matrix<PA_GF2>(ea, sig, reps, dim, cofactor, invert, inflate);

  case PA_zz_p_tag: 
    return new ThinStep2Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor, invert, inflate);

  default: return 0;
  }
}

template<class type> class ThinStep1Matrix : public MatMul1D_derived<type> 
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  shared_ptr<CubeSignature> sig;
  long dim;
  Mat<RX> A_deflated;

public:
  // constructor
  ThinStep1Matrix(const EncryptedArray& _ea, shared_ptr<CubeSignature> _sig,
              const Vec<long>& reps, long _dim, long cofactor)
    : base_ea(_ea), sig(_sig), dim(_dim)
  {
    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    const RXModulus G(ea.getG());
    long d = deg(G);

    long p = _ea.getAlMod().getZMStar().getP();
    long r = _ea.getAlMod().getR();

    long sz = sig->getDim(dim);
    assert(sz == reps.length());
    assert(dim == sig->getNumDims() - 1);
    assert(sig->getSize() == ea.size());

    // so sz == phi(m_last)/d, where d = deg(G) = order of p mod m

    Vec<RX> points(INIT_SIZE, sz);
    for (long j = 0; j < sz; j++) 
      points[j] = RX(reps[j]*cofactor, 1) % G;

    Mat<RX> AA(INIT_SIZE, sz*d, sz);
    for (long j = 0; j < sz; j++)
      AA[0][j] = 1;

    for (long i = 1; i < sz*d; i++)
      for (long j = 0; j < sz; j++)
	AA[i][j] = (AA[i-1][j] * points[j]) % G;

    Mat<mat_R> A;
    A.SetDims(sz, sz);
    for (long i = 0; i < sz; i++)
      for (long j = 0; j < sz; j++) {
	A[i][j].SetDims(d, d);
	for (long k = 0; k < d; k++)
	  VectorCopy(A[i][j][k], AA[i*d + k][j], d);
      }

    // if (invert) {
    // NOTE: this version is only used for the inverse matrix
    mat_R A1, A2;
    A1.SetDims(sz*d, sz*d);
    for (long i = 0; i < sz*d; i++)
      for (long j = 0; j < sz*d; j++)
	A1[i][j] = A[i/d][j/d][i%d][j%d];

    ppInvert(A2, A1, p, r);

    for (long i = 0; i < sz*d; i++)
      for (long j = 0; j < sz*d; j++)
	A[i/d][j/d][i%d][j%d] = A2[i][j];
    // }

    A_deflated.SetDims(sz, sz);
    vec_R v, w;
    v.SetLength(d);
    w.SetLength(d);

    RX h;  // set h = X^p mod G
    PowerXMod(h, p, G);

    Vec<R> trace_vec;
    trace_vec.SetLength(2*d-1);
    // set trace_vec[i] = trace(X^i mod G) 
    for (long i = 0; i < 2*d-1; i++) {
      RX trace_val;
      TraceMap(trace_val, (RX(i, 1) % G), d, G, h);
      assert(deg(trace_val) <= 0);
      trace_vec[i] = ConstTerm(trace_val);
    }

    Mat<R> trace_mat; 
    trace_mat.SetDims(d, d);
    // set trace_mat[i][j] = trace(X^{i+j} mod G)
    for (long i = 0; i < d; i++)
      for (long j = 0; j < d; j++)
         trace_mat[i][j] = trace_vec[i+j];

    Mat<R> trace_mat_inv;
    RelaxedInv(trace_mat_inv, trace_mat);

    for (long i = 0 ; i < sz; i++)
      for (long j = 0; j < sz; j++) {
         for (long i1 = 0; i1 < d; i1++)
            v[i1] = A[i][j][i1][0];
         mul(w, v, trace_mat_inv);
         conv(A_deflated[i][j], w);
      }
  } // constructor

  bool get(RX& out, long i, long j, long k) const override {
    out = A_deflated[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static MatMul1D*
buildThinStep1Matrix(const EncryptedArray& ea, shared_ptr<CubeSignature> sig,
                 const Vec<long>& reps, long dim, long cofactor)
{
  switch (ea.getTag()) {
  case PA_GF2_tag: 
    return new ThinStep1Matrix<PA_GF2>(ea, sig, reps, dim, cofactor);

  case PA_zz_p_tag: 
    return new ThinStep1Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor);

  default: return 0;
  }
}
//! \endcond

//================== thinRecrypt.h =================//

//FIXME
//#include "NumbTh.h"

class  PAlgebraMod;
class  EncryptedArray;
class  ThinEvalMap;
class  FHEcontext;
class  FHEPubKey;


//! @class RecryptData
//! @brief A structure to hold recryption-related data inside the FHEcontext
class ThinRecryptData {
public:
  //! default Hamming weight of recryption key
  static const long defSkHwt=100;

  //! Some data members that are only used for I/O
  Vec<long> mvec;     //! partition of m into co-prime factors
  long hwt;           //! Hamming weight of recryption secret-key
  bool conservative;  //! flag for choosing more conservatice parameters

  //! skey encrypted wrt space p^{e-e'+r}
  long e, ePrime;

  //! Hamming weight of recryption secret key
  long skHwt;

  //! an optimization parameter
  double alpha;

  //! for plaintext space p^{e-e'+r}
  PAlgebraMod *alMod;

  //! for plaintext space p^{e-e'+r}
  EncryptedArray *ea;

  bool build_cache;


  //! linear maps
  ThinEvalMap *coeffToSlot, *slotToCoeff;

  ThinRecryptData() {
    hwt=0; conservative=false; e=ePrime=0; alpha=0.0;
    alMod=NULL; ea=NULL; coeffToSlot=NULL; slotToCoeff=NULL; 
    build_cache = false;
  }
  ~ThinRecryptData();

  //! Initialize the recryption data in the context
  void init(const FHEcontext& context, const Vec<long>& mvec_,
            long t=0/*min Hwt for sk*/, 
            bool consFlag=false,
            bool build_cache=false,
            bool minimal=false);

  bool operator==(const ThinRecryptData& other) const;
  bool operator!=(const ThinRecryptData& other) const {
    return !(operator==(other));
  }
};

//================== thinRecrypt.cpp ===============//


ThinRecryptData::~ThinRecryptData()
{
  if (alMod!=NULL)     delete alMod;
  if (ea!=NULL)        delete ea;
  if (coeffToSlot!=NULL)  delete coeffToSlot;
  if (slotToCoeff!=NULL) delete slotToCoeff;
}

bool ThinRecryptData::operator==(const ThinRecryptData& other) const
{
  if (mvec != other.mvec) return false;
  if (hwt != other.hwt) return false;
  if (conservative != other.conservative) return false;

  return true;
}

// FIXME: these are copied from recryption.cpp.  Should be
// consolidated.

// Make every entry of vec divisible by p^e by adding/subtracting multiples
// of p^r and q, while keeping the added multiples small.
template<class VecInt>
long makeDivisible(VecInt& vec, long p2e, long p2r, long q, double alpha);

static inline double pow(long a, long b) {return pow(double(a), double(b));}



/** Computing the recryption parameters
 *
 * To get the smallest possible value of e-e', the params need to satisfy:
 *  (p^e +1)/4 =>
 *       max { (t+1)( 1+ (alpha/2)*(p^e/p^{ceil(log_p(t+2))}) ) + noise      }
 *           { (t+1)( 1+ ((1-alpha)/2)*(p^e/p^{ceil(log_p(t+2))}) +p^r/2) +1 },
 *
 * where noise is taken to be twice the mod-switching additive term, namely
 * noise = p^r *sqrt((t+1)*phi(m)/3). Denoting rho=(t+1)/p^{ceil(log_p(t+2))}
 * (and ignoring fome +1 terms), this is equivalent to:
 *
 *   p^e > max { 4(t+noise)/(1-2*alpha*rho), 2(t+1)p^r/(1-2(1-alpha)rho) }.
 *
 * We first compute the optimal value for alpha (which must be in [0,1]),
 * that makes the two terms in the max{...} as close as possible, and
 * then compute the smallest value of e satisfying this constraint.
 *
 * If this value is too big then we try again with e-e' one larger,
 * which means that rho is a factor of p smaller.
 */

// Some convenience functions
static double lowerBound1(long p, long r, long ePrime, long t,
			  double alpha, double noise)
{
  return (t+1)*(1+ alpha*pow(p,r+ePrime-1)/2)+noise;
}
static double lowerBound2(long p, long r, long ePrime, long t, double alpha)
{
  return (t+1)*(1+ (1-alpha)*pow(p,r+ePrime-1)/2 + pow(p,r)/2)+1;
}

static void setAlphaE(double& alpha, long& e, double rho, double gamma,
		      double noise, double logp, long p2r, long t)
{
  alpha = (1 +gamma*(2*rho-1))/(2*rho*(1+gamma));
  if (alpha<0) alpha=0;
  else if (alpha>1) alpha=1;

  if (alpha<1) {
    double ratio = 4*(t+noise)/(1-2*alpha*rho);
    e = floor(1+ log(ratio)/logp);
  }
  else
    e = floor(1+ log(2*(t+1)*p2r)/logp);
}

// The main method
void ThinRecryptData::init(const FHEcontext& context, const Vec<long>& mvec_,
		       long t, bool consFlag, bool build_cache_, bool minimal)
{
  if (alMod != NULL) { // were we called for a second time?
    cerr << "@Warning: multiple calls to RecryptData::init\n";
    return;
  }
  assert(computeProd(mvec_) == (long)context.zMStar.getM()); // sanity check

  // Record the arguments to this function
  mvec = mvec_;
  conservative = consFlag;
  build_cache = build_cache_;

  if (t <= 0) t = defSkHwt+1; // recryption key Hwt
  hwt = t;
  long p = context.zMStar.getP();
  long phim = context.zMStar.getPhiM();
  long r = context.alMod.getR();
  long p2r = context.alMod.getPPowR();
  double logp = log((double)p);

  double noise = p2r * sqrt((t+1)*phim/3.0);
  double gamma = 2*(t+noise)/((t+1)*p2r); // ratio between numerators

  long logT = ceil(log((double)(t+2))/logp); // ceil(log_p(t+2))
  double rho = (t+1)/pow(p,logT);

  if (!conservative) {   // try alpha, e with this "aggresive" setting
    setAlphaE(alpha, e, rho, gamma, noise, logp, p2r, t);
    ePrime = e -r +1 -logT;

    // If e is too large, try again with rho/p instead of rho
    long bound = (1L << (context.bitsPerLevel-1)); // halfSizePrime/2
    if (pow(p,e) > bound) { // try the conservative setting instead
      cerr << "* p^e="<<pow(p,e)<<" is too big (bound="<<bound<<")\n";
      conservative = true;
    }
  }
  if (conservative) { // set alpha, e with a "conservative" rho/p
    setAlphaE(alpha, e, rho/p, gamma, noise, logp, p2r, t);
    ePrime = e -r -logT;
  }

  // Compute highest key-Hamming-weight that still works (not more than 256)
  double qOver4 = (pow(p,e)+1)/4;
  for (t-=10; qOver4>=lowerBound2(p,r,ePrime,t,alpha)
	 &&  qOver4>=lowerBound1(p,r,ePrime,t,alpha,noise) && t<257; t++);
  skHwt = t-1;

  // First part of Bootstrapping works wrt plaintext space p^{r'}
  alMod = new PAlgebraMod(context.zMStar, e-ePrime+r);
  ea = new EncryptedArray(context, *alMod);
         // Polynomial defaults to F0, PAlgebraMod explicitly given

  coeffToSlot = new ThinEvalMap(*ea, minimal, mvec, true, build_cache);
  slotToCoeff = new ThinEvalMap(*context.ea, minimal, mvec, false, build_cache);
}


// Extract digits from thinly packed slots
void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime)
{
  FHE_TIMER_START;

  Ctxt unpacked(ctxt);
  unpacked.cleanUp();

  vector<Ctxt> scratch;


  // Step 2: extract the digits top-1,...,0 from the slots of unpacked[i]
  long p = ctxt.getContext().zMStar.getP();
  long p2r = power_long(p,r);
  long topHigh = botHigh + r-1;

#ifdef DEBUG_PRINTOUT
  cerr << "+ After unpack ";
  decryptAndPrint(cerr, unpacked, *dbgKey, *dbgEa, printFlag);
  cerr << "    extracting "<<(topHigh+1)<<" digits\n";
#endif

  if (p==2 && r>2)
    topHigh--; // For p==2 we sometime get a bit for free

  if (topHigh<=0) { // extracting LSB = no-op
    scratch.assign(1, unpacked);
  } else {          // extract digits topHigh...0, store them in scratch
    extractDigits(scratch, unpacked, topHigh+1);
  }

  // set upacked = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
  if (topHigh >= LONG(scratch.size())) {
    topHigh = scratch.size() -1;
    cerr << " @ suspect: not enough digits in extractDigitsPacked\n";
  }

  unpacked = scratch[topHigh];
  for (long j=topHigh-1; j>=botHigh; --j) {
    unpacked.multByP();
    unpacked += scratch[j];
  }
  if (p==2 && botHigh>0)   // For p==2, subtract also the previous bit
    unpacked += scratch[botHigh-1];
  unpacked.negate();

  if (r>ePrime) {          // Add in digits from the bottom part, if any
    long topLow = r-1 - ePrime;
    Ctxt tmp = scratch[topLow];
    for (long j=topLow-1; j>=0; --j) {
      tmp.multByP();
      tmp += scratch[j];
    }
    if (ePrime>0)
      tmp.multByP(ePrime); // multiply by p^e'
    unpacked += tmp;
  }
  unpacked.reducePtxtSpace(p2r); // Our plaintext space is now mod p^r

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before repack ";
  decryptAndPrint(cerr, unpacked[0], *dbgKey, *dbgEa, printFlag);
#endif

  ctxt = unpacked;

}


// Hack to get at private fields of public key
struct FHEPubKeyHack { // The public key
  const FHEcontext& context; // The context

  //! @var Ctxt pubEncrKey
  //! The public encryption key is an encryption of 0,
  //! relative to the first secret key
  Ctxt pubEncrKey;

  std::vector<long> skHwts; // The Hamming weight of the secret keys
  std::vector<KeySwitch> keySwitching; // The key-switching matrices

  // The keySwitchMap structure contains pointers to key-switching matrices
  // for re-linearizing automorphisms. The entry keySwitchMap[i][n] contains
  // the index j such that keySwitching[j] is the first matrix one needs to
  // use when re-linearizing s_i(X^n). 
  std::vector< std::vector<long> > keySwitchMap;

  NTL::Vec<int> KS_strategy; // NTL Vec's support I/O, which is
                             // more convenient

  // bootstrapping data

  long recryptKeyID; // index of the bootstrapping key
  Ctxt recryptEkey;  // the key itself, encrypted under key #0

};

#define PRINT_LEVELS

// bootstrap a ciphertext to reduce noise
void thinRecrypt(Ctxt &ctxt, const ThinRecryptData& rcData)
// FIXME:
// * rcData is supposed to be in the context, but for now we can just 
//   pass it in
// * recryptKeyID is stored in the public key, but we can get it for now
//   as the retval of genRecryptData
{
  FHE_TIMER_START;

  // FIXME: Some serious hacks to get access to private members
  const FHEPubKey& pk = ctxt.getPubKey();
  long recryptKeyID = ((FHEPubKeyHack&)pk).recryptKeyID;
  const Ctxt& recryptEkey = ((FHEPubKeyHack&)pk).recryptEkey;
  const std::vector<long>& skHwts =  ((FHEPubKeyHack&)pk).skHwts;



  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty()) return;

// FIXME: the following code refers to private members
// of Ctxt...
#if 0
  if (ctxt.parts.size()==1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i=0; i<poly.rep.length(); i++)
      poly[i] = to_ZZ( rem(poly[i],ptxtSpace) );
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return;
  }
#endif

  assert(recryptKeyID>=0); // check that we have bootstrapping data

  long p = ctxt.getContext().zMStar.getP();
  long r = ctxt.getContext().alMod.getR();
  long p2r = ctxt.getContext().alMod.getPPowR();

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  long e = rcData.e;
  long ePrime = rcData.ePrime;
  long p2ePrime = power_long(p,ePrime);
  long q = power_long(p,e)+1;
  assert(e>=r);

#ifdef DEBUG_PRINTOUT
  cerr << "reCrypt: p="<<p<<", r="<<r<<", e="<<e<<" ePrime="<<ePrime
       << ", q="<<q<<endl;
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assert(p2r % ptxtSpace == 0);

#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "init");
#endif

  // Move the slots to powerful-basis coefficients
  FHE_NTIMER_START(AAA_slotToCoeff);
  rcData.slotToCoeff->apply(ctxt);
  FHE_NTIMER_STOP(AAA_slotToCoeff);

#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after slotToCoeff");
#endif

  FHE_NTIMER_START(AAA_bootKeySwitch);

  // Make sure that this ciphertxt is in canonical form
  if (!ctxt.inCanonicalForm()) ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / ctxt.getContext().specialPrimes; // set minus
  if (s.card()>2) { // leave only bottom two primes
    long frst = s.first();
    long scnd = s.next(frst);
    IndexSet s2(frst,scnd);
    s.retain(s2); // retain only first two primes
  }
  ctxt.modDownToSet(s);

  // key-switch to the bootstrapping key
  ctxt.reLinearize(recryptKeyID);

  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  vector<ZZX> zzParts; // the mod-switched parts, in ZZX format
  double noise = ctxt.rawModSwitch(zzParts, q);
  noise = sqrt(noise);

  // Add multiples of p2r and q to make the zzParts divisible by p^{e'}
  long maxU=0;
  for (long i=0; i<(long)zzParts.size(); i++) {
    // make divisible by p^{e'}
    long newMax = makeDivisible(zzParts[i].rep, p2ePrime, p2r, q,
				rcData.alpha);
    zzParts[i].normalize();   // normalize after working directly on the rep
    if (maxU < newMax)  maxU = newMax;
  }

  // Check that the estimated noise is still low
  if (noise + maxU*p2r*(skHwts[recryptKeyID]+1) > q/2) 
    cerr << " * noise/q after makeDivisible = "
	 << ((noise + maxU*p2r*(skHwts[recryptKeyID]+1))/q) << endl;

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  // Multiply the post-processed cipehrtext by the encrypted sKey
#ifdef DEBUG_PRINTOUT
  cerr << "+ Before recryption ";
  decryptAndPrint(cerr, recryptEkey, *dbgKey, *dbgEa, printFlag);
#endif

  double p0size = to_double(coeffsL2Norm(zzParts[0]));
  double p1size = to_double(coeffsL2Norm(zzParts[1]));
  ctxt = recryptEkey;
  ctxt.multByConstant(zzParts[1], p1size*p1size);
  ctxt.addConstant(zzParts[0], p0size*p0size);

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif
  FHE_NTIMER_STOP(AAA_bootKeySwitch);

#ifdef PRINT_LEVELS
   CheckCtxt(ctxt, "after bootKeySwitch");
#endif

  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(AAA_coeffToSlot);
  rcData.coeffToSlot->apply(ctxt);
  FHE_NTIMER_STOP(AAA_coeffToSlot);


#ifdef PRINT_LEVELS
   CheckCtxt(ctxt, "after coeffToSlot");
#endif

#ifdef DEBUG_PRINTOUT
  cerr << "+ After linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  FHE_NTIMER_START(AAA_extractDigitsThin);
  extractDigitsThin(ctxt, e-ePrime, r, ePrime);
  FHE_NTIMER_STOP(AAA_extractDigitsThin);


#ifdef PRINT_LEVELS
   CheckCtxt(ctxt, "after extractDigitsThin");
#endif

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans2 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif
}

#if 0

//================== Test_ThinEvalMap.cpp ==================//

namespace std {} using namespace std;
namespace NTL {} using namespace NTL;

#include <NTL/BasicThreadPool.h>

//FIXME: #include "ThinEvalMap.h"
#include "EvalMap.h"

static bool dry = false; // a dry-run flag
static bool noPrint = false;

void  TestIt(long p, long r, long c, long _k, long w,
             long L, Vec<long>& mvec, 
             Vec<long>& gens, Vec<long>& ords, long useCache)
{
  if (lsize(mvec)<1) { // use default values
    mvec.SetLength(3); gens.SetLength(3); ords.SetLength(3);
    mvec[0] = 7;    mvec[1] = 3;    mvec[2] = 221;
    gens[0] = 3979; gens[1] = 3095; gens[2] = 3760;
    ords[0] = 6;    ords[1] = 2;    ords[2] = -8;
  }
  if (!noPrint)
    cout << "*** TestIt"
       << (dry? " (dry run):" : ":")
       << " p=" << p
       << ", r=" << r
       << ", c=" << c
       << ", k=" << _k
       << ", w=" << w
       << ", L=" << L
       << ", mvec=" << mvec << ", "
       << ", useCache = " << useCache
       << endl;

  setTimersOn();
  setDryRun(false); // Need to get a "real context" to test ThinEvalMap

  // mvec is supposed to include the prime-power factorization of m
  long nfactors = mvec.length();
  for (long i = 0; i < nfactors; i++)
    for (long j = i+1; j < nfactors; j++)
      assert(GCD(mvec[i], mvec[j]) == 1);

  // multiply all the prime powers to get m itself
  long m = computeProd(mvec);
  assert(GCD(p, m) == 1);

  // build a context with these generators and orders
  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);
  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, c);

  if (!noPrint) {
    context.zMStar.printout(); // print structure of Zm* /(p) to cout
    cout << endl;
  }
  long d = context.zMStar.getOrdP();
  long phim = context.zMStar.getPhiM();
  long nslots = phim/d;

  setDryRun(dry); // Now we can set the dry-run flag if desired

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need

  // GG defines the plaintext space Z_p[X]/GG(X)
  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  EncryptedArray ea(context, GG);

  zz_p::init(context.alMod.getPPowR());

  Vec<zz_p> val0;
  random(val0, nslots);

  vector<ZZX> val1;
  val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    val1[i] = conv<ZZX>(conv<ZZ>(rep(val0[i])));
  }

  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, val1);

  resetAllTimers();
  FHE_NTIMER_START(ALL);

  // Compute homomorphically the transformation that takes the
  // coefficients packed in the slots and produces the polynomial
  // corresponding to cube

  if (!noPrint) CheckCtxt(ctxt, "init");

  if (!noPrint) cout << "build ThinEvalMap\n";
  ThinEvalMap map(ea, /*minimal=*/false, mvec, 
    /*invert=*/false, /*build_cache=*/false); 
  // compute the transformation to apply

  if (!noPrint) cout << "apply ThinEvalMap\n";
  if (useCache) map.upgrade();
  map.apply(ctxt); // apply the transformation to ctxt
  if (!noPrint) CheckCtxt(ctxt, "ThinEvalMap");
  if (!noPrint) cout << "check results\n";

  if (!noPrint) cout << "build ThinEvalMap\n";
  ThinEvalMap imap(ea, /*minimal=*/false, mvec, 
    /*invert=*/true, /*build_cache=*/false); 
  // compute the transformation to apply
  if (!noPrint) cout << "apply ThinEvalMap\n";
  if (useCache) imap.upgrade();
  imap.apply(ctxt); // apply the transformation to ctxt
  if (!noPrint) {
    CheckCtxt(ctxt, "ThinEvalMap");
    cout << "check results\n";
  }

#if 1

  /* create dirty version of ctxt */
  Vec<zz_pX> dirty_val0;
  dirty_val0.SetLength(nslots);
  for (long i = 0; i < nslots; i++) {
    random(dirty_val0[i], d);
    SetCoeff(dirty_val0[i], 0, val0[i]);
  }
  
  vector<ZZX> dirty_val1;
  dirty_val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    dirty_val1[i] = conv<ZZX>(dirty_val0[i]);
  }

  Ctxt dirty_ctxt(publicKey);
  ea.encrypt(dirty_ctxt, publicKey, dirty_val1);


  EvalMap dirty_map(ea, /*minimal=*/false, mvec, 
    /*invert=*/false, /*build_cache=*/false); 

  dirty_map.apply(dirty_ctxt);
  imap.apply(dirty_ctxt);
#endif


  vector<ZZX> val2;
  ea.decrypt(ctxt, secretKey, val2);

  if (val1 == val2)
    cout << "ThinEvalMap: GOOD\n";
  else
    cout << "ThinEvalMap: BAD\n";

#if 1
  vector<ZZX> dirty_val2;
  ea.decrypt(dirty_ctxt, secretKey, dirty_val2);

  if (val1 == dirty_val2)
    cout << "ThinEvalMap: GOOD\n";
  else
    cout << "ThinEvalMap: BAD\n";
#endif


  FHE_NTIMER_STOP(ALL);

  if (!noPrint) {
    cout << "\n*********\n";
    printAllTimers();
    cout << endl;
  }
}


/* Usage: Test_EvalMap_x.exe [ name=value ]...
 *  p       plaintext base  [ default=2 ]
 *  r       lifting  [ default=1 ]
 *  c       number of columns in the key-switching matrices  [ default=2 ]
 *  k       security parameter  [ default=80 ]
 *  L       # of levels in the modulus chain  [ default=6 ]
 *  s       minimum number of slots  [ default=0 ]
 *  seed    PRG seed  [ default=0 ]
 *  mvec    use specified factorization of m
 *             e.g., mvec='[5 3 187]'
 *  gens    use specified vector of generators
 *             e.g., gens='[562 1871 751]'
 *  ords    use specified vector of orders
 *             e.g., ords='[4 2 -4]', negative means 'bad'
 */
int main(int argc, char *argv[])
{
  ArgMapping amap;

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  long c=2;
  amap.arg("c", c, "number of columns in the key-switching matrices");
  
  long k=80;
  amap.arg("k", k, "security parameter");

  long L=6;
  amap.arg("L", L, "# of levels in the modulus chain");

  long s=0;
  amap.arg("s", s, "minimum number of slots");

  long seed=0;
  amap.arg("seed", seed, "PRG seed");

  Vec<long> mvec;
  amap.arg("mvec", mvec, "use specified factorization of m", NULL);
  amap.note("e.g., mvec='[7 3 221]'");

  Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", NULL);
  amap.note("e.g., gens='[3979 3095 3760]'");

  Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", NULL);
  amap.note("e.g., ords='[6 2 -8]', negative means 'bad'");

  amap.arg("dry", dry, "a dry-run flag to check the noise");

  long nthreads=1;
  amap.arg("nthreads", nthreads, "number of threads");

  amap.arg("noPrint", noPrint, "suppress printouts");

  long useCache=0;
  amap.arg("useCache", useCache, "0: zzX cache, 2: DCRT cache");

  amap.parse(argc, argv);

  SetNumThreads(nthreads);

  SetSeed(conv<ZZ>(seed));
  TestIt(p, r, c, k, /*Key Hamming weight=*/64, L, mvec, gens, ords, useCache);
}

// ./Test_EvalMap_x mvec="[73 433]" gens="[18620 12995]" ords="[72 -6]"
#endif

#if 1

//===================== Test_ThinBootstrapping.cpp =======================//



static bool noPrint = false;

// #define DEBUG_PRINTOUT
#ifdef DEBUG_PRINTOUT
#include "debugging.h"
extern long printFlag;
#endif

static long mValues[][14] = { 
//{ p, phi(m),  m,    d, m1,  m2, m3,   g1,    g2,    g3,ord1,ord2,ord3, c_m}
  {  2,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 100},
  {  2,   600,  1023, 10, 11,  93,  0,   838,   584,    0, 10,  6,   0, 100}, // m=(3)*11*{31} m/phim(m)=1.7    C=24  D=2 E=1
  {  2,  1200,  1705, 20, 11, 155,  0,   156,   936,    0, 10,  6,   0, 100}, // m=(5)*11*{31} m/phim(m)=1.42   C=34  D=2 E=2
  {  2,  1728,  4095, 12,  7,  5, 117,  2341,  3277, 3641,  6,  4,   6, 100}, // m=(3^2)*5*7*{13} m/phim(m)=2.36 C=26 D=3 E=2
  {  2,  2304,  4641, 24,  7,  3, 221,  3979,  3095, 3760,  6,  2,  -8, 100}, // m=3*7*(13)*{17} :-( m/phim(m)=2.01 C=45 D=4 E=3
  {  2,  4096,  4369, 16, 17, 257,  0,   258,  4115,    0, 16,-16,   0, 100}, // m=17*(257) :-( m/phim(m)=1.06 C=61 D=3 E=4
  {  2, 12800, 17425, 40, 41, 425,  0,  5951,  8078,    0, 40, -8,   0, 100}, // m=(5^2)*{17}*41 m/phim(m)=1.36 C=93  D=3 E=3
  {  2, 15004, 15709, 22, 23, 683,  0,  4099, 13663,    0, 22, 31,   0, 100}, // m=23*(683) m/phim(m)=1.04      C=73  D=2 E=1
  {  2, 16384, 21845, 16, 17,   5,257,  8996, 17477, 21591, 16, 4, -16,1600}, // m=5*17*(257) :-( m/phim(m)=1.33 C=65 D=4 E=4
  {  2, 18000, 18631, 25, 31, 601,  0, 15627,  1334,    0, 30, 24,   0,  50}, // m=31*(601) m/phim(m)=1.03      C=77  D=2 E=0
  {  2, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29  C=84  D=2 E=2
  {  2, 21168, 27305, 28, 43, 635,  0, 10796, 26059,    0, 42, 18,   0, 100}, // m=(5)*43*{127} m/phim(m)=1.28  C=86  D=2 E=2
  {  2, 23040, 28679, 24, 17,  7, 241, 15184,  4098,28204, 16,  6, -10,1000}, // m=7*17*(241) m/phim(m)=1.24    C=63  D=4 E=3
  {  2, 24000, 31775, 20, 41, 775,  0,  6976, 24806,    0, 40, 30,   0, 100}, // m=(5^2)*{31}*41 m/phim(m)=1.32 C=88  D=2 E=2
  {  2, 26400, 27311, 55, 31, 881,  0, 21145,  1830,    0, 30, 16,   0, 100}, // m=31*(881) m/phim(m)=1.03      C=99  D=2 E=0
  {  2, 27000, 32767, 15, 31,  7, 151, 11628, 28087,25824, 30,  6, -10, 150},
  {  2, 31104, 35113, 36, 37, 949,  0, 16134,  8548,    0, 36, 24,   0, 400}, // m=(13)*37*{73} m/phim(m)=1.12  C=94  D=2 E=2
  {  2, 34848, 45655, 44, 23, 1985, 0, 33746, 27831,    0, 22, 36,   0, 100}, // m=(5)*23*{397} m/phim(m)=1.31  C=100 D=2 E=2
  {  2, 42336, 42799, 21, 127, 337, 0, 25276, 40133,    0,126, 16,   0,  20}, // m=127*(337) m/phim(m)=1.01     C=161 D=2 E=0
  {  2, 45360, 46063, 45, 73, 631,  0, 35337, 20222,    0, 72, 14,   0, 100}, // m=73*(631) m/phim(m)=1.01      C=129 D=2 E=0
  {  2, 46080, 53261, 24, 17, 13, 241, 43863, 28680,15913, 16, 12, -10, 100}, // m=13*17*(241) m/phim(m)=1.15   C=69  D=4 E=3
  {  2, 49500, 49981, 30, 151, 331, 0,  6952, 28540,    0,150, 11,   0, 100}, // m=151*(331) m/phim(m)=1        C=189 D=2 E=1
  {  2, 54000, 55831, 25, 31, 1801, 0, 19812, 50593,    0, 30, 72,   0, 100}, // m=31*(1801) m/phim(m)=1.03     C=125 D=2 E=0
  {  2, 60016, 60787, 22, 89, 683,  0,  2050, 58741,    0, 88, 31,   0, 200}, // m=89*(683) m/phim(m)=1.01      C=139 D=2 E=1

  {  7,    36,    57,  3,  3,  19,  0,    20,    40,    0,  2, -6,   0, 100}, // m=3*(19) :-( m/phim(m)=1.58 C=14 D=3 E=0

  { 17,    48,   105, 12,  3,  35,  0,    71,    76,    0,  2,  2,   0, 100}, // m=3*(5)*{7} m/phim(m)=2.18 C=14 D=2 E=2
  { 17,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36  C=22  D=3
  { 17, 18000, 21917, 30, 101, 217, 0,  5860,  5455,    0, 100, 6,   0, 100}, // m=(7)*{31}*101 m/phim(m)=1.21  C=134 D=2 
  { 17, 30000, 34441, 30, 101, 341, 0,  2729, 31715,    0, 100, 10,  0, 100}, // m=(11)*{31}*101 m/phim(m)=1.14 C=138 D=2
  { 17, 40000, 45551, 40, 101, 451, 0, 19394,  7677,    0, 100, 10,  0,2000}, // m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2
  { 17, 46656, 52429, 36, 109, 481, 0, 46658,  5778,    0, 108, 12,  0, 100}, // m=(13)*{37}*109 m/phim(m)=1.12 C=154 D=2
  { 17, 54208, 59363, 44, 23, 2581, 0, 25811,  5199,    0, 22, 56,   0, 100}, // m=23*(29)*{89} m/phim(m)=1.09  C=120 D=2
  { 17, 70000, 78881, 10, 101, 781, 0, 67167, 58581,    0, 100, 70,  0, 100}, // m=(11)*{71}*101 m/phim(m)=1.12 C=178 D=2

  {127,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36   C=22  D=3
  {127,  1200,  1925, 20,  11, 175, 0,  1751,   199,    0, 10, 6,    0, 100}, //  m=(5^2)*{7}*11 m/phim(m)=1.6   C=34 D=2
  {127,  2160,  2821, 30,  13, 217, 0,   652,   222,    0, 12, 6,    0, 100}, // m=(7)*13*{31} m/phim(m)=1.3     C=46 D=2
  {127, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29   C=84  D=2
  {127, 26112, 30277, 24, 17, 1781, 0, 14249, 10694,    0, 16, 68,   0, 100}, // m=(13)*17*{137} m/phim(m)=1.15  C=106 D=2
  {127, 31752, 32551, 14, 43,  757, 0,  7571, 28768,    0, 42, 54,   0, 100}, // m=43*(757) :-( m/phim(m)=1.02   C=161 D=3
  {127, 46656, 51319, 36, 37, 1387, 0, 48546, 24976,    0, 36, -36,  0, 200}, //m=(19)*37*{73}:-( m/phim(m)=1.09 C=141 D=3
  {127, 49392, 61103, 28, 43, 1421, 0,  1422, 14234,    0, 42, 42,   0,4000}, // m=(7^2)*{29}*43 m/phim(m)=1.23  C=110 D=2
  {127, 54400, 61787, 40, 41, 1507, 0, 30141, 46782,    0, 40, 34,   0, 100}, // m=(11)*41*{137} m/phim(m)=1.13  C=112 D=2
  {127, 72000, 77531, 30, 61, 1271, 0,  7627, 34344,    0, 60, 40,   0, 100}  // m=(31)*{41}*61 m/phim(m)=1.07   C=128 D=2
};
#define num_mValues (sizeof(mValues)/(14*sizeof(long)))

#define OUTER_REP (1)
#define INNER_REP (1)

static bool dry = false; // a dry-run flag

void TestIt(long idx, long p, long r, long L, long c, long B, long skHwt, bool cons=false, int build_cache=0)
{
  Vec<long> mvec;
  vector<long> gens;
  vector<long> ords;

  long phim = mValues[idx][1];
  long m = mValues[idx][2];
  assert(GCD(p, m) == 1);

  append(mvec, mValues[idx][4]);
  if (mValues[idx][5]>1) append(mvec, mValues[idx][5]);
  if (mValues[idx][6]>1) append(mvec, mValues[idx][6]);
  gens.push_back(mValues[idx][7]);
  if (mValues[idx][8]>1) gens.push_back(mValues[idx][8]);
  if (mValues[idx][9]>1) gens.push_back(mValues[idx][9]);
  ords.push_back(mValues[idx][10]);
  if (abs(mValues[idx][11])>1) ords.push_back(mValues[idx][11]);
  if (abs(mValues[idx][12])>1) ords.push_back(mValues[idx][12]);

  if (!noPrint) {
    cout << "*** TestIt";
    if (isDryRun()) cout << " (dry run)";
    cout << ": p=" << p
	 << ", r=" << r
	 << ", L=" << L
	 << ", B=" << B
	 << ", t=" << skHwt
	 << ", c=" << c
	 << ", m=" << m
	 << " (=" << mvec << "), gens="<<gens<<", ords="<<ords
	 << endl;
    cout << "Computing key-independent tables..." << std::flush;
  }
  setTimersOn();
  setDryRun(false); // Need to get a "real context" to test bootstrapping

  double t = -GetTime();
  FHEcontext context(m, p, r, gens, ords);
  context.bitsPerLevel = B;
  buildModChain(context, L, c,/*extraBits=*/7);

  // FIXME: The extraBits is an exceedingly ugly patch, used to bypass the
  //   issue that buildModChain must be called BEFORE the context is made
  //   bootstrappable (else the "powerful" basis is not initialized correctly.)
  //   This is a bug, the value 7 is sometimes the right one, but seriously??

  context.makeBootstrappable(mvec, /*t=*/0, cons, build_cache);

  ThinRecryptData rcData;

  FHE_NTIMER_START(AAA_rcDAta_init);
  rcData.init(context, mvec, 0, cons, build_cache);
  FHE_NTIMER_STOP(AAA_rcDAta_init);


  t += GetTime();

  if (skHwt>0) rcData.skHwt = skHwt;
  if (!noPrint) {
    cout << " done in "<<t<<" seconds\n";
    cout << "  e="    << context.rcData.e
	 << ", e'="   << context.rcData.ePrime
	 << ", alpha="<< context.rcData.alpha
	 << ", t="    << context.rcData.skHwt
	 << "\n  ";
    context.zMStar.printout();
  }
  setDryRun(dry); // Now we can set the dry-run flag if desired

  long nPrimes = context.numPrimes();
  IndexSet allPrimes(0,nPrimes-1);
  double bitsize = context.logOfProduct(allPrimes)/log(2.0);
  if (!noPrint)
    cout << "  "<<nPrimes<<" primes in chain, total bitsize="
	 << ceil(bitsize) << ", secparam="
	 << (7.2*phim/bitsize -110) << endl;

  long p2r = context.alMod.getPPowR();
  context.zMStar.set_cM(mValues[idx][13]/100.0);

  for (long numkey=0; numkey<OUTER_REP; numkey++) { // test with 3 keys

  t = -GetTime();
  if (!noPrint) cout << "Generating keys, " << std::flush;
  FHESecKey secretKey(context);
  FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64);      // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);
  if (!noPrint) cout << "computing key-dependent tables..." << std::flush;
  secretKey.genRecryptData();
  t += GetTime();
  if (!noPrint) cout << " done in "<<t<<" seconds\n";

  long d = context.zMStar.getOrdP();
  long phim = context.zMStar.getPhiM();
  long nslots = phim/d;

  // GG defines the plaintext space Z_p[X]/GG(X)
  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  EncryptedArray ea(context, GG);

  zz_p::init(p2r);
  Vec<zz_p> val0;
  random(val0, nslots);

  vector<ZZX> val1;
  val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    val1[i] = conv<ZZX>(conv<ZZ>(rep(val0[i])));
  }

  Ctxt c1(publicKey);
  ea.encrypt(c1, publicKey, val1);

  Ctxt c2(c1);
  CheckCtxt(c2, "before");

  thinRecrypt(c2, rcData);
  CheckCtxt(c2, "after");

  vector<ZZX> val2;
  ea.decrypt(c2, secretKey, val2);

  if (val1 == val2) 
    cerr << "GOOD\n";
  else
    cerr << "BAD\n";


  //cerr << convert<Vec<ZZX>>(val1) << "\n";

  }

  printAllTimers();
  
}

/********************************************************************
 ********************************************************************/
int main(int argc, char *argv[]) 
{
  ArgMapping amap;

  long p=2;
  long r=1;
  long c=3;
  long L=25;
  long B=23;
  long N=0;
  long t=0;
  bool cons=0;
  long nthreads=1;

  long seed=0;
  long useCache=1;

  amap.arg("p", p, "plaintext base");

  amap.arg("r", r,  "exponent");
  amap.note("p^r is the plaintext-space modulus");

  amap.arg("c", c, "number of columns in the key-switching matrices");
  amap.arg("L", L, "# of levels in the modulus chain");
  amap.arg("B", B, "# of bits per level");
  amap.arg("N", N, "lower-bound on phi(m)");
  amap.arg("t", t, "Hamming weight of recryption secret key", "heuristic");
  amap.arg("dry", dry, "dry=1 for a dry-run");
  amap.arg("cons", cons, "cons=1 for consevative settings (circuit deeper by 1)");
  amap.arg("nthreads", nthreads, "number of threads");
  amap.arg("seed", seed, "random number seed");
  amap.arg("noPrint", noPrint, "suppress printouts");
  amap.arg("useCache", useCache, "0: zzX cache, 1: DCRT cache");

  amap.parse(argc, argv);

  if (seed) 
    SetSeed(ZZ(seed));

  SetNumThreads(nthreads);


  if (B<=0) B=FHE_pSize;
  if (B>NTL_SP_NBITS/2) B = NTL_SP_NBITS/2;

  for (long i=0; i<(long)num_mValues; i++)
    if (mValues[i][0]==p && mValues[i][1]>=N) {
      TestIt(i,p,r,L,c,B,t,cons,useCache);
      break;
    }

  return 0;
}


#endif
