
#ifndef _HELIB_SPX_H_
#define _HELIB_SPX_H_
/* SPX.h - A MUX class that implements either GF2X or zz_pX: it keeps both
 *   objects, and decides in runitime to use one or the other as appropriate.
 *
 * Modifications:
 *
 * Created by Shai Halevi on May 27, 2015
 */
#include <iostream>
#include <NTL/version.h>
#include <NTL/ZZX.h>
#include <NTL/GF2X.h>
#include <NTL/lzz_pX.h>

// A helper class that maintains the NTL context (this is unnecessary
// as of NTL version 9.2.0)
#if (NTL_MAJOR_VERSION<9 || (NTL_MAJOR_VERSION==9 && NTL_MINOR_VERSION<2))
class SPmodulus { // single precision modulus
public:
  long p;                   // the modulus itself
  NTL::zz_pContext context; // NTL context, used if !gf2

  SPmodulus(): p(0) {}
  SPmodulus(long _p): p(_p), context(_p) {}
  SPmodulus(long _p, const NTL::zz_pContext& _c): p(_p), context(_c) {}

  bool null() const { return (p==0); }
  long modulus() const  { return p; }
  bool equals(const SPmodulus& other) const
  { return p==other.p; }
};
inline const NTL::zz_pContext& getContext(const SPmodulus& spm) { return spm.context; }
#else
typedef NTL::zz_pContext SPmodulus;
inline const NTL::zz_pContext& getContext(const SPmodulus& spm) { return spm; }
#endif
inline bool isGF2(const SPmodulus& spm) { return spm.modulus()==2; }



// macro for non-const unary operators (e.g., ++SPX)
#define sppolyUnary(fnc) SPX& fnc() {			\
    if (isGF2(modulus)) {NTL::fnc(gf2x);}			\
    else {NTL::zz_pPush push(getContext(modulus)); NTL::fnc(zzpx);}	\
    return *this;}

// macro for non-const binary operators that take a long (SPX += long)
#define sppolyBinaryLeft(fnc) SPX fnc(long y) {		\
  if (isGF2(modulus)) { NTL::fnc(gf2x,y); }			\
  else { NTL::zz_pPush push(getContext(modulus)); NTL::fnc(zzpx,y); }\
  return *this; }

// macro for non-const binay operators that take another SPX
// (SPX += SPX)
#define sppolyBinaryBoth(fnc) SPX fnc(const SPX& y) {	\
  ensureConsistency(*this,y);					\
  if (isGF2(modulus)) { NTL::fnc(gf2x,y.gf2x); }               \
  else { NTL::zz_pPush push(getContext(modulus)); NTL::fnc(zzpx,y.zzpx); }	\
  return *this;}

// The main class SPX, declerations follow GF2X, zz_pX
class SPX {
public:
  SPmodulus modulus; // not a reference, the object itself
  NTL::GF2X gf2x;
  NTL::zz_pX zzpx;

  // access methods
  const SPmodulus& getMod() const { return modulus; }

  // helper functions
  template <class T> void makeFromVec(const NTL::Vec<T>& v)
  {
    if (isGF2(modulus)) {
      NTL::Vec<NTL::GF2> v2;
      NTL::conv(v2, v);
      NTL::conv(gf2x, v2);
    } else {
      NTL::zz_pPush push(getContext(modulus));
      NTL::Vec<NTL::zz_p> vp;
      NTL::conv(vp, v);
      NTL::conv(zzpx, vp);
    }
  }

  long deg() const
  { return (isGF2(modulus)? NTL::deg(gf2x): NTL::deg(zzpx)); }

  template <class T> void getCoeffVec(NTL::Vec<T>& v, long n=-1) const
  {
    if (n<0) n=deg()+1;
    if (isGF2(modulus)) {
      NTL::Vec<NTL::GF2> v2 = NTL::VectorCopy(gf2x, n);
      NTL::conv(v, v2);
    } else {
      NTL::zz_pPush push(getContext(modulus));
      NTL::Vec<NTL::zz_p> vp = NTL::VectorCopy(zzpx, n);
      NTL::conv(v, vp);
    }
  }

  // internal constructors, useful for debugging
  explicit SPX(const NTL::GF2X& p2): modulus(2), gf2x(p2) {}
  SPX(const NTL::zz_pX& pp, const SPmodulus& mod): modulus(mod) {
    if (isGF2(mod)) { // need to convert to GF2X
      NTL::Vec<long> v;
      NTL::conv(v, pp.rep);
      makeFromVec(v);
    }
    else zzpx = pp;
  }

  void resetMod(const SPmodulus& mod)
  {
    if (modulus.null()) { // uninititalized object
      modulus = mod;
      return;
    }
    bool noteq = false;
    if (isGF2(mod)!=isGF2(modulus))          noteq = true;
    if (!isGF2(mod) && !mod.equals(modulus)) noteq = true;
    if (noteq) {
      gf2x.kill();
      zzpx.kill();
      modulus = mod;
    }
  }

  void setCoeff(long i, long c=1)
  {
    if (isGF2(modulus)) NTL::SetCoeff(gf2x, i, c);
    else                NTL::SetCoeff(zzpx, i, c);
  }

  // Checking consistency between arguments, T is either SPX or SPXModulus
  template<class T>
  static void ensureConsistency(const SPX& x, const T& y)
  { if (isGF2(x.getMod()) && isGF2(y.getMod())) return;
    if (isGF2(x.getMod()) || isGF2(y.getMod()) || !x.getMod().equals(y.getMod())) {
      throw std::invalid_argument("mismatched moduli in SPX arguments");
    }
  }
  template<class S, class T>
  static void ensureConsistency(const SPX& x, const S& y, const T& z)
  { ensureConsistency(x,y); ensureConsistency(x,z); }

  // Constructors

  SPX(){}                             // empty object
  explicit SPX(const SPmodulus& mod): modulus(mod) {} // equal to zero
  SPX(const NTL::Vec<long>& coefs, const SPmodulus& mod): // initialized
    modulus(mod) { makeFromVec<long>(coefs); }
  SPX(const NTL::ZZX& poly, const SPmodulus& mod):// initialized from ZZX
    modulus(mod) { makeFromVec<NTL::ZZ>(poly.rep); }

  // Convencience constructors: initialize to X^i or c*X^i
  SPX(NTL::INIT_MONO_TYPE, long i, const SPmodulus& mod):
    modulus(mod) { setCoeff(i); }
  SPX(NTL::INIT_MONO_TYPE, long i, long c, const SPmodulus& mod):
    modulus(mod) { setCoeff(i,c); }

  SPX(NTL::INIT_SIZE_TYPE, long n, const SPmodulus& mod):
    modulus(mod) { setCoeff(n-1,0); }
  // SPX(INIT_SIZE, n, mod) initializes to zero, but space is
  // pre-allocated for n coefficients

  // Implicit copy constructor, destructor, and assignment operator

  // f.SetMaxLength(n) pre-allocate spaces for n coefficients.  The
  // polynomial that f represents is unchanged.
  void SetMaxLength(long n)
  {
    if (modulus.null()) return; // uninitialized, cannot allocate space
    if (isGF2(modulus)) gf2x.SetMaxLength(n);
    else                zzpx.SetMaxLength(n);
  }

  // f.kill() sets f to 0 and frees all memory held by f.
  void kill()
  {
    if (modulus.null()) return; // uninitialized, nothing to do
    if (isGF2(modulus)) gf2x.kill();
    else                zzpx.kill();
  }

  sppolyBinaryLeft(operator+=) // SPX += long
  sppolyBinaryBoth(operator+=) // SPX += SPPoly

  template <class T> SPX operator+(T b) const
  { SPX tmp=*this; tmp += b; return tmp; }

  sppolyUnary(operator++)                // ++SPX
  void operator++(int) { operator++(); } // SPX++

  sppolyBinaryLeft(operator-=) // SPX -= long
  sppolyBinaryBoth(operator-=) // SPX -= SPPoly

  template <class T> SPX operator-(T b) const
  { SPX tmp=*this; tmp -= b; return tmp; }

  sppolyUnary(operator--)                // --SPX
  void operator--(int) { operator--(); } // SPX--

  void negate() {
    if (isGF2(modulus)) return; // nothing to do, over GF2 we have X == -X
    else {
      NTL::zz_pPush push(getContext(modulus));
      NTL::negate(zzpx, zzpx);
    }
  }
  SPX operator-() const
  { SPX tmp=*this; tmp.negate(); return tmp; }


  sppolyBinaryLeft(operator*=) // SPX *= long
  sppolyBinaryBoth(operator*=) // SPX *= SPPoly

  template <class T> SPX operator*(T b) const
  { SPX tmp=*this; tmp *= b; return tmp; }

  sppolyBinaryLeft(operator<<=) // SPX <<= long
  SPX operator<<(long n) const
  { SPX tmp=*this; tmp <<= n; return tmp; }

  sppolyBinaryLeft(operator>>=) // SPX >>= long
  SPX operator>>(long n) const
  { SPX tmp=*this; tmp >>= n; return tmp; }

  sppolyBinaryBoth(operator/=) // SPX /= SPX
  sppolyBinaryLeft(operator/=) // SPX /= long

  template <class T> SPX operator/(T b) const
  { SPX tmp=*this; tmp /= b; return tmp; }

  sppolyBinaryBoth(operator%=) // SPX %= SPX

  SPX operator%(const SPX& b) const
  { SPX tmp=*this; tmp %= b; return tmp; }

// SIZE INVARIANT: for any f in SPX, deg(f)+1 < 2^(NTL_BITS_PER_LONG-4).
};

typedef NTL::Vec<SPX> vec_SPX; // a conveneince vector decleration

// macro for function with signature void fnc(const SPX& x)
#define sppolyOnePoly(fnc) fnc(SPX& x) {	 \
    if (isGF2(x.getMod())) NTL::fnc(x.gf2x);	 \
    else { NTL::zz_pPush push(getContext(x.getMod())); \
           NTL::fnc(x.zzpx);}}

// macro for function with signature non-void fnc(const type& x)
// type is either SPX or SPXModulus
#define sppolyConstOnePoly(fnc, typ) fnc(const typ& x) { \
    if (isGF2(x.getMod())) return NTL::fnc(x.gf2x); \
    else { NTL::zz_pPush push(getContext(x.getMod()));      \
           return NTL::fnc(x.zzpx);}}

// macro for function with signature void fnc(SPX& x, const type& arg)
// type is either SPX or SPXModulus
#define sppolyTwoPoly(fnc,typ) fnc(SPX& x, typ& y) {\
    x.resetMod(y.getMod());				     \
    if (isGF2(x.getMod())) NTL::fnc(x.gf2x, y.gf2x);    \
    else { NTL::zz_pPush push(getContext(x.getMod()));	     \
           NTL::fnc(x.zzpx, y.zzpx);}}

// macro for function with signature
// void fnc(SPX& x, const SPX& a, const type& b)
// type is either SPX or SPXModulus
#define sppolyThreePoly(fnc,typ)\
  fnc(SPX& x, const SPX& a, const typ& b) {\
    SPX::ensureConsistency(a,b);\
    x.resetMod(a.getMod());	   \
    if (isGF2(x.getMod())) NTL::fnc(x.gf2x, a.gf2x, b.gf2x);\
    else { NTL::zz_pPush push(getContext(x.getMod()));	   \
           NTL::fnc(x.zzpx, a.zzpx, b.zzpx);}}

// macro for function with signature
// void fnc(SPX& x, const SPX& a, const type& b)
// type is either SPX or SPXModulus
#define sppolyFourPoly(fnc,typ)\
  fnc(SPX& x, const SPX& a, const SPX& b, const typ& f) {  \
    SPX::ensureConsistency(a,b,f);				    \
    x.resetMod(a.getMod());					    \
    if (isGF2(x.getMod()))					    \
      NTL::fnc(x.gf2x, a.gf2x, b.gf2x, f.gf2x);	    \
    else { NTL::zz_pPush push(getContext(x.getMod()));		    \
           NTL::fnc(x.zzpx, a.zzpx, b.zzpx, f.zzpx);}}

// macro for function with signature
// void fnc(SPX& x, const SPX& a, type b)
#define sppolyTwoPolyArg(fnc,typ) \
  fnc(SPX& x, const SPX& a, typ b) {		    \
    x.resetMod(a.getMod());					\
    if (isGF2(x.getMod())) NTL::fnc(x.gf2x, a.gf2x, b);\
    else { NTL::zz_pPush push(getContext(x.getMod()));	    \
           NTL::fnc(x.zzpx, a.zzpx, b);}}

// macro for function with signature
// void fnc(SPX& x, const SPX& g, const typ& F, long m)
#define sppolyThreePolyArg(fnc,typ)					\
  fnc(SPX& x, const SPX& g, const typ& F, long m) {\
    SPX::ensureConsistency(g,F);			      \
    x.resetMod(g.getMod());				      \
    if (isGF2(g.getMod())) NTL::fnc(x.gf2x, g.gf2x, F.gf2x, m);\
    else { NTL::zz_pPush push(getContext(g.getMod()));	      \
           NTL::fnc(x.zzpx, g.zzpx, F.zzpx, m);}}

/**************************************************************************\
                              Accessing coefficients

The degree of a polynomial f is obtained as deg(f),
where the zero polynomial, by definition, has degree -1.

A polynomial f is represented as a coefficient vector.
Coefficients may be accessed via coeff(f, i) to get the
coefficient of X^i in the polynomial f, and SetCoeff(f, i, a)
to set the coefficient of X^i in f to the scalar a.

\**************************************************************************/

inline long deg(const SPX& a) { return a.deg(); }
// return deg(a); deg(0) == -1.

inline long coeff(const SPX& x, long i)
{
  if (isGF2(x.getMod())) return NTL::conv<long>(NTL::coeff(x.gf2x, i));
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    return NTL::conv<long>(NTL::coeff(x.zzpx,i));
  }
}

// const long coeff(const SPX& a, long i);
// returns the coefficient of X^i, or zero if i not in range

inline long LeadCoeff(const SPX& a)
{ return coeff(a, deg(a)); }
// returns leading term of a, or zero if a == 0

const long ConstTerm(const SPX& a)
{ return coeff(a, 0); }
// returns constant term of a, or zero if a == 0

void SetCoeff(SPX& x, long i, long a)
{ x.setCoeff(i,a); }
// makes coefficient of X^i equal to a; error is raised if i < 0

void SetCoeff(SPX& x, long i)
{ x.setCoeff(i); }

// makes coefficient of X^i equal to 1;  error is raised if i < 0

inline void sppolyOnePoly(SetX)
// void SetX(SPX& x); // x is set to the monomial X

inline bool sppolyConstOnePoly(IsX,SPX)
// long IsX(const SPX& a); // test if x = X


/**************************************************************************\

                                  Comparison

\**************************************************************************/

bool operator==(const SPX& a, const SPX& b) {
  if (isGF2(a.getMod())!=isGF2(b.getMod())) return false;
  if (isGF2(a.getMod())) return (a.gf2x==b.gf2x);
  else return (a.getMod().equals(b.getMod()) && (a.zzpx==b.zzpx));
}
bool operator!=(const SPX& a, const SPX& b)
{ return !(a==b); }

inline bool sppolyConstOnePoly(IsZero,SPX)
// bool IsZero(const SPX& a); // test for 0

inline bool sppolyConstOnePoly(IsOne,SPX)
// bool IsOne(const SPX& a); // test for 1

// PROMOTIONS: operators ==, != promote {long, long} to SPX on (a, b)


/**************************************************************************\

                                   Addition

\**************************************************************************/

// operator notation:

inline SPX operator+(long a, const SPX& b)
{ return (b+a);}
//SPX operator-(long, const SPX& b); // not implemented

// procedural versions:

inline void add(SPX& x, const SPX& a, const SPX& b)// x = a + b
{ x = a; x+= b; }
inline void add(SPX& x, long a, const SPX& b) // x = a + b
{ x = b; x+= a; }
inline void add(SPX& x, const SPX& a, long b) // x = a + b
{ x = a; x+= b; }
inline void sub(SPX& x, const SPX& a, const SPX& b) // x = a - b
{ x=a; x-=b; }
// void sub(SPX& x, long a, const SPX& b); // x = a - b, not implemented
inline void sub(SPX& x, const SPX& a, long b) // x = a - b
{ x=a; x-=b; }
inline void negate(SPX& x, const SPX& a) // x = -a
{ x = a; x.negate(); }
//inline void negate(SPX& x, long a) // x = -a, not implemented


/**************************************************************************\

                               Multiplication

\**************************************************************************/

// operator notation:

inline SPX operator*(long a, const SPX& b)
{ return (b*a); }

// procedural versions:

inline void mul(SPX& x, const SPX& a, const SPX& b) // x = a * b
{ x=a; x*=b; }
inline void mul(SPX& x, long a, const SPX& b) // x = a * b
{ x=b; x*=a; }
inline void mul(SPX& x, const SPX& a, long b) // x = a * b
{ x=a; x*=b; }

inline void sppolyTwoPoly(sqr,const SPX)
//void sqr(SPX& x, const SPX& a); // x = a^2

inline SPX sqr(const SPX& a)
{ SPX x; sqr(x,a); return x; }

/**************************************************************************\

                               Shift Operations

LeftShift by n means multiplication by X^n
RightShift by n means division by X^n

A negative shift amount reverses the direction of the shift.

\**************************************************************************/

// procedural versions:

inline void LeftShift(SPX& x, const SPX& a, long n)
{ x=a; x <<= n; }
inline SPX LeftShift(const SPX& a, long n)
{ return a << n; }

inline void RightShift(SPX& x, const SPX& a, long n)
{ x=a; x >>= n; }
inline SPX RightShift(const SPX& a, long n)
{ return a >> n; }

/**************************************************************************\

                                  Division

\**************************************************************************/

// procedural versions:

// q = a/b, r = a%b (class T is either SPX or SPXModulus)
template<class T>
inline void DivRem(SPX& q, SPX& r, const SPX& a, const T& b)
{
  SPX::ensureConsistency(a,b);
  q.resetMod(a.getMod());
  r.resetMod(a.getMod());
  if (isGF2(a.getMod()))
    NTL::DivRem(q.gf2x, r.gf2x, a.gf2x, b.gf2x);
  else {
    NTL::zz_pPush push(getContext(a.getMod()));
    NTL::DivRem(q.zzpx, r.zzpx, a.zzpx, b.zzpx);
  }
}

// q = a/b
inline void sppolyThreePoly(div,SPX)
// void div(SPX& q, const SPX& a, const SPX& b);

// r = a%b
inline void sppolyThreePoly(rem,SPX)
// void rem(SPX& r, const SPX& a, const SPX& b);

// if b | a, sets q = a/b and returns 1; otherwise returns 0
inline bool divide(SPX& q, const SPX& a, const SPX& b)
{
  if (isGF2(a.getMod())!=isGF2(b.getMod())) return false;
  if (!isGF2(a.getMod()) && !a.getMod().equals(b.getMod())) return false;
  q.resetMod(a.getMod());
  if (isGF2(a.getMod())) return NTL::divide(q.gf2x, a.gf2x, b.gf2x);
  else {
    NTL::zz_pPush push(getContext(a.getMod()));
    return NTL::divide(q.zzpx, a.zzpx, b.zzpx);
  }
}

inline bool divide(const SPX& a, const SPX& b)
{
  if (isGF2(a.getMod())!=isGF2(b.getMod())) return false;
  if (!isGF2(a.getMod()) && !a.getMod().equals(b.getMod())) return false;
  if (isGF2(a.getMod())) return NTL::divide(a.gf2x, b.gf2x);
  else {
    NTL::zz_pPush push(getContext(a.getMod()));
    return NTL::divide(a.zzpx, b.zzpx);
  }
}
// if b | a returns 1; otherwise returns 0

/**************************************************************************\

                                   GCD's

\**************************************************************************/


inline void sppolyThreePoly(GCD,SPX)
//void GCD(SPX& x, const SPX& a, const SPX& b);

SPX GCD(const SPX& a, const SPX& b)
{ SPX x; GCD(x, a, b); return x; }
// x = GCD(a, b) (zero if a==b==0).

// d = gcd(a,b), a s + b t = d
inline void XGCD(SPX& d, SPX& s, SPX& t,
		 const SPX& a, const SPX& b)
{
  SPX::ensureConsistency(a,b);
  d.resetMod(a.getMod());
  s.resetMod(a.getMod());
  t.resetMod(a.getMod());
  if (isGF2(a.getMod()))
    NTL::XGCD(d.gf2x, s.gf2x, t.gf2x, a.gf2x, b.gf2x);
  else {
    NTL::zz_pPush push(getContext(a.getMod()));
    NTL::XGCD(d.zzpx, s.zzpx, t.zzpx, a.zzpx, b.zzpx);
  }
}

/**************************************************************************\

                                  Input/Output

I/O format:

   [a_0 a_1 ... a_n],

represents the polynomial a_0 + a_1*X + ... + a_n*X^n.

On output, all coefficients will be in [0,p-1], and a_n not zero
(the zero polynomial is [ ]).  On input, the coefficients may be
arbitrary integers which are reduced modulo p, and leading zeros
stripped.

x must have the right modulus set before calling these functions
\**************************************************************************/

inline std::istream& operator>>(std::istream& s, SPX& x)
{
  if (isGF2(x.getMod())) return s >> x.gf2x;
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    return s >> x.gf2x;
  }
}

inline std::ostream& operator<<(std::ostream& s, const SPX& x)
{
  if (isGF2(x.getMod())) return s << x.gf2x;
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    return s << x.zzpx;
  }
}

/**************************************************************************\

                              Some utility routines

\**************************************************************************/

// x = derivative of a
inline void sppolyTwoPoly(diff,const SPX)
// void diff(SPX& x, const SPX& a);

inline SPX diff(const SPX& a)
{ SPX x; diff(x,a); return x; }

inline void reverse(SPX& x, const SPX& a, long hi=-2)
{
  if (hi < -1) hi = deg(a);
  x.resetMod(a.getMod());
  if (isGF2(x.getMod())) NTL::reverse(x.gf2x, a.gf2x, hi);
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    NTL::reverse(x.zzpx, a.zzpx, hi);
  }
}

inline SPX reverse(const SPX& a, long hi=-2)
{ SPX x; reverse(x,a,hi); return x; }

// x = reverse of a[0]..a[hi] (hi >= -1);
// hi defaults to deg(a) in second version


/**************************************************************************\

                             Random Polynomials

\**************************************************************************/

// x = random polynomial of degree < n, with specified modulus
inline void random(SPX& x, long n, const SPmodulus& mod)
{
  x.resetMod(mod);
  if (isGF2(x.getMod())) NTL::random(x.gf2x, n);
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    NTL::random(x.zzpx, n);
  }
}

inline SPX random_SPX(long n, const SPmodulus& mod)
{ SPX x; random(x, n, mod); return x; }

/**************************************************************************\

                       Arithmetic mod X^n

Required: n >= 0; otherwise, an error is raised.

\**************************************************************************/

inline void sppolyTwoPolyArg(trunc,long)
//void trunc(SPX& x, const SPX& a, long n); // x = a % X^n

inline SPX trunc(const SPX& a, long n)
{ SPX x; trunc(x,a,n); return x; }

// x = a * b % X^n
inline void sppolyThreePolyArg(MulTrunc,SPX)
//void MulTrunc(SPX& x, const SPX& a, const SPX& b, long n);

SPX MulTrunc(const SPX& a, const SPX& b, long n)
{ SPX x; MulTrunc(x,a,b,n); return x; }

// x = a^2 % X^n
inline void sppolyTwoPolyArg(SqrTrunc,long)
//void SqrTrunc(SPX& x, const SPX& a, long n);

inline SPX SqrTrunc(const SPX& a, long n)
{ SPX x; SqrTrunc(x,a,n); return x; }

// computes x = a^{-1} % X^n.  Must have ConstTerm(a) invertible.

inline void sppolyTwoPolyArg(InvTrunc, long)
// void InvTrunc(SPX& x, const SPX& a, long n);

inline SPX InvTrunc(const SPX& a, long n)
{ SPX x; InvTrunc(x,a,n); return x; }

/**************************************************************************\

                Modular Arithmetic (without pre-conditioning)

Arithmetic mod f.

All inputs and outputs are polynomials of degree less than deg(f), and
deg(f) > 0.

NOTE: if you want to do many computations with a fixed f, use the
SPXModulus data structure and associated routines below for better
performance.

\**************************************************************************/

// x = (a * b) % f
inline void sppolyFourPoly(MulMod,SPX)
//void MulMod(SPX& x, const SPX& a, const SPX& b, const SPX& f);

inline SPX MulMod(const SPX& a, const SPX& b, const SPX& f)
{ SPX x; MulMod(x,a,b,f); return x; }


// x = a^2 % f
inline void sppolyThreePoly(SqrMod,SPX)
//void SqrMod(SPX& x, const SPX& a, const SPX& f);

inline SPX SqrMod(const SPX& a, const SPX& f)
{ SPX x; SqrMod(x,a,f); return x; }

// x = (a * X) mod f
inline void sppolyThreePoly(MulByXMod,SPX)
//void MulByXMod(SPX& x, const SPX& a, const SPX& f);

inline SPX MulByXMod(const SPX& a, const SPX& f)
{ SPX x; MulByXMod(x,a,f); return x; }

// x = a^{-1} % f, error is a is not invertible
inline void sppolyThreePoly(InvMod,SPX)
//void InvMod(SPX& x, const SPX& a, const SPX& f);

inline SPX InvMod(const SPX& a, const SPX& f)
{ SPX x; InvMod(x,a,f); return x; }

inline bool InvModStatus(SPX& x, const SPX& a, const SPX& f)
{
  SPX::ensureConsistency(a,f);
  x.resetMod(a.getMod());
  if (isGF2(a.getMod()))
    return NTL::InvModStatus(x.gf2x, a.gf2x, f.gf2x);
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    return NTL::InvModStatus(x.zzpx, a.zzpx, f.zzpx);
  }
}
// if (a, f) = 1, returns 0 and sets x = a^{-1} % f; otherwise,
// returns 1 and sets x = (a, f)


// for modular exponentiation, see below



/**************************************************************************\

                     Modular Arithmetic with Pre-Conditioning

If you need to do a lot of arithmetic modulo a fixed f, build
SPXModulus F for f.  This pre-computes information about f that
speeds up subsequent computations.

As an example, the following routine computes the product modulo f of a vector
of polynomials.

#include <NTL/SPX.h>

void product(SPX& x, const vec_SPX& v, const SPX& f)
{
   SPXModulus F(f);
   SPX res;
   res = 1;
   long i;
   for (i = 0; i < v.length(); i++)
      MulMod(res, res, v[i], F);
   x = res;
}


Note that automatic conversions are provided so that a SPX can
be used wherever a SPXModulus is required, and a SPXModulus
can be used wherever a SPX is required.

The SPXModulus routines optimize several important special cases:

  - f = X^n + X^k + 1, where k <= min((n+1)/2, n-NTL_BITS_PER_LONG)

  - f = X^n + X^{k_3} + X^{k_2} + X^{k_1} + 1,
      where k_3 <= min((n+1)/2, n-NTL_BITS_PER_LONG)

  - f = X^n + g, where deg(g) is small


\**************************************************************************/

class SPXModulus {
public:
  SPmodulus modulus; // not a reference, the object itself
  NTL::GF2XModulus gf2x;
  NTL::zz_pXModulus zzpx;
  SPXModulus() {} // initially in an unusable state

  // access methods
  const SPmodulus& getMod() const { return modulus; }

  SPXModulus(const SPX& f): // initialize with f, deg(f)>0
    modulus(f.getMod())
  {
    if (isGF2(f.getMod())) gf2x=f.gf2x;
    else                  zzpx=f.zzpx;
  }
  SPXModulus& operator=(const SPX& f) {
    modulus = f.getMod();
    if (isGF2(f.getMod())) gf2x=f.gf2x;
    else                  zzpx=f.zzpx;
    return *this;
  }
};

inline void build(SPXModulus& F, const SPX& f) { F=f; }
// pre-computes information about f and stores it in F; deg(f) > 0.
// Note that the declaration SPXModulus F(f) is equivalent to
// SPXModulus F; build(F, f).

// In the following, f refers to the polynomial f supplied to the
// build routine, and n = deg(f).

inline long sppolyConstOnePoly(deg,SPXModulus)
//long deg(const SPXModulus& F);  // return deg(f)

// x = (a * b) % f
inline void sppolyFourPoly(MulMod,SPXModulus)
//void MulMod(SPX& x, const SPX& a, const SPX& b, const SPXModulus& F);

inline SPX MulMod(const SPX& a, const SPX& b, const SPXModulus& f)
{ SPX x; MulMod(x,a,b,f); return x; }


// x = a^2 % f
inline void sppolyThreePoly(SqrMod,SPXModulus)
//void SqrMod(SPX& x, const SPX& a, const SPXModulus& F);

inline SPX SqrMod(const SPX& a, const SPXModulus& f)
{ SPX x; SqrMod(x,a,f); return x; }

// x = (a * X) mod f
inline void sppolyThreePoly(MulByXMod,SPXModulus)
//void MulByXMod(SPX& x, const SPX& a, const SPXModulus& F);

inline SPX MulByXMod(const SPX& a, const SPXModulus& f)
{ SPX x; MulByXMod(x,a,f); return x; }

/**
// x = a^{-1} % f, error is a is not invertible
inline void sppolyThreePoly(InvMod,SPXModulus)
//void InvMod(SPX& x, const SPX& a, const SPXModulus& f);

inline SPX InvMod(const SPX& a, const SPXModulus& f)
{ SPX x; InvMod(x,a,f); return x; }
**/

template <class T>
void PowerXMod(SPX& x, T e, const SPXModulus& F)
{
  x.resetMod(F.getMod());
  if (isGF2(x.getMod())) NTL::PowerXMod(x.gf2x, e, F.gf2x);
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    NTL::PowerXMod(x.zzpx, e, F.zzpx);
  }
}


template <class Z> // Z is either long or const ZZ&
void PowerMod(SPX& x, const SPX& a, Z e, const SPXModulus& F)
{
  SPX::ensureConsistency(a,F);
  x.resetMod(a.getMod());
  if (isGF2(a.getMod()))
    NTL::PowerMod(x.gf2x, a.gf2x, e, F.gf2x);
  else {
    NTL::zz_pPush push(getContext(a.getMod()));
    NTL::PowerMod(x.zzpx, a.zzpx, e, F.zzpx);
  }
}

template <class Z> // Z is either long or const ZZ&
SPX PowerMod(const SPX& a, Z e, const SPXModulus& F)
{ SPX x; PowerMod<Z>(x, a, e, F); return x; }


// q = a/b
inline void sppolyThreePoly(div,SPXModulus)
//void div(SPX& q, const SPX& a, const SPXModulus& F);

// r = a%b
inline void sppolyThreePoly(rem,SPXModulus)
//void rem(SPX& x, const SPX& a, const SPXModulus& F);


// x = a % f

//void DivRem(SPX& q, SPX& r, const SPX& a, const SPXModulus& F);
// q = a/f, r = a%f (defined above in a template)

// operator notation:

inline SPX operator/(const SPX& a, const SPXModulus& F)
{ SPX x; div(x,a,F); return x; }
inline SPX operator%(const SPX& a, const SPXModulus& F)
{ SPX x; rem(x,a,F); return x; }

inline SPX& operator/=(SPX& x, const SPXModulus& F)
{
  SPX::ensureConsistency(x,F);
  if (isGF2(x.getMod()))
    x.gf2x /= F.gf2x;
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    x.zzpx /= F.zzpx;
  }
  return x;
}

inline SPX& operator%=(SPX& x, const SPXModulus& F)
{
  SPX::ensureConsistency(x,F);
  if (isGF2(x.getMod()))
    x.gf2x %= F.gf2x;
  else {
    NTL::zz_pPush push(getContext(x.getMod()));
    x.zzpx %= F.zzpx;
  }
  return x;
}


/**************************************************************************\

                              Modular Composition

Modular composition is the problem of computing g(h) mod f for
polynomials f, g, and h.

The algorithm employed is that of Brent & Kung (Fast algorithms for
manipulating formal power series, JACM 25:581-595, 1978), which uses
O(n^{1/2}) modular polynomial multiplications, and O(n^2) scalar
operations.



\**************************************************************************/

// x = g(h) mod f; deg(h) < n
inline void sppolyFourPoly(CompMod,SPXModulus)
//void CompMod(SPX& x, const SPX& g, const SPX& h, const SPXModulus& F);

SPX CompMod(const SPX& g, const SPX& h, const SPXModulus& F)
{ SPX x; CompMod(x,g,h,F); return x; }

inline void Comp2Mod(SPX& x1, SPX& x2, const SPX& g1, const SPX& g2,
		     const SPX& h, const SPXModulus& F)
{
  SPX::ensureConsistency(g1,g2,F); SPX::ensureConsistency(h,F);
  x1.resetMod(g1.getMod());
  x2.resetMod(g2.getMod());
  if (isGF2(F.getMod()))
    NTL::Comp2Mod(x1.gf2x, x2.gf2x,
	     g1.gf2x, g2.gf2x, h.gf2x, F.gf2x);
  else {
    NTL::zz_pPush push(getContext(F.getMod()));
    NTL::Comp2Mod(x1.zzpx, x2.zzpx,
	     g1.zzpx, g2.zzpx, h.zzpx, F.zzpx);
  }
}
// xi = gi(h) mod f (i=1,2), deg(h) < n.

inline void Comp3Mod(SPX& x1, SPX& x2, SPX& x3,
		     const SPX& g1, const SPX& g2, const SPX& g3,
		     const SPX& h, const SPXModulus& F)
{
  SPX::ensureConsistency(g1,g2,F); SPX::ensureConsistency(g3,h,F);
  x1.resetMod(g1.getMod());
  x2.resetMod(g2.getMod());
  x3.resetMod(g3.getMod());
  if (isGF2(F.getMod()))
    NTL::Comp3Mod(x1.gf2x, x2.gf2x, x3.gf2x,
	     g1.gf2x, g2.gf2x, g3.gf2x, h.gf2x, F.gf2x);
  else {
    NTL::zz_pPush push(getContext(F.getMod()));
    NTL::Comp3Mod(x1.zzpx, x2.zzpx, x3.zzpx,
	     g1.zzpx, g2.zzpx, g3.zzpx, h.zzpx, F.zzpx);
  }
}

// xi = gi(h) mod f (i=1..3), deg(h) < n


/**************************************************************************\

                     Composition with Pre-Conditioning

If a single h is going to be used with many g's then you should build
a SPXArgument for h, and then use the compose routine below.  The
routine build computes and stores h, h^2, ..., h^m mod f.  After this
pre-computation, composing a polynomial of degree roughly n with h
takes n/m multiplies mod f, plus n^2 scalar multiplies.  Thus,
increasing m increases the space requirement and the pre-computation
time, but reduces the composition time.

\**************************************************************************/


struct SPXArgument {
  SPmodulus modulus;
  NTL::GF2XArgument gf2x;
  NTL::zz_pXArgument zzpx;

  // access methods
  const SPmodulus& getMod() const { return modulus; }

  void resetMod(const SPmodulus& mod)
  {
    if (modulus.null()) { // uninititalized object
      modulus = mod;
      return;
    }
    bool noteq = false;
    if (isGF2(mod)!=isGF2(modulus))          noteq = true;
    if (!isGF2(mod) && !mod.equals(modulus)) noteq = true;
    if (noteq) {
      gf2x.H.kill();
      zzpx.H.kill();
      modulus = mod;
    }
  }
};

inline void build(SPXArgument& H,
		  const SPX& h, const SPXModulus& F, long m)
{
  SPX::ensureConsistency(h,F);
  H.resetMod(h.getMod());
  if (isGF2(h.getMod()))
    NTL::build(H.gf2x, h.gf2x, F.gf2x, m);
  else {
    NTL::zz_pPush push(getContext(h.getMod()));
    NTL::build(H.zzpx, h.zzpx, F.zzpx, m);
  }
}
// Pre-Computes information about h.  m > 0, deg(h) < n

inline void CompMod(SPX& x, const SPX& g,
		    const SPXArgument& H, const SPXModulus& F)
{
  SPX::ensureConsistency(g,H,F);
  x.resetMod(F.getMod());
  if (isGF2(F.getMod()))
    NTL::CompMod(x.gf2x, g.gf2x, H.gf2x, F.gf2x);
  else {
    NTL::zz_pPush push(getContext(F.getMod()));
    NTL::CompMod(x.zzpx, g.zzpx, H.zzpx, F.zzpx);
  }
}

inline SPX CompMod(const SPX& g, const SPXArgument& H,
		      const SPXModulus& F)
{ SPX x; CompMod(x,g,H,F); return x; }

inline void setPolyArgBound(long B)
{
  extern long GF2XArgBound;
  extern long zz_pXArgBound;
  GF2XArgBound = zz_pXArgBound = B;
}
// Initially the bounds are 0.  If this is set to a value greater than
// zero, then composition routines will allocate a table of no than
// about SPXArgBound KB.  Setting this value affects all compose
// routines and the power projection and minimal polynomial routines
// below, and indirectly affects many routines in SPXFactoring.

/**************************************************************************\

                     Power Projection routines

\**************************************************************************/

inline void project(long& x, const NTL::Vec<long>& a, const SPX& b)
{
  if (isGF2(b.getMod())) {
    NTL::GF2 x2;
    NTL::vec_GF2 a2 = NTL::conv<NTL::vec_GF2>(a);
    NTL::project(x2, a2, b.gf2x);
    NTL::conv(x,x2);
  }
  else {
    NTL::zz_pPush push(getContext(b.getMod()));
    NTL::zz_p xp;
    NTL::vec_zz_p ap = NTL::conv<NTL::vec_zz_p>(a);
    NTL::project(xp, ap, b.zzpx);
    NTL::conv(x,xp);
  }
}
inline long project(const NTL::Vec<long>& a, const SPX& b)
{ long x; project(x,a,b); return x; }
// x = inner product of a with coefficient vector of b (conv needed)


// class T is either SPX or SPXArgument
template <class T>
inline void ProjectPowers(NTL::Vec<long>& x, const NTL::Vec<long>& a, long k,
			  const T& h, const SPXModulus& F)
{
  SPX::ensureConsistency(h,F);
  if (isGF2(F.getMod())) {
    NTL::vec_GF2 x2;
    NTL::vec_GF2 a2 = NTL::conv<NTL::vec_GF2>(a);
    NTL::ProjectPowers(x2, a2, k, h.gf2x, F.gf2x);
    NTL::conv(x,x2);
  }
  else {
    NTL::zz_pPush push(getContext(F.getMod()));
    NTL::vec_zz_p xp;
    NTL::vec_zz_p ap = NTL::conv<NTL::vec_zz_p>(a);
    NTL::ProjectPowers(xp, ap, k, h.zzpx, F.zzpx);
    NTL::conv(x,xp);
  }
}

template <class T>
inline NTL::Vec<long> ProjectPowers(const NTL::Vec<long>& a, long k,
                   const T& h, const SPXModulus& F)
{ NTL::Vec<long> x; ProjectPowers(x,a,k,h,F); return x; }

// Computes the vector

//   (project(a, 1), project(a, h), ..., project(a, h^{k-1} % f).

// Restriction: must have a.length <= deg(F) and deg(h) < deg(F).
// This operation is really the "transpose" of the modular composition
// operation.

/**************************************************************************\

                              Minimum Polynomials

All of these routines implement the algorithm from [Shoup, J. Symbolic
Comp. 17:371-391, 1994] and [Shoup, J. Symbolic Comp. 20:363-397,
1995], based on transposed modular composition and the
Berlekamp/Massey algorithm.

\**************************************************************************/


inline void
MinPolySeq(SPX& h, const NTL::Vec<long>& a, long m, const SPmodulus& mod)
{
  h.resetMod(mod);
  if (isGF2(h.getMod())) {
    NTL::vec_GF2 a2 = NTL::conv<NTL::vec_GF2>(a);
    NTL::MinPolySeq(h.gf2x, a2, m);
  }
  else {
    NTL::zz_pPush push(getContext(h.getMod()));
    NTL::vec_zz_p ap = NTL::conv<NTL::vec_zz_p>(a);
    NTL::MinPolySeq(h.zzpx, ap, m);
  }
}
// computes the minimum polynomial of a linealy generated sequence; m is a
// bound on the degree of the polynomial; required: a.length() >= 2*m

inline void sppolyThreePolyArg(ProbMinPolyMod,SPXModulus)
//void ProbMinPolyMod(SPX& h,const SPX& g,const SPXModulus& F,long m);

inline SPX ProbMinPolyMod(const SPX& g, const SPXModulus& F, long m)
{ SPX h; ProbMinPolyMod(h,g,F,m); return h; }

inline void sppolyThreePoly(ProbMinPolyMod,SPXModulus)
//void ProbMinPolyMod(SPX& h, const SPX& g, const SPXModulus& F);

inline SPX ProbMinPolyMod(const SPX& g, const SPXModulus& F)
{ SPX h; ProbMinPolyMod(h,g,F); return h; }

// computes the monic minimal polynomial of (g mod f).  m is a bound on
// the degree of the minimal polynomial; in the second version, this
// argument defaults to n.  The algorithm is probabilistic; it always
// returns a divisor of the minimal polynomial, possibly a proper divisor.

inline void sppolyThreePolyArg(MinPolyMod,SPXModulus)
//void MinPolyMod(SPX& h, const SPX& g, const SPXModulus& F, long m);

SPX MinPolyMod(const SPX& g, const SPXModulus& F, long m)
{ SPX h; MinPolyMod(h,g,F,m); return h; }

inline void sppolyThreePoly(MinPolyMod,SPXModulus)
//void MinPolyMod(SPX& h, const SPX& g, const SPXModulus& F);

inline SPX MinPolyMod(const SPX& g, const SPXModulus& F)
{ SPX h; MinPolyMod(h,g,F); return h; }


// same as above, but guarantees that result is correct

inline void sppolyThreePolyArg(IrredPolyMod,SPXModulus)
//void IrredPolyMod(SPX& h, const SPX& g, const SPXModulus& F, long m)

inline SPX IrredPolyMod(const SPX& g, const SPXModulus& F, long m)
{ SPX h; IrredPolyMod(h,g,F,m); return h; }

inline void sppolyThreePoly(IrredPolyMod,SPXModulus)
//void IrredPolyMod(SPX& h, const SPX& g, const SPXModulus& F);

inline SPX IrredPolyMod(const SPX& g, const SPXModulus& F)
{ SPX h; IrredPolyMod(h,g,F); return h; }


// same as above, but assumes that F is irreducible, or at least that
// the minimal poly of g is itself irreducible.  The algorithm is
// deterministic (and is always correct).


/**************************************************************************\

                                Traces

\**************************************************************************/

// class T can be either SPX or SPXModulus
template<class T>
void TraceMod(long& x, const SPX& a, const T& F)
{
  SPX::ensureConsistency(a,F);
  if (isGF2(a.getMod())) {
    NTL::GF2 x2;
    NTL::TraceMod(x2, a.gf2x, F.gf2x);
    NTL::conv(x,x2);
  } else {
    NTL::zz_pPush push(getContext(a.getMod()));
    NTL::zz_p xp;
    NTL::TraceMod(xp, a.zzpx, F.zzpx);
    NTL::conv(x,xp);
  }
}
template<class T>
long TraceMod(const SPX& a, const T& F)
{ long x; TraceMod(x,a,F); return x; }

//void TraceMod(long& x, const SPX& a, const SPX& f);
//long TraceMod(const SPX& a, const SPX& f);
// x = Trace(a mod f); deg(a) < deg(f)

inline void TraceVec(NTL::Vec<long>& S, const SPX& f)
{
  if (isGF2(f.getMod())) {
    NTL::vec_GF2 S2;
    TraceVec(S2, f.gf2x);
    NTL::conv(S,S2);
  }
  else {
    NTL::zz_pPush push(getContext(f.getMod()));
    NTL::vec_zz_p Sp;
    TraceVec(Sp, f.zzpx);
    NTL::conv(S,Sp);
  }
}

inline NTL::Vec<long> TraceVec(const SPX& f)
{ NTL::Vec<long> S; TraceVec(S,f); return S; }
// S[i] = Trace(X^i mod f), i = 0..deg(f)-1; 0 < deg(f)

// The above routines implement the asymptotically fast trace
// algorithm from [von zur Gathen and Shoup, Computational Complexity,
// 1992].


/**************************************************************************\
                           Miscellany
\**************************************************************************/

inline void sppolyOnePoly(clear)
//void clear(SPX& x); // x = 0

inline void sppolyOnePoly(set)
//void set(SPX& x);   // x = 1

inline void sppolyTwoPoly(swap, SPX)
//void swap(SPX& x, SPX& y);
// swap (via "pointer swapping" -- if possible)


/**************************************************************************\
                           Factoring
\**************************************************************************/
#include <NTL/GF2XFactoring.h>
#include <NTL/pair_GF2X_long.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/pair_zz_pX_long.h>

typedef NTL::Pair<SPX,long> pair_SPX_long;
typedef NTL::Vec<pair_SPX_long> vec_pair_SPX_long;
/*
void SquareFreeDecomp(vec_pair_SPX_long& u, const SPX& f);
inline vec_pair_SPX_long SquareFreeDecomp(const SPX& f)
{ vec_pair_SPX_long vp; SquareFreeDecomp(vp, f); return vp; }

// Performs square-free decomposition. f must be monic, if f = prod_i gi^ei
// then u is set to a list of pairs (gi,ei). The list is is increasing order
// of ei, with trivial terms (i.e., gi=1) deleted.


void DDF(vec_pair_SPX_long& factors, const SPX& f, long verbose=0);
inline vec_pair_SPX_long DDF(const SPX& f, long verbose=0)
{ vec_pair_SPX_long vp; DDF(vp, f, verbose); return vp; }

// This computes a distinct-degree factorization.  The input must be
// monic and square-free.  factors is set to a list of pairs (g, d),
// where g is the product of all irreducible factors of f of degree d.
// Only nontrivial pairs (i.e., g != 1) are included.
*/


inline void EDF(vec_SPX& factors, const SPX& f, long d, long verbose=0)
{
  SPX tmp(f.getMod());
  factors.kill(); // remove everything from factors (if any)
  if (isGF2(f.getMod())) {
    NTL::vec_GF2X fac2;
    NTL::EDF(fac2, f.gf2x, d, verbose);
    factors.SetLength(fac2.length(), tmp); // allocate space
    for (long i=0; i<factors.length(); i++)
      factors[i].gf2x = fac2[i];
  } else {
    NTL::zz_pPush push(getContext(f.getMod()));
    NTL::vec_zz_pX facp;
    NTL::EDF(facp, f.zzpx, PowerXMod(NTL::zz_p::modulus(),f.zzpx), d, verbose);
    factors.SetLength(facp.length(), tmp); // allocate space
    for (long i=0; i<factors.length(); i++)
      factors[i].zzpx = facp[i];
  }
}
inline vec_SPX EDF(const SPX& f, long d, long verbose=0)
{ vec_SPX factors;  EDF(factors, f, d, verbose); return factors; }

// Performs equal-degree factorization.  f is monic, square-free,
// and all irreducible factors have same degree.  d = degree of
// irreducible factors of f

inline void SFCanZass(vec_SPX& factors, const SPX& f, long verbose=0)
{
  SPX tmp(f.getMod());
  factors.kill(); // remove everything from factors (if any)
  if (isGF2(f.getMod())) {
    NTL::vec_GF2X fac2;
    NTL::SFCanZass(fac2, f.gf2x, verbose);
    factors.SetLength(fac2.length(), tmp); // allocate space
    for (long i=0; i<factors.length(); i++)
      factors[i].gf2x = fac2[i];
  } else {
    NTL::zz_pPush push(getContext(f.getMod()));
    NTL::vec_zz_pX facp;
    NTL::SFCanZass(facp, f.zzpx, verbose);
    factors.SetLength(facp.length(), tmp); // allocate space
    for (long i=0; i<factors.length(); i++)
      factors[i].zzpx = facp[i];
  }
}
inline vec_SPX SFCanZass(const SPX& f, long verbose=0)
{ vec_SPX factors; SFCanZass(factors, f, verbose); return factors; }

// Assumes f is monic and square-free.  returns list of factors of f.


inline void CanZass(vec_pair_SPX_long& factors, const SPX& f, long verbose=0)
{
  pair_SPX_long tmp = NTL::Pair<SPX,long>(SPX(f.getMod()),0);
  factors.kill(); // remove everything from factors (if any)
  if (isGF2(f.getMod())) {
    NTL::vec_pair_GF2X_long fac2;
    NTL::CanZass(fac2, f.gf2x, verbose);
    factors.SetLength(fac2.length(), tmp); // allocate space
    for (long i=0; i<factors.length(); i++) {
      factors[i].a.gf2x = fac2[i].a;
      factors[i].b = fac2[i].b;
    }
  } else {
    NTL::zz_pPush push(getContext(f.getMod()));
    NTL::vec_pair_zz_pX_long facp;
    NTL::CanZass(facp, f.zzpx, verbose);
    factors.SetLength(facp.length(), tmp); // allocate space
    for (long i=0; i<factors.length(); i++) {
      factors[i].a.zzpx = facp[i].a;
      factors[i].b = facp[i].b;
    }
  }
}
inline vec_pair_SPX_long CanZass(const SPX& f, long verbose=0)
{ vec_pair_SPX_long factors; CanZass(factors, f, verbose); return factors; }

// returns a list of factors, with multiplicities.  f must be monic.
// Calls SquareFreeDecomp and SFCanZass.


inline void mul(SPX& f, const vec_pair_SPX_long& v, const SPmodulus& mod)
{
  f.resetMod(mod);
  if (isGF2(f.getMod())) {
    NTL::vec_pair_GF2X_long v2;
    v2.SetLength(v.length());
    for (long i=0; i<v2.length(); i++) {
      v2[i].a = v[i].a.gf2x;
      v2[i].b = v[i].b;
    }
    NTL::mul(f.gf2x, v2);
  } else {
    NTL::zz_pPush push(getContext(f.getMod()));
    NTL::vec_pair_zz_pX_long vp;
    vp.SetLength(v.length());
    for (long i=0; i<vp.length(); i++) {
      vp[i].a = v[i].a.zzpx;
      vp[i].b = v[i].b;
    }
    NTL::mul(f.zzpx, vp);
  }
}
inline SPX mul(const vec_pair_SPX_long& v, const SPmodulus& mod)
{ SPX f; mul(f,v,mod); return f; }

// multiplies polynomials, with multiplicities


/**************************************************************************\
                            Irreducible Polynomials
\**************************************************************************/

inline bool IterIrredTest(const SPX& f)
{
  if (isGF2(f.getMod()))
    return NTL::IterIrredTest(f.gf2x);
  else {
    NTL::zz_pPush push(getContext(f.getMod()));
    return NTL::IterIrredTest(f.zzpx);
  }
}

// performs an iterative deterministic irreducibility test, based on
// DDF.  Fast on average (when f has a small factor).

/**
inline void BuildSparseIrred(SPX& f, long n, const SPmodulus& mod)
{
  f.resetMod(mod);
  if (isGF2(mod))
    NTL::BuildSparseIrred(f.gf2x, n);
  else {
    NTL::zz_pPush push(getContext(f.getMod()));
    NTL::BuildSparseIrred(f.zzpx, n);
  }
}
inline SPX BuildSparseIrred_SPX(long n, const SPmodulus& mod)
{ SPX f; BuildSparseIrred(f,n,mod); return f; }
**/

// Builds a monic irreducible polynomial of degree n.
// If there is an irreducible trinomial X^n + X^k + 1,
// then the one with minimal k is chosen.
// Otherwise, if there is an irreducible pentanomial
// X^n + X^k3 + X^k2 + x^k1 + 1, then the one with minimal
// k3 is chosen (minimizing first k2 and then k1).
// Otherwise, if there is niether an irreducible trinomial
// or pentanomial, the routine result from BuildIrred (see below)
// is chosen---this is probably only of academic interest,
// since it a reasonable, but unproved, conjecture that they
// exist for every n > 1.

// For n <= 2048, the polynomial is constructed
// by table lookup in a pre-computed table.

// The SPXModulus data structure and routines (and indirectly GF2E)
// are optimized to deal with the output from BuildSparseIrred.

inline void BuildIrred(SPX& f, long n, const SPmodulus& mod)
{
  f.resetMod(mod);
  if (isGF2(mod))
    NTL::BuildIrred(f.gf2x, n);
  else {
    NTL::zz_pPush push(getContext(f.getMod()));
    NTL::BuildIrred(f.zzpx, n);
  }
}

inline SPX BuildIrred_SPX(long n, const SPmodulus& mod)
{ SPX f; BuildIrred(f,n,mod); return f; }

// Build a monic irreducible poly of degree n.  The polynomial
// constructed is "canonical" in the sense that it is of the form
// f=X^n + g, where the bits of g are the those of the smallest
// non-negative integer that make f irreducible.

// The SPXModulus data structure and routines (and indirectly GF2E)
// are optimized to deal with the output from BuildIrred.

// Note that the output from BuildSparseIrred will generally yield
// a "better" representation (in terms of efficiency) for
// GF(2^n) than the output from BuildIrred.



inline void BuildRandomIrred(SPX& f, const SPX& g)
{
  f.resetMod(g.getMod());
  if (isGF2(g.getMod()))
    NTL::BuildRandomIrred(f.gf2x, g.gf2x);
  else {
    NTL::zz_pPush push(getContext(g.getMod()));
    NTL::BuildRandomIrred(f.zzpx, g.zzpx);
  }
}

inline SPX BuildRandomIrred(const SPX& g)
{ SPX f; BuildRandomIrred(f, g); return f; }

// g is a monic irreducible polynomial.  Constructs a random monic
// irreducible polynomial f of the same degree.

#endif // #ifdef _HELIB_SPX_H_
