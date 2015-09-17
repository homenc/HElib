

#ifndef NTL_ZZ__H
#define NTL_ZZ__H




/********************************************************

   LIP INTERFACE 

   The class ZZ implements signed, arbitrary length integers.

**********************************************************/


#include <NTL/lip.h>
#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/SmartPtr.h>
#include <NTL/sp_arith.h>

NTL_OPEN_NNS


class ZZ_p; // forward declaration
class ZZX;

class ZZ {
public:
typedef ZZ_p residue_type;
typedef ZZX poly_type;



class Deleter {
public:
   static void apply(NTL_verylong& p) { NTL_zfree(&p); }
};

WrappedPtr<NTL_verylong_body, Deleter> rep;
// This is currently public for "emergency" situations
// May be private in future versions.

ZZ() { }


explicit ZZ(long a) { *this = a; }

ZZ(INIT_SIZE_TYPE, long k) 
// initial value is 0, but space is pre-allocated so that numbers
// x with x.size() <= k can be stored without re-allocation.
// Call with ZZ(INIT_SIZE, k).
// The purpose for the INIT_SIZE argument is to prevent automatic
// type conversion from long to ZZ, which would be tempting, but wrong.


{
   NTL_zsetlength(&rep, k); 
}

ZZ(const ZZ& a) 
// initial value is a.

{
   NTL_zcopy(a.rep, &rep);
}


ZZ(INIT_VAL_TYPE, long a)  { NTL_zintoz(a, &rep); }
ZZ(INIT_VAL_TYPE, int a)  { NTL_zintoz(a, &rep); }

ZZ(INIT_VAL_TYPE, unsigned long a)  { NTL_zuintoz(a, &rep); }
ZZ(INIT_VAL_TYPE, unsigned int a)  { NTL_zuintoz((unsigned long) a, &rep); }

inline ZZ(INIT_VAL_TYPE, const char *);
inline ZZ(INIT_VAL_TYPE, float);
inline ZZ(INIT_VAL_TYPE, double);


ZZ& operator=(const ZZ& a) { NTL_zcopy(a.rep, &rep); return *this; }

ZZ& operator=(long a) { NTL_zintoz(a, &rep); return *this; }


void kill()
// force the space held by this ZZ to be released.
// The value then becomes 0.

{ rep.kill(); }


void swap(ZZ& x)
{ NTL_zswap(&rep, &x.rep); }

void SetSize(long k)
// pre-allocates space for k-digit numbers (base 2^NTL_ZZ_NBITS);  
// does not change the value.

{ NTL_zsetlength(&rep, k); }

long size() const
// returns the number of (NTL_ZZ_NBIT-bit) digits of |a|; the size of 0 is 0.
   { return NTL_zsize(rep); }

long null() const
// test of rep is null
   { return !rep; }

long MaxAlloc() const
// returns max allocation request, possibly rounded up a bit...
   { return NTL_zmaxalloc(rep); }


long SinglePrecision() const
   { return NTL_zsptest(rep); }

// tests if less than NTL_SP_BOUND in absolute value

long WideSinglePrecision() const
   { return NTL_zwsptest(rep); }

// tests if less than NTL_WSP_BOUND in absolute value

static const ZZ& zero();


ZZ(ZZ& x, INIT_TRANS_TYPE) { rep.swap(x.rep); }
// used to cheaply hand off memory management of return value,
// without copying, assuming compiler implements the
// "return value optimization".  This is probably obsolete by
// now, as modern compilers can and should optimize
// the copy constructor in the situations where this is used.
// This should only be used for simple, local variables
// that are not be subject to special memory management.


// mainly for internal consumption by ZZWatcher

void KillBig() { if (MaxAlloc() > NTL_RELEASE_THRESH) kill(); }

};


class ZZWatcher {
public:
   ZZ& watched;
   explicit
   ZZWatcher(ZZ& _watched) : watched(_watched) {}

   ~ZZWatcher() { watched.KillBig(); }
};

#define NTL_ZZRegister(x) NTL_THREAD_LOCAL static ZZ x; ZZWatcher _WATCHER__ ## x(x)





const ZZ& ZZ_expo(long e);


inline void clear(ZZ& x)
// x = 0

   { NTL_zzero(&x.rep); }

inline void set(ZZ& x)
// x = 1

   { NTL_zone(&x.rep); }


inline void swap(ZZ& x, ZZ& y)
// swap the values of x and y (swaps pointers only)

   { x.swap(y); }


inline double log(const ZZ& a)
   { return NTL_zlog(a.rep); }




/**********************************************************

   Conversion routines.

***********************************************************/



inline void conv(ZZ& x, const ZZ& a) { x = a; }
inline ZZ to_ZZ(const ZZ& a) { return a; }

inline void conv(ZZ& x, long a) { NTL_zintoz(a, &x.rep); }
inline ZZ to_ZZ(long a) { return ZZ(INIT_VAL, a); }

inline void conv(ZZ& x, int a) { NTL_zintoz(long(a), &x.rep); }
inline ZZ to_ZZ(int a) { return ZZ(INIT_VAL, a); }

inline void conv(ZZ& x, unsigned long a) { NTL_zuintoz(a, &x.rep); }
inline ZZ to_ZZ(unsigned long a) { return ZZ(INIT_VAL, a); }

inline void conv(ZZ& x, unsigned int a) { NTL_zuintoz((unsigned long)(a), &x.rep); }
inline ZZ to_ZZ(unsigned int a) { return ZZ(INIT_VAL, a); }

void conv(ZZ& x, const char *s);
inline ZZ::ZZ(INIT_VAL_TYPE, const char *s)  { conv(*this, s); }
inline ZZ to_ZZ(const char *s) { return ZZ(INIT_VAL, s); }

inline void conv(ZZ& x, double a) { NTL_zdoubtoz(a, &x.rep); }
inline ZZ::ZZ(INIT_VAL_TYPE, double a)  { conv(*this, a); }
inline ZZ to_ZZ(double a) { return ZZ(INIT_VAL, a); }

inline void conv(ZZ& x, float a) { NTL_zdoubtoz(double(a), &x.rep); }
inline ZZ::ZZ(INIT_VAL_TYPE, float a)  { conv(*this, a); }
inline ZZ to_ZZ(float a) { return ZZ(INIT_VAL, a); }

inline void conv(long& x, const ZZ& a) { x = NTL_ztoint(a.rep); }
inline long to_long(const ZZ& a)  { return NTL_ztoint(a.rep); }

inline void conv(int& x, const ZZ& a) 
   { unsigned int res = (unsigned int) NTL_ztouint(a.rep); 
     x = NTL_UINT_TO_INT(res); }

inline int to_int(const ZZ& a)  
   { unsigned int res = (unsigned int) NTL_ztouint(a.rep); 
     return NTL_UINT_TO_INT(res); }

inline void conv(unsigned long& x, const ZZ& a) { x = NTL_ztouint(a.rep); }
inline unsigned long to_ulong(const ZZ& a)  { return NTL_ztouint(a.rep); }

inline void conv(unsigned int& x, const ZZ& a) 
   { x = (unsigned int)(NTL_ztouint(a.rep)); }
inline unsigned int to_uint(const ZZ& a)  
   { return (unsigned int)(NTL_ztouint(a.rep)); }

inline void conv(double& x, const ZZ& a) { x = NTL_zdoub(a.rep); }
inline double to_double(const ZZ& a) { return NTL_zdoub(a.rep); }

inline void conv(float& x, const ZZ& a) { x = float(NTL_zdoub(a.rep)); }
inline float to_float(const ZZ& a) { return float(NTL_zdoub(a.rep)); }

inline void ZZFromBytes(ZZ& x, const unsigned char *p, long n)
   { NTL_zfrombytes(&x.rep, p, n); }

inline ZZ ZZFromBytes(const unsigned char *p, long n)
   { ZZ x; ZZFromBytes(x, p, n); NTL_OPT_RETURN(ZZ, x); }

inline void BytesFromZZ(unsigned char *p, const ZZ& a, long n)
   { NTL_zbytesfromz(p, a.rep, n); }




// ****** comparisons


inline long sign(const ZZ& a)
// returns the sign of a (-1, 0, or 1).

   { return NTL_zsign(a.rep); }


inline long compare(const ZZ& a, const ZZ& b)
// returns the sign of a-b (-1, 0, or 1).

{
   return NTL_zcompare(a.rep, b.rep);
}

inline long IsZero(const ZZ& a)
// zero test

   { return NTL_ziszero(a.rep); }


inline long IsOne(const ZZ& a)
   { return NTL_zisone(a.rep); }
// test for 1
   

/* the usual comparison operators */

inline long operator==(const ZZ& a, const ZZ& b)
  { return NTL_zcompare(a.rep, b.rep) == 0; }
inline long operator!=(const ZZ& a, const ZZ& b)
  { return NTL_zcompare(a.rep, b.rep) != 0; }
inline long operator<(const ZZ& a, const ZZ& b)
  { return NTL_zcompare(a.rep, b.rep) < 0; }
inline long operator>(const ZZ& a, const ZZ& b)
  { return NTL_zcompare(a.rep, b.rep) > 0; }
inline long operator<=(const ZZ& a, const ZZ& b)
  { return NTL_zcompare(a.rep, b.rep) <= 0; }
inline long operator>=(const ZZ& a, const ZZ& b)
  { return NTL_zcompare(a.rep, b.rep) >= 0; }

/* single-precision versions of the above */

inline long compare(const ZZ& a, long b) { return NTL_zscompare(a.rep, b); }
inline long compare(long a, const ZZ& b) { return -NTL_zscompare(b.rep, a); }

inline long operator==(const ZZ& a, long b) { return NTL_zscompare(a.rep, b) == 0; }
inline long operator!=(const ZZ& a, long b) { return NTL_zscompare(a.rep, b) != 0; }
inline long operator<(const ZZ& a, long b) { return NTL_zscompare(a.rep, b) < 0; }
inline long operator>(const ZZ& a, long b) { return NTL_zscompare(a.rep, b) > 0; }
inline long operator<=(const ZZ& a, long b) { return NTL_zscompare(a.rep, b) <= 0; }
inline long operator>=(const ZZ& a, long b) { return NTL_zscompare(a.rep, b) >= 0; }


inline long operator==(long a, const ZZ& b) { return b == a; }
inline long operator!=(long a, const ZZ& b) { return b != a; }
inline long operator<(long a, const ZZ& b) { return b > a; }
inline long operator>(long a, const ZZ& b) { return b < a; }
inline long operator<=(long a, const ZZ& b) { return b >= a; }
inline long operator>=(long a, const ZZ& b) { return b <= a; }

/**************************************************

                 Addition

**************************************************/


inline void add(ZZ& x, const ZZ& a, const ZZ& b)
// x = a + b

   { NTL_zadd(a.rep, b.rep, &x.rep); }

inline void sub(ZZ& x, const ZZ& a, const ZZ& b)
// x = a - b

   { NTL_zsub(a.rep, b.rep, &x.rep); }

inline void SubPos(ZZ& x, const ZZ& a, const ZZ& b)
// x = a - b;  assumes a >= b >= 0.

   { NTL_zsubpos(a.rep, b.rep, &x.rep); }

inline void negate(ZZ& x, const ZZ& a)
// x = -a

   { NTL_zcopy(a.rep, &x.rep); NTL_znegate(&x.rep); }

inline void abs(ZZ& x, const ZZ& a)
// x = |a|
{ NTL_zcopy(a.rep, &x.rep); NTL_zabs(&x.rep); }


/* single-precision versions of the above */

inline void add(ZZ& x, const ZZ& a, long b)
   { NTL_zsadd(a.rep, b, &x.rep); }

inline void add(ZZ& x, long a, const ZZ& b) { add(x, b, a); }


void sub(ZZ& x, const ZZ& a, long b);
void sub(ZZ& x, long a, const ZZ& b);

/* operator/function notation */

inline ZZ operator+(const ZZ& a, const ZZ& b) 
  { ZZ x; add(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator+(const ZZ& a, long b) 
  { ZZ x; add(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator+(long  a, const ZZ& b) 
  { ZZ x; add(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(const ZZ& a, const ZZ& b) 
  { ZZ x; sub(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(const ZZ& a, long b) 
  { ZZ x; sub(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(long  a, const ZZ& b) 
  { ZZ x; sub(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(const ZZ& a)
  { ZZ x; negate(x, a); NTL_OPT_RETURN(ZZ, x); }

inline ZZ abs(const ZZ& a)
  { ZZ x; abs(x, a); NTL_OPT_RETURN(ZZ, x); }

/* op= notation */

inline ZZ& operator+=(ZZ& x, const ZZ& a)
  { add(x, x, a); return x; }

inline ZZ& operator+=(ZZ& x, long a)
  { add(x, x, a); return x; }

inline ZZ& operator-=(ZZ& x, const ZZ& a)
  { sub(x, x, a); return x; }

inline ZZ& operator-=(ZZ& x, long a)
  { sub(x, x, a); return x; }

/* inc/dec */

inline ZZ& operator++(ZZ& x) { add(x, x, 1); return x; }

inline void operator++(ZZ& x, int) { add(x, x, 1); }

inline ZZ& operator--(ZZ& x) { add(x, x, -1); return x; }

inline void operator--(ZZ& x, int) { add(x, x, -1); }



/*******************************************************

                 Multiplication.

********************************************************/


inline void mul(ZZ& x, const ZZ& a, const ZZ& b)
// x = a * b

   { NTL_zmul(a.rep, b.rep, &x.rep); }


inline void sqr(ZZ& x, const ZZ& a)
// x = a*a

   { NTL_zsq(a.rep, &x.rep); }

inline ZZ sqr(const ZZ& a)
   { ZZ x; sqr(x, a); NTL_OPT_RETURN(ZZ, x); }


/* single-precision versions */

inline void mul(ZZ& x, const ZZ& a, long b)
   { NTL_zsmul(a.rep, b, &x.rep); }

inline void mul(ZZ& x, long a, const ZZ& b)
    { mul(x, b, a); }

/* operator notation */

inline ZZ operator*(const ZZ& a, const ZZ& b)
  { ZZ x; mul(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator*(const ZZ& a, long b)
  { ZZ x; mul(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator*(long a, const ZZ& b)
  { ZZ x; mul(x, a, b); NTL_OPT_RETURN(ZZ, x); }

/* op= notation */

inline ZZ& operator*=(ZZ& x, const ZZ& a)
  { mul(x, x, a); return x; }

inline ZZ& operator*=(ZZ& x, long a)
  { mul(x, x, a); return x; }

// x += a*b

inline void
MulAddTo(ZZ& x, const ZZ& a, long b)
{
   NTL_zsaddmul(a.rep, b, &x.rep);
}

inline void
MulAddTo(ZZ& x, const ZZ& a, const ZZ& b)
{
   NTL_zaddmul(a.rep, b.rep, &x.rep);
}

// x -= a*b

inline void
MulSubFrom(ZZ& x, const ZZ& a, long b)
{
   NTL_zssubmul(a.rep, b, &x.rep);
}

inline void
MulSubFrom(ZZ& x, const ZZ& a, const ZZ& b)
{
   NTL_zsubmul(a.rep, b.rep, &x.rep);
}



// Special routines for implementing CRT in ZZ_pX arithmetic
// These are verbose, but fairly boilerplate



class ZZ_CRTStructAdapter;
class ZZ_RemStructAdapter;

class ZZ_TmpVecAdapter {
public:
   UniquePtr<_ntl_tmp_vec> rep;

   inline void fetch(const ZZ_CRTStructAdapter&);
   inline void fetch(ZZ_CRTStructAdapter&);
   inline void fetch(const ZZ_RemStructAdapter&);
};


class ZZ_CRTStructAdapter {
public:
   UniquePtr<_ntl_crt_struct> rep;

   void init(long n, const ZZ& p, long (*primes)(long))
   {
      rep.reset(_ntl_crt_struct_build(n, p.rep, primes));
   }

   void insert(long i, const ZZ& m)
   {
       rep->insert(i, m.rep);
   }

   void eval(ZZ& t, const long *a, ZZ_TmpVecAdapter& tmp_vec) const
   {
      rep->eval(&t.rep, a, tmp_vec.rep.get());
   }

   bool special() const
   { 
      return rep->special();
   }
};


class ZZ_RemStructAdapter {
public:
   UniquePtr<_ntl_rem_struct> rep;

   void init(long n, const ZZ& p, long (*primes)(long))
   {
      rep.reset(_ntl_rem_struct_build(n, p.rep, primes));
   }

   void eval(long *x, const ZZ& a, ZZ_TmpVecAdapter& tmp_vec) const
   {
      rep->eval(x, a.rep, tmp_vec.rep.get());
   }
};


inline void ZZ_TmpVecAdapter::fetch(const ZZ_CRTStructAdapter& crt_struct)
{
   rep.reset(crt_struct.rep->fetch()); 
}

inline void ZZ_TmpVecAdapter::fetch(ZZ_CRTStructAdapter& crt_struct)
{
   rep.reset(crt_struct.rep->extract()); // EXTRACT!!
}


inline void ZZ_TmpVecAdapter::fetch(const ZZ_RemStructAdapter& rem_struct)
{
   rep.reset(rem_struct.rep->fetch()); 
}




/*******************************************************

                    Division

*******************************************************/


inline void DivRem(ZZ& q, ZZ& r, const ZZ& a, const ZZ& b)
// q = [a/b], r = a - b*q
// |r| < |b|, and if r != 0, sign(r) = sign(b)

   { NTL_zdiv(a.rep, b.rep, &q.rep, &r.rep); }



inline void div(ZZ& q, const ZZ& a, const ZZ& b)
// q = a/b

   { NTL_zdiv(a.rep, b.rep, &q.rep, 0); }

inline void rem(ZZ& r, const ZZ& a, const ZZ& b)
// r = a%b

   { NTL_zmod(a.rep, b.rep, &r.rep); }


inline void QuickRem(ZZ& r, const ZZ& b)
// r = r%b
// assumes b > 0 and r >=0
// division is performed in place and may cause r to be re-allocated.

   { NTL_zquickmod(&r.rep, b.rep); }

long divide(ZZ& q, const ZZ& a, const ZZ& b);
// if b | a, sets q = a/b and returns 1; otherwise returns 0.

long divide(const ZZ& a, const ZZ& b);
// if b | a, returns 1; otherwise returns 0.


/* non-standard single-precision versions */

inline long DivRem(ZZ& q, const ZZ& a, long b)
   { return NTL_zsdiv(a.rep, b, &q.rep); } 

inline long rem(const ZZ& a, long b)
   { return NTL_zsmod(a.rep, b); }


/* single precision versions */

inline void div(ZZ& q, const ZZ& a, long b)
   { (void) NTL_zsdiv(a.rep, b, &q.rep); }


long divide(ZZ& q, const ZZ& a, long b);
// if b | a, sets q = a/b and returns 1; otherwise returns 0.

long divide(const ZZ& a, long b);
// if b | a, returns 1; otherwise returns 0.


inline ZZ operator/(const ZZ& a, const ZZ& b)
   { ZZ x; div(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator/(const ZZ& a, long b)
   { ZZ x; div(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator%(const ZZ& a, const ZZ& b)
   { ZZ x; rem(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline long operator%(const ZZ& a, long b)
   { return rem(a, b); }

inline ZZ& operator/=(ZZ& x, const ZZ& b)
   { div(x, x, b); return x; } 

inline ZZ& operator/=(ZZ& x, long b)
   { div(x, x, b); return x; } 

inline ZZ& operator%=(ZZ& x, const ZZ& b)
   { rem(x, x, b); return x; } 


/**********************************************************

                        GCD's

***********************************************************/


inline void GCD(ZZ& d, const ZZ& a, const ZZ& b)
// d = gcd(a, b)

   { NTL_zgcd(a.rep, b.rep, &d.rep); }

inline ZZ GCD(const ZZ& a, const ZZ& b)
   { ZZ x; GCD(x, a, b); NTL_OPT_RETURN(ZZ, x); }


inline void XGCD(ZZ& d, ZZ& s, ZZ& t, const ZZ& a, const ZZ& b)
//  d = gcd(a, b) = a*s + b*t;

   { NTL_zexteucl(a.rep, &s.rep, b.rep, &t.rep, &d.rep); }

// single-precision versions
long GCD(long a, long b);

void XGCD(long& d, long& s, long& t, long a, long b);







/************************************************************

                      Bit Operations

*************************************************************/


inline void LeftShift(ZZ& x, const ZZ& a, long k)
// x = (a << k), k < 0 => RightShift

   { NTL_zlshift(a.rep, k, &x.rep); }

inline ZZ LeftShift(const ZZ& a, long k)
   { ZZ x; LeftShift(x, a, k); NTL_OPT_RETURN(ZZ, x); }


inline void RightShift(ZZ& x, const ZZ& a, long k)
// x = (a >> k), k < 0 => LeftShift

   { NTL_zrshift(a.rep, k, &x.rep); }

inline ZZ RightShift(const ZZ& a, long k)
   { ZZ x; RightShift(x, a, k); NTL_OPT_RETURN(ZZ, x); }

#ifndef NTL_TRANSITION

inline ZZ operator>>(const ZZ& a, long n)
   { ZZ x; RightShift(x, a, n); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator<<(const ZZ& a, long n)
   { ZZ x; LeftShift(x, a, n); NTL_OPT_RETURN(ZZ, x); }

inline ZZ& operator<<=(ZZ& x, long n)
   { LeftShift(x, x, n); return x; }

inline ZZ& operator>>=(ZZ& x, long n)
   { RightShift(x, x, n); return x; }

#endif


inline long MakeOdd(ZZ& x)
// removes factors of 2 from x, returns the number of 2's removed
// returns 0 if x == 0

   { return NTL_zmakeodd(&x.rep); }

inline long NumTwos(const ZZ& x)
// returns max e such that 2^e divides x if x != 0, and returns 0 if x == 0.

   { return NTL_znumtwos(x.rep); }


inline long IsOdd(const ZZ& a)
// returns 1 if a is odd, otherwise 0

   { return NTL_zodd(a.rep); }


inline long NumBits(const ZZ& a)
// returns the number of bits in |a|; NumBits(0) = 0
   { return NTL_z2log(a.rep); }



inline long bit(const ZZ& a, long k)
// returns bit k of a, 0 being the low-order bit

   { return NTL_zbit(a.rep, k); }

#ifndef NTL_GMP_LIP

// only defined for the "classic" long integer package, for backward
// compatability.

inline long digit(const ZZ& a, long k)
   { return NTL_zdigit(a.rep, k); }

#endif

// returns k-th digit of |a|, 0 being the low-order digit.

inline void trunc(ZZ& x, const ZZ& a, long k)
// puts k low order bits of |a| into x

   { NTL_zlowbits(a.rep, k, &x.rep); }

inline ZZ trunc_ZZ(const ZZ& a, long k)
   { ZZ x; trunc(x, a, k); NTL_OPT_RETURN(ZZ, x); }

inline long trunc_long(const ZZ& a, long k)
// returns k low order bits of |a|

   { return NTL_zslowbits(a.rep, k); }

inline long SetBit(ZZ& x, long p)
// returns original value of p-th bit of |a|, and replaces
// p-th bit of a by 1 if it was zero;
// error if p < 0 

   { return NTL_zsetbit(&x.rep, p); }

inline long SwitchBit(ZZ& x, long p)
// returns original value of p-th bit of |a|, and switches
// the value of p-th bit of a;
// p starts counting at 0;
//   error if p < 0

   { return NTL_zswitchbit(&x.rep, p); }

inline long weight(long a)
// returns Hamming weight of |a|

   { return NTL_zweights(a); }

inline long weight(const ZZ& a)
// returns Hamming weight of |a|

   { return NTL_zweight(a.rep); }

inline void bit_and(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| AND |b|

   { NTL_zand(a.rep, b.rep, &x.rep); }

void bit_and(ZZ& x, const ZZ& a, long b);
inline void bit_and(ZZ& x, long a, const ZZ& b)
   { bit_and(x, b, a); }


inline void bit_or(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| OR |b|

   { NTL_zor(a.rep, b.rep, &x.rep); }

void bit_or(ZZ& x, const ZZ& a, long b);
inline void bit_or(ZZ& x, long a, const ZZ& b)
   { bit_or(x, b, a); }

inline void bit_xor(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| XOR |b|

   { NTL_zxor(a.rep, b.rep, &x.rep); }

void bit_xor(ZZ& x, const ZZ& a, long b);
inline void bit_xor(ZZ& x, long a, const ZZ& b)
   { bit_xor(x, b, a); }


inline ZZ operator&(const ZZ& a, const ZZ& b)
   { ZZ x; bit_and(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator&(const ZZ& a, long b)
   { ZZ x; bit_and(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator&(long a, const ZZ& b)
   { ZZ x; bit_and(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator|(const ZZ& a, const ZZ& b)
   { ZZ x; bit_or(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator|(const ZZ& a, long b)
   { ZZ x; bit_or(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator|(long a, const ZZ& b)
   { ZZ x; bit_or(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator^(const ZZ& a, const ZZ& b)
   { ZZ x; bit_xor(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator^(const ZZ& a, long b)
   { ZZ x; bit_xor(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator^(long a, const ZZ& b)
   { ZZ x; bit_xor(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ& operator&=(ZZ& x, const ZZ& b) 
   { bit_and(x, x, b); return x; }

inline ZZ& operator&=(ZZ& x, long b) 
   { bit_and(x, x, b); return x; }

inline ZZ& operator|=(ZZ& x, const ZZ& b) 
   { bit_or(x, x, b); return x; }

inline ZZ& operator|=(ZZ& x, long b) 
   { bit_or(x, x, b); return x; }

inline ZZ& operator^=(ZZ& x, const ZZ& b) 
   { bit_xor(x, x, b); return x; }

inline ZZ& operator^=(ZZ& x, long b) 
   { bit_xor(x, x, b); return x; }



long NumBits(long a);

long bit(long a, long k);

long NextPowerOfTwo(long m);
// returns least nonnegative k such that 2^k >= m

inline
long NumBytes(const ZZ& a)
   { return (NumBits(a)+7)/8; }

inline
long NumBytes(long a)
   { return (NumBits(a)+7)/8; }



/***********************************************************

          Some specialized routines

************************************************************/


inline long ZZ_BlockConstructAlloc(ZZ& x, long d, long n)
   { return NTL_zblock_construct_alloc(&x.rep, d, n); }

inline void ZZ_BlockConstructSet(ZZ& x, ZZ& y, long i)
   { NTL_zblock_construct_set(x.rep, &y.rep, i); }

inline long ZZ_BlockDestroy(ZZ& x)
   { return NTL_zblock_destroy(x.rep);  }

inline long ZZ_storage(long d)
   { return NTL_zblock_storage(d); }

inline long ZZ_RoundCorrection(const ZZ& a, long k, long residual)
   { return NTL_zround_correction(a.rep, k, residual); }


/***********************************************************

                  Psuedo-random Numbers

************************************************************/


// ================ NEW PRG STUFF =================


// Low-level key-derivation 


void DeriveKey(unsigned char *hash, long hlen,  
               const unsigned char *data, long dlen);



// Low-level chacha stuff

#define NTL_PRG_KEYLEN (32)

class RandomByteStream {
private:
   _ntl_uint32 state[16];
   unsigned char buf[64];
   long pos;

   void do_get(unsigned char *res, long n); 

public:
   explicit
   RandomByteStream(const unsigned char *key);

   // No default constructor 
   // default copy and assignment

   void get(unsigned char *res, long n) 
   {
      // optimize short reads
      if (n >= 0 && n <= 64-pos) {
         long i;
         for (i = 0; i < n; i++) {
            res[i] = buf[pos+i];
         }
         pos += n;
      }
      else {
         do_get(res, n);
      }
   }

};

   


RandomByteStream& GetCurrentRandomByteStream();
// get reference to the current random by stream --
// if SetSeed has not ben called, it is called with
// a default value (which should be unique to each
// process/thread


void SetSeed(const ZZ& s);
void SetSeed(const unsigned char *data, long dlen);
void SetSeed(const RandomByteStream& s);
// initialize random number generator
// in the first two version, a PRG key is derived from
// the data using DeriveKey.


void RandomBnd(ZZ& x, const ZZ& n);
// x = "random number" in the range 0..n-1, or 0  if n <= 0

inline ZZ RandomBnd(const ZZ& n)
   { ZZ x; RandomBnd(x, n); NTL_OPT_RETURN(ZZ, x); }


void RandomLen(ZZ& x, long NumBits);
// x = "random number" with precisely NumBits bits.


inline ZZ RandomLen_ZZ(long NumBits)
   { ZZ x; RandomLen(x, NumBits); NTL_OPT_RETURN(ZZ, x); }


void RandomBits(ZZ& x, long NumBits);
// x = "random number", 0 <= x < 2^NumBits 

inline ZZ RandomBits_ZZ(long NumBits)
   { ZZ x; RandomBits(x, NumBits); NTL_OPT_RETURN(ZZ, x); }


// single-precision version of the above

long RandomBnd(long n);

long RandomLen_long(long l);

long RandomBits_long(long l);

unsigned long RandomWord();
unsigned long RandomBits_ulong(long l);




/**********************************************************

             Incremental Chinese Remaindering

***********************************************************/

long CRT(ZZ& a, ZZ& p, const ZZ& A, const ZZ& P);
long CRT(ZZ& a, ZZ& p, long A, long P);
// 0 <= A < P, (p, P) = 1;
// computes b such that b = a mod p, b = A mod p,
//   and -p*P/2 < b <= p*P/2;
// sets a = b, p = p*P, and returns 1 if a's value
//   has changed, otherwise 0

inline long CRTInRange(const ZZ& gg, const ZZ& aa)
   { return NTL_zcrtinrange(gg.rep, aa.rep); }

// an auxilliary routine used by newer CRT routines to maintain
// backward compatability.

// test if a > 0 and -a/2 < g <= a/2
// this is "hand crafted" so as not too waste too much time
// in the CRT routines.



/**********************************************************

                  Rational Reconstruction

***********************************************************/

inline
long ReconstructRational(ZZ& a, ZZ& b, const ZZ& u, const ZZ& m, 
                         const ZZ& a_bound, const ZZ& b_bound)
{
   return NTL_zxxratrecon(u.rep, m.rep, a_bound.rep, b_bound.rep, &a.rep, &b.rep);

}




/************************************************************

                      Primality Testing 

*************************************************************/


void GenPrime(ZZ& n, long l, long err = 80);
inline ZZ GenPrime_ZZ(long l, long err = 80) 
{ ZZ x; GenPrime(x, l, err); NTL_OPT_RETURN(ZZ, x); }

long GenPrime_long(long l, long err = 80);
// This generates a random prime n of length l so that the
// probability of erroneously returning a composite is bounded by 2^(-err).

void GenGermainPrime(ZZ& n, long l, long err = 80);
inline ZZ GenGermainPrime_ZZ(long l, long err = 80) 
{ ZZ x; GenGermainPrime(x, l, err); NTL_OPT_RETURN(ZZ, x); }

long GenGermainPrime_long(long l, long err = 80);
// This generates a random prime n of length l so that the


long ProbPrime(const ZZ& n, long NumTrials = 10);
// tests if n is prime;  performs a little trial division,
// followed by a single-precision MillerWitness test, followed by
// up to NumTrials general MillerWitness tests.

long MillerWitness(const ZZ& n, const ZZ& w);
// Tests if w is a witness to primality a la Miller.
// Assumption: n is odd and positive, 0 <= w < n.

void RandomPrime(ZZ& n, long l, long NumTrials=10);
// n =  random l-bit prime

inline ZZ RandomPrime_ZZ(long l, long NumTrials=10)
   { ZZ x; RandomPrime(x, l, NumTrials); NTL_OPT_RETURN(ZZ, x); }

void NextPrime(ZZ& n, const ZZ& m, long NumTrials=10);
// n = smallest prime >= m.

inline ZZ NextPrime(const ZZ& m, long NumTrials=10)
   { ZZ x; NextPrime(x, m, NumTrials); NTL_OPT_RETURN(ZZ, x); }

// single-precision versions

long ProbPrime(long n, long NumTrials = 10);


long RandomPrime_long(long l, long NumTrials=10);

long NextPrime(long l, long NumTrials=10);


/************************************************************

                      Exponentiation

*************************************************************/

inline void power(ZZ& x, const ZZ& a, long e)
   { NTL_zexp(a.rep, e, &x.rep); }

inline ZZ power(const ZZ& a, long e)
   { ZZ x; power(x, a, e); NTL_OPT_RETURN(ZZ, x); }

inline void power(ZZ& x, long a, long e)
   {  NTL_zexps(a, e, &x.rep); }

inline ZZ power_ZZ(long a, long e)
   { ZZ x; power(x, a, e); NTL_OPT_RETURN(ZZ, x); }

long power_long(long a, long e); 

void power2(ZZ& x, long e);

inline ZZ power2_ZZ(long e)
   { ZZ x; power2(x, e); NTL_OPT_RETURN(ZZ, x); }





/*************************************************************

                       Square Roots

**************************************************************/




inline void SqrRoot(ZZ& x, const ZZ& a)
// x = [a^{1/2}], a >= 0

{
   NTL_zsqrt(a.rep, &x.rep);
}

inline ZZ SqrRoot(const ZZ& a)
   { ZZ x; SqrRoot(x, a); NTL_OPT_RETURN(ZZ, x); }


inline long SqrRoot(long a) { return NTL_zsqrts(a); }
// single-precision version



/***************************************************************

                      Modular Arithmetic

***************************************************************/

// The following routines perform arithmetic mod n, n positive.
// All args (other than exponents) are assumed to be in the range 0..n-1.



inline void AddMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a+b)%n
   { NTL_zaddmod(a.rep, b.rep, n.rep, &x.rep); }
   

inline ZZ AddMod(const ZZ& a, const ZZ& b, const ZZ& n)
   { ZZ x; AddMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void SubMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a-b)%n

   { NTL_zsubmod(a.rep, b.rep, n.rep, &x.rep); }

inline ZZ SubMod(const ZZ& a, const ZZ& b, const ZZ& n)
   { ZZ x; SubMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void NegateMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = -a % n

   { NTL_zsubmod(0, a.rep, n.rep, &x.rep); }

inline ZZ NegateMod(const ZZ& a, const ZZ& n)
   { ZZ x; NegateMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }

void AddMod(ZZ& x, const ZZ& a, long b, const ZZ& n);
inline ZZ AddMod(const ZZ& a, long b, const ZZ& n)
   { ZZ x; AddMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void AddMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
   { AddMod(x, b, a, n); }
inline ZZ AddMod(long a, const ZZ& b, const ZZ& n)
   { ZZ x; AddMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

void SubMod(ZZ& x, const ZZ& a, long b, const ZZ& n);
inline ZZ SubMod(const ZZ& a, long b, const ZZ& n)
   { ZZ x; SubMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

void SubMod(ZZ& x, long a, const ZZ& b, const ZZ& n);
inline ZZ SubMod(long a, const ZZ& b, const ZZ& n)
   { ZZ x; SubMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void MulMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a*b)%n

   { NTL_zmulmod(a.rep, b.rep, n.rep, &x.rep); }

inline ZZ MulMod(const ZZ& a, const ZZ& b, const ZZ& n)
   { ZZ x; MulMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void MulMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
// x = (a*b)%n

   { NTL_zsmulmod(a.rep, b, n.rep, &x.rep); }

inline ZZ MulMod(const ZZ& a, long b, const ZZ& n)
   { ZZ x; MulMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void MulMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
   { MulMod(x, b, a, n); }

inline ZZ MulMod(long a, const ZZ& b, const ZZ& n)
   { ZZ x; MulMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }


inline void SqrMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = a^2 % n

   { NTL_zsqmod(a.rep, n.rep, &x.rep); }

inline ZZ SqrMod(const ZZ& a, const ZZ& n)
   {  ZZ x; SqrMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }

void InvMod(ZZ& x, const ZZ& a, const ZZ& n);
// defined in ZZ.c in terms of InvModStatus

inline ZZ InvMod(const ZZ& a, const ZZ& n)
   {  ZZ x; InvMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }


inline long InvModStatus(ZZ& x, const ZZ& a, const ZZ& n)
// if gcd(a,b) = 1, then ReturnValue = 0, x = a^{-1} mod n
// otherwise, ReturnValue = 1, x = gcd(a, n)

  { return NTL_zinv(a.rep, n.rep, &x.rep); }


void PowerMod(ZZ& x, const ZZ& a, const ZZ& e, const ZZ& n);
// defined in ZZ.c in terms of LowLevelPowerMod

inline void LowLevelPowerMod(ZZ& x, const ZZ& a, const ZZ& e, const ZZ& n)
   { NTL_zpowermod(a.rep, e.rep, n.rep, &x.rep); }

inline ZZ PowerMod(const ZZ& a, const ZZ& e, const ZZ& n)
   { ZZ x; PowerMod(x, a, e, n); NTL_OPT_RETURN(ZZ, x); }

inline void PowerMod(ZZ& x, const ZZ& a, long e, const ZZ& n)
   { PowerMod(x, a, ZZ_expo(e), n); }

inline ZZ PowerMod(const ZZ& a, long e, const ZZ& n)
   { ZZ x; PowerMod(x, a, e, n); NTL_OPT_RETURN(ZZ, x); }






/*************************************************************

   Jacobi symbol and modular squre roots

**************************************************************/


long Jacobi(const ZZ& a, const ZZ& n);
//  compute Jacobi symbol of a and n;
//  assumes 0 <= a < n, n odd

void SqrRootMod(ZZ& x, const ZZ& a, const ZZ& n);
//  computes square root of a mod n;
//  assumes n is an odd prime, and that a is a square mod n

inline ZZ SqrRootMod(const ZZ& a, const ZZ& n)
   { ZZ x; SqrRootMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }




/*************************************************************


                    Small Prime Generation


*************************************************************/


// primes are generated in sequence, starting at 2, 
// and up until (2*NTL_PRIME_BND+1)^2, which is less than NTL_SP_BOUND.

#if (NTL_SP_NBITS > 30)
#define NTL_PRIME_BND ((1L << 14) - 1)
#else
#define NTL_PRIME_BND ((1L << (NTL_SP_NBITS/2-1)) - 1)
#endif


class PrimeSeq {

const char *movesieve;
Vec<char> movesieve_mem;
long pindex;
long pshift;
long exhausted;

public:

PrimeSeq();

long next();
// returns next prime in the sequence.
// returns 0 if list of small primes is exhausted.

void reset(long b);
// resets generator so that the next prime in the sequence
// is the smallest prime >= b.

private:

PrimeSeq(const PrimeSeq&);        // disabled
void operator=(const PrimeSeq&);  // disabled

// auxilliary routines

void start();
void shift(long);

};




/**************************************************************

                      Input/Output

***************************************************************/

NTL_SNS istream& operator>>(NTL_SNS istream& s, ZZ& x);  
NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const ZZ& a); 




// Some additional SP arithmetic routines, not defined in sp_arith.h


long InvMod(long a, long n);
// computes a^{-1} mod n.  Error is raised if undefined.

long PowerMod(long a, long e, long n);
// computes a^e mod n, e >= 0


// Error handling

#ifdef NTL_EXCEPTIONS

class InvModErrorObject : public ArithmeticErrorObject {
private:
   SmartPtr<ZZ> a_ptr;
   SmartPtr<ZZ> n_ptr;
public:
   InvModErrorObject(const char *s, const ZZ& a, const ZZ& n)
      : ArithmeticErrorObject(s) , a_ptr(MakeSmart<ZZ>(a)),
        n_ptr(MakeSmart<ZZ>(n)) { }

   const ZZ& get_a() const { return *a_ptr; }
   const ZZ& get_n() const { return *n_ptr; }
};

#else

// We need this alt definition to keep pre-C++11 
// compilers happy (NTL_EXCEPTIONS should only be used
// with C++11 compilers).

class InvModErrorObject : public ArithmeticErrorObject {
public:
   InvModErrorObject(const char *s, const ZZ& a, const ZZ& n)
      : ArithmeticErrorObject(s) { }

   const ZZ& get_a() const { return ZZ::zero(); }
   const ZZ& get_n() const { return ZZ::zero(); }
};

#endif


void InvModError(const char *s, const ZZ& a, const ZZ& n); 

NTL_CLOSE_NNS


#endif

