
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/Lazy.h>
#include <NTL/fileio.h>

#include <cstring>



NTL_START_IMPL





const ZZ& ZZ::zero()
{
   NTL_THREAD_LOCAL static ZZ z;
   return z;
}


const ZZ& ZZ_expo(long e)
{
   NTL_THREAD_LOCAL static ZZ expo_helper;
   conv(expo_helper, e);
   return expo_helper;
}



void AddMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   NTL_ZZRegister(B);
   conv(B, b);
   AddMod(x, a, B, n);
}


void SubMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   NTL_ZZRegister(B);
   conv(B, b);
   SubMod(x, a, B, n);
}

void SubMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
{
   NTL_ZZRegister(A);
   conv(A, a);
   SubMod(x, A, b, n);
}



// ****** input and output

NTL_THREAD_LOCAL static long iodigits = 0;
NTL_THREAD_LOCAL static long ioradix = 0;

// iodigits is the greatest integer such that 10^{iodigits} < NTL_WSP_BOUND
// ioradix = 10^{iodigits}

static void InitZZIO()
{
   long x;

   x = (NTL_WSP_BOUND-1)/10;
   iodigits = 0;
   ioradix = 1;

   while (x) {
      x = x / 10;
      iodigits++;
      ioradix = ioradix * 10;
   }

   if (iodigits <= 0) TerminalError("problem with I/O");
}


istream& operator>>(istream& s, ZZ& x)
{
   long c;
   long cval;
   long sign;
   long ndigits;
   long acc;
   NTL_ZZRegister(a);

   if (!s) NTL_INPUT_ERROR(s, "bad ZZ input");

   if (!iodigits) InitZZIO();

   a = 0;

   SkipWhiteSpace(s);
   c = s.peek();

   if (c == '-') {
      sign = -1;
      s.get();
      c = s.peek();
   }
   else
      sign = 1;

   cval = CharToIntVal(c);

   if (cval < 0 || cval > 9) NTL_INPUT_ERROR(s, "bad ZZ input");

   ndigits = 0;
   acc = 0;
   while (cval >= 0 && cval <= 9) {
      acc = acc*10 + cval;
      ndigits++;

      if (ndigits == iodigits) {
         mul(a, a, ioradix);
         add(a, a, acc);
         ndigits = 0;
         acc = 0;
      }

      s.get();
      c = s.peek();
      cval = CharToIntVal(c);
   }

   if (ndigits != 0) {
      long mpy = 1;
      while (ndigits > 0) {
         mpy = mpy * 10;
         ndigits--;
      }

      mul(a, a, mpy);
      add(a, a, acc);
   }

   if (sign == -1)
      negate(a, a);

   x = a;
   return s;
}


// The class _ZZ_local_stack should be defined in an empty namespace,
// but since I don't want to rely on namespaces, we just give it a funny 
// name to avoid accidental name clashes.

struct _ZZ_local_stack {
   long top;
   Vec<long> data;

   _ZZ_local_stack() { top = -1; }

   long pop() { return data[top--]; }
   long empty() { return (top == -1); }
   void push(long x);
};

void _ZZ_local_stack::push(long x)
{
   if (top+1 >= data.length()) 
      data.SetLength(max(32, long(1.414*data.length())));

   top++;
   data[top] = x;
}


static
void PrintDigits(ostream& s, long d, long justify)
{
   NTL_THREAD_LOCAL static Vec<char> buf(INIT_SIZE, iodigits);

   long i = 0;

   while (d) {
      buf[i] = IntValToChar(d % 10);
      d = d / 10;
      i++;
   }

   if (justify) {
      long j = iodigits - i;
      while (j > 0) {
         s << "0";
         j--;
      }
   }

   while (i > 0) {
      i--;
      s << buf[i];
   }
}
      

   

ostream& operator<<(ostream& s, const ZZ& a)
{
   ZZ b;
   _ZZ_local_stack S;
   long r;
   long k;

   if (!iodigits) InitZZIO();

   b = a;

   k = sign(b);

   if (k == 0) {
      s << "0";
      return s;
   }

   if (k < 0) {
      s << "-";
      negate(b, b);
   }

   do {
      r = DivRem(b, b, ioradix);
      S.push(r);
   } while (!IsZero(b));

   r = S.pop();
   PrintDigits(s, r, 0);

   while (!S.empty()) {
      r = S.pop();
      PrintDigits(s, r, 1);
   }
      
   return s;
}



long GCD(long a, long b)
{
   long u, v, t, x;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) ResourceError("GCD: integer overflow");
      a = -a;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) ResourceError("GCD: integer overflow");
      b = -b;
   }


   if (b==0)
      x = a;
   else {
      u = a;
      v = b;
      do {
         t = u % v;
         u = v; 
         v = t;
      } while (v != 0);

      x = u;
   }

   return x;
}

         

void XGCD(long& d, long& s, long& t, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) ResourceError("XGCD: integer overflow");
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) ResourceError("XGCD: integer overflow");
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d = u;
   s = u1;
   t = v1;
}
   

long InvMod(long a, long n)
{
   long d, s, t;

   XGCD(d, s, t, a, n);
   if (d != 1) {
      InvModError("InvMod: inverse undefined");
   }
   if (s < 0)
      return s + n;
   else
      return s;
}


long PowerMod(long a, long ee, long n)
{
   long x, y;

   unsigned long e;

   if (ee < 0)
      e = - ((unsigned long) ee);
   else
      e = ee;

   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = MulMod(x, y, n);
      y = MulMod(y, y, n);
      e = e >> 1;
   }

   if (ee < 0) x = InvMod(x, n);

   return x;
}

long ProbPrime(long n, long NumTests)
{
   long m, x, y, z;
   long i, j, k;

   if (n <= 1) return 0;


   if (n == 2) return 1;
   if (n % 2 == 0) return 0;

   if (n == 3) return 1;
   if (n % 3 == 0) return 0;

   if (n == 5) return 1;
   if (n % 5 == 0) return 0;

   if (n == 7) return 1;
   if (n % 7 == 0) return 0;

   if (n >= NTL_SP_BOUND) {
      return ProbPrime(to_ZZ(n), NumTests);
   }

   m = n - 1;
   k = 0;
   while((m & 1) == 0) {
      m = m >> 1;
      k++;
   }

   // n - 1 == 2^k * m, m odd

   for (i = 0; i < NumTests; i++) {
      do {
         x = RandomBnd(n);
      } while (x == 0);
      // x == 0 is not a useful candidtae for a witness!


      if (x == 0) continue;
      z = PowerMod(x, m, n);
      if (z == 1) continue;
   
      j = 0;
      do {
         y = z;
         z = MulMod(y, y, n);
         j++;
      } while (j != k && z != 1);

      if (z != 1 || y !=  n-1) return 0;
   }

   return 1;
}


long MillerWitness(const ZZ& n, const ZZ& x)
{
   ZZ m, y, z;
   long j, k;

   if (x == 0) return 0;

   add(m, n, -1);
   k = MakeOdd(m);
   // n - 1 == 2^k * m, m odd

   PowerMod(z, x, m, n);
   if (z == 1) return 0;

   j = 0;
   do {
      y = z;
      SqrMod(z, y, n);
      j++;
   } while (j != k && z != 1);

   if (z != 1) return 1;
   add(y, y, 1);
   if (y != n) return 1;
   return 0;
}


// ComputePrimeBound computes a reasonable bound for trial
// division in the Miller-Rabin test.
// It is computed a bit on the "low" side, since being a bit
// low doesn't hurt much, but being too high can hurt a lot.

static
long ComputePrimeBound(long bn)
{
   long wn = (bn+NTL_ZZ_NBITS-1)/NTL_ZZ_NBITS;

   long fn;

   if (wn <= 36)
      fn = wn/4 + 1;
   else
      fn = long(1.67*sqrt(double(wn)));

   long prime_bnd;

   if (NumBits(bn) + NumBits(fn) > NTL_SP_NBITS)
      prime_bnd = NTL_SP_BOUND;
   else
      prime_bnd = bn*fn;

   return prime_bnd;
}


long ProbPrime(const ZZ& n, long NumTrials)
{
   if (n <= 1) return 0;

   if (n.SinglePrecision()) {
      return ProbPrime(to_long(n), NumTrials);
   }


   long prime_bnd = ComputePrimeBound(NumBits(n));


   PrimeSeq s;
   long p;

   p = s.next();
   while (p && p < prime_bnd) {
      if (rem(n, p) == 0)
         return 0;

      p = s.next();
   }

   ZZ W;
   W = 2;

   // first try W == 2....the exponentiation
   // algorithm runs slightly faster in this case

   if (MillerWitness(n, W))
      return 0;


   long i;

   for (i = 0; i < NumTrials; i++) {
      do {
         RandomBnd(W, n);
      } while (W == 0);
      // W == 0 is not a useful candidate for a witness!

      if (MillerWitness(n, W)) 
         return 0;
   }

   return 1;
}


void RandomPrime(ZZ& n, long l, long NumTrials)
{
   if (l <= 1)
      LogicError("RandomPrime: l out of range");

   if (l == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }

   do {
      RandomLen(n, l);
      if (!IsOdd(n)) add(n, n, 1);
   } while (!ProbPrime(n, NumTrials));
}

void NextPrime(ZZ& n, const ZZ& m, long NumTrials)
{
   ZZ x;

   if (m <= 2) {
      n = 2;
      return;
   }

   x = m;

   while (!ProbPrime(x, NumTrials))
      add(x, x, 1);

   n = x;
}

long NextPrime(long m, long NumTrials)
{
   long x;

   if (m <= 2) 
      return 2;

   x = m;

   while (x < NTL_SP_BOUND && !ProbPrime(x, NumTrials))
      x++;

   if (x >= NTL_SP_BOUND)
      ResourceError("NextPrime: no more primes");

   return x;
}



long NextPowerOfTwo(long m)
{
   long k; 
   unsigned long n, um;

   if (m < 0) return 0;

   um = m;
   n = 1;
   k = 0;

   while (n < um) {
      n = n << 1;
      k++;
   }

   if (k >= NTL_BITS_PER_LONG-1)
      ResourceError("NextPowerOfTwo: overflow");

   return k;
}



long NumBits(long a)
{
   unsigned long aa;
   if (a < 0) 
      aa = - ((unsigned long) a);
   else
      aa = a;

   long k = 0;
   while (aa) {
      k++;
      aa = aa >> 1;
   }

   return k;
}


long bit(long a, long k)
{
   unsigned long aa;
   if (a < 0)
      aa = - ((unsigned long) a);
   else
      aa = a;

   if (k < 0 || k >= NTL_BITS_PER_LONG) 
      return 0;
   else
      return long((aa >> k) & 1);
}



long divide(ZZ& q, const ZZ& a, const ZZ& b)
{
   NTL_ZZRegister(qq);
   NTL_ZZRegister(r);

   if (IsZero(b)) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }


   if (IsOne(b)) {
      q = a;
      return 1;
   }

   DivRem(qq, r, a, b);
   if (!IsZero(r)) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, const ZZ& b)
{
   NTL_ZZRegister(r);

   if (IsZero(b)) return IsZero(a);
   if (IsOne(b)) return 1;

   rem(r, a, b);
   return IsZero(r);
}

long divide(ZZ& q, const ZZ& a, long b)
{
   NTL_ZZRegister(qq);

   if (!b) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }

   if (b == 1) {
      q = a;
      return 1;
   }

   long r = DivRem(qq, a, b);
   if (r) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, long b)
{
   if (!b) return IsZero(a);
   if (b == 1) {
      return 1;
   }

   long r = rem(a,  b);
   return (r == 0);
}


void InvMod(ZZ& x, const ZZ& a, const ZZ& n)
{
   // NOTE: the underlying LIP routines write to the first argument,
   // even if inverse is undefined

   NTL_ZZRegister(xx);
   if (InvModStatus(xx, a, n)) 
      InvModError("InvMod: inverse undefined", a, n);
   x = xx;
}

void PowerMod(ZZ& x, const ZZ& a, const ZZ& e, const ZZ& n)
{
   // NOTE: this ensures that all modular inverses are computed
   // in the routine InvMod above, rather than the LIP-internal
   // modular inverse routine
   if (e < 0) {
      ZZ a_inv;
      ZZ e_neg;

      InvMod(a_inv, a, n);
      negate(e_neg, e);
      LowLevelPowerMod(x, a_inv, e_neg, n);
   }
   else
      LowLevelPowerMod(x, a, e, n); 
}
   
#ifdef NTL_EXCEPTIONS

void InvModError(const char *s, const ZZ& a, const ZZ& n)
{
   throw InvModErrorObject(s, a, n); 
}

#else

void InvModError(const char *s, const ZZ& a, const ZZ& n)
{
   TerminalError(s);
}


#endif

long RandomPrime_long(long l, long NumTrials)
{
   if (l <= 1 || l >= NTL_BITS_PER_LONG)
      ResourceError("RandomPrime: length out of range");

   long n;
   do {
      n = RandomLen_long(l);
   } while (!ProbPrime(n, NumTrials));

   return n;
}


static Lazy< Vec<char> > lowsieve_storage;
// This is a GLOBAL VARIABLE


PrimeSeq::PrimeSeq()
{
   movesieve = 0;
   pshift = -1;
   pindex = -1;
   exhausted = 0;
}


long PrimeSeq::next()
{
   if (exhausted) {
      return 0;
   }

   if (pshift < 0) {
      shift(0);
      return 2;
   }

   for (;;) {
      const char *p = movesieve;
      long i = pindex;

      while ((++i) < NTL_PRIME_BND) {
         if (p[i]) {
            pindex = i;
            return pshift + 2 * i + 3;
         }
      }

      long newshift = pshift + 2*NTL_PRIME_BND;

      if (newshift > 2 * NTL_PRIME_BND * (2 * NTL_PRIME_BND + 1)) {
         /* end of the road */
         exhausted = 1;
         return 0;
      }

      shift(newshift);
   }
}

void PrimeSeq::shift(long newshift)
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibound;
   char *p;

   if (!lowsieve_storage.built())
      start();

   const char *lowsieve = lowsieve_storage->elts();


   if (newshift < 0) {
      pshift = -1;
   }
   else if (newshift == 0) {
      pshift = 0;
      movesieve = lowsieve;
   } 
   else if (newshift != pshift) {
      if (movesieve_mem.length() == 0) {
         movesieve_mem.SetLength(NTL_PRIME_BND);
      }

      pshift = newshift;
      movesieve = p = movesieve_mem.elts();
      for (i = 0; i < NTL_PRIME_BND; i++)
         p[i] = 1;

      jstep = 3;
      ibound = pshift + 2 * NTL_PRIME_BND + 1;
      for (i = 0; jstep * jstep <= ibound; i++) {
         if (lowsieve[i]) {
            if (!((jstart = (pshift + 2) / jstep + 1) & 1))
               jstart++;
            if (jstart <= jstep)
               jstart = jstep;
            jstart = (jstart * jstep - pshift - 3) / 2;
            for (j = jstart; j < NTL_PRIME_BND; j += jstep)
               p[j] = 0;
         }
         jstep += 2;
      }
   }

   pindex = -1;
   exhausted = 0;
}


void PrimeSeq::start()
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibnd;
   char *p;

   do {
      Lazy< Vec<char> >::Builder builder(lowsieve_storage);
      if (!builder()) break;

      UniquePtr< Vec<char> > ptr;
      ptr.make();
      ptr->SetLength(NTL_PRIME_BND);

      p = ptr->elts();

      for (i = 0; i < NTL_PRIME_BND; i++)
         p[i] = 1;
         
      jstep = 1;
      jstart = -1;
      ibnd = (SqrRoot(2 * NTL_PRIME_BND + 1) - 3) / 2;
      for (i = 0; i <= ibnd; i++) {
         jstart += 2 * ((jstep += 2) - 1);
         if (p[i])
            for (j = jstart; j < NTL_PRIME_BND; j += jstep)
               p[j] = 0;
      }

      builder.move(ptr);
   } while (0);

}

void PrimeSeq::reset(long b)
{
   if (b > (2*NTL_PRIME_BND+1)*(2*NTL_PRIME_BND+1)) {
      exhausted = 1;
      return;
   }

   if (b <= 2) {
      shift(-1);
      return;
   }

   if ((b & 1) == 0) b++;

   shift(((b-3) / (2*NTL_PRIME_BND))* (2*NTL_PRIME_BND));
   pindex = (b - pshift - 3)/2 - 1;
}
 
long Jacobi(const ZZ& aa, const ZZ& nn)
{
   ZZ a, n;
   long t, k;
   long d;

   a = aa;
   n = nn;
   t = 1;

   while (a != 0) {
      k = MakeOdd(a);
      d = trunc_long(n, 3);
      if ((k & 1) && (d == 3 || d == 5)) t = -t;

      if (trunc_long(a, 2) == 3 && (d & 3) == 3) t = -t;
      swap(a, n);
      rem(a, a, n);
   }

   if (n == 1)
      return t;
   else
      return 0;
}


void SqrRootMod(ZZ& x, const ZZ& aa, const ZZ& nn)
{
   if (aa == 0 || aa == 1) {
      x = aa;
      return;
   }

   // at this point, we must have nn >= 5

   if (trunc_long(nn, 2) == 3) {  // special case, n = 3 (mod 4)
      ZZ n, a, e, z;

      n = nn;
      a  = aa;

      add(e, n, 1);
      RightShift(e, e, 2);

      PowerMod(z, a, e, n);
      x = z;

      return;
   }

   ZZ n, m;
   int h, nlen;

   n = nn;
   nlen = NumBits(n);

   sub(m, n, 1);
   h = MakeOdd(m);  // h >= 2


   if (nlen > 50 && h < SqrRoot(nlen)) {
      long i, j;
      ZZ a, b, a_inv, c, r, m1, d;

      a = aa;
      InvMod(a_inv, a, n);

      if (h == 2) 
         b = 2;
      else {
         do {
            RandomBnd(b, n);
         } while (Jacobi(b, n) != -1);
      }


      PowerMod(c, b, m, n);
      
      add(m1, m, 1);
      RightShift(m1, m1, 1);
      PowerMod(r, a, m1, n);

      for (i = h-2; i >= 0; i--) {
         SqrMod(d, r, n);
         MulMod(d, d, a_inv, n);
         for (j = 0; j < i; j++)
            SqrMod(d, d, n);
         if (!IsOne(d))
            MulMod(r, r, c, n);
         SqrMod(c, c, n);
      } 

      x = r;
      return;
   } 





   long i, k;
   ZZ ma, t, u, v, e;
   ZZ t1, t2, t3, t4;

   n = nn;
   NegateMod(ma, aa, n);

   // find t such that t^2 - 4*a is not a square

   MulMod(t1, ma, 4, n);
   do {
      RandomBnd(t, n);
      SqrMod(t2, t, n);
      AddMod(t2, t2, t1, n);
   } while (Jacobi(t2, n) != -1);

   // compute u*X + v = X^{(n+1)/2} mod f, where f = X^2 - t*X + a

   add(e, n, 1);
   RightShift(e, e, 1);

   u = 0;
   v = 1;

   k = NumBits(e);

   for (i = k - 1; i >= 0; i--) {
      add(t2, u, v);
      sqr(t3, t2);  // t3 = (u+v)^2
      sqr(t1, u);
      sqr(t2, v);
      sub(t3, t3, t1);
      sub(t3, t3, t2); // t1 = u^2, t2 = v^2, t3 = 2*u*v
      rem(t1, t1, n);
      mul(t4, t1, t);
      add(t4, t4, t3);
      rem(u, t4, n);

      mul(t4, t1, ma);
      add(t4, t4, t2);
      rem(v, t4, n);
      
      if (bit(e, i)) {
         MulMod(t1, u, t, n);
         AddMod(t1, t1, v, n);
         MulMod(v, u, ma, n);
         u = t1;
      }

   }

   x = v;
}



// Chinese Remaindering.
//
// This version in new to v3.7, and is significantly
// simpler and faster than the previous version.
//
// This function takes as input g, a, G, p,
// such that a > 0, 0 <= G < p, and gcd(a, p) = 1.
// It computes a' = a*p and g' such that 
//   * g' = g (mod a);
//   * g' = G (mod p);
//   * -a'/2 < g' <= a'/2.
// It then sets g := g' and a := a', and returns 1 iff g has changed.
//
// Under normal use, the input value g satisfies -a/2 < g <= a/2;
// however, this was not documented or enforced in earlier versions,
// so to maintain backward compatability, no restrictions are placed
// on g.  This routine runs faster, though, if -a/2 < g <= a/2,
// and the first thing the routine does is to make this condition
// hold.
//
// Also, under normal use, both a and p are odd;  however, the routine
// will still work even if this is not so.
//
// The routine is based on the following simple fact.
//
// Let -a/2 < g <= a/2, and let h satisfy
//   * g + a h = G (mod p);
//   * -p/2 < h <= p/2.
// Further, if p = 2*h and g > 0, set
//   g' := g - a h;
// otherwise, set
//   g' := g + a h.
// Then g' so defined satisfies the above requirements.
//
// It is trivial to see that g's satisfies the congruence conditions.
// The only thing is to check that the "balancing" condition
// -a'/2 < g' <= a'/2 also holds.


long CRT(ZZ& gg, ZZ& a, long G, long p)
{
   if (p >= NTL_SP_BOUND) {
      ZZ GG, pp;
      conv(GG, G);
      conv(pp, p);
      return CRT(gg, a, GG, pp);
   }

   long modified = 0;

   NTL_ZZRegister(g);

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   long p1;
   p1 = p >> 1;

   long a_inv;
   a_inv = rem(a, p);
   a_inv = InvMod(a_inv, p);

   long h;
   h = rem(g, p);
   h = SubMod(G, h, p);
   h = MulMod(h, a_inv, p);
   if (h > p1)
      h = h - p;

   if (h != 0) {
      modified = 1;

      if (!(p & 1) && g > 0 && (h == p1))
         MulSubFrom(g, a, h);
      else
         MulAddTo(g, a, h);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}

long CRT(ZZ& gg, ZZ& a, const ZZ& G, const ZZ& p)
{
   long modified = 0;

   ZZ g;

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   ZZ p1;
   RightShift(p1, p, 1);

   ZZ a_inv;
   rem(a_inv, a, p);
   InvMod(a_inv, a_inv, p);

   ZZ h;
   rem(h, g, p);
   SubMod(h, G, h, p);
   MulMod(h, h, a_inv, p);
   if (h > p1)
      sub(h, h, p);

   if (h != 0) {
      modified = 1;
      ZZ ah;
      mul(ah, a, h);

      if (!IsOdd(p) && g > 0 &&  (h == p1))
         sub(g, g, ah);
      else
         add(g, g, ah);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}



void sub(ZZ& x, const ZZ& a, long b)
{
   NTL_ZZRegister(B);
   conv(B, b);
   sub(x, a, B);
}

void sub(ZZ& x, long a, const ZZ& b)
{
   NTL_ZZRegister(A);
   conv(A, a);
   sub(x, A, b);
}


void power2(ZZ& x, long e)
{
   if (e < 0) ArithmeticError("power2: negative exponent");
   set(x);
   LeftShift(x, x, e);
}

   
void conv(ZZ& x, const char *s)
{
   long c;
   long cval;
   long sign;
   long ndigits;
   long acc;
   long i = 0;

   NTL_ZZRegister(a);

   if (!s) InputError("bad ZZ input");

   if (!iodigits) InitZZIO();

   a = 0;

   c = s[i];
   while (IsWhiteSpace(c)) {
      i++;
      c = s[i];
   }

   if (c == '-') {
      sign = -1;
      i++;
      c = s[i];
   }
   else
      sign = 1;

   cval = CharToIntVal(c);
   if (cval < 0 || cval > 9) InputError("bad ZZ input");

   ndigits = 0;
   acc = 0;
   while (cval >= 0 && cval <= 9) {
      acc = acc*10 + cval;
      ndigits++;

      if (ndigits == iodigits) {
         mul(a, a, ioradix);
         add(a, a, acc);
         ndigits = 0;
         acc = 0;
      }

      i++;
      c = s[i];
      cval = CharToIntVal(c);
   }

   if (ndigits != 0) {
      long mpy = 1;
      while (ndigits > 0) {
         mpy = mpy * 10;
         ndigits--;
      }

      mul(a, a, mpy);
      add(a, a, acc);
   }

   if (sign == -1)
      negate(a, a);

   x = a;
}



void bit_and(ZZ& x, const ZZ& a, long b)
{
   NTL_ZZRegister(B);
   conv(B, b);
   bit_and(x, a, B);
}

void bit_or(ZZ& x, const ZZ& a, long b)
{
   NTL_ZZRegister(B);
   conv(B, b);
   bit_or(x, a, B);
}

void bit_xor(ZZ& x, const ZZ& a, long b)
{
   NTL_ZZRegister(B);
   conv(B, b);
   bit_xor(x, a, B);
}


long power_long(long a, long e)
{
   if (e < 0) ArithmeticError("power_long: negative exponent");

   if (e == 0) return 1;

   if (a == 1) return 1;
   if (a == -1) {
      if (e & 1)
         return -1;
      else
         return 1;
   }

   // no overflow check --- result is computed correctly
   // modulo word size

   unsigned long res = 1;
   unsigned long aa = a;
   long i;

   for (i = 0; i < e; i++)
      res *= aa;

   return to_long(res);
}



// ======================= new PRG stuff ======================




#if (NTL_BITS_PER_INT32 == 32)
#define INT32MASK(x) (x)
#else
#define INT32MASK(x) ((x) & _ntl_uint32(0xffffffff))
#endif



// SHA256 code adapted from an implementauin by Brad Conte.
// The following is from his original source files.
/*********************************************************************
* Filename:   sha256.c
* Author:     Brad Conte (brad AT bradconte.com)
* Copyright:
* Disclaimer: This code is presented "as is" without any guarantees.
* Details:    Implementation of the SHA-256 hashing algorithm.
              SHA-256 is one of the three algorithms in the SHA2
              specification. The others, SHA-384 and SHA-512, are not
              offered in this implementation.
              Algorithm specification can be found here:
               * http://csrc.nist.gov/publications/fips/fips180-2/fips180-2withchangenotice.pdf
              This implementation uses little endian byte order.
*********************************************************************/




#define SHA256_BLOCKSIZE (64)
#define SHA256_HASHSIZE  (32)

// DBL_INT_ADD treats two unsigned ints a and b as one 64-bit integer and adds c to it
static inline
void DBL_INT_ADD(_ntl_uint32& a, _ntl_uint32& b, _ntl_uint32 c)
{
   _ntl_uint32 aa = INT32MASK(a);
   if (aa > INT32MASK(_ntl_uint32(0xffffffff) - c)) b++;
   a = aa + c;
}

#define ROTLEFT(a,b) (((a) << (b)) | (INT32MASK(a) >> (32-(b))))
#define ROTRIGHT(a,b) ((INT32MASK(a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ (INT32MASK(x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ (INT32MASK(x) >> 10))

struct SHA256_CTX {
   unsigned char data[64];
   _ntl_uint32 datalen;
   _ntl_uint32 bitlen[2];
   _ntl_uint32 state[8];
};

static const _ntl_uint32 sha256_const[64] = {
   0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
   0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
   0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
   0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
   0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
   0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
   0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
   0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};


static
void sha256_transform(SHA256_CTX& ctx, unsigned char *data)
{  
   _ntl_uint32 a,b,c,d,e,f,g,h,i,j,t1,t2,m[64];
      
   for (i=0,j=0; i < 16; ++i, j += 4)
      m[i] = (data[j] << 24) | (data[j+1] << 16) | (data[j+2] << 8) | (data[j+3]);
   for ( ; i < 64; ++i)
      m[i] = SIG1(m[i-2]) + m[i-7] + SIG0(m[i-15]) + m[i-16];

   a = ctx.state[0];
   b = ctx.state[1];
   c = ctx.state[2];
   d = ctx.state[3];
   e = ctx.state[4];
   f = ctx.state[5];
   g = ctx.state[6];
   h = ctx.state[7];
   
   for (i = 0; i < 64; ++i) {
      t1 = h + EP1(e) + CH(e,f,g) + sha256_const[i] + m[i];
      t2 = EP0(a) + MAJ(a,b,c);
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
   }
   
   ctx.state[0] += a;
   ctx.state[1] += b;
   ctx.state[2] += c;
   ctx.state[3] += d;
   ctx.state[4] += e;
   ctx.state[5] += f;
   ctx.state[6] += g;
   ctx.state[7] += h;
}  

static
void sha256_init(SHA256_CTX& ctx)
{  
   ctx.datalen = 0; 
   ctx.bitlen[0] = 0; 
   ctx.bitlen[1] = 0; 
   ctx.state[0] = 0x6a09e667;
   ctx.state[1] = 0xbb67ae85;
   ctx.state[2] = 0x3c6ef372;
   ctx.state[3] = 0xa54ff53a;
   ctx.state[4] = 0x510e527f;
   ctx.state[5] = 0x9b05688c;
   ctx.state[6] = 0x1f83d9ab;
   ctx.state[7] = 0x5be0cd19;
}

static
void sha256_update(SHA256_CTX& ctx, const unsigned char *data, _ntl_uint32 len)
{  
   _ntl_uint32 t,i;
   
   for (i=0; i < len; ++i) { 
      ctx.data[ctx.datalen] = data[i]; 
      ctx.datalen++; 
      if (ctx.datalen == 64) { 
         sha256_transform(ctx,ctx.data);
         DBL_INT_ADD(ctx.bitlen[0],ctx.bitlen[1],512); 
         ctx.datalen = 0; 
      }  
   }  
}  

static
void sha256_final(SHA256_CTX& ctx, unsigned char *hash, 
                  long hlen=SHA256_HASHSIZE)
{  
   _ntl_uint32 i, j; 
   
   i = ctx.datalen; 
   
   // Pad whatever data is left in the buffer. 
   if (ctx.datalen < 56) { 
      ctx.data[i++] = 0x80; 
      while (i < 56) 
         ctx.data[i++] = 0x00; 
   }  
   else { 
      ctx.data[i++] = 0x80; 
      while (i < 64) 
         ctx.data[i++] = 0x00; 
      sha256_transform(ctx,ctx.data);
      memset(ctx.data,0,56); 
   }  
   
   // Append to the padding the total message's length in bits and transform. 
   DBL_INT_ADD(ctx.bitlen[0],ctx.bitlen[1],ctx.datalen * 8);

   ctx.data[63] = ctx.bitlen[0]; 
   ctx.data[62] = ctx.bitlen[0] >> 8; 
   ctx.data[61] = ctx.bitlen[0] >> 16; 
   ctx.data[60] = ctx.bitlen[0] >> 24; 
   ctx.data[59] = ctx.bitlen[1]; 
   ctx.data[58] = ctx.bitlen[1] >> 8; 
   ctx.data[57] = ctx.bitlen[1] >> 16;  
   ctx.data[56] = ctx.bitlen[1] >> 24; 
   sha256_transform(ctx,ctx.data);
   
   for (i = 0; i < 8; i++) {
      _ntl_uint32 w = ctx.state[i];
      for (j = 0; j < 4; j++) {
         if (hlen <= 0) break;
         hash[4*i + j] = w >> (24-j*8); 
         hlen--;
      }
   }

}  



static
void sha256(const unsigned char *data, long dlen, unsigned char *hash, 
            long hlen=SHA256_HASHSIZE)
{
   if (dlen < 0) dlen = 0;
   if (hlen < 0) hlen = 0;

   SHA256_CTX ctx;
   sha256_init(ctx);

   const long BLKSIZE = 4096;

   long i;
   for (i = 0; i <= dlen-BLKSIZE; i += BLKSIZE) 
      sha256_update(ctx, data + i, BLKSIZE);

   if (i < dlen)
      sha256_update(ctx, data + i, dlen - i);

   sha256_final(ctx, hash, hlen);
}


static
void hmac_sha256(const unsigned char *key, long klen, 
                 const unsigned char *data, long dlen,
                 unsigned char *hash, long hlen=SHA256_HASHSIZE)
{
   if (klen < 0) klen = 0;
   if (dlen < 0) dlen = 0;
   if (hlen < 0) hlen = 0;

   unsigned char K[SHA256_BLOCKSIZE];
   unsigned char tmp[SHA256_HASHSIZE];

   long i;

   if (klen <= SHA256_BLOCKSIZE) {
      for (i = 0; i < klen; i++)
         K[i] = key[i];
      for (i = klen; i < SHA256_BLOCKSIZE; i++) 
         K[i] = 0;
   }
   else {
      sha256(key, klen, K, SHA256_BLOCKSIZE); 
      for (i = SHA256_HASHSIZE; i < SHA256_BLOCKSIZE; i++)
         K[i] = 0;
   }

   for (i = 0; i < SHA256_BLOCKSIZE; i++)
      K[i] ^= 0x36;

   SHA256_CTX ctx;
   sha256_init(ctx);
   sha256_update(ctx, K, SHA256_BLOCKSIZE);
   sha256_update(ctx, data, dlen);
   sha256_final(ctx, tmp);

   for (i = 0; i < SHA256_BLOCKSIZE; i++)
      K[i] ^= (0x36 ^ 0x5C);

   sha256_init(ctx);
   sha256_update(ctx, K, SHA256_BLOCKSIZE);
   sha256_update(ctx, tmp, SHA256_HASHSIZE);
   sha256_final(ctx, hash, hlen);
}


// This key derivation uses HMAC with a zero key to derive
// an intermediate key K from the data, and then uses HMAC
// as a PRF in counter mode with key K to derive the final key

void DeriveKey(unsigned char *hash, long hlen,  
               const unsigned char *data, long dlen)
{
   if (dlen < 0) LogicError("DeriveKey: bad args");
   if (hlen < 0) LogicError("DeriveKey: bad args");

   long i, j;


   unsigned char K[SHA256_HASHSIZE];
   hmac_sha256(0, 0, data, dlen, K); 

   // initialize 64-bit counter to zero
   unsigned char counter[8];
   for (j = 0; j < 8; j++) counter[j] = 0;

   for (i = 0; i <= hlen-SHA256_HASHSIZE; i += SHA256_HASHSIZE) {
      hmac_sha256(K, SHA256_HASHSIZE, counter, 8, hash+i); 

      // increment counter
      for (j = 0; j < 8; j++) {
         counter[j]++;
         if (counter[j] != 0) break; 
      }
   }

   if (i < hlen) 
      hmac_sha256(K, SHA256_HASHSIZE, counter, 8, hash+i, hlen-i);
}




// ******************** ChaCha20 stuff ***********************

static const _ntl_uint32 chacha_const[4] = 
   { 0x61707865, 0x3320646e, 0x79622d32, 0x6b206574 };


#define LE(p) (((_ntl_uint32)((p)[0])) + ((_ntl_uint32)((p)[1]) << 8) + \
    ((_ntl_uint32)((p)[2]) << 16) + ((_ntl_uint32)((p)[3]) << 24))

#define FROMLE(p, x) (p)[0] = (x), (p)[1] = ((x) >> 8), \
   (p)[2] = ((x) >> 16), (p)[3] = ((x) >> 24)


#define QUARTERROUND(x, a, b, c, d) \
    x[a] += x[b], x[d] = ROTLEFT(x[d] ^ x[a], 16), \
    x[c] += x[d], x[b] = ROTLEFT(x[b] ^ x[c], 12), \
    x[a] += x[b], x[d] = ROTLEFT(x[d] ^ x[a], 8), \
    x[c] += x[d], x[b] = ROTLEFT(x[b] ^ x[c], 7)


static
void salsa20_core(_ntl_uint32* data)
{
   long i;

   for (i = 0; i < 10; i++) {
      QUARTERROUND(data, 0, 4, 8, 12);
      QUARTERROUND(data, 1, 5, 9, 13);
      QUARTERROUND(data, 2, 6, 10, 14);
      QUARTERROUND(data, 3, 7, 11, 15);
      QUARTERROUND(data, 0, 5, 10, 15);
      QUARTERROUND(data, 1, 6, 11, 12);
      QUARTERROUND(data, 2, 7, 8, 13);
      QUARTERROUND(data, 3, 4, 9, 14);
   }
}


// key K must be exactly 32 bytes
static
void salsa20_init(_ntl_uint32 *state, const unsigned char *K)  
{
   long i;

   for (i = 0; i < 4; i++)
      state[i] = chacha_const[i];

   for (i = 4; i < 12; i++)
      state[i] = LE(K + 4*(i-4));

   for (i = 12; i < 16; i++)
      state[i] = 0;
}



// state and data are of length 16
static
void salsa20_apply(_ntl_uint32 *state, _ntl_uint32 *data)
{
   long i;

   for (i = 0; i < 16; i++) data[i] = state[i];

   salsa20_core(data);

   for (i = 0; i < 16; i++) data[i] += state[i];

   for (i = 12; i < 16; i++) {
      state[i]++;
      state[i] = INT32MASK(state[i]);
      if (state[i] != 0) break;
   }
}


// state is 16 words, data is 64 bytes
static
void salsa20_apply(_ntl_uint32 *state, unsigned char *data)
{
   _ntl_uint32 wdata[16];
   salsa20_apply(state, wdata);

   long i;
   for (i = 0; i < 16; i++)
      FROMLE(data + 4*i, wdata[i]);

   // FIXME: could use memcpy for above if everything 
   // is right
}



RandomByteStream::RandomByteStream(const unsigned char *key)
{
   salsa20_init(state, key);
   pos = 64;
}


void RandomByteStream::do_get(unsigned char *res, long n)
{
   if (n < 0) LogicError("RandomByteStream::get: bad args");

   long i, j;

   if (n <= 64-pos) {
      for (i = 0; i < n; i++) res[i] = buf[pos+i];
      pos += n;
      return;
   }

   // read remainder of buffer
   for (i = 0; i < 64-pos; i++) res[i] = buf[pos+i];
   n -= 64-pos;
   res += 64-pos;

   _ntl_uint32 wdata[16];

   // read 64-byte chunks
   for (i = 0; i <= n-64; i += 64) {
      salsa20_apply(state, wdata);
      for (j = 0; j < 16; j++)
         FROMLE(res + i + 4*j, wdata[j]);
   }

   if (i < n) { 
      salsa20_apply(state, wdata);

      for (j = 0; j < 16; j++)
         FROMLE(buf + 4*j, wdata[j]);

      pos = n-i;
      for (j = 0; j < pos; j++)
         res[i+j] = buf[j];
   }
   else {
      pos = 64;
   }
}


NTL_THREAD_LOCAL static UniquePtr<RandomByteStream> CurrentRandomByteStream;


void SetSeed(const RandomByteStream& s)
{
   if (!CurrentRandomByteStream)
      CurrentRandomByteStream.make(s);
   else
      *CurrentRandomByteStream = s;
}


void SetSeed(const unsigned char *data, long dlen)
{
   if (dlen < 0) LogicError("SetSeed: bad args");

   Vec<unsigned char> key;
   key.SetLength(NTL_PRG_KEYLEN);
   DeriveKey(key.elts(), NTL_PRG_KEYLEN, data, dlen);
 
   SetSeed(RandomByteStream(key.elts()));
}

void SetSeed(const ZZ& seed)
{
   long nb = NumBytes(seed);

   Vec<unsigned char> buf;
   buf.SetLength(nb);

   BytesFromZZ(buf.elts(), seed, nb);

   SetSeed(buf.elts(), nb);
}


static
void InitRandomByteStream()
{
   const string& id = UniqueID();
   SetSeed((const unsigned char *) id.c_str(), id.length());
}

static inline
RandomByteStream& LocalGetCurrentRandomByteStream()
{
   if (!CurrentRandomByteStream) InitRandomByteStream();
   return *CurrentRandomByteStream;
}

RandomByteStream& GetCurrentRandomByteStream()
{
   return LocalGetCurrentRandomByteStream();
}





static 
void ran_bytes(unsigned char *bytes, long n)
{
   RandomByteStream& stream = LocalGetCurrentRandomByteStream();
   stream.get(bytes, n);
}


unsigned long RandomWord()
{
   unsigned char buf[NTL_BITS_PER_LONG/8];
   long i;
   unsigned long res;

   ran_bytes(buf, NTL_BITS_PER_LONG/8);

   res = 0;
   for (i = NTL_BITS_PER_LONG/8 - 1; i >= 0; i--) {
      res = res << 8;
      res = res | buf[i];
   }

   return res;
}

long RandomBits_long(long l)
{
   if (l <= 0) return 0;
   if (l >= NTL_BITS_PER_LONG) 
      ResourceError("RandomBits: length too big");

   unsigned char buf[NTL_BITS_PER_LONG/8];
   unsigned long res;
   long i;

   long nb = (l+7)/8;
   ran_bytes(buf, nb);

   res = 0;
   for (i = nb - 1; i >= 0; i--) {
      res = res << 8;
      res = res | buf[i];
   }

   return long(res & ((1UL << l)-1UL)); 
}

unsigned long RandomBits_ulong(long l)
{
   if (l <= 0) return 0;
   if (l > NTL_BITS_PER_LONG) 
      ResourceError("RandomBits: length too big");

   unsigned char buf[NTL_BITS_PER_LONG/8];
   unsigned long res;
   long i;

   long nb = (l+7)/8;
   ran_bytes(buf, nb);

   res = 0;
   for (i = nb - 1; i >= 0; i--) {
      res = res << 8;
      res = res | buf[i];
   }

   if (l < NTL_BITS_PER_LONG)
      res = res & ((1UL << l)-1UL);

   return res;
}

long RandomLen_long(long l)
{
   if (l <= 0) return 0;
   if (l == 1) return 1;
   if (l >= NTL_BITS_PER_LONG) 
      ResourceError("RandomLen: length too big");

   return RandomBits_long(l-1) + (1L << (l-1)); 
}


void RandomBits(ZZ& x, long l)
{
   if (l <= 0) {
      x = 0;
      return;
   }

   if (NTL_OVERFLOW(l, 1, 0))
      ResourceError("RandomBits: length too big");

   long nb = (l+7)/8;

   NTL_THREAD_LOCAL static Vec<unsigned char> buf_mem;
   Vec<unsigned char>::Watcher watch_buf_mem(buf_mem);

   buf_mem.SetLength(nb);
   unsigned char *buf = buf_mem.elts();

   ran_bytes(buf, nb);

   NTL_ZZRegister(res);

   ZZFromBytes(res, buf, nb);
   trunc(res, res, l);

   x = res;
}


void RandomLen(ZZ& x, long l)
{
   if (l <= 0) {
      x = 0;
      return;
   }

   if (l == 1) {
      x = 1;
      return;
   }

   if (NTL_OVERFLOW(l, 1, 0))
      ResourceError("RandomLen: length too big");

   // pre-allocate space to avoid two allocations
   long nw = (l + NTL_ZZ_NBITS - 1)/NTL_ZZ_NBITS;
   x.SetSize(nw);

   RandomBits(x, l-1);
   SetBit(x, l-1);
}



void RandomBnd(ZZ& x, const ZZ& bnd)
{
   if (bnd <= 1) {
      x = 0;
      return;
   }

   NTL_ZZRegister(tmp);
   sub(tmp, bnd, 1);
   
   long k = NumBits(tmp);

   do {
      RandomBits(tmp, k);
   } while (tmp >= bnd);

   x = tmp;
}

long RandomBnd(long bnd)
{
   if (bnd <= 1) return 0;

   long k = NumBits(bnd-1);

   long tmp;

   do {
      tmp = RandomBits_long(k);
   } while (tmp >= bnd);

   return tmp;
}




// More prime generation stuff...

static
double Log2(double x)
{
   NTL_THREAD_LOCAL static double log2 = log(2.0);
   return log(x)/log2;
}

// Define p(k,t) to be the conditional probability that a random, odd, k-bit 
// number is composite, given that it passes t iterations of the 
// Miller-Rabin test.
// This routine returns 0 or 1, and if it returns 1 then
// p(k,t) <= 2^{-n}.
// This basically encodes the estimates of Damgard, Landrock, and Pomerance;
// it uses floating point arithmetic, but is coded in such a way
// that its results should be correct, assuming that the log function
// is computed with reasonable precision.
// 
// It is assumed that k >= 3 and t >= 1; if this does not hold,
// then 0 is returned.

static
long ErrBoundTest(long kk, long tt, long nn)

{
   const double fudge = (1.0 + 1024.0/NTL_FDOUBLE_PRECISION);
   const double log2_3 = Log2(3.0);
   const double log2_7 = Log2(7.0);
   const double log2_20 = Log2(20.0);

   double k = kk;
   double t = tt;
   double n = nn;

   if (k < 3 || t < 1) return 0;
   if (n < 1) return 1;

   // the following test is largely academic
   if (9*t > NTL_FDOUBLE_PRECISION) LogicError("ErrBoundTest: t too big");

   double log2_k = Log2(k);

   if ((n + log2_k)*fudge <= 2*t)
      return 1;

   if ((2*log2_k + 4.0 + n)*fudge <= 2*sqrt(k))
      return 2;

   if ((t == 2 && k >= 88) || (3 <= t && 9*t <= k && k >= 21)) {
      if ((1.5*log2_k + t + 4.0 + n)*fudge <= 0.5*Log2(t) + 2*(sqrt(t*k)))
         return 3;
   }

   if (k <= 9*t && 4*t <= k && k >= 21) {
      if ( ((log2_3 + log2_7 + log2_k + n)*fudge <= log2_20 + 5*t)  &&
           ((log2_3 + (15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t) &&
           ((2*log2_3 + 2 + log2_k + n)*fudge <= k/4 + 3*t) )
         return 4; 
   }

   if (4*t >= k && k >= 21) {
      if (((15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t)
         return 5;
   }

   return 0;
}


void GenPrime(ZZ& n, long k, long err)
{
   if (k <= 1) LogicError("GenPrime: bad length");

   if (k > (1L << 20)) ResourceError("GenPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }


   long t;

   t = 1;
   while (!ErrBoundTest(k, t, err))
      t++;

   RandomPrime(n, k, t);
}


long GenPrime_long(long k, long err)
{
   if (k <= 1) LogicError("GenPrime: bad length");

   if (k >= NTL_BITS_PER_LONG) ResourceError("GenPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         return 3;
      else
         return 2;
   }

   long t;

   t = 1;
   while (!ErrBoundTest(k, t, err))
      t++;

   return RandomPrime_long(k, t);
}


void GenGermainPrime(ZZ& n, long k, long err)
{
   if (k <= 1) LogicError("GenGermainPrime: bad length");

   if (k > (1L << 20)) ResourceError("GenGermainPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }


   long prime_bnd = ComputePrimeBound(k);

   if (NumBits(prime_bnd) >= k/2)
      prime_bnd = (1L << (k/2-1));


   ZZ two;
   two = 2;

   ZZ n1;

   
   PrimeSeq s;

   ZZ iter;
   iter = 0;


   for (;;) {
      iter++;

      RandomLen(n, k);
      if (!IsOdd(n)) add(n, n, 1);

      s.reset(3);
      long p;

      long sieve_passed = 1;

      p = s.next();
      while (p && p < prime_bnd) {
         long r = rem(n, p);

         if (r == 0) {
            sieve_passed = 0;
            break;
         }

         // test if 2*r + 1 = 0 (mod p)
         if (r == p-r-1) {
            sieve_passed = 0;
            break;
         }

         p = s.next();
      }

      if (!sieve_passed) continue;


      if (MillerWitness(n, two)) continue;

      // n1 = 2*n+1
      mul(n1, n, 2);
      add(n1, n1, 1);


      if (MillerWitness(n1, two)) continue;

      // now do t M-R iterations...just to make sure
 
      // First compute the appropriate number of M-R iterations, t
      // The following computes t such that 
      //       p(k,t)*8/k <= 2^{-err}/(5*iter^{1.25})
      // which suffices to get an overall error probability of 2^{-err}.
      // Note that this method has the advantage of not requiring 
      // any assumptions on the density of Germain primes.

      long err1 = max(1, err + 7 + (5*NumBits(iter) + 3)/4 - NumBits(k));
      long t;
      t = 1;
      while (!ErrBoundTest(k, t, err1))
         t++;

      ZZ W;
      long MR_passed = 1;

      long i;
      for (i = 1; i <= t; i++) {
         do {
            RandomBnd(W, n);
         } while (W == 0);
         // W == 0 is not a useful candidate witness!

         if (MillerWitness(n, W)) {
            MR_passed = 0;
            break;
         }
      }

      if (MR_passed) break;
   }
}

long GenGermainPrime_long(long k, long err)
{
   if (k >= NTL_BITS_PER_LONG-1)
      ResourceError("GenGermainPrime_long: length too long");

   ZZ n;
   GenGermainPrime(n, k, err);
   return to_long(n);
}


NTL_END_IMPL
