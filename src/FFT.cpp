#include "FFT.h"


static inline 
Vec<long> *get_brc_mem()
{
   NTL_TLS_LOCAL_INIT(Vec< Vec<long> >, brc_mem_vec, (INIT_SIZE, NTL_FFTMaxRoot+1));
   return brc_mem_vec.elts();
}


static
long RevInc(long a, long k)
{
   long j, m;

   j = k; 
   m = 1L << (k-1);

   while (j && (m & a)) {
      a ^= m;
      m >>= 1;
      j--;
   }
   if (j) a ^= m;
   return a;
}


static
void BitReverseCopy(unsigned long * NTL_RESTRICT A, const long * NTL_RESTRICT a, long k)
{
   Vec<long> *brc_mem = get_brc_mem();

   long n = 1L << k;
   long* NTL_RESTRICT rev;
   long i, j;

   rev = brc_mem[k].elts();
   if (!rev) {
      brc_mem[k].SetLength(n);
      rev = brc_mem[k].elts();
      for (i = 0, j = 0; i < n; i++, j = RevInc(j, k))
         rev[i] = j;
   }

   for (i = 0; i < n; i++)
      A[rev[i]] = a[i];
}




// FFT with  lazy multiplication

#if (defined(NTL_LONGLONG_SP_MULMOD))


#if (NTL_BITS_PER_LONG >= NTL_SP_NBITS+4) 

static inline unsigned long 
sp_NormalizedLazyPrepMulModPreconWithRem(unsigned long& rres, long b, long n, unsigned long ninv)
{
   unsigned long H = cast_unsigned(b);
   unsigned long Q = ll_mul_hi(H << 4, ninv);
   unsigned long L = cast_unsigned(b) << (NTL_SP_NBITS+2);
   long r = L - Q*cast_unsigned(n);  // r in [0..2*n)

   r = sp_CorrectExcessQuo(Q, r, n);
   rres = r;
   return Q; // NOTE: not shifted
}

static inline unsigned long 
sp_NormalizedLazyPrepMulModPrecon(long b, long n, unsigned long ninv)
{
   unsigned long H = cast_unsigned(b);
   unsigned long Q = ll_mul_hi(H << 4, ninv);
   unsigned long L = cast_unsigned(b) << (NTL_SP_NBITS+2);
   long r = L - Q*cast_unsigned(n);  // r in [0..2*n)

   Q += 1L + sp_SignMask(r-n);
   return Q; // NOTE: not shifted
}


#else

// NTL_BITS_PER_LONG == NTL_SP_NBITS+2
static inline unsigned long 
sp_NormalizedLazyPrepMulModPreconWithRem(unsigned long& rres, long b, long n, unsigned long ninv)
{
   unsigned long H = cast_unsigned(b) << 2;
   unsigned long Q = ll_mul_hi(H, (ninv << 1)) + H;
   unsigned long rr = -Q*cast_unsigned(n);  // r in [0..3*n)

   long r = sp_CorrectExcessQuo(Q, rr, n);
   r = sp_CorrectExcessQuo(Q, r, n);
   rres = r;
   return Q;  // NOTE: not shifted
}

static inline unsigned long 
sp_NormalizedLazyPrepMulModPrecon(long b, long n, unsigned long ninv)
{
   unsigned long H = cast_unsigned(b) << 2;
   unsigned long Q = ll_mul_hi(H, (ninv << 1)) + H;
   unsigned long rr = -Q*cast_unsigned(n);  // r in [0..3*n)
   Q += 2L + sp_SignMask(rr-n) + sp_SignMask(rr-2*n);
   return Q; // NOTE: not shifted
}


#endif


static inline unsigned long
LazyPrepMulModPrecon(long b, long n, sp_inverse ninv)
{
   return sp_NormalizedLazyPrepMulModPrecon(b << ninv.shamt, n << ninv.shamt, ninv.inv) << (NTL_BITS_PER_LONG-NTL_SP_NBITS-2);
}


static inline unsigned long
LazyPrepMulModPreconWithRem(unsigned long& rres, long b, long n, sp_inverse ninv)
{
   unsigned long qq, rr;
   qq = sp_NormalizedLazyPrepMulModPreconWithRem(rr, b << ninv.shamt, n << ninv.shamt, ninv.inv); 
   rres = rr >> ninv.shamt;
   return qq << (NTL_BITS_PER_LONG-NTL_SP_NBITS-2);
}








#elif (NTL_BITS_PER_LONG - NTL_SP_NBITS >= 4 && NTL_WIDE_DOUBLE_PRECISION - NTL_SP_NBITS >= 4)


// slightly faster functions, which should kick in on x86-64, where 
//    NTL_BITS_PER_LONG == 64
//    NTL_SP_NBITS == 60 (another reason for holding this back to 60 bits)
//    NTL_WIDE_DOUBLE_PRECISION == 64

// DIRT: if the relative error in floating point calcuations (muls and reciprocals)
//   is <= epsilon, the relative error in the calculations is <= 3*epsilon +
//   O(epsilon^2), and we require that this relative error is at most
//   2^{-(NTL_SP_NBITS+2)}, so it should be pretty safe as long as
//   epsilon is at most, or not much geater than, 2^{-NTL_WIDE_DOUBLE_PRECISION}.

static inline 
unsigned long LazyPrepMulModPrecon(long b, long n, wide_double ninv)
{
   long q = (long) ( (((wide_double) b) * wide_double(4*NTL_SP_BOUND)) * ninv ); 

   unsigned long rr = (cast_unsigned(b) << (NTL_SP_NBITS+2)) 
                       - cast_unsigned(q)*cast_unsigned(n);

   q += sp_SignMask(rr) + sp_SignMask(rr-n) + 1L;

   return cast_unsigned(q) << (NTL_BITS_PER_LONG - NTL_SP_NBITS - 2);
}

static inline 
unsigned long LazyPrepMulModPreconWithRem(unsigned long& rres, long b, long n, wide_double ninv)
{
   long q = (long) ( (((wide_double) b) * wide_double(4*NTL_SP_BOUND)) * ninv ); 

   unsigned long rr = (cast_unsigned(b) << (NTL_SP_NBITS+2)) 
                       - cast_unsigned(q)*cast_unsigned(n);

   long r = sp_CorrectDeficitQuo(q, rr, n);
   r = sp_CorrectExcessQuo(q, r, n);

   unsigned long qres = cast_unsigned(q) << (NTL_BITS_PER_LONG - NTL_SP_NBITS - 2);
   rres = r;
   return qres;
}

#else


static inline 
unsigned long LazyPrepMulModPrecon(long b, long n, wide_double ninv)
{
   long q = (long) ( (((wide_double) b) * wide_double(NTL_SP_BOUND)) * ninv ); 

   unsigned long rr = (cast_unsigned(b) << (NTL_SP_NBITS)) 
                       - cast_unsigned(q)*cast_unsigned(n);

   long r = sp_CorrectDeficitQuo(q, rr, n);
   r = sp_CorrectExcessQuo(q, r, n);

   unsigned long qq = q;

   qq = 2*qq;
   r = 2*r;
   r = sp_CorrectExcessQuo(qq, r, n);

   qq = 2*qq;
   r = 2*r;
   qq += sp_SignMask(r-n) + 1L;

   return qq << (NTL_BITS_PER_LONG - NTL_SP_NBITS - 2);
}





static inline 
unsigned long LazyPrepMulModPreconWithRem(unsigned long& rres, long b, long n, wide_double ninv)
{
   long q = (long) ( (((wide_double) b) * wide_double(NTL_SP_BOUND)) * ninv ); 

   unsigned long rr = (cast_unsigned(b) << (NTL_SP_NBITS)) 
                       - cast_unsigned(q)*cast_unsigned(n);

   long r = sp_CorrectDeficitQuo(q, rr, n);
   r = sp_CorrectExcessQuo(q, r, n);

   unsigned long qq = q;

   qq = 2*qq;
   r = 2*r;
   r = sp_CorrectExcessQuo(qq, r, n);

   qq = 2*qq;
   r = 2*r;
   r = sp_CorrectExcessQuo(qq, r, n);

   rres = r;
   return qq << (NTL_BITS_PER_LONG - NTL_SP_NBITS - 2);
}

#endif



static inline
unsigned long LazyMulModPreconQuo(unsigned long a, unsigned long b, 
                                  unsigned long n, unsigned long bninv)
{
   unsigned long q = ll_mul_hi(a, bninv);
   unsigned long r = a*b - q*n;
   q += sp_SignMask(r-n) + 1L;
   return q << (NTL_BITS_PER_LONG - NTL_SP_NBITS - 2);
}


static inline 
unsigned long LazyMulModPrecon(unsigned long a, unsigned long b, 
                               unsigned long n, unsigned long bninv)
{
   unsigned long q = ll_mul_hi(a, bninv);
   unsigned long res = a*b - q*n;
   return res;
}


static inline 
unsigned long LazyReduce1(unsigned long a, long q)
{
  return sp_CorrectExcess(long(a), q);
}

static inline 
unsigned long LazyReduce2(unsigned long a, long q)
{
  return sp_CorrectExcess(a, 2*q);
}





static
void LazyPrecompFFTMultipliers(long k, long q, mulmod_t qinv, const long *root, const FFTMultipliers& tab)
{
   if (k < 1) LogicError("LazyPrecompFFTMultipliers: bad input");

   do { // NOTE: thread safe lazy init
      FFTMultipliers::Builder bld(tab, k+1);
      long amt = bld.amt();
      if (!amt) break;

      long first = k+1-amt;
      // initialize entries first..k


      for (long s = first; s <= k; s++) {
         UniquePtr<FFTVectorPair> item;

         if (s == 0) {
            bld.move(item); // position 0 not used
            continue;
         }

         if (s == 1) {
            item.make();
            item->wtab_precomp.SetLength(1);
            item->wqinvtab_precomp.SetLength(1);
            item->wtab_precomp[0] = 1;
            item->wqinvtab_precomp[0] = LazyPrepMulModPrecon(1, q, qinv);
            bld.move(item);
            continue;
         }

         item.make();
         item->wtab_precomp.SetLength(1L << (s-1));
         item->wqinvtab_precomp.SetLength(1L << (s-1));

         long m = 1L << s;
         long m_half = 1L << (s-1);
         long m_fourth = 1L << (s-2);

         const long *wtab_last = tab[s-1]->wtab_precomp.elts();
         const mulmod_precon_t *wqinvtab_last = tab[s-1]->wqinvtab_precomp.elts();

         long *wtab = item->wtab_precomp.elts();
         mulmod_precon_t *wqinvtab = item->wqinvtab_precomp.elts();

         for (long i = 0; i < m_fourth; i++) {
            wtab[i] = wtab_last[i];
            wqinvtab[i] = wqinvtab_last[i];
         } 

         long w = root[s];
         mulmod_precon_t wqinv = LazyPrepMulModPrecon(w, q, qinv);

         // prepare wtab...

         if (s == 2) {
            wtab[1] = LazyReduce1(LazyMulModPrecon(wtab[0], w, q, wqinv), q);
            wqinvtab[1] = LazyPrepMulModPrecon(wtab[1], q, qinv);
         }
         else {
            // some software pipelining
            long i, j;

            i = m_half-1; j = m_fourth-1;
            wtab[i-1] = wtab[j];
            wqinvtab[i-1] = wqinvtab[j];
            wtab[i] = LazyReduce1(LazyMulModPrecon(wtab[i-1], w, q, wqinv), q);

            i -= 2; j --;

            for (; i >= 0; i -= 2, j --) {
               long wp2 = wtab[i+2];
               long wm1 = wtab[j];
               wqinvtab[i+2] = LazyPrepMulModPrecon(wp2, q, qinv);
               wtab[i-1] = wm1;
               wqinvtab[i-1] = wqinvtab[j];
               wtab[i] = LazyReduce1(LazyMulModPrecon(wm1, w, q, wqinv), q);
            }

            wqinvtab[1] = LazyPrepMulModPrecon(wtab[1], q, qinv);
         }

         bld.move(item);
      }
   } while (0);
}





void InitAltFFTPrimeInfo(AltFFTPrimeInfo& alt_info, const FFTPrimeInfo& info, long k)
{
   // start by copying...we can't use operator= because there is
   // no operator= for bigtab (which is a UniquePtr).
   // We don't copy zz_p_context either, as we don't need it here.

   if (k < 1 || alt_info.initialized) Error("bad call to InitAltFFTPrimeInfo");

   alt_info.q = info.q; 
   alt_info.qinv = info.qinv;
   alt_info.qrecip = info.qrecip;
   alt_info.RootTable[0] = info.RootTable[0];
   alt_info.RootTable[1] = info.RootTable[1];
   alt_info.TwoInvTable = info.TwoInvTable;
   alt_info.TwoInvPreconTable = info.TwoInvPreconTable;

   // here we make bigtab
   alt_info.bigtab.make();

   for (long dir = 0; dir < 2; dir++) {
      long q = alt_info.q;
      const long *root = alt_info.RootTable[dir].elts();
      mulmod_t qinv = alt_info.qinv;
      const FFTMultipliers& tab = alt_info.bigtab->MulTab[dir];

      LazyPrecompFFTMultipliers(k, q, qinv, root, tab);
   }

   alt_info.initialized = true;
}


void AltFFT(long* A, const long* a, long k, const AltFFTPrimeInfo& info, long dir)

// performs a 2^k-point convolution modulo q

{
   if (!info.initialized) Error("AltFFTPrimeInfo not initialized");

   long q = info.q;
   const long *root = info.RootTable[dir].elts();
   mulmod_t qinv = info.qinv;
   const FFTMultipliers& tab = info.bigtab->MulTab[dir];

   if (k <= 1) {
      if (k == 0) {
	 A[0] = a[0];
	 return;
      }
      if (k == 1) {
	 long a0 = AddMod(a[0], a[1], q);
	 long a1 = SubMod(a[0], a[1], q);
         A[0] = a0;
         A[1] = a1;
	 return;
      }
   }

   // assume k > 1

   if (k >= tab.length()) Error("cannot grow this FFT table to requested length");

   NTL_TLS_LOCAL(Vec<unsigned long>, AA_store);
   AA_store.SetLength(1L << k);
   unsigned long *AA = AA_store.elts();


   long n = 1L << k;

   BitReverseCopy(AA, a, k);


   /* we work with redundant representations, in the range [0, 4q) */



   long s, m, m_half, m_fourth, i, j; 
   unsigned long t, u, t1, u1;


   // s = 1
   for (i = 0; i < n; i += 2) {
      t = AA[i + 1];
      u = AA[i];
      AA[i] = u + t;
      AA[i+1] = u - t + q;
   }


   // s = 2
   {
      const long * NTL_RESTRICT wtab = tab[2]->wtab_precomp.elts();
      const mulmod_precon_t * NTL_RESTRICT wqinvtab = tab[2]->wqinvtab_precomp.elts();

      const long w1 = wtab[1];
      const mulmod_precon_t wqi1 = wqinvtab[1];

      for (i = 0; i < n; i += 4) {

         unsigned long * NTL_RESTRICT AA0 = &AA[i];
         unsigned long * NTL_RESTRICT AA1 = &AA[i + 2];

         {
            const unsigned long a11 = AA1[0];
            const unsigned long a01 = AA0[0];

            const unsigned long tt1 = a11;
            const unsigned long uu1 = a01;
            const unsigned long b01 = uu1 + tt1; 
            const unsigned long b11 = uu1 - tt1 + 2*q;

            AA0[0] = b01;
            AA1[0] = b11;
         }
         {
            const unsigned long a11 = AA1[1];
            const unsigned long a01 = AA0[1];

            const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
            const unsigned long uu1 = a01;
            const unsigned long b01 = uu1 + tt1; 
            const unsigned long b11 = uu1 - tt1 + 2*q;

            AA0[1] = b01;
            AA1[1] = b11;
         }
      }
   }


   //  s = 3..k

   for (s = 3; s <= k; s++) {
      m = 1L << s;
      m_half = 1L << (s-1);
      m_fourth = 1L << (s-2);

      const long* NTL_RESTRICT wtab = tab[s]->wtab_precomp.elts();
      const mulmod_precon_t * NTL_RESTRICT wqinvtab = tab[s]->wqinvtab_precomp.elts();

      for (i = 0; i < n; i += m) {

         unsigned long * NTL_RESTRICT AA0 = &AA[i];
         unsigned long * NTL_RESTRICT AA1 = &AA[i + m_half];

         // a little loop unrolling: this gives the best code

         for (j = 0; j < m_half; j += 4) {
            {
               const long w1 = wtab[j+0];
               const mulmod_precon_t wqi1 = wqinvtab[j+0];
               const unsigned long a11 = AA1[j+0];
               const unsigned long a01 = AA0[j+0];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce2(a01, q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + 2*q;

               AA0[j+0] = b01;
               AA1[j+0] = b11;
            }
            {
               const long w1 = wtab[j+1];
               const mulmod_precon_t wqi1 = wqinvtab[j+1];
               const unsigned long a11 = AA1[j+1];
               const unsigned long a01 = AA0[j+1];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce2(a01, q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + 2*q;

               AA0[j+1] = b01;
               AA1[j+1] = b11;
            }
            {
               const long w1 = wtab[j+2];
               const mulmod_precon_t wqi1 = wqinvtab[j+2];
               const unsigned long a11 = AA1[j+2];
               const unsigned long a01 = AA0[j+2];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce2(a01, q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + 2*q;

               AA0[j+2] = b01;
               AA1[j+2] = b11;
            }
            {
               const long w1 = wtab[j+3];
               const mulmod_precon_t wqi1 = wqinvtab[j+3];
               const unsigned long a11 = AA1[j+3];
               const unsigned long a01 = AA0[j+3];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce2(a01, q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + 2*q;

               AA0[j+3] = b01;
               AA1[j+3] = b11;
            }
         }
      }
   }

   /* need to reduce redundant representations */

   for (i = 0; i < n; i++) {
      unsigned long tmp = LazyReduce2(AA[i], q);
      A[i] = LazyReduce1(tmp, q);
   }
}
