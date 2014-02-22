/* Copyright IBM Corporation 2012 All rights reserved.
 */
#include "NTL/zz_pX.h"
#include "FHE.h"
#include "timing.h"

#include <complex>



void rotate1D(DoubleCRT& d, long i, long amt, 
              const vector< vector< GF2X > > & maskTable)

// rotate d in dimension i by amt

{
  const FHEcontext& context = d.getContext();
  const PAlgebra& al = context.zMstar;
  //  const PAlgebraModTwo& al2 = context.modTwo;

  long ngens = al.numOfGens();
  //  long nslots = al.NSlots();

  assert(i >= 0 && i < ngens);
  long ord = al.OrderOf(i);

  amt = amt % ord;
  if (amt < 0) amt += ord;

  if (amt == 0) return;

  if (al.SameOrd(i)) {
    // "native" rotation
    long val = PowerMod(al.ZmStarGen(i), amt, al.M());
    d.automorph(val);
  }
  else {
    // more expensive "non-native" rotation
    assert(maskTable[i].size() > 0);
    long val = PowerMod(al.ZmStarGen(i), amt, al.M());
    //    long ival = InvMod(val, al.M());
    long ival = PowerMod(al.ZmStarGen(i), amt-ord, al.M());

    ZZX m1 = to_ZZX(maskTable[i][ord-amt]);

    DoubleCRT d1(d);

    d1 *= m1;
    d -= d1;
    d.automorph(val);
    d1.automorph(ival);
    d += d1;
  }
}

long coordinate(const PAlgebra& al, long i, long k)

// returns ith coord of index k
{
  return al.dLog(al.ith_rep(k))[i];
}

void rotate(DoubleCRT& d, long amt, 
            const vector< vector< GF2X > > & maskTable)

// rotate d by amt

{
  const FHEcontext& context = d.getContext();
  const PAlgebra& al = context.zMstar;
  const PAlgebraModTwo& al2 = context.modTwo;

  long ngens = al.numOfGens();
  long nslots = al.NSlots();

  amt = amt % nslots;
  if (amt < 0) amt += nslots;

  if (amt == 0) return;

  long i, v;

  i = ngens-1;
  v = coordinate(al, i, amt);
  rotate1D(d, i, v, maskTable);

  if (i == 0) return;

  GF2X mask = maskTable[i][v];

  DoubleCRT tmp(context, d.getIndexSet());

  for (i--; i >= 0; i--) {
    v = coordinate(al, i, amt);

    tmp.SetZero();
    tmp += d;
    tmp *= to_ZZX(mask);
    d -= tmp;
    rotate1D(d, i, v+1, maskTable);
    rotate1D(tmp, i, v, maskTable); 
    d += tmp;

    if (i > 0) {
      mask = ((mask * (maskTable[i][v] - maskTable[i][v+1])) % al2.PhimXMod())
             + maskTable[i][v+1];
    }
  }
}

//evaluate f at e^{2 pi i/m}, returning a complex number

complex<double> evalPoly(const ZZX& f, long i, long m)
{
  complex<double> t(0.0, 2*M_PI*i/((double) m));
  complex<double> x = exp(t);

  complex<double> res = 0.0;
  for (long j = deg(f); j >= 0; j--)
    res = res*x + to_double(coeff(f, j));

  return res;
}

// return sum_{i in Z_m^*} f(x^i) conj(f(x^i)),
// where x = e^{2 pi i/m}, and m defined by context

double canonEmbedSum(const ZZX& f, const FHEcontext& context)
{
  complex<double> val;
  double nval, res;

  long m = context.zMstar.M();

  res = 0.0;

  for (long i = 0; i < m; i++) {
    if (context.zMstar.inZmStar(i)) {
      val = evalPoly(f, i, m);
      nval = norm(val);
      res += nval;
    }
  }

  return res;
}

// return max_{i in Z_m^*} f(x^i) conj(f(x^i)),
// where x = e^{2 pi i/m}, and m defined by context

double canonEmbedMax(const ZZX& f, const FHEcontext& context)
{
  complex<double> val;
  double nval, res;

  long m = context.zMstar.M();

  res = 0.0;

  for (long i = 0; i < m; i++) {
    if (context.zMstar.inZmStar(i)) {
      val = evalPoly(f, i, m);
      nval = norm(val);
      if (nval > res) res = nval; 
    }
  }

  return res;
}

#define sigma 3.2  /* 3.2 is a nice round number */
#define pSize 18   /* empirical average size of small primes */
#define p0Size 18  /* size of first 1-2 small primes */

void checkCiphertext(const Ctxt& ctxt, const ZZX& ptxt, const FHESecKey& sk);
double c_m = 0.0;

int main(int argc, char *argv[]) 
{
  if (argc<2) {
    cout << "\nUsage: " << argv[0] << " L [c=2 w=64 k=80 d=1]" << endl;
    cout << "  L is the number of levels\n";
    cout << "  optional c is number of columns in the key-switching matrices (default=2)\n";
    cout << "  optional w is Hamming weight of the secret key (default=64)\n";
    cout << "  optional k is the security parameter (default=80)\n";
    cout << "  optional d specifies GF(2^d) arithmetic (default=1, must be <=16)\n";
    //    cout << "  k is the security parameter\n";
    //    cout << "  m determines the ring mod Phi_m(X)" << endl;
    cout << endl;
    exit(0);
  }
  cout.unsetf(ios::floatfield);
  cout.precision(4);

  long L = atoi(argv[1]);
  long c = 2;
  long w = 64;
  long k = 80;
  long d = 1;
  if (argc>2) c = atoi(argv[2]);
  if (argc>3) w = atoi(argv[3]);
  if (argc>4) k = atoi(argv[4]);
  if (argc>5) d = atoi(argv[5]);

  if (d>16) Error("d cannot be larger than 16\n");

  cout << "\nTesting FHE with parameters L="<<L
       << ", c="<<c<<", w="<<w<<", k="<<k<<", d="<<d<< endl;

  // get a lower-bound on the parameter N=phi(m):
  // 1. Empirically, we use ~20-bit small primes in the modulus chain (the main
  //    constraints is that 2m must divide p-1 for every prime p). The first
  //    prime is larger, a 40-bit prime. (If this is a 32-bit machine then we
  //    use two 20-bit primes instead.)
  // 2. With L levels, the largest modulus for "fresh ciphertexts" has size
  //          q0 ~ p0 * p^{L} ~ 2^{40+20L}
  // 3. We break each ciphertext into upto c digits, do each digit is as large
  //    as    D=2^{(40+20L)/c}
  // 4. The added noise variance term from the key-switching operation is
  //    c*N*sigma^2*D^2, and this must be mod-switched down to w*N (so it is
  //    on part with the added noise from modulus-switching). Hence the ratio
  //    P that we use for mod-switching must satisfy c*N*sigma^2*D^2/P^2<w*N,
  //    or    P > sqrt(c/w) * sigma * 2^{(40+20L)/c}
  // 5. With this extra P factor, the key-switching matrices are defined
  //    relative to a modulus of size
  //          Q0 = q0*P ~ sqrt{c/w} sigma 2^{(40+20L)(1+1/c)}
  // 6. To get k-bit security we need N>log(Q0/sigma)(k+110)/7.2, i.e. roughly
  //          N > (40+20L)(1+1/c)(k+110) / 7.2

  long ptxtSpace = 2;
  double cc = 1.0+(1.0/(double)c);
  long N = (long) ceil((pSize*L+p0Size)*cc*(k+110)/7.2);
  cout << "  bounding phi(m) > " << N << endl;

#if 0  // A small m for debugging purposes
  long m = 15;
#else
  // pre-computed values of [phi(m),m,d]
  long ms[][4] = {
    //phi(m)  m  ord(2) c_m*1000
    { 1176,  1247, 28,  3736},
    { 1936,  2047, 11,  3870},
    { 2880,  3133, 24,  3254},
    { 4096,  4369, 16,  3422},
    { 5292,  5461, 14,  4160},
    { 5760,  8435, 24,  8935},
    { 8190,  8191, 13,  1273},
    {10584, 16383, 14,  8358},
    {10752, 11441, 48,  3607},
    {12000, 13981, 20,  2467},
    {11520, 15665, 24, 14916},
    {14112, 18415, 28, 11278},
    {15004, 15709, 22,  3867},
    {15360, 20485, 24, 12767},
 // {16384, 21845, 16, 12798},
    {17208 ,21931, 24, 18387},
    {18000, 18631, 25,  4208},
    {18816, 24295, 28, 16360},
    {19200, 21607, 40, 35633},
    {21168, 27305, 28, 15407},
    {23040, 23377, 48,  5292},
    {24576, 24929, 48,  5612},
    {27000, 32767, 15, 20021},
    {31104, 31609, 71,  5149},
    {42336, 42799, 21,  5952},
    {46080, 53261, 24, 33409},
    {49140, 57337, 39,  2608},
    {51840, 59527, 72, 21128},
    {61680, 61681, 40,  1273},
    {65536, 65537, 32,  1273},
    {75264, 82603, 56, 36484},
    {84672, 92837, 56, 38520}
  };

#if 0

  for (long i = 0; i < 25; i++) {
    long m = ms[i][1];
    PAlgebra alg(m);
    alg.printout();
    cout << "\n";
    // compute phi(m) directly
    long phim = 0;
    for (long j = 0; j < m; j++)
      if (GCD(j, m) == 1) phim++;

    if (phim != alg.phiM()) cout << "ERROR\n";
  }

  exit(0);

#endif

  // find the first m satisfying phi(m)>=N and d | ord(2) in Z_m^*
  long m = 0;
  for (unsigned i=0; i<sizeof(ms)/sizeof(long[3]); i++) 
    if (ms[i][0]>=N && (ms[i][2] % d) == 0) {
      m = ms[i][1];
      c_m = 0.001 * (double) ms[i][3];
      break;
    }
  if (m==0) Error("Cannot support this L,d combination");
#endif
  m = 771;
  FHEcontext context(m);
#if 0
  context.stdev = to_xdouble(0.5); // very low error
#endif
  activeContext = &context; // Mark this as the "current" context

  context.zMstar.printout();
  cout << endl;

  // Set the modulus chain

#if 1
  // The first 1-2 primes of total p0size bits
  #if (NTL_SP_NBITS > p0Size)
    AddPrimesByNumber(context, 1, 1UL<<p0Size); // add a single prime
  #else
    AddPrimesByNumber(context, 2, 1UL<<(p0Size/2)); // add two primes
  #endif
#endif

  // The next L primes, as small as possible
  AddPrimesByNumber(context, L);

  ZZ productOfCtxtPrimes = context.productOfPrimes(context.ctxtPrimes);
  double productSize = context.logOfProduct(context.ctxtPrimes);

  // might as well test that the answer is roughly correct
  cout << "  context.logOfProduct(...)-log(context.productOfPrimes(...)) = "
       << productSize-log(productOfCtxtPrimes) << endl;

  // calculate the size of the digits

  context.digits.resize(c);
  IndexSet s1;
#if 0
  for (long i=0; i<c-1; i++) context.digits[i] = IndexSet(i,i);
  context.digits[c-1] = context.ctxtPrimes / IndexSet(0,c-2);
  AddPrimesByNumber(context, 2, 1, true);
#else
  double sizeSoFar = 0.0;
  double maxDigitSize = 0.0;
  if (c>1) {   // break ciphetext into a few digits
    double dsize = productSize/c;  // initial estimate
    double target = dsize-(pSize/3.0);
    long idx = context.ctxtPrimes.first();
    for (long i=0; i<c-1; i++) { // compute next digit
      IndexSet s;
      while (idx <= context.ctxtPrimes.last() && sizeSoFar < target) {
        s.insert(idx);
	sizeSoFar += log((double)context.ithPrime(idx));
	idx = context.ctxtPrimes.next(idx);
      }
      context.digits[i] = s;
      s1.insert(s);
      double thisDigitSize = context.logOfProduct(s);
      if (maxDigitSize < thisDigitSize) maxDigitSize = thisDigitSize;
      cout << "  digit #"<<i+1<< " " <<s << ": size " << thisDigitSize << endl;
      target += dsize;
    }
    IndexSet s = context.ctxtPrimes / s1; // all the remaining primes
    context.digits[c-1] = s;
    double thisDigitSize = context.logOfProduct(s);
    if (maxDigitSize < thisDigitSize) maxDigitSize = thisDigitSize;
    cout << "  digit #"<<c<< " " <<s << ": size " << thisDigitSize << endl;
  }
  else { 
    maxDigitSize = context.logOfProduct(context.ctxtPrimes);
    context.digits[0] = context.ctxtPrimes;
  }

  // Add primes to the chain for the P factor of key-switching
  double sizeOfSpecialPrimes 
    = maxDigitSize + log(c/(double)w)/2 + log(context.stdev *2);

  AddPrimesBySize(context, sizeOfSpecialPrimes, true);
#endif

  cout << "* ctxtPrimes: " << context.ctxtPrimes 
       << ", log(q0)=" << context.logOfProduct(context.ctxtPrimes) << endl;
  cout << "* specialPrimes: " << context.specialPrimes
       << ", log(P)=" << context.logOfProduct(context.specialPrimes) << endl;

  for (long i=0; i<context.numPrimes(); i++) {
    cout << "  modulus #" << i << " " << context.ithPrime(i) << endl;
  }
  cout << endl;

  setTimersOn();
  const ZZX& PhimX = context.zMstar.PhimX(); // The polynomial Phi_m(X)
  long phim = context.zMstar.phiM();         // The integer phi(m)
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;

#if 0 // Debug mode: use sk=1,2
  DoubleCRT newSk(to_ZZX(2), context);
  long id1 = secretKey.ImportSecKey(newSk, 64, ptxtSpace);
  newSk -= 1;
  long id2 = secretKey.ImportSecKey(newSk, 64, ptxtSpace);
#else
  long id1 = secretKey.GenSecKey(w,ptxtSpace); // A Hamming-weight-w secret key
  long id2 = secretKey.GenSecKey(w,ptxtSpace); // A second Hamming-weight-w secret key
#endif

  ZZX zero = to_ZZX(0);
//  Ctxt zeroCtxt(publicKey);

  /******************************************************************/
  /**                      TESTS BEGIN HERE                       ***/
  /******************************************************************/

  cout << "ptxtSpace = " << ptxtSpace << endl;

  GF2X G;          // G is the AES polynomial, G(X)= X^8 +X^4 +X^3 +X +1
  SetCoeff(G,8); SetCoeff(G,4); SetCoeff(G,3); SetCoeff(G,1); SetCoeff(G,0);
  GF2X X;
  SetX(X);
  
#if 1
  // code for rotations...

  {
    GF2X::HexOutput = 1;
    
    const PAlgebra& al = context.zMstar;
    const PAlgebraModTwo& al2 = context.modTwo;

    long ngens = al.numOfGens();
    long nslots = al.NSlots();

    vector< vector< GF2X > > maskTable;

    maskTable.resize(ngens);
    for (long i = 0; i < ngens; i++) {
      if (i==0 && al.SameOrd(i)) continue;
      long ord = al.OrderOf(i);
      maskTable[i].resize(ord+1);
      for (long j = 0; j <= ord; j++) {
        // initialize the mask that is 1 whenever
        // the ith coordinate is at least j

        vector<GF2X> maps, alphas;

        al2.mapToSlots(maps, G); // Change G to X to get bits in the slots
        alphas.resize(nslots);

        for (long k = 0; k < nslots; k++) 
          if (coordinate(al, i, k) >= j)
            alphas[k] = 1;

       GF2X ptxt;
       al2.embedInSlots(ptxt, alphas, maps);

       maskTable[i][j] = ptxt;
      }
    }

  vector<GF2X> maps;
  al2.mapToSlots(maps, G);

  vector<GF2X> alphas(nslots);
  for (long i=0; i < nslots; i++) 
    random(alphas[i], 8); // random degree-7 polynomial mod 2

  for (long amt = 0; amt < 20; amt++) {

     cout << ".";

     GF2X ptxt;
     al2.embedInSlots(ptxt, alphas, maps);

     DoubleCRT pp(context);
     pp = to_ZZX(ptxt);

     rotate(pp, amt, maskTable);

     GF2X ptxt1 = to_GF2X(to_ZZX(pp));

     vector<GF2X> betas;
     al2.decodePlaintext(betas, ptxt1, G, maps);

     for (long i = 0; i < nslots; i++) {
       if (alphas[i] != betas[(i+amt)%nslots]) {
          cout << "oops\n";
          return 0;
       }
     }
   }

   cout << "\n";

#if 0
  long ord0 = al.OrderOf(0);

  for (long i = 0; i < nslots; i++) {
    cout << alphas[i] << " ";
    if ((i+1) % (nslots/ord0) == 0) cout << "\n";
  }
 
  cout << "\n\n";
  cout << betas.size() << "\n";

  for (long i = 0; i < nslots; i++) {
    cout << betas[i] << " ";
    if ((i+1) % (nslots/ord0) == 0) cout << "\n";
  }
#endif

  return 0;
  }

#endif

  // an initial sanity check on noise estimates,
  // comparing the estimated variance to the actual average
  cout << "pk:"; checkCiphertext(publicKey.pubEncrKey, zero, secretKey);

  ZZX ptxt[6]; // first four are plaintext, last two are constants
  std::vector<Ctxt> ctxt(4, Ctxt(publicKey));

  // Initialize the plaintext and constants to random 0-1 polynomials
  for (size_t j=0; j<6; j++) {
    ptxt[j].rep.SetLength(phim);
    for (long i = 0; i < phim; i++)
      ptxt[j].rep[i] = RandomBnd(ptxtSpace);
    ptxt[j].normalize();

    if (j<4) { 
      publicKey.Encrypt(ctxt[j], ptxt[j], ptxtSpace);
      cout << "c"<<j<<":"; checkCiphertext(ctxt[j], ptxt[j], secretKey);
    }
  }

  // perform upto 2L levels of computation, each level computing:
  //    1. c0 += c1
  //    2. c1 *= c2            // L1' = max(L1,L2)+1
  //    3. c1.reLinearlize
  //    4. c2 *= p4
  //    5. c2.automorph(k)     // k is the first generator of Zm^* /(2)
  //    6. c2.reLinearlize
  //    7. c3 += p5
  //    8. c3 *= c0            // L3' = max(L3,L0,L1)+1
  //    9. c2 *= c3            // L2' = max(L2,L0+1,L1+1,L3+1)+1
  //   10. c0 *= c0            // L0' = max(L0,L1)+1
  //   11. c0.reLinearlize
  //   12. c2.reLinearlize
  //   13. c3.reLinearlize
  //
  // The levels of the four ciphertexts behave as follows:
  // 0, 0, 0, 0  =>  1, 1, 2, 1  =>  2, 3, 3, 2
  //             =>  4, 4, 5, 4  =>  5, 6, 6, 5
  //             =>  7, 7, 8, 7  =>  8,,9, 9, 10  => [...]
  //
  // We perform the same operations on the plaintext, and after each operation
  // we check that decryption still works, and print the curretn modulus and
  // noise estimate. We stop when we get the first decryption error, or when
  // we reach 2L levels (which really should not happen).

  zz_pContext zzpc;
  zz_p::init(ptxtSpace);
  zzpc.save();
  const zz_pXModulus F = to_zz_pX(PhimX);
  long g = context.zMstar.ZmStarGen(0); // the first generator in Zm*
  zz_pX x2g(g, 1);
  zz_pX p2;

  // generate a key-switching matrix from s(X^g) to s(X)
  secretKey.GenKeySWmatrix(/*powerOfS= */  1,
			   /*powerOfX= */  g,
			   0, 0,
			   /*ptxtSpace=*/  ptxtSpace);

  // generate a key-switching matrix from s^2 to s
  secretKey.GenKeySWmatrix(/*powerOfS= */  2,
			   /*powerOfX= */  1,
			   0, 0,
			   /*ptxtSpace=*/  ptxtSpace);

  // generate a key-switching matrix from s^3 to s
  secretKey.GenKeySWmatrix(/*powerOfS= */  3,
			   /*powerOfX= */  1,
			   0, 0,
			   /*ptxtSpace=*/  ptxtSpace);

  for (long lvl=0; lvl<2*L; lvl++) {
    cout << "=======================================================\n";
    ctxt[0] += ctxt[1];
    ptxt[0] += ptxt[1];
    PolyRed(ptxt[0], ptxtSpace, true);
    cout << "c0+=c1:  "; checkCiphertext(ctxt[0], ptxt[0], secretKey);

    ctxt[1].multiplyBy(ctxt[2]);
    ptxt[1] = (ptxt[1] * ptxt[2]) % PhimX;
    PolyRed(ptxt[1], ptxtSpace, true);
    cout << "c1*=c2:  "; checkCiphertext(ctxt[1], ptxt[1], secretKey);

    ctxt[2].multByConstant(ptxt[4]);
    ptxt[2] = (ptxt[2] * ptxt[4]) % PhimX;
    PolyRed(ptxt[2], ptxtSpace, true);
    cout <<  "c2*=p4:  "; checkCiphertext(ctxt[2], ptxt[2], secretKey);

    ctxt[2] >>= g;
    zzpc.restore();
    p2 = to_zz_pX(ptxt[2]);
    CompMod(p2, p2, x2g, F);
    ptxt[2] = to_ZZX(p2);
    cout << "c2>>="<<g<<":"; checkCiphertext(ctxt[2], ptxt[2], secretKey);

    ctxt[2].reLinearize();
    cout << "c2.relin:"; checkCiphertext(ctxt[2], ptxt[2], secretKey);

    ctxt[3].addConstant(ptxt[5]);
    ptxt[3] += ptxt[5];
    PolyRed(ptxt[3], ptxtSpace, true);
    cout << "c3+=p5:  "; checkCiphertext(ctxt[3], ptxt[3], secretKey);

    ctxt[3].multiplyBy(ctxt[0]);
    ptxt[3] = (ptxt[3] * ptxt[0]) % PhimX;
    PolyRed(ptxt[3], ptxtSpace, true);
    cout << "c3*=c0:  ";    checkCiphertext(ctxt[3], ptxt[3], secretKey);

    ctxt[0].square();
    ptxt[0] = (ptxt[0] * ptxt[0]) % PhimX;
    PolyRed(ptxt[0], ptxtSpace, true);
    cout << "c0*=c0:  ";    checkCiphertext(ctxt[0], ptxt[0], secretKey);

    ctxt[2].multiplyBy(ctxt[3]);
    ptxt[2] = (ptxt[2] * ptxt[3]) % PhimX;
    PolyRed(ptxt[2], ptxtSpace, true);
    cout << "c2*=c3:  ";    checkCiphertext(ctxt[2], ptxt[2], secretKey);
  }
  /******************************************************************/
  /**                       TESTS END HERE                        ***/
  /******************************************************************/
  cout << endl;
  return 0;
}

void checkCiphertext(const Ctxt& ctxt, const ZZX& ptxt, const FHESecKey& sk)
{
  const FHEcontext& context = ctxt.getContext();
  /*
  IndexSet base = baseSetOf(ctxt);
  double addedNoise = log(ctxt.modSwitchAddedNoiseVar());
  Ctxt tmp = ctxt;
  tmp.modDownToSet(base);
  double totalNoise = log(tmp.getNoiseVar());
  cout << "   @@@ log(added-noise)="<<addedNoise
       << ", log(total-noise)="<<totalNoise<<endl;
  */
  cout << " ln(q)="<< context.logOfProduct(ctxt.getPrimeSet())
       << ", ln(nVar)/2="<< log(ctxt.getNoiseVar())/2;
  //       << ", ln(nMag)="<< log(ctxt.getNoiseMag());

  ZZX res;
  //  sk.Decrypt(res, ctxt);
  ZZX f;
  sk.Decrypt(res, ctxt, f);
  cout << ", ln(mxPtxtCoef)=" << log(largestCoeff(f));

  // ensure we reduce the same way on both
  PolyRed((ZZX&)res,res,ctxt.getPtxtSpace(),true);
  PolyRed((ZZX&)ptxt,ptxt,ctxt.getPtxtSpace(),true);
  if (res != ptxt) {
    cout << ", failed\n";
    for (long i=0; i<=deg(ptxt); i++) if (coeff(res,i)!=coeff(ptxt,i)) {
	cout << "first mismatch in coeff "<<i<<": "
	     << coeff(res,i)<<"!="<<coeff(ptxt,i)<<"\n";
	break;
      }

    cout << "Timing information:\n";
    printAllTimers();
    cout << "\n";
    exit(0);
  }
  else cout << ", succeeded\n";
}

/***********************************************************************
 ***************************** UNUSED CODE *****************************
 ***********************************************************************/

#if 0
void adjustLevelForMult(Ctxt& c1, const char name1[], const ZZX& p1,
			Ctxt& c2, const char name2[], const ZZX& p2, 
			const FHESecKey& sk)
{
  const FHEcontext& context = c1.getContext();

  // The highest possible level for this multiplication is the
  // intersection of the two primeSets, without the special primes.
  IndexSet primes = c1.getPrimeSet() & c2.getPrimeSet();
  primes.remove(context.specialPrimes);
  assert (!empty(primes));

  //  double phim = (double) context.zMstar.phiM();
  //  double factor = c_m*sqrt(log(phim))*4;

  xdouble n1,n2,d1,d2;
  xdouble dvar1 = c1.modSwitchAddedNoiseVar();
  xdouble dvar2 = c2.modSwitchAddedNoiseVar();
  //  xdouble dmag1 = c1.modSwitchAddedNoiseMag(c_m);
  //  xdouble dmag2 = c2.modSwitchAddedNoiseMag(c_m);
  //  cout << " ** log(dvar1)=" << log(dvar1) 
  //       << ", log(dvar2)=" << log(dvar2) <<endl;

  double logF1, logF2;
  xdouble n1var, n2var, modSize; // n1mag, n2mag,
  // init to large number
  xdouble noiseVarRatio=xexp(2*(context.logOfProduct(context.ctxtPrimes)
				+ context.logOfProduct(context.specialPrimes)));
  //  xdouble noiseMagRatio=noiseVarRatio;

  // Find the level that minimizes the noise-to-modulus ratio
  bool oneLevelMore = false;
  for (IndexSet levelDown = primes; 
       !empty(levelDown); levelDown.remove(levelDown.last())) {

    // compute noise variane/magnitude after mod-switchign to this level
    logF1 = context.logOfProduct(c1.getPrimeSet() / levelDown);
    n1var = c1.getNoiseVar()/xexp(2*logF1);

    logF2 = context.logOfProduct(c2.getPrimeSet() / levelDown);
    n2var = c2.getNoiseVar()/xexp(2*logF2);

    // compute modulus/noise ratio at this level
    modSize = xexp(context.logOfProduct(levelDown));
    xdouble nextNoiseVarRatio = sqrt((n1var+dvar1)*(n2var+dvar2))/modSize;

    if (nextNoiseVarRatio < 2*noiseVarRatio || oneLevelMore) {
      noiseVarRatio = nextNoiseVarRatio;
      primes = levelDown;  // record the currently best prime set
      n1 = n1var; d1=dvar1; n2 = n2var; d2=dvar2;
    }
    oneLevelMore = (n1var > dvar1 || n2var > dvar2);
  }

  if (primes < c1.getPrimeSet()) {
    cout << "          ** " << c1.getPrimeSet()<<"=>"<<primes << endl;
    cout << "             n1var="<<n1<<", d1var="<<d1<<endl;;
    c1.modDownToSet(primes);
    cout << name1 << ".mDown:"; checkCiphertext(c1, p1, sk);
  }
  if (primes < c2.getPrimeSet()) {
    cout << "          ** " << c2.getPrimeSet()<<"=>"<<primes << endl;
    cout << "             n2var="<<n2<<", d2var="<<d2<<endl;;
    c2.modDownToSet(primes);
    cout << name2 << ".mDown:"; checkCiphertext(c2, p2, sk);
  }
}
#endif

#if 0
  ZZX ptxt, ptxt2, res;
  Ctxt ctxt(publicKey), ctxt2(publicKey);

  ptxt.rep.SetLength(phim);
  for (long i = 0; i < phim; i++)
    ptxt.rep[i] = RandomBnd(2);
  ptxt.normalize();

  // Basis test for encryption and decryption
  publicKey.Encrypt(ctxt, ptxt, ptxtSpace);
  ctxt2 = ctxt;                  // test also the Ctxt assignment method
  secretKey.Decrypt(res, ctxt2,ptxtSpace);

  cout << "test  #1: " << ((ptxt == res) ? "PASS" : "FAIL") << " Enc/Dec\n";

  // Some tests here to check the mod-switch up and down
  s1 = ctxt2.getPrimeSet();
  IndexSet s = s1 | context.specialPrimes;
  ctxt2.modUpToSet(s);
  secretKey.Decrypt(res, ctxt2);

  cout << "test  #2: " <<  ((ptxt == res) ? "PASS" : "FAIL") 
       << " Enc/Mod-UP/Dec\n";

  ctxt2.modDownToSet(s1);
  secretKey.Decrypt(res, ctxt2);

  cout << "test  #3: " <<  ((ptxt == res) ? "PASS" : "FAIL") 
       << " Enc/Mod-UP/mod-DOWN/Dec\n";
  
  // generate a key-switching matrix from key id2 to key id1
  secretKey.GenKeySWmatrix(/*powerOfS= */  1,
			   /*powerOfX= */  1,
			   /*fromKeyID=*/ id2,
			   /*toKeyID=  */ id1,
			   /*ptxtSpace=*/  ptxtSpace);

  secretKey.Encrypt(ctxt, ptxt, ptxtSpace, id2);
  secretKey.Decrypt(res, ctxt);
  cout << "test  #4: " 
       << ((ptxt == res) ? "PASS" : "FAIL") 
       << " skEnc/Dec (key-ID "<<id2<<")\n";

  ctxt.reLinearize(id1); // switch from key-id2 to key-id1
  secretKey.Decrypt(res, ctxt);

  cout << "test  #5: "
       << ((ptxt == res) ? "PASS" : "FAIL") 
       << " Enc/key-switch/dec\n";

  // check basic encryption equations

  ZZX c0, c1, sk1;
  ZZX tmp1, tmp2, tmp3;

  publicKey.Encrypt(ctxt, ptxt);

  ctxt[0].toPoly(c0);
  ctxt[1].toPoly(c1);
  secretKey.sKeys[id1].toPoly(sk1);

  tmp1 = (c0 + c1*sk1) % PhimX;
  PolyRed(tmp1, productOfCtxtPrimes);
  PolyRed(tmp1, ptxtSpace, true);

  cout << "test  #6: " <<
          ((ptxt == tmp1) ? "PASS" : "FAIL") 
       << " verify decryption equations\n";

  // Test automorphism
  long pOfX = 3;
  while (!context.zMstar.inZmStar(pOfX)) pOfX += 2;

  ctxt.automorph(pOfX);
  secretKey.Decrypt(res, ctxt);

  ModComp(ptxt2, ptxt, ZZX(pOfX, 1), PhimX);
  PolyRed(ptxt2, ptxtSpace, true);

  cout << "test  #7: " << ((ptxt2 == res) ? "PASS" : "FAIL") 
       << " Enc/automorph("<<pOfX<<")/Dec\n";

  // generate a key-switching matrix from s(X^pOfX) to s(X)
  secretKey.GenKeySWmatrix(/*powerOfS= */  1,
			   /*powerOfX= */  pOfX,
			   /*fromKeyID=*/ id1,
			   /*toKeyID=  */ id1,
			   /*ptxtSpace=*/  ptxtSpace);
  ctxt.reLinearize(id1); 
  secretKey.Decrypt(res, ctxt);

  cout << "test  #8: " << ((ptxt2 == res) ? "PASS" : "FAIL") 
       << " Enc/automorph("<<pOfX<<")/relin/Dec\n";

  // Test add/mult by constant

  ptxt2.rep.SetLength(phim);
  for (long i = 0; i < phim; i++)
    ptxt2.rep[i] = RandomBnd(2);
  ptxt2.normalize();

  publicKey.Encrypt(ctxt, ptxt);
  ctxt.addConstant(ptxt2);
  secretKey.Decrypt(res, ctxt);

  ptxt += ptxt2;
  PolyRed(ptxt, ptxtSpace, true);

  cout << "test  #9: " << ((ptxt == res) ? "PASS" : "FAIL") 
       << " Enc/addConst/Dec\n";

  ctxt.multByConstant(ptxt2);
  secretKey.Decrypt(res, ctxt);

  ptxt = (ptxt * ptxt2) % PhimX;
  PolyRed(ptxt, ptxtSpace, true);

  cout << "test #10: " << ((ptxt == res) ? "PASS" : "FAIL") 
       << " Enc/multByConst/Dec\n";

  // Test addition

  publicKey.Encrypt(ctxt, ptxt);
  publicKey.Encrypt(ctxt2, ptxt2);

  ctxt2 += ctxt;
  ptxt2 += ptxt;
  PolyRed(ptxt2, ptxtSpace, true);
  secretKey.Decrypt(res, ctxt2);

  cout << "test #11: " << ((ptxt2 == res) ? "PASS" : "FAIL") 
       << " Enc/add/Dec\n";

  // Test multiplication

  publicKey.Encrypt(ctxt, ptxt);
  publicKey.Encrypt(ctxt2, ptxt2);

  ctxt2 *= ctxt;
  ptxt2 = (ptxt * ptxt2) % PhimX;
  PolyRed(ptxt2, ptxtSpace, true);
  secretKey.Decrypt(res, ctxt2);

  cout << "test #12: " << ((ptxt2 == res) ? "PASS" : "FAIL") 
       << " Enc/mul/Dec\n";

  // generate a key-switching matrix from s^2 to s (keyID=id1) 
  secretKey.GenKeySWmatrix(/*powerOfS= */  2,
			   /*powerOfX= */  1,
			   /*fromKeyID=*/ id1,
			   /*toKeyID=  */ id1,
			   /*ptxtSpace=*/  ptxtSpace);
  ctxt2.reLinearize(id1); 

  secretKey.Decrypt(res, ctxt2);

  cout << "test #13: " << ((ptxt2 == res) ? "PASS" : "FAIL") 
       << " Enc/mul/relin/Dec\n";

#endif

#if 0

  // Generate a key-switching matrix from s^2(X^1) to s with plaintext space 2
secretKey.GenKeySWmatrix(2,1,0,0,ptxtSpace);

  // For each generator gi of Zm*/<2> of order oi, generate a key-switching
  // matrix from s(X^{gi^j}) to s for j=1,2,3,...,oi-1. In addition, if the
  // order of gi in Zm* is larger than oi, then generate a key-switching
  // matrix from s(X^{gi^j}) to s for for j=-1,-2,...,1-oi

  long nGens = context.zMstar.numOfGens();
  for (long i=0; i<nGens; i++) {
    long g = context.zMstar.ZmStarGen(i);
    long ord = context.zMstar.OrderOf(i);
    long gj = 1;
    for (long j=1; j<ord; j++) {
      gj = MulMod(gj, g, m); // g^j mod m
      // Generate a key-switch matrix from s(X^gj) to s with plaintext space 2
      secretKey.GenKeySWmatrix(1,gj,0,0,ptxtSpace);
      cout << ".";
    }
    if (!context.zMstar.SameOrd(i)){ // order in Zm* not the same as in Zm*/<2>
      g = InvMod(g,m);
      gj = 1;
      for (long j=1; j<ord; j++) {
	gj = MulMod(gj, g, m); // g^{-j} mod m
	// Generate a key-switch matrix from s(X^gj) to s with plaintext space 2
	secretKey.GenKeySWmatrix(1,gj,0,0,ptxtSpace);
	cout << ".";
      }
    }
  }
  cout << "  generated "<< secretKey.keySwitching.size() 
       << " key-switching matrices\n";

  /******************************************************************/
  /******************************************************************/
  DoubleCRT sss;
  sss = secretKey.sKeys[id1];
  sss.removePrimes( sss.getIndexSet() / context.ctxtPrimes );

  DoubleCRT tcrt;
  DoubleCRT tcrt1;

/*
  tcrt = ctxt[2];
  tcrt *= sss;
  tcrt += ctxt[1];
  tcrt *= sss;
  tcrt += ctxt[0];
*/

  DoubleCRT testvec[3];

  tcrt = ctxt[0];
  testvec[0] = tcrt;

  tcrt1 = sss;
  tcrt1 *= ctxt[1];
  tcrt += tcrt1;
  testvec[1] = tcrt;


  tcrt1 = sss;
  tcrt1.Exp(2);
  tcrt1 *= ctxt[2];
  tcrt += tcrt1;
  testvec[2] = tcrt;


  tcrt.toPoly(tmp1);
  PolyRed(tmp1, ptxtSpace, true);


  cout << "test #10: " <<
          ((ptxt2 == tmp1) ? "PASS" : "FAIL") 
       << " Enc/mul/decrypt by hand\n";

  secretKey.Decrypt(tcrt1, ctxt, testvec);

  cout << "test #11: " <<
          ((tcrt1 == tcrt) ? "PASS" : "FAIL") 
       << " Enc/mul/dcrt decrypt\n";

  sss.Exp(2);
  sss.toPoly(tmp1);
  PolyRed(tmp1, productOfCtxtPrimes);

  tmp2 = (sk1 * sk1) % PhimX;
  PolyRed(tmp2, productOfCtxtPrimes);

  cout << "test #12: " <<
          ((tmp1 == tmp2) ? "PASS" : "FAIL") 
       << " Exp test\n";

  tmp1 = 0;
  for (size_t i = 0; i < ctxt.size(); i++) {
    ctxt[i].toPoly(tmp2);
    for (long j = 0; j < ctxt[i].skHandle.getPowerOfS(); j++) {
      tmp2 *= sk1;
      tmp2 %= PhimX;
      PolyRed(tmp2, productOfCtxtPrimes);
    }

    tmp1 += tmp2;
    PolyRed(tmp1, productOfCtxtPrimes);
  }

  PolyRed(tmp1, ptxtSpace, true);

  cout << "test #13: " <<
          ((ptxt2 == tmp1) ? "PASS" : "FAIL") 
       << " Enc/mul/verify equations\n";
#endif

#if 0
  {
    long ITER = 10;

    cout << "fresh ciphertexts\n";
    for (long i = 0; i < ITER; i++) {
      ZZX ptxt = to_ZZX(0);
      ZZX res, f;
      Ctxt ctxt(publicKey);

      publicKey.Encrypt(ctxt, ptxt, ptxtSpace);
      secretKey.Decrypt(res, ctxt, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt.getNoiseVar();
      cout << vbar << " " << est << " " << (vbar - est)/est << "\n";
    }

    cout << "product ciphertexts\n";
    for (long i = 0; i < ITER; i++) {
      ZZX ptxt = to_ZZX(0);
      ZZX res, f;
      Ctxt ctxt(publicKey), ctxt2(publicKey);

      publicKey.Encrypt(ctxt, ptxt, ptxtSpace);
      publicKey.Encrypt(ctxt2, ptxt, ptxtSpace);
      ctxt *= ctxt2;
      secretKey.Decrypt(res, ctxt, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt.getNoiseVar();
      cout << vbar << " " << est << " " << (vbar - est)/est << "\n";
    }

    cout << "product ciphertexts modded down\n";
    for (long i = 0; i < ITER; i++) {
      ZZX ptxt = to_ZZX(0);
      ZZX res, f;
      Ctxt ctxt(publicKey), ctxt2(publicKey);

      publicKey.Encrypt(ctxt, ptxt, ptxtSpace);
      publicKey.Encrypt(ctxt2, ptxt, ptxtSpace);

      IndexSet s = ctxt.getPrimeSet();
      s.remove(s.last());
      s.remove(s.last());
      s.remove(s.last());
      ctxt.modDownToSet(s);
      ctxt2.modDownToSet(s);

      ctxt *= ctxt2;
      secretKey.Decrypt(res, ctxt, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt.getNoiseVar();
      cout << vbar << " " << est << " " << (vbar - est)/est << "\n";
    }

    cout << "squared ciphertexts modded down\n";
    for (long i = 0; i < ITER; i++) {
      ZZX ptxt = to_ZZX(0);
      ZZX res, f;
      Ctxt ctxt(publicKey), ctxt2(publicKey);

      publicKey.Encrypt(ctxt, ptxt, ptxtSpace);
      publicKey.Encrypt(ctxt2, ptxt, ptxtSpace);

      IndexSet s = ctxt.getPrimeSet();
      s.remove(s.last());
      s.remove(s.last());
      s.remove(s.last());
      ctxt.modDownToSet(s);
      ctxt2.modDownToSet(s);

      ctxt *= ctxt;
      secretKey.Decrypt(res, ctxt, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt.getNoiseVar();
      cout << vbar << " " << est << " " << (vbar - est)/est << "\n";
    }

    cout << "some test...\n";
    for (long i = 0; i < ITER; i++) {
      ZZX ptxt = to_ZZX(0);
      ZZX res, f;
      Ctxt ctxt1(publicKey), ctxt2(publicKey), 
           ctxt3(publicKey), ctxt4(publicKey);
      
      cout << "***\n";

      publicKey.Encrypt(ctxt1, ptxt, ptxtSpace);
      publicKey.Encrypt(ctxt2, ptxt, ptxtSpace);
      publicKey.Encrypt(ctxt3, ptxt, ptxtSpace);
      publicKey.Encrypt(ctxt4, ptxt, ptxtSpace);

      IndexSet s = ctxt1.getPrimeSet();
{
      secretKey.Decrypt(res, ctxt1, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt1.getNoiseVar();
      cout << vbar << " " << est << " " << vbar/est << "\n";
}


      s.remove(s.last());
      ctxt1.modDownToSet(s);
      ctxt2.modDownToSet(s);
      ctxt3.modDownToSet(s);
      ctxt4.modDownToSet(s);
      ctxt4.modDownToSet(s);
{
      secretKey.Decrypt(res, ctxt1, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt1.getNoiseVar();
      cout << vbar << " " << est << " " << vbar/est << "\n";
}

      ctxt1 *= ctxt2;
      // ctxt3 *= ctxt4;

{
      secretKey.Decrypt(res, ctxt1, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt1.getNoiseVar();
      cout << vbar << " " << est << " " << vbar/est << "\n";
}

      s.remove(s.last());
      s.remove(s.last());
      ctxt1.modDownToSet(s);
      ctxt3.modDownToSet(s);
{
      secretKey.Decrypt(res, ctxt1, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt1.getNoiseVar();
      cout << vbar << " " << est << " " << vbar/est << "\n";
}

      ctxt1 *= ctxt3;

{
      secretKey.Decrypt(res, ctxt1, f);
      double vbar = canonEmbedSum(f, context)/context.zMstar.phiM();
      xdouble est = ctxt1.getNoiseVar();
      cout << vbar << " " << est << " " << vbar/est << "\n";
}
    }
    exit(0);
  }
#endif
