#include "SPX.h"
#include <NTL/lzz_pXFactoring.h>

using namespace std;
using namespace NTL;

#if ((NTL_MAJOR_VERSION<9) || ((NTL_MAJOR_VERSION==9)&&(NTL_MINOR_VERSION<2)))
#warning "NTL < 9.2.0, using SPmodulus wrapper"
#else
#warning "NTL >= 9.2.0, using zz_pContext directly"
#endif

static void zpxFromSPX(zz_pX& zpx, const SPX& spp)
{
  Vec<long> v;
  spp.getCoeffVec(v);  // convert to zz_pX
  conv(zpx.rep, v);
  zpx.normalize();
}

static bool testGCD(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

   SPX b, c, ss, tt;
   zz_pX a_p, b_p, c_p, ss_p, tt_p;

  // test polynomials of various degrees, from 32 to 8192
  for (long n = 32; n <= 8192; n <<= 2) {
    random(a, n, mod);   // get random polynomial
    zpxFromSPX(a_p,a);// convert to zz_pX
    random(b_p, n);      // get random zz_pX
    b = SPX(b_p,mod); // convert to SPX

    GCD(c_p, a_p, b_p);  // work on zz_pX
    GCD(c, a, b);        // work on SPX
    if (c != SPX(c_p,mod)) {
      cerr << "**** SPX Test GCD FAILED!\n";
      return false;
    }

    XGCD(c_p, ss_p, tt_p, a_p, b_p);
    XGCD(c, ss, tt, a, b);
    if (c != SPX(c_p,mod) ||
	ss != SPX(ss_p,mod) || tt != SPX(tt_p,mod)) {
      cerr << "**** SPX Test XGCD FAILED!\n";
      return false;
    }
  }
  return true;
}

static bool testArith(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

   SPX b, c, d;
   zz_pX a_p, b_p, c_p, d_p;

  // test polynomials of various degrees, from 32 to 8192
  for (long n = 32; n <= 8192; n <<= 2) {
    random(a, n, mod);   // get random polynomial
    zpxFromSPX(a_p,a);// convert to zz_pX
    random(b_p, n);      // get random zz_pX
    b = SPX(b_p,mod); // convert to SPX

    c   = (a  +b);
    c_p = (a_p+b_p);
    if (c != SPX(c_p,mod)) {
      cerr << "  ** SPX testArith (1) FAILED!\n";
      return false;
    }
    //    cout << "." << std::flush;

    DivRem(c,  d,  c,  b);  // c = c/b, d=c%b
    DivRem(c_p,d_p,c_p,b_p);
    if (c != SPX(c_p,mod) || d != SPX(d_p,mod)) {
      cerr << "  ** SPX testArith (2) FAILED!\n";
      return false;
    }
    //    cout << "," << std::flush;

    d  = (c << 5) - (b >> 2);
    d_p= (c_p<<5) - (b_p>>2);
    if (d != SPX(d_p,mod)) {
      cerr << "  ** SPX testArith (3) FAILED!\n";
      return false;
    }
    //    cout << ";" << std::flush;

    --d;   ++c;
    --d_p; ++c_p;
    if (c != SPX(c_p,mod) || d != SPX(d_p,mod)) {
      cerr << "  ** SPX testArith (4) FAILED!\n";
      return false;
    }
    //    cout << ";" << std::flush;
  }
  return true;
}

static bool testCoeffs(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

   zz_pX a_p;

  // test polynomials of various degrees, from 32 to 8192
  for (long n = 32; n <= 8192; n <<= 2) {
    random(a, n, mod);   // get random polynomial
    zpxFromSPX(a_p,a);// convert to zz_pX

    if (LeadCoeff(a) != conv<long>(LeadCoeff(a_p))
	|| ConstTerm(a) != conv<long>(ConstTerm(a_p))) {
      cerr << "**** SPX testCoeffs (1) FAILED!\n";
      return false;
    }

    SetCoeff(a, 10, 3);
    SetCoeff(a_p, 10, 3);
    if (a != SPX(a_p,mod)) {
      cerr << "**** SPX testCoeffs (2) FAILED!\n";
      return false;
    }

    SetX(a);
    SetX(a_p);
    if (!IsX(a) || a != SPX(a_p,mod))  {
      cerr << "**** SPX testCoeffs (3) FAILED!\n";
      return false;
    }

    clear(a);
    clear(a_p);
    if (!IsZero(a) || a != SPX(a_p,mod))  {
      cerr << "**** SPX testCoeffs (4) FAILED!\n";
      return false;
    }

    set(a);
    set(a_p);
    if (!IsOne(a) || a != SPX(a_p,mod))  {
      cerr << "**** SPX testCoeffs (5) FAILED!\n";
      return false;
    }
  }
  return true;
}

static bool testMisc(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

  SPX b;
  zz_pX a_p, b_p;

  // test polynomials of various degrees, from 32 to 8192
  for (long n = 32; n <= 8192; n <<= 2) {
    random(a, n, mod);   // get random polynomial
    zpxFromSPX(a_p,a);// convert to zz_pX

    diff(b,  a);
    diff(b_p,a_p);
    if (b != SPX(b_p,mod)) {
      cerr << "**** SPX testMisc (1) FAILED!\n";
      return false;
    }
    reverse(b,  a);
    reverse(b_p,a_p);
    if (b != SPX(b_p,mod)) {
      cerr << "**** SPX testMisc (2) FAILED!\n";
      return false;
    }
    reverse(b,  a,  10);
    reverse(b_p,a_p,10);
    if (b != SPX(b_p,mod)) {
      cerr << "**** SPX testMisc (3) FAILED!\n";
      return false;
    }
  }
  return true;
}


static bool testTrunc(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

  SPX b, c;
  zz_pX a_p, b_p, c_p;

  // test polynomials of various degrees, from 32 to 8192
  for (long n = 32; n <= 8192; n <<= 2) {
    random(a, n, mod);   // get random polynomial
    zpxFromSPX(a_p,a);// convert to zz_pX
    random(b_p, n);      // get random zz_pX
    b = SPX(b_p,mod); // convert to SPX

    MulTrunc(c,  a,  b,  32);    MulTrunc(c_p,a_p,b_p,32);
    SqrTrunc(c,  c,  48);        SqrTrunc(c_p,c_p,48);
    if (c != SPX(c_p,mod)) {
      cerr << "**** SPX testTrunc (1) FAILED!\n";
      return false;
    }

    if (coeff(c,0)==0) {
      SetCoeff(c,0);  SetCoeff(c_p,0);
    }
    InvTrunc(c,  c,  40);        InvTrunc(c_p,c_p,40);
    if (c != SPX(c_p,mod)) {
      cerr << "**** SPX testTrunc (2) FAILED!\n";
      return false;
    }
  }
  return true;
}

static bool testMod(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

  SPX b, c;
  zz_pX a_p, b_p, c_p, F_p;

  // test polynomials of various degrees, from 32 to 512
  for (long n = 32; n <= 512; n <<= 2) {
    BuildIrred(F_p,n+1);  // build an irreducible polynomial
    SPX F(F_p,mod);    // convert to SPX

    random(a, n, mod);   // get random polynomial
    zpxFromSPX(a_p,a);// convert to zz_pX
    random(b_p, n);      // get random zz_pX
    b = SPX(b_p,mod); // convert to SPX

    MulMod(c,  a,  b,  F);    MulMod(c_p,a_p,b_p,F_p);
    SqrMod(c,  c,  F);        SqrMod(c_p,c_p,F_p);
    if (c != SPX(c_p,mod)) {
      cerr << "**** SPX testMod (1) FAILED!\n";
      return false;
    }

    MulByXMod(c, a, F); MulByXMod(c_p,a_p,F_p);
    InvMod(c, c, F);    InvMod(c_p,c_p,F_p);
    if (c != SPX(c_p,mod)) {
      cerr << "**** SPX testMod (2) FAILED!\n";
      return false;
    }
  }
  return true;
}

static bool testModPre(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

  SPX b, c;
  zz_pX a_p, b_p, c_p;

  // test polynomials of various degrees, from 32 to 512
  for (long n = 32; n <= 512; n <<= 2) {
    BuildIrred(c_p,n+1);  // build an irreducible polynomial
    zz_pXModulus F_p=c_p; // set it as modulus
    c = SPX(c_p,mod);  // convert to SPX
    SPXModulus F=c;

    random(a, n, mod);   // get random polynomial
    zpxFromSPX(a_p,a);// convert to zz_pX
    random(b_p, n);      // get random zz_pX
    b = SPX(b_p,mod); // convert to SPX

    MulMod(c,  a,  b,  F);    MulMod(c_p,a_p,b_p,F_p);
    SqrMod(c,  c,  F);        SqrMod(c_p,c_p,F_p);
    if (c != SPX(c_p,mod)) {
      cerr << "**** SPX testModPre (1) FAILED!\n";
      return false;
    }

    PowerXMod(b, n+5, F); PowerXMod(b_p,n+5,F_p);
    MulByXMod(c, b, F);   MulByXMod(c_p,b_p,F_p);
    if (c != SPX(c_p,mod)) {
      cerr << "**** SPX testModPre (2) FAILED!\n";
      return false;
    }

    PowerMod(b, c, n/2, F); PowerMod(b_p, c_p, n/2, F_p);
    if (b != SPX(b_p,mod)) {
      cerr << "**** SPX testModPre (3) FAILED!\n";
      return false;
    }

    random(a_p, 2*n);    // get random zz_pX
    a = SPX(a_p,mod); // convert to SPX
    DivRem(b, c, a, F);  DivRem(b_p,c_p,a_p,F_p);
    if (b != SPX(b_p,mod) || c != SPX(c_p,mod)) {
      cerr << "**** SPX testModPre (4) FAILED!\n";
      return false;
    }

    b = c = a;
    b_p = c_p = a_p;
    b /= F; b_p /= F_p;
    c %= F; c_p %= F_p;
    if (b != SPX(b_p,mod) || c != SPX(c_p,mod)) {
      cerr << "**** SPX testModPre (5) FAILED!\n";
      return false;
    }
  }
  return true;
}


static bool testComp(SPX& a1)
{
  const SPmodulus& mod = a1.getMod();
  getContext(mod).restore();

  SPX a2, a3, b1, b2, b3, c;
  zz_pX a1_p,a2_p,a3_p,b1_p,b2_p,b3_p,c_p;

  // test polynomials of various degrees, from 32 to 512
  for (long n = 32; n <= 512; n <<= 2) {
    BuildIrred(c_p,n+1);  // build an irreducible polynomial
    zz_pXModulus F_p=c_p; // set it as modulus
    c = SPX(c_p,mod);  // convert to SPX
    SPXModulus F=c;

    random(b1_p, n);      // get random zz_pX
    b1 = SPX(b1_p,mod); // convert to SPX

    random(b2_p, n);      // get random zz_pX
    b2 = SPX(b2_p,mod); // convert to SPX

    random(b3_p, n);      // get random zz_pX
    b3 = SPX(b3_p,mod); // convert to SPX

    random(c_p, 3);      // get random zz_pX
    c = SPX(c_p,mod); // convert to SPX

    CompMod(a1, b1, c, F); CompMod(a1_p,b1_p,c_p,F_p);
    if (a1 != SPX(a1_p,mod)) {
      cerr << "**** SPX testComp (1) FAILED!\n";
      return false;
    }

    Comp2Mod(a1, a2, b1, b2, c, F); Comp2Mod(a1_p,a2_p,b1_p,b2_p,c_p,F_p);
    if (a1 != SPX(a1_p,mod) || a2 != SPX(a2_p,mod)) {
      cerr << "**** SPX testComp (2) FAILED!\n";
      return false;
    }

    Comp3Mod(a1,  a2,  a3,  b1,  b2,  b3,  c,  F);
    Comp3Mod(a1_p,a2_p,a3_p,b1_p,b2_p,b3_p,c_p,F_p);
    if (a1 != SPX(a1_p,mod)
	|| a2 != SPX(a2_p,mod) || a3 != SPX(a3_p,mod)) {
      cerr << "**** SPX testComp (3) FAILED!\n";
      return false;
    }

    SPXArgument H;  build(H,  c,  F, 4);
    zz_pXArgument H_p; build(H_p,c_p,F_p,4);

    CompMod(a1, b1, H, F); CompMod(a1_p,b1_p,H_p,F_p);
    if (a1 != SPX(a1_p,mod)) {
      cerr << "**** SPX testComp (4) FAILED!\n";
      return false;
    }
  }
  return true;
}

static bool testProj(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

  Vec<long> v;
  SPX f;

  Vec<zz_p> v_p;
  zz_pX a_p,f_p;

  // test polynomials of various degrees, from 32 to 512
  for (long n = 32; n <= 512; n <<= 2) {
    BuildIrred(f_p,n+1);  // build an irreducible polynomial
    f = SPX(f_p,mod);  // convert to SPX

    random(a_p, n);      // get random zz_pX
    a = SPX(a_p,mod); // convert to SPX

    VectorCopy(v_p,a_p,deg(a_p)+1); // get coefficient vector
    a.getCoeffVec<long>(v);

    random(a_p, n);      // get random zz_pX
    a = SPX(a_p,mod); // convert to SPX

    {long x = project(v, a);
    zz_p x_p = project(v_p, a_p);
    if (x != conv<long>(x_p)) {
      cerr << "**** SPX testProj (1) FAILED! ("
	   << x << "!=" << conv<long>(x_p) << ")\n";
      return false;
    }}
    Vec<long> x;
    Vec<zz_p> x_p;
    ProjectPowers(x,  v,  5, a, f);
    ProjectPowers(x_p,v_p,5,a_p,f_p);
    if (x != conv< Vec<long> >(x_p)) {
      cerr << "**** SPX testProj (2) FAILED!\n";
      return false;
    }

    SPXModulus F=f;
    ProjectPowers(x,  v,  5, a, F);
    if (x != conv< Vec<long> >(x_p)) {
      cerr << "**** SPX testProj (3) FAILED!\n";
      return false;
    }
  }
  return true;
}

/*
static bool testMinp(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();
  // ...
}
*/

static bool testTrace(SPX& a)
{
  const SPmodulus& mod = a.getMod();
  getContext(mod).restore();

  SPX f;
  zz_pX ap, fp;

  // test polynomials of various degrees, from 32 to 512
  for (long n = 32; n <= 512; n <<= 2) {
    BuildIrred(fp,n+1);  // build an irreducible polynomial
    f = SPX(fp,mod);  // convert to SPX
    SPXModulus F=f;

    random(ap, n);      // get random zz_pX
    a = SPX(ap,mod); // convert to SPX

    long t = TraceMod(a, f);
    zz_p tp= TraceMod(ap,fp);
    if (t != conv<long>(tp)) {
      cerr << "**** SPX testTrace (1) FAILED!\n";
      return false;
    }

    t = TraceMod(a, F);
    if (t != conv<long>(tp)) {
      cerr << "**** SPX testTrace (2) FAILED!\n";
      return false;
    }

    Vec<long> tv = TraceVec(f);
    Vec<zz_p> tvp= TraceVec(fp);
    if (tv != conv<vec_long>(tvp)) {
      cerr << "**** SPX testTrace (3) FAILED!\n";
      return false;
    }
  }
  return true;
}

int main()
{
  for (long p=11; p>1; p =(p-1)/2) { // test p=11, p=5, and p=2
    SPmodulus mod(p);
    SPX a(mod);
    getContext(mod).restore();

    cout << "\np="<<p<<endl;
    if (testGCD(a))   cout << "  testGCD\tPASS\n" << std::flush;
    if (testArith(a)) cout << "  testArith\tPASS\n" << std::flush;
    if (testCoeffs(a)) cout<< "  testCoeffs\tPASS\n" << std::flush;
    if (testMisc(a))  cout << "  testMisc\tPASS\n" << std::flush;
    if (testTrunc(a)) cout << "  testTrunc\tPASS\n" << std::flush;
    if (testMod(a))   cout << "  testMod\tPASS\n" << std::flush;
    if (testModPre(a)) cout<< "  testModPre\tPASS\n" << std::flush;
    if (testComp(a))  cout << "  testComp\tPASS\n" << std::flush;
    if (testProj(a))  cout << "  testProj\tPASS\n" << std::flush;
    //    if (testMinp(a))  cout << "  testMinp\tPASS\n" << std::flush;
    if (testTrace(a)) cout << "  testTrace\tPASS\n" << std::flush;
  }
}
