#include <helib/helib.h>
#include <helib/debugging.h>

#include <NTL/BasicThreadPool.h>

#include "gtest/gtest.h"
#include "test_common.h"

// only run if HEXL has been linked.
//#ifdef USE_INTEL_HEXL
//#include "../src/intelExt.h" // Private header
//
namespace {

// FIXME Copied from GTestGeneral. Should really have common functionality.
::testing::AssertionResult ciphertextMatches(const helib::EncryptedArray& ea,
                                             const helib::SecKey& sk,
                                             const helib::PtxtArray& p,
                                             const helib::Ctxt& c)
{
  helib::PtxtArray pp(ea);
  pp.decrypt(c, sk);
  if (pp == p) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure()
           << "Ciphertext does not match plaintext:" << std::endl
           << "p = " << p << std::endl
           << "pp = " << pp << std::endl;
  }
}

struct Parameters
{
  Parameters(unsigned m,
             unsigned p,
             unsigned r,
             unsigned bits,
             const std::vector<long>& gens = {},
             const std::vector<long>& ords = {}) :
      m(m), p(p), r(r), bits(bits), gens(gens), ords(ords){};

  const unsigned m;
  const unsigned p;
  const unsigned r;
  const unsigned bits;
  const std::vector<long> gens;
  const std::vector<long> ords;

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << ", "
              << "gens = " << helib::vecToStr(params.gens) << ", "
              << "ords = " << helib::vecToStr(params.ords) << ", "
              << "bits = " << params.bits << "}";
  }
};

class TestHEXL_BGV : public ::testing::TestWithParam<Parameters>
{
protected:
  const unsigned long m;
  const unsigned long p;
  const unsigned long r;
  const unsigned long bits;
  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestHEXL_BGV() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      context(helib::ContextBuilder<helib::BGV>()
                  .m(m)
                  .p(p)
                  .r(r)
                  .bits(bits)
                  .build()),
      secretKey(context),
      publicKey((secretKey.GenSecKey(),
                 helib::addSome1DMatrices(secretKey),
                 secretKey)),
      ea(context.getEA())
  {}

  virtual void SetUp() override
  {
    NTL::SetNumThreads(1);
    if (helib_test::verbose) {
      ea.getPAlgebra().printout();
      std::cout << "r = " << context.getAlMod().getR() << std::endl;
      std::cout << "ctxtPrimes=" << context.getCtxtPrimes()
                << ", specialPrimes=" << context.getSpecialPrimes() << "\n"
                << std::endl;
    }

    helib::setupDebugGlobals(&secretKey, context.shareEA());
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

TEST(TestHEXL, hexlInUse)
{
  long N = 8; // This is phi_m
  long modulus = 769;
  long root = 2;

  std::vector<long> args = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  auto expected_outputs = args;

  intel::FFTFwd(args.data(), args.data(), N, modulus, root);
  intel::FFTRev1(args.data(), args.data(), N, modulus, root);

  EXPECT_TRUE(std::equal(args.begin(), args.end(), expected_outputs.begin()));
}

TEST_P(TestHEXL_BGV, multiplyTwoCtxts)
{
  helib::PtxtArray p0(ea), p1(ea);
  p0.random();
  p1.random();

  helib::Ctxt c0(publicKey), c1(publicKey);
  p0.encrypt(c0); 
  p1.encrypt(c1); 

  p0 *= p1;
  c0 *= c1;

  EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0));
}

TEST_P(TestHEXL_BGV, encryptDecrypt)
{
    NTL::SetNumThreads(1);
  helib::PtxtArray p0(ea);
 
  //std::vector<long> v0 = {1,2,3};
  //std::vector<long> v0 = {1,1,1,1,1,1,1,1};
  std::vector<long> v0 = {0,5};

  p0.load(v0);

  helib::Ctxt c0(publicKey);
  p0.encrypt(c0); 

  EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0));
}

TEST(TestHEXL_BGV, CModulusFFT)
{
  NTL::SetNumThreads(1);
 
  helib::PAlgebra zms(8, 769);
  helib::Cmodulus cmod(zms, 1135873, 0);

  NTL::ZZX poly;
  poly.SetLength(2);
  poly[0] = 0;
  poly[1] = 5;
  std::cout << "TEST: input " << poly << std::endl;

  NTL::vec_long transformed;
  cmod.FFT(transformed, poly);
  std::cout << "TEST: transformed " << transformed << std::endl;

  NTL::zz_pX inverse;
  cmod.iFFT(inverse, transformed);
  std::cout << "TEST: inverse " << inverse << std::endl;

  NTL::ZZX inverse_conv = NTL::conv<NTL::ZZX>(inverse);
  EXPECT_EQ(inverse_conv, poly);
}


// clang-format off
INSTANTIATE_TEST_SUITE_P(typicalParameters, TestHEXL_BGV, ::testing::Values(
    //Parameters(16, 3, 1, 300) // m power of 2.
    Parameters(8, 769, 1, 50) // m power of 2.
    ));
// clang-format on

} // namespace

//#else
//
//namespace {
//
//TEST(TestHEXL, noTestRequired) {}
//
//} // namespace
//
//#endif // USE_INTEL_HEXL
