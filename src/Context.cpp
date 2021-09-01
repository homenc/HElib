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
#include <cstring>
#include <algorithm>
#include <optional>

#include <json.hpp>
using json = ::nlohmann::json;

#include <helib/Context.h>
#include <helib/EvalMap.h>
#include <helib/powerful.h>
#include <helib/sample.h>
#include <helib/EncryptedArray.h>
#include <helib/PolyModRing.h>
#include <helib/fhe_stats.h>

#include "macro.h"
#include "PrimeGenerator.h"
#include "binio.h"
#include "io.h"

namespace helib {

double lweEstimateSecurity(int n, double log2AlphaInv, int hwt)
{
  if (hwt < 0 || (hwt > 0 && hwt < MIN_SK_HWT)) {
    return 0;
  }

  // clang-format off
  constexpr double hwgts[] =
      {120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450};
  constexpr double slopes[] =
      {2.4, 2.67, 2.83, 3.0, 3.1, 3.3, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55};
  constexpr double cnstrms[] =
      {19, 13, 10, 6, 3, 1, -3, -4, -5, -7, -10, -12};
  // clang-format on

  constexpr size_t numWghts = sizeof(hwgts) / sizeof(hwgts[0]);

  const size_t idx = (hwt - 120) / 30; // index into the array above
  double slope = 0, consterm = 0;
  if (hwt == 0) { // dense keys
    slope = 3.8;
    consterm = -20;
  } else if (idx < numWghts - 1) {
    // estimate prms on a line from prms[i] to prms[i+1]
    // how far into this interval
    double a = double(hwt - hwgts[idx]) / (hwgts[idx + 1] - hwgts[idx]);
    slope = slopes[idx] + a * (slopes[idx + 1] - slopes[idx]);
    consterm = cnstrms[idx] + a * (cnstrms[idx + 1] - cnstrms[idx]);
  } else {
    // Use the params corresponding to largest weight (450 above)
    slope = slopes[numWghts - 1];
    consterm = cnstrms[numWghts - 1];
  }

  double x = n / log2AlphaInv;
  double ret = slope * x + consterm;

  return ret < 0.0 ? 0.0 : ret; // If ret is negative then return 0.0
}

// Useful params objects (POD) to simplify calls between ContextBuilder and
// Context.
struct Context::ModChainParams
{
  long bits;
  long c;
  bool bootstrappableFlag;
  long skHwt;
  long resolution;
  long bitsInSpecialPrimes;
  double stdev;
  double scale;
};

struct Context::BootStrapParams
{
  NTL::Vec<long> mvec;
  bool buildCacheFlag;
  bool thickFlag;
};

struct Context::SerializableContent
{
  long p;
  long r;
  long m;
  std::vector<long> gens;
  std::vector<long> ords;
  NTL::xdouble stdev;
  double scale;
  IndexSet smallPrimes;
  IndexSet specialPrimes;
  std::vector<long> qs;
  std::vector<IndexSet> digits;
  long hwt_param;
  long e_param;
  long ePrime_param;
  NTL::Vec<long> mvec;
  bool build_cache;
  bool alsoThick;
};

long FindM(long k,
           long nBits,
           long c,
           long p,
           long d,
           long s,
           long chosen_m,
           bool verbose)
{
  // get a lower-bound on the parameter N=phi(m):
  // 1. Each level in the modulus chain corresponds to pSize=p2Size/2
  //    bits (where we have one prime of this size, and all the others are of
  //    size p2Size).
  //    When using DoubleCRT, we need 2m to divide q-1 for every prime q.
  // 2. With nBits of ctxt primes,
  //    the largest modulus for "fresh ciphertexts" has size
  //          Q0 ~ 2^{nBits}
  // 3. We break each ciphertext into upto c digits, do each digit is as large
  //    as    D=2^{nBits/c}
  // 4. The added noise variance term from the key-switching operation is
  //    c*N*sigma^2*D^2, and this must be mod-switched down to w*N (so it is
  //    on par with the added noise from modulus-switching). Hence the ratio
  //    P that we use for mod-switching must satisfy c*N*sigma^2*D^2/P^2<w*N,
  //    or    P > sqrt(c/w) * sigma * 2^{(L+1)*pSize/c}
  // 5. With this extra P factor, the key-switching matrices are defined
  //    relative to a modulus of size
  //          Q0 = q0*P ~ sqrt{c/w} sigma 2^{nBits*(1+1/c)}
  // 6. To get k-bit security we need N>log(Q0/sigma)(k+110)/7.2, i.e. roughly
  //          N > nBits*(1+1/c)(k+110) / 7.2

  // Compute a bound on m, and make sure that it is not too large
  double cc = 1.0 + (1.0 / static_cast<double>(c));

  double dN = ceil(nBits * cc * (k + 110) / 7.2);
  // FIXME: the bound for dN is not conservative enough...
  // this should be re-worked.

  long N = NTL_SP_BOUND;
  if (N > dN)
    N = dN;
  else {
    std::stringstream ss;
    ss << "Cannot support a bound of " << dN;
    throw RuntimeError(ss.str());
  }

  long m = 0;
  size_t i = 0;

  // find the first m satisfying phi(m)>=N and d | ord(p) in Z_m^*
  // and phi(m)/ord(p) >= s
  if (chosen_m) {
    if (NTL::GCD(p, chosen_m) == 1) {
      long ordP = multOrd(p, chosen_m);
      if (d == 0 || ordP % d == 0) {
        // chosen_m is OK
        m = chosen_m;
      }
    }
  } else if (p == 2) { // use pre-computed table, divisors of 2^n-1 for some n's

    // clang-format off
    static long ms[][4] = {
        // pre-computed values of [phi(m),m,d]
        // phi(m),     m, ord(2),c_m*1000 (not used anymore)
          {  1176,  1247,     28,  3736}, // gens=5(42)
          {  2880,  3133,     24,  3254}, // gens=6(60), 7(!2)
          {  4050,  4051,     50,     0}, // gens=130(81)
          {  4096,  4369,     16,  3422}, // gens=129(16),3(!16)
          {  4704,  4859,     28,     0}, // gens=7(42),3(!4)
          {  5292,  5461,     14,  4160}, // gens=3(126),509(3)
          {  5760,  8435,     24,  8935}, // gens=58(60),1686(2),11(!2)
          {  7500,  7781,     50,     0}, // gens=353(30),3(!5)
          {  8190,  8191,     13,  1273}, // gens=39(630)
          {  9900, 10261,     30,     0}, // gens=3(330)
          { 10752, 11441,     48,  3607}, // gens=7(112),5(!2)
          { 10800, 11023,     45,     0}, // gens=270(24),2264(2),3(!5)
          { 12000, 13981,     20,  2467}, // gens=10(30),23(10),3(!2)
          { 11520, 15665,     24, 14916}, // gens=6(60),177(4),7(!2)
          { 14112, 14351,     18,     0}, // gens=7(126),3(!4)
          { 15004, 15709,     22,  3867}, // gens=5(682)
          { 18000, 18631,     25,  4208}, // gens=17(120),1177(6)
          { 15360, 20485,     24, 12767}, // gens=6(80),242(4),7(2)
          { 16384, 21845,     16, 12798}, // gens=129(16),273(4),3(!16)
          { 17280, 21931,     24, 18387}, // gens=6(60),467(6),11(!2)
          { 19200, 21607,     40, 35633}, // gens=13(120),2789(2),3(!2)
          { 21168, 27305,     28, 15407}, // gens=6(126),781(6)
          { 23040, 23377,     48,  5292}, // gens=35(240),5(!2)
          { 23310, 23311,     45,     0}, // gens=489(518)
          { 24576, 24929,     48,  5612}, // gens=12(256),5(!2)
          { 27000, 32767,     15, 20021}, // gens=3(150),873(6),6945(2)
          { 31104, 31609,     72,  5149}, // gens=11(216),5(!2)
          { 43690, 43691,     34,     0}, // gens=69(1285)
          { 49500, 49981,     30,     0}, // gens=3(1650)
          { 46080, 53261,     24, 33409}, // gens=3(240),242(4),7(!2)
          { 54000, 55831,     25,     0}, // gens=22(360),3529(6)
          { 49140, 57337,     39,  2608}, // gens=39(630),40956(2)
          { 51840, 59527,     72, 21128}, // gens=58(60),1912(6),7(!2)
          { 61680, 61681,     40,  1273}, // gens=33(771),17(!2)
          { 65536, 65537,     32,  1273}, // gens=2(32),3(!2048)
          { 75264, 82603,     56, 36484}, // gens=3(336),24294(2),7(!2)
          { 84672, 92837,     56, 38520}  // gens=18(126),1886(6),3(!2)
    };
    // clang-format on
    for (i = 0; i < sizeof(ms) / sizeof(long[4]); i++) {
      if (ms[i][0] < N || NTL::GCD(p, ms[i][1]) != 1)
        continue;
      long ordP = multOrd(p, ms[i][1]);
      long nSlots = ms[i][0] / ordP;
      if (d != 0 && ordP % d != 0)
        continue;
      if (nSlots < s)
        continue;

      m = ms[i][1];
      break;
    }
  }

  // If m is not set yet, just set it close to N. This may be a lousy
  // choice of m for this p, since you will get a small number of slots.

  if (m == 0) {
    // search only for odd values of m, to make phi(m) a little closer to m
    for (long candidate = N | 1; candidate < 10 * N; candidate += 2) {
      if (NTL::GCD(p, candidate) != 1)
        continue;

      long ordP = multOrd(p, candidate); // the multiplicative order of p mod m
      if (d > 1 && ordP % d != 0)
        continue;
      if (ordP > 100)
        continue; // order too big, we will get very few slots

      long n = phi_N(candidate); // compute phi(m)
      if (n < N)
        continue; // phi(m) too small

      m = candidate; // all tests passed, return this value of m
      break;
    }
  }

  if (verbose) {
    std::cerr << "*** Bound N=" << N << ", choosing m=" << m
              << ", phi(m)=" << phi_N(m) << std::endl;
  }

  return m;
}

// A global variable, pointing to the "current" context
Context* activeContext = nullptr;

void Context::productOfPrimes(NTL::ZZ& p, const IndexSet& s) const
{
  p = 1;
  for (long i : s)
    p *= ithPrime(i);
}

bool Context::operator==(const Context& other) const
{
  if (&other == this)
    return true;

  if (zMStar != other.zMStar)
    return false;
  if (alMod != other.alMod)
    return false;

  if (moduli.size() != other.moduli.size())
    return false;
  for (size_t i = 0; i < moduli.size(); i++) {
    const Cmodulus& m1 = moduli[i];
    const Cmodulus& m2 = other.moduli[i];
    if (m1.getQ() != m2.getQ())
      return false;
  }

  if (smallPrimes != other.smallPrimes)
    return false;
  if (ctxtPrimes != other.ctxtPrimes)
    return false;
  if (specialPrimes != other.specialPrimes)
    return false;

  if (digits.size() != other.digits.size())
    return false;
  for (size_t i = 0; i < digits.size(); i++)
    if (digits[i] != other.digits[i])
      return false;

  if (stdev != other.getStdev())
    return false;

  if (scale != other.scale)
    return false;

  if (hwt_param != other.hwt_param)
    return false;

  if (rcData != other.rcData)
    return false;

  return true;
}

void Context::writeTo(std::ostream& str) const
{
  SerializeHeader<Context>().writeTo(str);

  writeEyeCatcher(str, EyeCatcher::CONTEXT_BEGIN);

  write_raw_int(str, this->zMStar.getP());
  write_raw_int(str, this->alMod.getR());
  write_raw_int(str, this->zMStar.getM());

  write_raw_int(str, this->zMStar.numOfGens());

  // There aren't simple getters to get the gens and ords vectors
  for (long i = 0; i < this->zMStar.numOfGens(); i++) {
    write_raw_int(str, this->zMStar.ZmStarGen(i));
  }

  write_raw_int(str, this->zMStar.numOfGens());

  // Copying the way it is done in ASCII IO.
  // Bad dimensions are represented as a negated ord
  for (long i = 0; i < this->zMStar.numOfGens(); i++) {
    if (this->zMStar.SameOrd(i))
      write_raw_int(str, this->zMStar.OrderOf(i));
    else
      write_raw_int(str, -this->zMStar.OrderOf(i));
  }

  // standard-deviation
  write_raw_xdouble(str, this->stdev);

  // scale
  write_raw_double(str, this->scale);

  // the "small" index
  this->smallPrimes.writeTo(str);

  // the "special" index
  this->specialPrimes.writeTo(str);

  // output the primes in the chain
  write_raw_int(str, this->moduli.size());

  for (const auto& modulo : this->moduli) {
    write_raw_int(str, modulo.getQ());
  }

  // output the digits
  write_raw_int(str, this->digits.size());

  for (const auto& digit : this->digits) {
    digit.writeTo(str);
  }

  write_raw_int(str, this->hwt_param);
  write_raw_int(str, this->e_param);
  write_raw_int(str, this->ePrime_param);

  write_ntl_vec_long(str, this->rcData.mvec);
  write_raw_int(str, static_cast<long>(this->rcData.build_cache));
  write_raw_int(str, static_cast<long>(this->rcData.alsoThick));

  writeEyeCatcher(str, EyeCatcher::CONTEXT_END);
}

Context::SerializableContent Context::readParamsFrom(std::istream& str)
{
  const auto header = SerializeHeader<Context>::readFrom(str);
  assertEq<IOError>(header.version,
                    Binio::VERSION_0_0_1_0,
                    "Header: version " + header.versionString() +
                        " not supported");

  bool eyeCatcherFound = readEyeCatcher(str, EyeCatcher::CONTEXT_BEGIN);
  assertTrue<IOError>(eyeCatcherFound,
                      "Could not find pre-context eye catcher");

  Context::SerializableContent context_params;
  context_params.p = read_raw_int(str);
  context_params.r = read_raw_int(str);
  context_params.m = read_raw_int(str);

  // Number of gens and ords saved in front of vectors
  read_raw_vector(str, context_params.gens);
  read_raw_vector(str, context_params.ords);

  // Get the standard deviation
  context_params.stdev = read_raw_xdouble(str);

  // Get the scale
  context_params.scale = read_raw_double(str);

  context_params.smallPrimes = IndexSet::readFrom(str);
  context_params.specialPrimes = IndexSet::readFrom(str);

  read_raw_vector<long>(str, context_params.qs);

  long nDigits = read_raw_int(str);
  context_params.digits.reserve(nDigits);

  for (long i = 0; i < (long)nDigits; i++) {
    context_params.digits.emplace_back(IndexSet::readFrom(str));
  }
  context_params.hwt_param = read_raw_int(str);
  context_params.e_param = read_raw_int(str);
  context_params.ePrime_param = read_raw_int(str);

  // Read in the partition of m into co-prime factors (if bootstrappable)
  read_ntl_vec_long(str, context_params.mvec);

  context_params.build_cache = read_raw_int(str);
  context_params.alsoThick = read_raw_int(str);

  eyeCatcherFound = readEyeCatcher(str, EyeCatcher::CONTEXT_END);
  assertTrue<IOError>(eyeCatcherFound,
                      "Could not find post-context eye catcher");

  return context_params;
}

Context Context::readFrom(std::istream& str)
{
  return Context(readParamsFrom(str));
}

Context* Context::readPtrFrom(std::istream& str)
{
  return new Context(readParamsFrom(str));
}

Context::SerializableContent Context::readParamsFromJSON(
    const JsonWrapper& jwrap)
{
  auto body = [&]() {
    const json j = fromTypedJson<Context>(unwrap(jwrap));

    Context::SerializableContent content;

    // This way stops ordering inconsistencies
    content.m = j.at("m");
    content.p = j.at("p");
    content.r = j.at("r");
    content.gens = j.at("gens").get<std::vector<long>>();
    content.ords = j.at("ords").get<std::vector<long>>();
    content.stdev = j.at("stdev").get<NTL::xdouble>();
    content.scale = j.at("scale");
    content.smallPrimes = IndexSet::readFromJSON(wrap(j.at("smallPrimes")));
    content.specialPrimes = IndexSet::readFromJSON(wrap(j.at("specialPrimes")));
    content.qs = j.at("qs").get<std::vector<long>>();
    content.digits = readVectorFromJSON<IndexSet>(j.at("digits"));
    content.hwt_param = j.at("hwt_param");
    content.e_param = j.at("e_param");
    content.ePrime_param = j.at("ePrime_param");
    content.mvec = j.at("mvec");
    content.build_cache = j.at("build_cache");
    content.alsoThick = j.at("alsoThick");

    return content;
  };

  return executeRedirectJsonError<Context::SerializableContent>(body);
}

Context Context::readFromJSON(std::istream& is)
{
  return executeRedirectJsonError<Context>([&]() {
    json j;
    is >> j;
    return Context::readFromJSON(wrap(j));
  });
}

Context Context::readFromJSON(const JsonWrapper& jw)
{
  return Context(readParamsFromJSON(jw));
}

Context* Context::readPtrFromJSON(std::istream& is)
{
  return executeRedirectJsonError<Context*>([&]() {
    json j;
    is >> j;
    return new Context(readParamsFromJSON(wrap(j)));
  });
}

JsonWrapper Context::writeToJSON() const
{
  std::function<JsonWrapper()> body = [this]() {
    std::vector<long> gens(this->zMStar.numOfGens());
    // There aren't simple getters to get the gens and ords vectors
    for (long i = 0; i < this->zMStar.numOfGens(); i++) {
      gens[i] = this->zMStar.ZmStarGen(i);
    }

    std::vector<long> ords(this->zMStar.numOfGens());
    // Bad dimensions are represented as a negated ord
    for (long i = 0; i < this->zMStar.numOfGens(); i++) {
      if (this->zMStar.SameOrd(i))
        ords[i] = this->zMStar.OrderOf(i);
      else
        ords[i] = -this->zMStar.OrderOf(i);
    }

    // output the primes in the chain
    std::vector<long> qs;
    qs.reserve(this->moduli.size());
    for (const auto& modulo : this->moduli) {
      qs.emplace_back(modulo.getQ());
    }

    // m
    // p
    // r
    // gens
    // ords
    // stdev
    // scale
    json j = {{"m", this->zMStar.getM()},
              {"p", this->zMStar.getP()},
              {"r", this->alMod.getR()},
              {"gens", gens},
              {"ords", ords},
              {"stdev", this->stdev},
              {"scale", this->scale},
              {"smallPrimes", unwrap(this->smallPrimes.writeToJSON())},
              {"specialPrimes", unwrap(this->specialPrimes.writeToJSON())},
              {"qs", qs},
              {"digits", writeVectorToJSON(this->digits)},
              {"hwt_param", this->hwt_param},
              {"e_param", this->e_param},
              {"ePrime_param", this->ePrime_param},
              {"mvec", this->rcData.mvec},
              {"build_cache", this->rcData.build_cache},
              {"alsoThick", this->rcData.alsoThick}};
    return wrap(toTypedJson<Context>(j));
  };

  return executeRedirectJsonError<JsonWrapper>(body);
} // namespace helib

void Context::writeToJSON(std::ostream& os) const
{
  // We need to wrap as this->writeToJSON() returns a JsonWrapper, so os << js
  // may throw a json-related exception
  return executeRedirectJsonError<void>(
      [&, this]() { os << this->writeToJSON(); });
}

std::ostream& operator<<(std::ostream& os, const Context& context)
{
  context.writeToJSON(os);
  return os;
}

NTL::ZZX getG(const EncryptedArray& ea)
{
  NTL::ZZX G;
  switch (ea.getTag()) {
  case PA_GF2_tag:
    G = NTL::conv<NTL::ZZX>(ea.getDerived(PA_GF2()).getG());
    break;
  case PA_zz_p_tag:
    convert(G, ea.getDerived(PA_zz_p()).getG());
    break;
  case PA_cx_tag:
    throw LogicError("Cannot get polynomial modulus G when scheme is CKKS");
    break;
  default:
    throw LogicError("No valid tag found in EncryptedArray");
  }
  return G;
}

// Constructors must ensure that alMod points to zMStar, and
// rcEA (if set) points to rcAlmod which points to zMStar
Context::Context(unsigned long m,
                 unsigned long p,
                 unsigned long r,
                 const std::vector<long>& gens,
                 const std::vector<long>& ords) :
    zMStar(m, p, gens, ords),
    alMod(zMStar, r),

    // VJS-FIXME: I'm not sure this makes sense.
    // This constrictor was provided mainly for bootstrapping.
    // Most BGV applications will *not* use fully packed slots
    ea(std::make_shared<EncryptedArray>(*this, alMod)),

    // ea(std::make_shared<EncryptedArray>(*this)),
    // This constructor uses the polynomial X, which
    // corresponds to thinly packed slots for BGV.
    // Unfortunately, this doesn't work. The problem is that
    // the bootstrapping rountines, even thin bootstrapping,
    // require G=F0 here.  What we really need to do is
    // put an EA wit G=F0 in rcData, separate from the
    // default EA.  They could, in fact, be the same if default
    // EA has G=F0 (if we allow the default EA G-value to be set
    // by the user.

    pwfl_converter(nullptr),
    stdev(3.2),
    scale(10.0)
{
  // NOTE: pwfl_converter will be set in buildModChain (or endBuildModChain),
  // after the prime chain has been built, as it depends on the primeChain

  if (!isCKKS()) {
    slotRing =
        std::make_shared<PolyModRing>(zMStar.getP(), alMod.getR(), getG(*ea));
  }
}

void Context::printout(std::ostream& out) const
{
  ea->getPAlgebra().printout(out);
  out << "r = " << alMod.getR() << "\n"
      << "nslots = " << ea->size() << "\n"
      << "hwt = " << hwt_param << "\n"
      << "ctxtPrimes = " << ctxtPrimes << "\n"
      << "specialPrimes = " << specialPrimes << "\n"
      << "number of bits = " << bitSizeOfQ() << "\n\n"
      << "security level = " << securityLevel() << std::endl;
}

Context::Context(long m,
                 long p,
                 long r,
                 const std::vector<long>& gens,
                 const std::vector<long>& ords,
                 const std::optional<Context::ModChainParams>& mparams,
                 const std::optional<Context::BootStrapParams>& bparams) :
    Context(m, p, r, gens, ords)
{
  if (mparams) {
    this->stdev = mparams->stdev;
    this->scale = mparams->scale;

    this->buildModChain(mparams->bits,
                        mparams->c,
                        mparams->bootstrappableFlag,
                        mparams->skHwt,
                        mparams->resolution,
                        mparams->bitsInSpecialPrimes);

    if (mparams->bootstrappableFlag && bparams) {
      this->enableBootStrapping(bparams->mvec,
                                bparams->buildCacheFlag,
                                bparams->thickFlag);
    }
  }
}

Context::Context(const SerializableContent& content) :
    Context(content.m, content.p, content.r, content.gens, content.ords)
{
  this->stdev = content.stdev;
  this->scale = content.scale;
  this->digits = content.digits;
  this->hwt_param = content.hwt_param;
  this->e_param = content.e_param;
  this->ePrime_param = content.ePrime_param;

  for (long i = 0; i < lsize(content.qs); i++) {
    long q = content.qs[i];

    this->moduli.emplace_back(this->zMStar, q, 0);

    // FIXME: Consider serializing all 3 sets and setting them directly.
    if (content.smallPrimes.contains(i))
      this->smallPrimes.insert(i); // small prime
    else if (content.specialPrimes.contains(i))
      this->specialPrimes.insert(i); // special prime
    else
      this->ctxtPrimes.insert(i); // ciphertext prime
  }

  endBuildModChain();

  // Read in the partition of m into co-prime factors (if bootstrappable)
  if (content.mvec.length() > 0) {
    // VJS-FIXME: what about the build_cache and alsoThick params?
    this->enableBootStrapping(content.mvec,
                              content.build_cache,
                              content.alsoThick);
  }
}

static void CheckPrimes(const Context& context,
                        const IndexSet& s,
                        const char* name)
{
  for (long i : s) {
    NTL::zz_pPush push;
    context.ithModulus(i).restoreModulus();
    if (!NTL::zz_p::IsFFTPrime()) {
      Warning(__func__ + std::string(": non-FFT prime in ") + name);
    }
  }
}

// Add small primes to get target resolution
// FIXME: there is some black magic here.
// we need to better document the strategy.
void Context::addSmallPrimes(long resolution, long cpSize)
{
  // cpSize is the size of the ciphertext primes
  // Sanity-checks, cpSize \in [0.9*HELIB_SP_NBITS, HELIB_SP_NBITS]
  assertTrue(cpSize >= 30, "cpSize is too small (minimum is 30)");
  assertInRange(cpSize * 10,
                9l * HELIB_SP_NBITS,
                10l * HELIB_SP_NBITS,
                "cpSize not in [0.9*HELIB_SP_NBITS, HELIB_SP_NBITS]",
                true);

  long m = getM();
  if (m <= 0 || m > (1 << 20)) // sanity checks
    throw RuntimeError("addSmallPrimes: m undefined or larger than 2^20");
  // NOTE: Below we are ensured that 16m*log(m) << NTL_SP_BOUND

  if (resolution < 1 || resolution > 10) // set to default of 3-bit resolution
    resolution = 3;

  std::vector<long> sizes;
  long smallest; // size of the smallest of the smallPrimes
  // We need at least two of this size, maybe three

  if (cpSize >= 54)
    smallest = divc(2 * cpSize, 3);
  else if (cpSize >= 45)
    smallest = divc(7 * cpSize, 10);
  else { // Make the smallest ones at least 22-bit primes
    smallest = divc(11 * cpSize, 15);
    sizes.push_back(smallest); // need three of them
  }
  sizes.push_back(smallest);
  sizes.push_back(smallest);

  // This ensures we can express everything to given resolution.

  // use sizes cpSize-r, cpSize-2r, cpSize-4r,... down to the sizes above
  for (long delta = resolution; cpSize - delta > smallest; delta *= 2)
    sizes.push_back(cpSize - delta);

  // This helps to minimize the number of small primes needed
  // to express any particular resolution.
  // This could be removed...need to experiment.

  // Special cases: add also cpSize-3*resolution,
  // and for resolution=1 also cpSize-11
  if (cpSize - 3 * resolution > smallest)
    sizes.push_back(cpSize - 3 * resolution);
  if (resolution == 1 && cpSize - 11 > smallest)
    sizes.push_back(cpSize - 11);

  std::sort(sizes.begin(), sizes.end()); // order by size

  long last_sz = 0;
  std::unique_ptr<PrimeGenerator> gen;
  for (long sz : sizes) {
    if (sz != last_sz)
      gen.reset(new PrimeGenerator(sz, m));
    long q = gen->next();
    addSmallPrime(q);
    last_sz = sz;
  }
}

void Context::addCtxtPrime(long q)
{
  assertFalse(inChain(q), "Prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back(Cmodulus(zMStar, q, 0));
  ctxtPrimes.insert(i);
}

void Context::addSpecialPrime(long q)
{
  assertFalse(inChain(q), "Special prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back(Cmodulus(zMStar, q, 0));
  specialPrimes.insert(i);
}

// Determine the target size of the ctxtPrimes. The target size is
// set at 2^n, where n is at most HELIB_SP_NBITS and at least
// ceil(0.9*HELIB_SP_NBITS), so that we don't overshoot nBits by too
// much.
// The reason that we do not allow to go below 0.9*HELIB_SP_NBITS is
// that we need some of the smallPrimes to be sufficiently smaller
// than the ctxtPrimes, and still we need these smallPrimes to have
// m'th roots of unity.
static long ctxtPrimeSize(long nBits)
{
  double bit_loss =
      -std::log1p(-1.0 / double(1L << PrimeGenerator::B)) / std::log(2.0);
  // std::cerr << "*** bit_loss=" << bit_loss;

  // How many primes of size HELIB_SP_NBITS it takes to get to nBits
  double maxPsize = HELIB_SP_NBITS - bit_loss;
  // primes of length len are guaranteed to be at least (1-1/2^B)*2^len,

  long nPrimes = long(ceil(nBits / maxPsize));
  // this is sufficiently many primes

  // now we want to trim the size to avoid unnecessary overshooting
  // so we decrease targetSize, while guaranteeing that
  // nPrimes primes of length targetSize multiply out to
  // at least nBits bits.

  long targetSize = HELIB_SP_NBITS;
  while (10 * (targetSize - 1) >= 9 * HELIB_SP_NBITS &&
         (targetSize - 1) >= 30 &&
         ((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    targetSize--;

  if (((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    Warning(__func__ + std::string(": non-optimal targetSize"));

  return targetSize;
}

void Context::addCtxtPrimes(long nBits, long targetSize)
{
  // We add enough primes of size targetSize until their product is
  // at least 2^{nBits}

  // Sanity-checks, targetSize \in [0.9*HELIB_SP_NBITS, HELIB_SP_NBITS]
  assertTrue(targetSize >= 30,
             "Target prime is too small (minimum size is 30)");
  assertInRange(targetSize * 10,
                9l * HELIB_SP_NBITS,
                10l * HELIB_SP_NBITS,
                "targetSize not in [0.9*HELIB_SP_NBITS, HELIB_SP_NBITS]",
                true);
  const PAlgebra& palg = getZMStar();
  long m = palg.getM();

  PrimeGenerator gen(targetSize, m);
  double bitlen = 0; // how many bits we already have
  while (bitlen < nBits - 0.5) {
    long q = gen.next(); // generate the next prime
    addCtxtPrime(q);     // add it to the list
    bitlen += std::log2(q);
  }

  // std::cerr << "*** ctxtPrimes excess: " << (bitlen - nBits) << "\n";
  HELIB_STATS_UPDATE("excess-ctxtPrimes", bitlen - nBits);
}

void Context::addSpecialPrimes(long nDgts,
                               bool willBeBootstrappable,
                               long bitsInSpecialPrimes)
{
  const PAlgebra& palg = getZMStar();
  long p = std::abs(palg.getP()); // for CKKS, palg.getP() == -1
  long m = palg.getM();
  long phim = palg.getPhiM();
  long p2r = isCKKS() ? 1 : getAlMod().getPPowR();

  long p2e = p2r;
  if (willBeBootstrappable && !isCKKS()) {
    // bigger p^e for bootstrapping
    long e, ePrime;
    RecryptData::setAE(e, ePrime, *this);
    p2e *= NTL::power_long(p, e - ePrime);

    // initialize e and ePrime parameters in the context
    this->e_param = e;
    this->ePrime_param = ePrime;
  }

  long nCtxtPrimes = getCtxtPrimes().card();
  if (nDgts > nCtxtPrimes)
    nDgts = nCtxtPrimes; // sanity checks
  if (nDgts <= 0)
    nDgts = 1;

  digits.resize(nDgts); // allocate space

  if (nDgts > 1) { // we break ciphertext into a few digits when key-switching
    // NOTE: The code below assumes that all the ctxtPrimes have roughly the
    // same size

    IndexSet remaining = getCtxtPrimes();
    for (long dgt = 0; dgt < nDgts - 1; dgt++) {
      long digitCard = divc(remaining.card(), nDgts - dgt);
      // ceiling(#-of-remaining-primes, #-or-remaining-digits)

      for (long i : remaining) {
        digits[dgt].insert(i);
        if (digits[dgt].card() >= digitCard)
          break;
      }
      remaining.remove(digits[dgt]); // update the remaining set
    }
    // The last digit has everything else
    if (empty(remaining)) { // sanity check, use one less digit
      nDgts--;
      digits.resize(nDgts);
    } else
      digits[nDgts - 1] = remaining;
  } else { // only one digit
    digits[0] = getCtxtPrimes();
  }

  double maxDigitLog = 0.0;
  for (auto& digit : digits) {
    double size = logOfProduct(digit);
    if (size > maxDigitLog)
      maxDigitLog = size;
  }

  // Add special primes to the chain for the P factor of key-switching
  double nBits;

  if (bitsInSpecialPrimes)
    nBits = bitsInSpecialPrimes;
  else {
#if 0
    nBits = (maxDigitLog + std::log(nDgts) + NTL::log(stdev * 2) +
             std::log(p2e)) /
            std::log(2.0);
    // FIXME: Victor says: the above calculation does not make much sense to me
#else
    double h;
    if (getHwt() == 0)
      h = phim / 2.0;
    else
      h = getHwt();

    double log_phim = std::log(phim);
    if (log_phim < 1)
      log_phim = 1;

    if (isCKKS()) {
      // This is based on a smaller noise estimate so as
      // to better protect precision...this is based on
      // a noise level equal to the mod switch added noise.
      // Note that the relin_CKKS_adjust function in Ctxt.cpp
      // depends on this estimate.
      nBits = (maxDigitLog + NTL::log(getStdev()) + std::log(nDgts) -
               0.5 * std::log(h)) /
              std::log(2.0);
    } else if (palg.getPow2()) {
      nBits = (maxDigitLog + std::log(p2e) + NTL::log(getStdev()) +
               0.5 * std::log(12.0) + std::log(nDgts) -
               0.5 * std::log(log_phim) - 2 * std::log(p) - std::log(h)) /
              std::log(2.0);
    } else {
      nBits =
          (maxDigitLog + std::log(m) + std::log(p2e) + NTL::log(getStdev()) +
           0.5 * std::log(12.0) + std::log(nDgts) - 0.5 * log_phim -
           0.5 * std::log(log_phim) - 2 * std::log(p) - std::log(h)) /
          std::log(2.0);
    }

    // Both of the above over-estimate nBits by a factor of
    // log2(scale). That should provide a sufficient safety margin.
    // See design document

#endif
  }

  if (nBits < 1)
    nBits = 1;

  double bit_loss =
      -std::log1p(-1.0 / double(1L << PrimeGenerator::B)) / std::log(2.0);

  // How many primes of size HELIB_SP_NBITS it takes to get to nBits
  double maxPsize = HELIB_SP_NBITS - bit_loss;
  // primes of length len are guaranteed to be at least (1-1/2^B)*2^len,

  long nPrimes = long(ceil(nBits / maxPsize));
  // this is sufficiently many prime

  // now we want to trim the size to avoid unnecessary overshooting
  // so we decrease targetSize, while guaranteeing that
  // nPrimes primes of length targetSize multiply out to
  // at least nBits bits.

  long targetSize = HELIB_SP_NBITS;
  while ((targetSize - 1) >= 0.55 * HELIB_SP_NBITS && (targetSize - 1) >= 30 &&
         ((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    targetSize--;

  if (((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    Warning(__func__ + std::string(": non-optimal targetSize"));

  PrimeGenerator gen(targetSize, m);

  while (nPrimes > 0) {
    long q = gen.next();

    if (inChain(q))
      continue;
    // nbits could equal NTL_SP_BITS or the size of one
    // of the small primes, so we have to check for duplicates here...
    // this is not the most efficient way to do this,
    // but it doesn't make sense to optimize this any further

    addSpecialPrime(q);
    nPrimes--;
  }

  // std::cerr << "*** specialPrimes excess: " <<
  // (logOfProduct(specialPrimes)/std::log(2.0) - nBits) <<
  // "\n";
  HELIB_STATS_UPDATE("excess-specialPrimes",
                     logOfProduct(getSpecialPrimes()) / std::log(2.0) - nBits);
}

void Context::buildModChain(long nBits,
                            long nDgts,
                            bool willBeBootstrappable,
                            long skHwt,
                            long resolution,
                            long bitsInSpecialPrimes)
{
  // Cannot build modulus chain with nBits < 0
  assertTrue<InvalidArgument>(nBits > 0,
                              "Cannot initialise modulus chain with nBits < 1");

  assertTrue(skHwt >= 0, "invalid skHwt parameter");

  // ignore for CKKS
  if (isCKKS())
    willBeBootstrappable = false;

  if (skHwt == 0) {
    // default skHwt: if bootstrapping, set to BOOT_DFLT_SK_HWT
    if (willBeBootstrappable)
      skHwt = BOOT_DFLT_SK_HWT;
  }

  // initialize hwt param in context
  hwt_param = skHwt;

  long pSize = ctxtPrimeSize(nBits);
  addSmallPrimes(resolution, pSize);
  addCtxtPrimes(nBits, pSize);
  addSpecialPrimes(nDgts, willBeBootstrappable, bitsInSpecialPrimes);

  CheckPrimes(*this, smallPrimes, "smallPrimes");
  CheckPrimes(*this, ctxtPrimes, "ctxtPrimes");
  CheckPrimes(*this, specialPrimes, "specialPrimes");

  endBuildModChain();
}

void Context::addSmallPrime(long q)
{
  assertFalse(inChain(q), "Small prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back(Cmodulus(zMStar, q, 0));
  smallPrimes.insert(i);
}

void Context::endBuildModChain()
{
  setModSizeTable();
  long m = getM();
  std::vector<long> mvec;
  pp_factorize(mvec, m);
  NTL::Vec<long> mmvec;
  convert(mmvec, mvec);
  pwfl_converter = std::make_shared<PowerfulDCRT>(*this, mmvec);
}

// Helper for the build and buildPtr methods
template <typename SCHEME>
const std::pair<std::optional<Context::ModChainParams>,
                std::optional<Context::BootStrapParams>>
ContextBuilder<SCHEME>::makeParamsArgs() const
{
  const auto mparams =
      buildModChainFlag_
          ? std::make_optional<Context::ModChainParams>({bits_,
                                                         c_,
                                                         bootstrappableFlag_,
                                                         skHwt_,
                                                         resolution_,
                                                         bitsInSpecialPrimes_,
                                                         stdev_,
                                                         scale_})
          : std::nullopt;

  const auto bparams = bootstrappableFlag_
                           ? std::make_optional<Context::BootStrapParams>(
                                 {mvec_, buildCacheFlag_, thickFlag_})
                           : std::nullopt;

  return {mparams, bparams};
}

template <typename SCHEME>
Context ContextBuilder<SCHEME>::build() const
{
  auto [mparams, bparams] = makeParamsArgs();
  return Context(m_, p_, r_, gens_, ords_, mparams, bparams);
}

template <typename SCHEME>
Context* ContextBuilder<SCHEME>::buildPtr() const
{
  auto [mparams, bparams] = makeParamsArgs();
  return new Context(m_, p_, r_, gens_, ords_, mparams, bparams);
}

// Essentially serialization of params.
template <>
std::ostream& operator<<<BGV>(std::ostream& os, const ContextBuilder<BGV>& cb)
{
  const json j = {{"scheme", "bgv"},
                  {"m", cb.m_},
                  {"p", cb.p_},
                  {"r", cb.r_},
                  {"c", cb.c_},
                  {"gens", cb.gens_},
                  {"ords", cb.ords_},
                  {"buildModChainFlag", cb.buildModChainFlag_},
                  {"bits", cb.bits_},
                  {"skHwt", cb.skHwt_},
                  {"resolution", cb.resolution_},
                  {"bitsInSpecialPrimes", cb.bitsInSpecialPrimes_},
                  {"bootstrappableFlag", cb.bootstrappableFlag_},
                  {"mvec", cb.mvec_},
                  {"buildCacheFlag", cb.buildCacheFlag_},
                  {"thickFlag", cb.thickFlag_}};
  os << toTypedJson<ContextBuilder<BGV>>(j);
  return os;
}

template <>
std::ostream& operator<<<CKKS>(std::ostream& os, const ContextBuilder<CKKS>& cb)
{
  const json j = {{"scheme", "ckks"},
                  {"m", cb.m_},
                  {"precision", cb.r_},
                  {"c", cb.c_},
                  {"gens", cb.gens_},
                  {"ords", cb.ords_},
                  {"buildModChainFlag", cb.buildModChainFlag_},
                  {"bits", cb.bits_},
                  {"skHwt", cb.skHwt_},
                  {"resolution", cb.resolution_},
                  {"bitsInSpecialPrimes", cb.bitsInSpecialPrimes_}};
  os << toTypedJson<ContextBuilder<CKKS>>(j);
  return os;
}

template class ContextBuilder<BGV>;
template class ContextBuilder<CKKS>;

} // namespace helib
