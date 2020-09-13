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
#include <helib/Context.h>
#include <helib/EvalMap.h>
#include <helib/powerful.h>
#include <helib/binio.h>
#include <helib/sample.h>
#include <helib/EncryptedArray.h>
#include <helib/PolyModRing.h>

namespace helib {

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
  double cc = 1.0 + (1.0 / (double)c);

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

  if (stdev != other.stdev)
    return false;

  if (scale != other.scale)
    return false;

  if (rcData != other.rcData)
    return false;
  return true;
}

void writeContextBaseBinary(std::ostream& str, const Context& context)
{
  writeEyeCatcher(str, BINIO_EYE_CONTEXTBASE_BEGIN);

  write_raw_int(str, context.zMStar.getP());
  write_raw_int(str, context.alMod.getR());
  write_raw_int(str, context.zMStar.getM());

  write_raw_int(str, context.zMStar.numOfGens());

  // There aren't simple getters to get the gens and ords vectors
  for (long i = 0; i < context.zMStar.numOfGens(); i++) {
    write_raw_int(str, context.zMStar.ZmStarGen(i));
  }

  write_raw_int(str, context.zMStar.numOfGens());

  // Copying the way it is done in ASCII IO.
  // Bad dimensions are represented as a negated ord
  for (long i = 0; i < context.zMStar.numOfGens(); i++) {
    if (context.zMStar.SameOrd(i))
      write_raw_int(str, context.zMStar.OrderOf(i));
    else
      write_raw_int(str, -context.zMStar.OrderOf(i));
  }

  writeEyeCatcher(str, BINIO_EYE_CONTEXTBASE_END);
}

void readContextBaseBinary(std::istream& str,
                           unsigned long& m,
                           unsigned long& p,
                           unsigned long& r,
                           std::vector<long>& gens,
                           std::vector<long>& ords)
{
  int eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_CONTEXTBASE_BEGIN);
  assertEq(eyeCatcherFound, 0, "Could not find pre-context-base eye catcher");

  p = read_raw_int(str);
  r = read_raw_int(str);
  m = read_raw_int(str);

  // Number of gens and ords saved in front of vectors
  read_raw_vector(str, gens);
  read_raw_vector(str, ords);

  eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_CONTEXTBASE_END);
  assertEq(eyeCatcherFound, 0, "Could not find post-context-base eye catcher");
}

std::unique_ptr<Context> buildContextFromBinary(std::istream& str)
{
  unsigned long m, p, r;
  std::vector<long> gens, ords;
  readContextBaseBinary(str, m, p, r, gens, ords);
  return std::unique_ptr<Context>(new Context(m, p, r, gens, ords));
}

void writeContextBinary(std::ostream& str, const Context& context)
{

  writeEyeCatcher(str, BINIO_EYE_CONTEXT_BEGIN);

  // standard-deviation
  write_raw_xdouble(str, context.stdev);

  // scale
  write_raw_double(str, context.scale);

  // the "small" index
  write_raw_int(str, context.smallPrimes.card());
  for (long tmp : context.smallPrimes) {
    ;
    write_raw_int(str, tmp);
  }

  // the "special" index
  write_raw_int(str, context.specialPrimes.card());
  for (long tmp : context.specialPrimes) {
    write_raw_int(str, tmp);
  }

  // output the primes in the chain
  write_raw_int(str, context.moduli.size());

  for (long i = 0; i < (long)context.moduli.size(); i++) {
    write_raw_int(str, context.moduli[i].getQ());
  }

  // output the digits
  write_raw_int(str, context.digits.size());

  for (long i = 0; i < (long)context.digits.size(); i++) {
    write_raw_int(str, context.digits[i].card());
    for (long tmp : context.digits[i]) {
      write_raw_int(str, tmp);
    }
  }

  write_ntl_vec_long(str, context.rcData.mvec);

  write_raw_int(str, context.rcData.skHwt);

  writeEyeCatcher(str, BINIO_EYE_CONTEXT_END);
}

void readContextBinary(std::istream& str, Context& context)
{
  int eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_CONTEXT_BEGIN);
  assertEq(eyeCatcherFound, 0, "Could not find pre-context eye catcher");

  // Get the standard deviation
  context.stdev = read_raw_xdouble(str);

  // Get the scale
  context.scale = read_raw_double(str);

  IndexSet smallPrimes;
  long smallPrimes_sz = read_raw_int(str);
  for (long tmp, i = 0; i < smallPrimes_sz; i++) {
    tmp = read_raw_int(str);
    smallPrimes.insert(tmp);
  }

  IndexSet specialPrimes;
  long specialPrimes_sz = read_raw_int(str);
  for (long tmp, i = 0; i < specialPrimes_sz; i++) {
    tmp = read_raw_int(str);
    specialPrimes.insert(tmp);
  }

  context.moduli.clear();
  context.smallPrimes.clear();
  context.specialPrimes.clear();
  context.ctxtPrimes.clear();

  long nPrimes = read_raw_int(str);

  for (long p, i = 0; i < nPrimes; i++) {
    p = read_raw_int(str);

    context.moduli.push_back(Cmodulus(context.zMStar, p, 0));

    if (smallPrimes.contains(i))
      context.smallPrimes.insert(i); // small prime
    else if (specialPrimes.contains(i))
      context.specialPrimes.insert(i); // special prime
    else
      context.ctxtPrimes.insert(i); // ciphertext prime
  }

  long nDigits = read_raw_int(str);

  context.digits.resize(nDigits);
  for (long i = 0; i < (long)context.digits.size(); i++) {
    long sizeOfS = read_raw_int(str);

    for (long tmp, n = 0; n < sizeOfS; n++) {
      tmp = read_raw_int(str);
      context.digits[i].insert(tmp);
    }
  }

  endBuildModChain(context);

  // Read in the partition of m into co-prime factors (if bootstrappable)
  NTL::Vec<long> mv;
  read_ntl_vec_long(str, mv);

  long t = read_raw_int(str);

  if (mv.length() > 0) {
    context.makeBootstrappable(mv, t);
  }

  eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_CONTEXT_END);
  assertEq(eyeCatcherFound, 0, "Could not find post-context eye catcher");
}

void writeContextBase(std::ostream& str, const Context& context)
{
  str << "[" << context.zMStar.getM() << " " << context.zMStar.getP() << " "
      << context.alMod.getR() << " [";
  for (long i = 0; i < (long)context.zMStar.numOfGens(); i++) {
    str << context.zMStar.ZmStarGen(i)
        << ((i == (long)context.zMStar.numOfGens() - 1) ? "]" : " ");
  }
  str << " [";
  for (long i = 0; i < (long)context.zMStar.numOfGens(); i++) {
    long ord = context.zMStar.OrderOf(i);
    if (context.zMStar.SameOrd(i))
      str << ord;
    else
      str << (-ord);
    if (i < (long)context.zMStar.numOfGens() - 1)
      str << ' ';
  }
  str << "]]";
}

std::ostream& operator<<(std::ostream& str, const Context& context)
{
  str << "[\n";

  // standard-deviation
  str << context.stdev << "\n";

  // scale
  str << context.scale << "\n";

  // the "small" index
  str << context.smallPrimes << "\n ";

  // the "special" index
  str << context.specialPrimes << "\n ";

  // output the primes in the chain
  str << context.moduli.size() << "\n";
  for (long i = 0; i < (long)context.moduli.size(); i++)
    str << context.moduli[i].getQ() << " ";
  str << "\n ";

  // output the digits
  str << context.digits.size() << "\n";
  for (long i = 0; i < (long)context.digits.size(); i++)
    str << context.digits[i] << " ";

  str << "\n";

  str << context.rcData.mvec;
  str << " " << context.rcData.skHwt;
  str << " " << context.rcData.build_cache;

  str << "]\n";

  return str;
}

void readContextBase(std::istream& str,
                     unsigned long& m,
                     unsigned long& p,
                     unsigned long& r,
                     std::vector<long>& gens,
                     std::vector<long>& ords)
{
  // Advance str beyond first '['
  seekPastChar(str, '['); // this function is defined in NumbTh.cpp

  str >> m >> p >> r;
  str >> gens;
  str >> ords;

  seekPastChar(str, ']');
}
std::unique_ptr<Context> buildContextFromAscii(std::istream& str)
{
  unsigned long m, p, r;
  std::vector<long> gens, ords;
  readContextBase(str, m, p, r, gens, ords);
  return std::unique_ptr<Context>(new Context(m, p, r, gens, ords));
}

std::istream& operator>>(std::istream& str, Context& context)
{
  seekPastChar(str, '['); // this function is defined in NumbTh.cpp

  // Get the standard deviation
  str >> context.stdev;

  // Get the scale
  str >> context.scale;

  IndexSet smallPrimes;
  str >> smallPrimes;

  IndexSet specialPrimes;
  str >> specialPrimes;

  context.moduli.clear();
  context.smallPrimes.clear();
  context.specialPrimes.clear();
  context.ctxtPrimes.clear();

  long nPrimes;
  str >> nPrimes;
  for (long i = 0; i < nPrimes; i++) {
    long p;
    str >> p;

    context.moduli.push_back(Cmodulus(context.zMStar, p, 0));

    if (smallPrimes.contains(i))
      context.smallPrimes.insert(i); // small prime
    else if (specialPrimes.contains(i))
      context.specialPrimes.insert(i); // special prime
    else
      context.ctxtPrimes.insert(i); // ciphertext prime
  }

  // read in the partition to digits
  long nDigits;
  str >> nDigits;
  context.digits.resize(nDigits);
  for (long i = 0; i < (long)context.digits.size(); i++)
    str >> context.digits[i];

  endBuildModChain(context);

  // Read in the partition of m into co-prime factors (if bootstrappable)
  NTL::Vec<long> mv;
  long t;
  int build_cache;
  str >> mv;
  str >> t;
  str >> build_cache;
  if (mv.length() > 0) {
    context.makeBootstrappable(mv, t, build_cache);
  }
  seekPastChar(str, ']');
  return str;
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
    ea(std::make_shared<EncryptedArray>(*this, alMod)),
    pwfl_converter(nullptr),
    stdev(3.2),
    scale(10.0)
{
  // NOTE: pwfl_converter will be set in buildModChain (or endBuildModChain),
  // after the prime chain has been built, as it depends on the primeChain

  if (this->alMod.getTag() != PA_cx_tag) {
    slotRing =
        std::make_shared<PolyModRing>(zMStar.getP(), alMod.getR(), getG(*ea));
  }
}

void Context::printout(std::ostream& out) const
{
  ea->getPAlgebra().printout(out);
  out << "r = " << alMod.getR() << "\n"
      << "nslots = " << ea->size() << "\n"
      << "ctxtPrimes = " << ctxtPrimes << "\n"
      << "specialPrimes = " << specialPrimes << "\n"
      << "number of bits = " << bitSizeOfQ() << "\n\n"
      << "security level = " << securityLevel() << std::endl;
}

} // namespace helib
