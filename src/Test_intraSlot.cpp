/* Copyright (C) 2012-2019 IBM Corp.
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
// testPacking.cxx - testing uppack/repack functionality
#include <helib/intraSlot.h>
#include <helib/ArgMap.h>

NTL_CLIENT
using namespace helib;

bool verbose = false;

int main(int argc, char **argv)
{
  ArgMap amap;
  long p=2;
  amap.arg("p", p, "plaintext base");
  long n=2;
  amap.arg("n", n, "number of packed ciphertexts");
  long r=1;
  amap.arg("r", r,  "lifting");
  long L=10;
  amap.arg("L", L, "# of levels in the modulus chain");
  long m=91;
  amap.arg("m", m, "use specified value as modulus", nullptr);
  long seed=0;
  amap.arg("seed", seed, "PRG seed");
  amap.arg("verbose", verbose, "suppress printouts");
  amap.parse(argc, argv);

  SetSeed(ZZ(seed));

  Context context(m, p, r);
  if (verbose)
    context.zMStar.printout();
  buildModChain(context, L, 3);

  ZZX G = context.alMod.getFactorsOverZZ()[0];
  EncryptedArray ea(context, G);

  SecKey secretKey(context);
  const PubKey& publicKey = secretKey;
  secretKey.GenSecKey(); // A +-1/0 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);

  long d = ea.getDegree(); // size of each slot

  vector<Ctxt> unpacked(d*n -1, Ctxt(publicKey));

  // generate (almost) d*n ciphertexts, with only integrs in the slots
  std::vector<PlaintextArray> p1(lsize(unpacked), PlaintextArray(ea));
  for (long i=0; i<lsize(unpacked); i++) {
    vector<long> slots;
    ea.random(slots);
    encode(ea, p1[i] ,slots);
    ea.encrypt(unpacked[i], publicKey, p1[i]);
  }

  // Pack (almost) d*n ciphetexts into only n of them
  std::vector<Ctxt> ct(n, Ctxt(publicKey));
  repack(CtPtrs_vectorCt(ct), CtPtrs_vectorCt(unpacked), ea);

  // Unpack them back
  vector<zzX> unpackSlotEncoding;
  buildUnpackSlotEncoding(unpackSlotEncoding, ea);
  unpack(CtPtrs_vectorCt(unpacked), CtPtrs_vectorCt(ct), ea, unpackSlotEncoding);

  PlaintextArray p2(ea);
  for (long i=0; i<lsize(unpacked); i++) {
    ea.decrypt(unpacked[i], secretKey, p2);

    if (!equals(ea, p1[i], p2)) {
      cout << "BAD";
      if (verbose)
        cout << "p2["<<i<<"]="<<p2 << endl;
      exit(0);
    }
  }
  cout << "GOOD\n";
}
