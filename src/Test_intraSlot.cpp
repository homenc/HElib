/* Copyright (C) 2012-2017 IBM Corp.
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
#include "intraSlot.h"

int main(int argc, char **argv)
{
  ArgMapping amap;
  long p=2;
  amap.arg("p", p, "plaintext base");
  long n=2;
  amap.arg("n", n, "number of packed ciphertexts");
  long r=1;
  amap.arg("r", r,  "lifting");
  long L=10;
  amap.arg("L", L, "# of levels in the modulus chain");
  long m=91;
  amap.arg("m", m, "use specified value as modulus", NULL);
  long seed=0;
  amap.arg("seed", seed, "PRG seed");
  amap.parse(argc, argv);

  SetSeed(ZZ(seed));

  FHEcontext context(m, p, r);
  context.zMStar.printout();
  buildModChain(context, L, 3);

  ZZX G = context.alMod.getFactorsOverZZ()[0];
  EncryptedArray ea(context, G);

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(100); // A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);

  long d = ea.getDegree(); // size of each slot
  cout << "packing/unpaking "<<n<<"<-->"<<(n*d -1)<<" ciphertexts\n";

  vector<Ctxt> unpacked(d*n -1, Ctxt(publicKey));

  // generate (almost) d*n ciphertexts, with only integrs in the slots
  std::vector<NewPlaintextArray> p1(lsize(unpacked), NewPlaintextArray(ea));
  for (long i=0; i<lsize(unpacked); i++) {
    vector<long> slots;
    ea.random(slots);
    encode(ea, p1[i] ,slots);
    //    cout << "p1["<<i<<"]="<<p1[i] << endl;
    ea.encrypt(unpacked[i], publicKey, p1[i]);
  }

  // Pack (almost) d*n ciphetexts into only n of them
  std::vector<Ctxt> ct(n, Ctxt(publicKey));
  repack(CtPtrs_vectorCt(ct), CtPtrs_vectorCt(unpacked), ea);

  // Unpack them back
  vector<zzX> unpackSlotEncoding;
  buildUnpackSlotEncoding(unpackSlotEncoding, ea);
  unpack(CtPtrs_vectorCt(unpacked), CtPtrs_vectorCt(ct), ea, unpackSlotEncoding);

  NewPlaintextArray p2(ea);
  for (long i=0; i<lsize(unpacked); i++) {
    ea.decrypt(unpacked[i], secretKey, p2);

    if (!equals(ea, p1[i], p2)) {
      cout << "BAD, ";
      cout << "p2["<<i<<"]="<<p2 << endl;
      exit(0);
    }
  }
  cout << "Good\n";
}
