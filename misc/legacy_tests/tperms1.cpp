/* Copyright (C) 2020 IBM Corp.
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

/* Test_Permutations.cpp - Applying plaintext permutation to encrypted vector
 */
#include <NTL/ZZ.h>
NTL_CLIENT

#include <helib/NumbTh.h>
#include <helib/timing.h>
#include <helib/permutations.h>
#include <helib/EncryptedArray.h>
#include <helib/ArgMap.h>

using namespace helib;

static bool noPrint = true;

void testCtxt(long m, long p, long depthBound, long L, long r)
{
  if (!noPrint)
    std::cout << "@testCtxt(m="<<m<<",p="<<p<<",depth="<<depthBound<< ",r="<<r<<")\n";

  Context context(m,p,r);
  buildModChain(context, /*nLevels=*/L);

  // Generate a sk/pk pair
  SecKey secretKey(context);
  secretKey.GenSecKey(); // A +-1/0 secret key
  const PubKey& publicKey = secretKey;


  long n = context.getNSlots();

  std::cout << "n=" << n << "\n";

  PermIndepPrecomp pip(context, depthBound);

  std::cout << "depth=" << pip.getDepth() << "\n";;
  std::cout << "cost=" << pip.getCost() << "\n";;

  Permut pi;
  randomPerm(pi, n);

  PermPrecomp pp(pip, pi);

  addSome1DMatrices(secretKey);
  // addMatrices4Network(secretKey, net); 
  // I can't remember if this is significantly better than 
  // addSome1DMatrices.  I mean, I *think* we always stuck to the
  // convention of only having KS matrices for powers of generators.
  // I guess I can play around and see...

  Ctxt ctxt(publicKey);
  PtxtArray v(context);
  v.random();


  if (p < 0)
    v.encrypt(ctxt, n); // CKKS encryption
  else
    v.encrypt(ctxt);    // BGV encryption

  pp.apply(ctxt);
  pp.apply(v);

  PtxtArray w(context);
  w.decrypt(ctxt, secretKey);


  if (w == Approx(v))
    std::cout << "GOOD\n";
  else
    std::cout << "BAD\n";
}



int main(int argc, char *argv[])
{
  long p = 2;
  long r = 0;
  long m = 4369;
  long L = 1000;
  long depth = 5;

  noPrint=0;

  ArgMap amap;
  amap.arg("p", p);
  amap.arg("r", r);
  amap.arg("m", m);
  amap.arg("L", L);
  amap.arg("depth", depth);
  amap.arg("noPrint", noPrint);
  amap.parse(argc, argv);

  if (p < 0 && r == 0) r = 20; // CKKS default
  if (p > 0 && r == 0) r = 1;  // BGV default


  testCtxt(m,p,depth,L,r);
}
