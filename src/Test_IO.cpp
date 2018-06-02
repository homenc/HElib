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

/* Test_IO.cpp - Testing the I/O of the important classes of the library
 * (context, keys, ciphertexts).
 */
#include <fstream>
#include <unistd.h>

#include <NTL/ZZX.h>
#include <NTL/vector.h>

#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"

#define N_TESTS 3
static long ms[N_TESTS][10] = {
  //nSlots  m   phi(m) ord(2)
  //  {   2,    7,    6,    3,   0,0,0,0,0,0},
  {   6,   31,   30,    5,   0,0,0,0,0,0},
  {   6,  127,  126,    7,  127,  1,  108,  24,   6,-3}, // gens=108(6), 24(!3)
  {  60, 1023,  600,   10,   11, 93,  838, 584,  10, 6}, // gens=838(10),584(6)
  //  {  378,  5461,  5292, 14}, // gens=3(126),509(3)
  //  {  630,  8191,  8190, 13}, // gens=39(630)
  //  {  600, 13981, 12000, 20}, // gens=10(30),23(10),3(!2)
  //  {  682, 15709, 15004, 22} // gens=5(682)
};

void checkCiphertext(const Ctxt& ctxt, const ZZX& ptxt, const FHESecKey& sk);

// Testing the I/O of the important classes of the library
// (context, keys, ciphertexts).
int main(int argc, char *argv[])
{
  ArgMapping amap;

  long r=1;
  long p=2;
  long c = 2;
  long w = 64;
  long L = 5;
  long mm=0;
  amap.arg("p", p, "plaintext base");
  amap.arg("r", r,  "lifting");
  amap.arg("c", c, "number of columns in the key-switching matrices");
  amap.arg("m", mm, "cyclotomic index","{31,127,1023}");
  amap.parse(argc, argv);

  bool useTable = (mm==0 && p==2);
  long ptxtSpace = power_long(p,r);
  long numTests = useTable? N_TESTS : 1;

  std::unique_ptr<FHEcontext> contexts[numTests];
  std::unique_ptr<FHESecKey> sKeys[numTests];
  std::unique_ptr<Ctxt> ctxts[numTests];
  std::unique_ptr<EncryptedArray> eas[numTests];
  vector<ZZX> ptxts[numTests];

  // first loop: generate stuff and write it to cout

  // open file for writing
  {fstream keyFile("iotest.txt", fstream::out|fstream::trunc);
   assert(keyFile.is_open());
  for (long i=0; i<numTests; i++) {
    long m = (mm==0)? ms[i][1] : mm;

    cout << "Testing IO: m="<<m<<", p^r="<<p<<"^"<<r<<endl;

    Vec<long> mvec(INIT_SIZE,2);
    mvec[0] = ms[i][4];  mvec[1] = ms[i][5];
    vector<long> gens(2);
    gens[0] = ms[i][6];  gens[1] = ms[i][7];
    vector<long> ords(2);
    ords[0] = ms[i][8];  ords[1] = ms[i][9];

    if (useTable && gens[0]>0)
      contexts[i].reset(new FHEcontext(m, p, r, gens, ords));
    else
      contexts[i].reset(new FHEcontext(m, p, r));
    contexts[i]->zMStar.printout();

    buildModChain(*contexts[i], L, c);  // Set the modulus chain
    if (mm==0 && m==1023) contexts[i]->makeBootstrappable(mvec);

    // Output the FHEcontext to file
    writeContextBase(keyFile, *contexts[i]);
    writeContextBase(cout, *contexts[i]);
    keyFile << *contexts[i] << endl;

    sKeys[i].reset(new FHESecKey(*contexts[i]));
    const FHEPubKey& publicKey = *sKeys[i];
    sKeys[i]->GenSecKey(w,ptxtSpace); // A Hamming-weight-w secret key
    addSome1DMatrices(*sKeys[i]);// compute key-switching matrices that we need
    eas[i].reset(new EncryptedArray(*contexts[i]));

    long nslots = eas[i]->size();

    // Output the secret key to file, twice. Below we will have two copies
    // of most things.
    keyFile << *sKeys[i] << endl;;
    keyFile << *sKeys[i] << endl;;

    vector<ZZX> b;
    long p2r = eas[i]->getContext().alMod.getPPowR();
    ZZX poly = RandPoly(0,to_ZZ(p2r)); // choose a random constant polynomial
    eas[i]->decode(ptxts[i], poly);

    ctxts[i].reset(new Ctxt(publicKey));
    eas[i]->encrypt(*ctxts[i], publicKey, ptxts[i]);
    eas[i]->decrypt(*ctxts[i], *sKeys[i], b);
    assert(ptxts[i].size() == b.size());
    for (long j = 0; j < nslots; j++) assert (ptxts[i][j] == b[j]);

    // output the plaintext
    keyFile << "[ ";
    for (long j = 0; j < nslots; j++) keyFile << ptxts[i][j] << " ";
    keyFile << "]\n";

    eas[i]->encode(poly,ptxts[i]);
    keyFile << poly << endl;

    // Output the ciphertext to file
    keyFile << *ctxts[i] << endl;
    keyFile << *ctxts[i] << endl;
    cerr << "okay " << i << endl<< endl;
  }
  keyFile.close();}
  cerr << "so far, so good\n\n";

  // second loop: read from input and repeat the computation

  // open file for read
  {fstream keyFile("iotest.txt", fstream::in);
  for (long i=0; i<numTests; i++) {

    // Read context from file
    unsigned long m1, p1, r1;
    vector<long> gens, ords;
    readContextBase(keyFile, m1, p1, r1, gens, ords);
    FHEcontext tmpContext(m1, p1, r1, gens, ords);
    keyFile >> tmpContext;
    assert (*contexts[i] == tmpContext);
    cerr << i << ": context matches input\n";

    // We define some things below wrt *contexts[i], not tmpContext.
    // This is because the various operator== methods check equality of
    // references, not equality of the referenced FHEcontext objects.
    FHEcontext& context = *contexts[i];
    FHESecKey secretKey(context);
    FHESecKey secretKey2(tmpContext);
    const FHEPubKey& publicKey = secretKey;
    const FHEPubKey& publicKey2 = secretKey2;

    keyFile >> secretKey;
    keyFile >> secretKey2;
    assert(secretKey == *sKeys[i]);
    cerr << "   secret key matches input\n";

    EncryptedArray ea(context);
    EncryptedArray ea2(tmpContext);

    long nslots = ea.size();

    // Read the plaintext from file
    vector<ZZX> a;
    a.resize(nslots);
    assert(nslots == (long)ptxts[i].size());
    seekPastChar(keyFile, '['); // defined in NumbTh.cpp
    for (long j = 0; j < nslots; j++) {
      keyFile >> a[j];
      assert(a[j] == ptxts[i][j]);
    }
    seekPastChar(keyFile, ']');
    cerr << "   ptxt matches input\n";

    // Read the encoded plaintext from file
    ZZX poly1, poly2;
    keyFile >> poly1;
    eas[i]->encode(poly2,a);
    assert(poly1 == poly2);
    cerr << "   eas[i].encode(a)==poly1 okay\n";

    ea.encode(poly2,a);
    assert(poly1 == poly2);
    cerr << "   ea.encode(a)==poly1 okay\n";

    ea2.encode(poly2,a);
    assert(poly1 == poly2);
    cerr << "   ea2.encode(a)==poly1 okay\n";

    eas[i]->decode(a,poly1);
    assert(nslots == (long)a.size());
    for (long j = 0; j < nslots; j++) assert(a[j] == ptxts[i][j]);
    cerr << "   eas[i].decode(poly1)==ptxts[i] okay\n";

    ea.decode(a,poly1);
    assert(nslots == (long)a.size());
    for (long j = 0; j < nslots; j++) assert(a[j] == ptxts[i][j]);
    cerr << "   ea.decode(poly1)==ptxts[i] okay\n";

    ea2.decode(a,poly1);
    assert(nslots == (long)a.size());
    for (long j = 0; j < nslots; j++) assert(a[j] == ptxts[i][j]);
    cerr << "   ea2.decode(poly1)==ptxts[i] okay\n";

    // Read ciperhtext from file
    Ctxt ctxt(publicKey);
    Ctxt ctxt2(publicKey2);
    keyFile >> ctxt;
    keyFile >> ctxt2;
    assert(ctxts[i]->equalsTo(ctxt,/*comparePkeys=*/false));
    cerr << "   ctxt matches input\n";

    sKeys[i]->Decrypt(poly2,*ctxts[i]);
    assert(poly1 == poly2);
    cerr << "   sKeys[i]->decrypt(*ctxts[i]) == poly1 okay\n";

    secretKey.Decrypt(poly2,*ctxts[i]);
    assert(poly1 == poly2);
    cerr << "   secretKey.decrypt(*ctxts[i]) == poly1 okay\n";

    secretKey.Decrypt(poly2,ctxt);
    assert(poly1 == poly2);
    cerr << "   secretKey.decrypt(ctxt) == poly1 okay\n";

    secretKey2.Decrypt(poly2,ctxt2);
    assert(poly1 == poly2);
    cerr << "   secretKey2.decrypt(ctxt2) == poly1 okay\n";

    eas[i]->decrypt(ctxt, *sKeys[i], a);
    assert(nslots == (long)a.size());
    for (long j = 0; j < nslots; j++) assert(a[j] == ptxts[i][j]);
    cerr << "   eas[i].decrypt(ctxt, *sKeys[i])==ptxts[i] okay\n";

    ea.decrypt(ctxt, secretKey, a);
    assert(nslots == (long)a.size());
    for (long j = 0; j < nslots; j++) assert(a[j] == ptxts[i][j]);
    cerr << "   ea.decrypt(ctxt, secretKey)==ptxts[i] okay\n";

    ea2.decrypt(ctxt2, secretKey2, a);
    assert(nslots == (long)a.size());
    for (long j = 0; j < nslots; j++) assert(a[j] == ptxts[i][j]);
    cerr << "   ea2.decrypt(ctxt2, secretKey2)==ptxts[i] okay\n";

    cerr << "test "<<i<<" okay\n\n";
  }}
  unlink("iotest.txt"); // clean up before exiting
}
