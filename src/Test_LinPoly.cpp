/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>
#include <cassert>
#include <cstdio>

void  TestIt(long m, long p, long r, long d)
{
  cout << "\n\n******** TestIt: m=" << m 
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, /*L=*/3, /*c=*/2);

  // context.lazy = false;
  if (context.lazy)
    cout << "LAZY REDUCTIONS\n";
  else
    cout << "NON-LAZY REDUCTIONS\n";

  context.zMStar.printout();
  cout << endl;

  cout << "generating keys and key-switching matrices... " << std::flush;
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64);// A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  cout << "done\n";

  ZZX G;
  if (d == 0) {
    G = context.alMod.getFactorsOverZZ()[0];
    d = deg(G);
  }
  else
    G = makeIrredPoly(p, d);
  cout << "G = " << G << "\n";
  cout << "computing masks and tables for rotation... " << std::flush;
  EncryptedArray ea(context, G);
  cout << "done\n";

  long nslots = ea.size();

  PlaintextArray p0(ea);
  PlaintextArray pp0(ea);

  Ctxt c0(publicKey);

  cout << "\nTest #1: pply the same linear transformation to all slots\n";
  {
  vector<ZZX> LM(d); // LM selects even coefficients
  for (long j = 0; j < d; j++) 
    if (j % 2 == 0) LM[j] = ZZX(j, 1);

  // "building" the linearized-polynomial coefficients
  vector<ZZX> C;
  ea.buildLinPolyCoeffs(C, LM);

  p0.random();  
  ea.encrypt(c0, publicKey, p0);
  applyLinPoly1(ea, c0, C);
  ea.decrypt(c0, secretKey, pp0);

  // FIXME: Put a code here to check the result
  //  p0.print(cout); cout << "\n";
  //  pp0.print(cout); cout << "\n";
  }

  cout << "\nTest #2: Apply different transformations to the different slots\n";
  {
  vector< vector<ZZX> > LM(nslots); 
  // LM[i] rotates the coefficients in the i'th slot by (i % d)
  for (long i = 0; i < nslots; i++) {
    LM[i].resize(d);
    for (long j = 0; j < d; j++)  {
      long jj = (i+j) % d;
      LM[i][j] = ZZX(jj, 1);
    }
  }

  // "building" the linearized-polynomial coefficients
  vector< vector<ZZX> > C(nslots);
  for (long i = 0; i < nslots; i++)
    ea.buildLinPolyCoeffs(C[i], LM[i]);

  p0.random();
  ea.encrypt(c0, publicKey, p0);
  applyLinPolyMany(ea, c0, C); // apply the linearized polynomials
  ea.decrypt(c0, secretKey, pp0);

  // FIXME: Put a code here to check the result
  }

  cout << "\nTest #3: Testing low-level (cached) implementation\n";
  {
  vector< vector<ZZX> > LM(nslots); 
  // LM[i] adds coefficients (i % d) and (i+1 % d) in the i'th slot
  for (long i = 0; i < nslots; i++) {
    LM[i].resize(d);
    for (long j = 0; j < d; j++)  {
      if ( j == (i % d) || j == ((i+1)%d) )
	LM[i][j] = conv<ZZX>(1L);
    }
  }

  // "building" the linearized-polynomial coefficients
  vector< vector<ZZX> > C(nslots);
  for (long i = 0; i < nslots; i++)
    ea.buildLinPolyCoeffs(C[i], LM[i]);

  // "encoding" the linearized-polynomial coefficients
  vector<ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long i = 0; i < nslots; i++) v[i] = LM[i][j];
    ea.encode(encodedC[j], v);
  }

  p0.random();  
  ea.encrypt(c0, publicKey, p0);
  applyLinPolyLL(ea, c0, encodedC); // apply the linearized polynomials
  ea.decrypt(c0, secretKey, pp0);

  // FIXME: Put a code here to check the result
  }
}


int main(int argc, char *argv[]) 
{
  ArgMapping amap;

  long m=91;
  amap.arg("m", m, "use specified value as modulus");

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  long d=1;
  amap.arg("d", d, "degree of the field extension");
  amap.note("d == 0 => factors[0] defines extension");

  amap.parse(argc, argv);

  long repeat = 2;
  setTimersOn();
  for (long repeat_cnt = 0; repeat_cnt < repeat; repeat_cnt++) {
    TestIt(m, p, r, d);
  }

}
