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


template<class type> 
class RunningSumMatrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

public:
  ~RunningSumMatrix() { cerr << "destructor: running sum matrix\n"; }

  RunningSumMatrix(const EncryptedArray& _ea) : ea(_ea) { }

  virtual const EncryptedArray& getEA() const {
    return ea;
  }

  virtual void get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if (j >= i)
      out = 1;
    else
      out = 0;
  }
};


PlaintextMatrixBaseInterface *
buildRunningSumMatrix(const EncryptedArray& ea)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RunningSumMatrix<PA_GF2>(ea);
    }

    case PA_zz_p_tag: {
      return new RunningSumMatrix<PA_zz_p>(ea);
    }

    default: return 0;
  }
}


template<class type> 
class TotalSumMatrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

public:
  ~TotalSumMatrix() { cerr << "destructor: total sum matrix\n"; }

  TotalSumMatrix(const EncryptedArray& _ea) : ea(_ea) { }

  virtual const EncryptedArray& getEA() const {
    return ea;
  }

  virtual void get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    out = 1;
  }
};


PlaintextMatrixBaseInterface *
buildTotalSumMatrix(const EncryptedArray& ea)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new TotalSumMatrix<PA_GF2>(ea);
    }

    case PA_zz_p_tag: {
      return new TotalSumMatrix<PA_zz_p>(ea);
    }

    default: return 0;
  }
}


template<class type> 
class RandomMatrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< RX > > data;

public:
  ~RandomMatrix() { cerr << "destructor: random matrix\n"; }

  RandomMatrix(const EncryptedArray& _ea) : ea(_ea) { 
    long n = ea.size();
    long d = ea.getDegree();

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.resize(n);
    for (long i = 0; i < n; i++) {
      data[i].resize(n);
      for (long j = 0; j < n; j++)
        random(data[i][j], d);
    }
  }

  virtual const EncryptedArray& getEA() const {
    return ea;
  }

  virtual void get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    out = data[i][j];
  }
};


PlaintextMatrixBaseInterface *
buildRandomMatrix(const EncryptedArray& ea)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix<PA_GF2>(ea);
    }

    case PA_zz_p_tag: {
      return new RandomMatrix<PA_zz_p>(ea);
    }

    default: return 0;
  }
}



void  TestIt(long R, long p, long r, long d, long c, long k, long w, 
               long L, long m)
{
  cerr << "*** TestIt: R=" << R 
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", c=" << c
       << ", k=" << k
       << ", w=" << w
       << ", L=" << L
       << ", m=" << m
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, L, c);

  context.zMStar.printout();
  cerr << endl;

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key


  ZZX G;

  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  cerr << "G = " << G << "\n";
  cerr << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  cerr << "done\n";

  printAllTimers();
  resetAllTimers();

  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cerr << "done\n";

  // choose a random plaintext square matrix
  //  PlaintextMatrixBaseInterface *ptr = buildTotalSumMatrix(ea);
  PlaintextMatrixBaseInterface *ptr = buildRandomMatrix(ea);

  // choose a random plaintext vector
  PlaintextArray v(ea);

/* v.encode(1);  
   v.print(cout); cout << "\n";
   v.alt_mul(*ptr);
*/
  v.random();
  // v.print(cout); cout << "\n";

  // encrypt the random vector
  Ctxt ctxt(publicKey);
  startFHEtimer("ea.encrypt");
  ea.encrypt(ctxt, publicKey, v);
  stopFHEtimer("ea.encrypt");

  v.mat_mul(*ptr);         // multiply the plaintext vector
  startFHEtimer("ea.mat_mul");
  ea.mat_mul(ctxt, *ptr);  // multiply the ciphertext vector
  stopFHEtimer("ea.mat_mul");
  //  totalSums(ea, ctxt);  // multiply the ciphertext vector

  PlaintextArray v1(ea);
  //  v1 = v;
  startFHEtimer("ea.decrypt");
  ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector
  stopFHEtimer("ea.decrypt");

  if (v.equals(v1))        // check that we've got the right answer
    cout << "Nice!!\n";
  else
    cout << "Grrr...\n";

  // v.print(cout); cout << "\n";
  // v1.print(cout); cout << "\n";
}




void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=1 p=2 k=80'\n\n";
  cerr << "  R is the number of rounds\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==1]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chai [default=4*R]\n";
  cerr << "  s is the minimum number of slots [default=4]\n";
  cerr << "  m defined the cyclotomic polynomial Phi_m(X)\n";
  exit(0);
}


int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["R"] = "1";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["d"] = "1";
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["L"] = "0";
  argmap["s"] = "0";
  argmap["m"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long R = atoi(argmap["R"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  //  long z = atoi(argmap["z"]);
  long L = atoi(argmap["L"]);
  if (L==0) { // determine L based on R,r
    if (r==1) L = 2*R+2;
    else      L = 4*R;
  }
  long s = atoi(argmap["s"]);
  long chosen_m = atoi(argmap["m"]);

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels

  long m = FindM(k, L, c, p, d, s, chosen_m, true);

  //  setTimersOn();
  TestIt(R, p, r, d, c, k, w, L, m);

  cerr << endl;
  printAllTimers();
  cerr << endl;

}

