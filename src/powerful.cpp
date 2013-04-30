/* a prelimary test program for playing around
 * with Peikert's "powerful" basis.
 */

#include "NumbTh.h"

#include <cassert>

using namespace std;
using namespace NTL;


long computePhi(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e - 1) * (p-1);
}

long computePow(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e);
}


long computeProd(const Vec<long>& vec)
{
  long prod = 1;
  long k = vec.length();
  for (long i = 0; i < k; i++)
    prod = prod * vec[i];
  return prod;
}

long computeProd(const Vec< Pair<long, long> >& vec)
{
  long prod = 1;
  long k = vec.length();
  for (long i = 0; i < k; i++) {
    prod = prod * computePow(vec[i]);
  }
  return prod;
}

void computePhiVec(Vec<long>& phiVec, 
                   const Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  phiVec.SetLength(k);

  for (long i = 0; i < k; i++) 
    phiVec[i] = computePhi(factors[i]);
}

void computeDivVec(Vec<long>& divVec,
                   const Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  divVec.SetLength(k);

  long m = computeProd(factors);

  for (long i = 0; i < k; i++)
    divVec[i] = m/computePow(factors[i]);
}
        
void mapIndexToPowerful(Vec<long>& pow, long j, const Vec<long>& phiVec)
// this maps an index j in [phi(m)] to a vector
// representing the powerful basis coordinates

{
  long k = phiVec.length();
  long phim = computeProd(phiVec);
  assert(j >= 0 && j < phim);

  pow.SetLength(k);

  for (long i = k-1; i >= 0; i--) {
    pow[i] = j % phiVec[i];
    j = (j - pow[i])/phiVec[i];
  }
}


void mapPowerfulToPoly(ZZX& poly, 
                       const Vec<long>& pow, 
                       const Vec<long>& divVec,
                       long m,
                       const ZZX& phimX)
{
  long k = pow.length();
  assert(divVec.length() == k);

  long j = 0;
  for (long i = 0; i < k; i++)
    j += pow[i] * divVec[i];

  j %= m;

  ZZX f = ZZX(j, 1);

  poly = f % phimX;
}

void usage()
{
  cerr << "bad args\n";
  exit(0);
}

int main(int argc, char *argv[])
{
  argmap_t argmap;

  argmap["q"] = "101";
  argmap["m"] = "100";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage();

  long q = atoi(argmap["q"]);
  long m = atoi(argmap["m"]);

  cout << "q=" << q << "\n";
  cout << "m=" << m << "\n";

  Vec< Pair<long, long> > factors;

  factorize(factors, m);

  cout << factors << "\n";

  Vec<long> phiVec;
  computePhiVec(phiVec, factors);
  cout << phiVec << "\n";

  long phim = computeProd(phiVec);
  cout << phim << "\n";

  Vec<long> divVec;
  computeDivVec(divVec, factors);
  cout << divVec << "\n";


  ZZX phimX = Cyclotomic(m);


  for (long j = 0; j < phim; j++) {
    Vec<long> pow;
    mapIndexToPowerful(pow, j, phiVec);
    ZZX poly;
    mapPowerfulToPoly(poly, pow, divVec, m, phimX);
    cout << pow << "  " << poly << "\n";
  }

}

