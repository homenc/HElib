#include <fstream>
#include <unistd.h>

#include <NTL/ZZX.h>
#include <NTL/vector.h>

#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"

int main(int argc, char *argv[])
{
  
  ArgMapping amap;

  long m=7;
  long r=1;
  long p=2;
  long c=2;
  long w=64;
  long L=5;
 
  amap.arg("m", m, "order of cyclotomic polynomial");
  amap.arg("p", p, "plaintext base");
  amap.arg("r", r, "lifting");
  amap.arg("c", c, "number of columns in the key-switching matrices");
  amap.parse(argc, argv);

  
// Write ASCII and bin files.
{    
  FHEcontext* context;
 
  ofstream asciiFile("iotest_ascii1.txt");
  ofstream binFile("iotest_bin.bin", ios::binary);
  assert(asciiFile.is_open());  

  context = new FHEcontext(m, p, r);
  buildModChain(*context, L, c);  // Set the modulus chain

  FHESecKey* secKey = new FHESecKey(*context);
  FHEPubKey* pubKey = secKey;
  secKey->GenSecKey(w);

  // ASCII 
  cout << "Writing ASCII1 file." << endl;
  writeContextBase(asciiFile, *context);
  asciiFile << *context << endl << endl;
  asciiFile << *pubKey << endl << endl;

  // Bin
  cout << "Writing Binary file." << endl;
  writeContextBaseBinary(binFile, *context);
  writeContextBinary(binFile, *context);
  writePubKeyBinary(binFile, *pubKey);

  asciiFile.close();
  binFile.close();

  delete context; 

} 

// Read in bin files and write out ASCII.
{
  
  ifstream inFile("iotest_bin.bin", ios::binary);
  ofstream outFile("iotest_ascii2.txt");
  
  // Read in context,
  FHEcontext* context;

  cout << "Reading Binary Base." << endl;
  readContextBaseBinary(inFile, context);  
  cout << "Reading Binary Context." << endl;
  readContextBinary(inFile, *context);  

  // Read in PubKey.
  FHEPubKey* pubKey = new FHEPubKey(*context);
  readPubKeyBinary(inFile, *pubKey);
 
  // ASCII 
  cout << "Writing ASCII2 file." << endl;
  writeContextBase(outFile, *context);
  outFile << *context << endl << endl;
  outFile << *pubKey << endl << endl;

  inFile.close();
  outFile.close();

  delete context; 

}

  return 0;
}


