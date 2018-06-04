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
#include <cstring>
#include <fstream>
#include <unistd.h>

#include <NTL/ZZX.h>
#include <NTL/vector.h>

#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"

bool isLittleEndian()
{
    int i=1;
    return static_cast<bool>(*reinterpret_cast<char *>(&i));
}

// Compare two binary files, return 0 if they are equal, -1 if they
// have different length, and 1 if they have the same length but
// different content.
long compareFiles(string filename1, string filename2)
{

  ifstream file1(filename1);
  ifstream file2(filename2);
 
  if(!file1.is_open() && !file2.is_open()){
    cerr << "Could not open one of the following files:" << endl;
    cerr << filename1 << " and/or " << filename2 << endl;
    exit(EXIT_FAILURE);
  }
 
  fstream::pos_type file1size, file2size;

  file1size = file1.seekg(0, ifstream::end).tellg();
  file1.seekg(0, ifstream::beg);

  file2size = file2.seekg(0, ifstream::end).tellg();
  file2.seekg(0, ifstream::beg);

  // Quick compare sizes.
  if(file1size != file2size){ 
    file1.close();
    file2.close();
    cerr << "Files "<<filename1<<" and "<<filename2<<" not the same size :("<< endl;
    return -1;
  }

  // Now compare byte blocks at a time.
  const size_t BLOCKSIZE = 4096; // 4 kB

  char buffer1[BLOCKSIZE];
  char buffer2[BLOCKSIZE];
  size_t curBlckSz = 0;

  for(size_t i=file1size, cnt=0; i > 0; i-=curBlckSz, cnt++) {

    curBlckSz = std::min(BLOCKSIZE, i);

    file1.read(buffer1, curBlckSz);
    file2.read(buffer2, curBlckSz);

    if(memcmp(buffer1, buffer2, curBlckSz)) {
      cerr << "Block "
           <<cnt<<" (block size: "<<BLOCKSIZE<<" bytes) "
           <<cnt<<" does not match :("<<endl;
      return 1; 
    }
  }  
  return 0; // Files are the same!
}


int main(int argc, char *argv[])
{ 
  cout << "\n*** TEST BINARY IO " << endl << endl;
   
  ArgMapping amap;

  long m=7;
  long r=1;
  long p=2;
  long c=2;
  long w=64;
  long L=5;
  long cleanup=1;
 
  amap.arg("m", m, "order of cyclotomic polynomial");
  amap.arg("p", p, "plaintext base");
  amap.arg("r", r, "lifting");
  amap.arg("c", c, "number of columns in the key-switching matrices");
  amap.arg("L", L, "number of levels wanted");
  amap.arg("cleanup", cleanup, "cleanup files created");
  amap.parse(argc, argv);

  const char* asciiFile1 = "misc/iotest_ascii1.txt"; 
  const char* asciiFile2 = "misc/iotest_ascii2.txt"; 
  const char* binFile1 = "misc/iotest_bin.bin"; 

  { // 1. Write ASCII and bin files.

    cout << "Test to write out ASCII and Binary Files.\n";    
 
    ofstream asciiFile(asciiFile1);
    ofstream binFile(binFile1, ios::binary);
    assert(asciiFile.is_open());  

    std::unique_ptr<FHEcontext> context(new FHEcontext(m, p, r));
    buildModChain(*context, L, c);  // Set the modulus chain

    context->zMStar.printout(); // Printout context params
    cout << "\tSecurity Level: " << context->securityLevel() << endl;

    std::unique_ptr<FHESecKey> secKey(new FHESecKey(*context));
    FHEPubKey* pubKey = (FHEPubKey*) secKey.get();
    secKey->GenSecKey(w);
    addSome1DMatrices(*secKey);
    addFrbMatrices(*secKey);

    // ASCII 
    cout << "\tWriting ASCII1 file " << asciiFile1 << endl;
    writeContextBase(asciiFile, *context);
    asciiFile << *context << endl << endl;
    asciiFile << *pubKey << endl << endl;
    asciiFile << *secKey << endl << endl;

    // Bin
    cout << "\tWriting Binary file " << binFile1<< endl;
    writeContextBaseBinary(binFile, *context);
    writeContextBinary(binFile, *context);
    writePubKeyBinary(binFile, *pubKey);
    writeSecKeyBinary(binFile, *secKey);

    asciiFile.close();
    binFile.close();

    cout << "Test successful.\n\n";
  }
  { // 2. Read in bin files and write out ASCII.
    cout << "Test to read binary file and write it out as ASCII" << endl;
  
    ifstream inFile(binFile1, ios::binary);
    ofstream outFile(asciiFile2);
  
    // Read in context,
    //  cout << "\tReading Binary Base." << endl;
    std::unique_ptr<FHEcontext> context = buildContextFromBinary(inFile);  
    //  cout << "\tReading Binary Context." << endl;
    readContextBinary(inFile, *context);  

    // Read in SecKey and PubKey.
    // Got to insert pubKey into seckey obj first.
    std::unique_ptr<FHESecKey> secKey(new FHESecKey(*context));
    FHEPubKey* pubKey = (FHEPubKey*) secKey.get();
  
    readPubKeyBinary(inFile, *pubKey);
    readSecKeyBinary(inFile, *secKey);
 
    // ASCII 
    cout << "\tWriting ASCII2 file." << endl;
    writeContextBase(outFile, *context);
    outFile << *context << endl << endl;
    outFile << *pubKey << endl << endl;
    outFile << *secKey << endl << endl;

    inFile.close();
    outFile.close();

    cout << "Test successful.\n\n";
  }
  { // 3. Compare byte-wise the two ASCII files
    cout << "Comparing the two ASCII files\n"; 
  
    long differ = compareFiles(asciiFile1, asciiFile2); 

    if(differ != 0){
      cout << "\tFAIL - Files differ. Return Code: " << differ << endl;
      cout << "Test failed.\n";
      exit(EXIT_FAILURE);
    } else {
      cout << "\tSUCCESS - Files are identical.\n";
    }

    cout << "Test successful.\n\n";
  }
  { // 4. Read in binary and perform operation.
    cout << "Test reading in Binary files and performing an operation between two ctxts\n";  

    ifstream inFile(binFile1, ios::binary);

    // Read in context,
    //  cout << "Reading Binary Base." << endl;
    std::unique_ptr<FHEcontext> context = buildContextFromBinary(inFile);
    //  cout << "Reading Binary Context." << endl;
    readContextBinary(inFile, *context);  

    // Read in PubKey.
    std::unique_ptr<FHESecKey> secKey(new FHESecKey(*context));
    FHEPubKey* pubKey = (FHEPubKey*) secKey.get();
    readPubKeyBinary(inFile, *pubKey);
    readSecKeyBinary(inFile, *secKey);
    inFile.close(); 

    // Get the ea
    const EncryptedArray& ea = *context->ea;
 
    // Setup some ptxts and ctxts.
    Ctxt c1(*pubKey), c2(*pubKey);
    NewPlaintextArray p1(ea),  p2(ea);

    random(ea, p1);
    random(ea, p2);

    ea.encrypt(c1, *pubKey, p1);
    ea.encrypt(c2, *pubKey, p2);

    // Operation multiply and add.
    mul(ea, p1, p2);
    //  random(ea, p1);
    c1 *= c2;

    // Decrypt and Compare.
    NewPlaintextArray pp1(ea);
    ea.decrypt(c1, *secKey, pp1);     
	
    if(!equals(ea, p1, pp1)){
      cout << "\tFAIL - polynomials are not equal\n";
      cout << "Test failed.\n";
      exit(EXIT_FAILURE);
    }
    cout << "\tGOOD - polynomials are equal\n";

    if(cleanup) {
      cout << "Clean up. Deleting created files." << endl;
      if(unlink(asciiFile1)) cerr<< "Delete of "<<asciiFile1<<" failed."<<endl; 
      if(unlink(asciiFile2)) cerr<< "Delete of "<<asciiFile2<<" failed."<<endl; 
      if(unlink(binFile1)) cerr << "Delete of "<<binFile1<<" failed."<<endl;  
    }

    cout << "Test successful.\n\n";
  }
  { // 5. Read in binary from opposite little endian and print ASCII and compare.
    cout << "Test to read binary file and write it out as ASCII" << endl;
 
    string otherEndianFileIn
      = isLittleEndian()? "misc/iotest_binBE.bin" : "misc/iotest_binLE.bin";
    string otherEndianASCII
      = isLittleEndian()? "misc/iotest_asciiBE.txt" : "misc/iotest_asciiLE.txt";
    cout << "\tSample file used: " << otherEndianFileIn << endl;

    string otherEndianFileOut = "misc/iotest_ascii3.txt";  
    ifstream inFile(otherEndianFileIn, ios::binary);
    ofstream outFile(otherEndianFileOut);

    if(!inFile.is_open()) {
      cout << "File " << otherEndianFileIn 
           << " could not be opened.\n";
      cout << "Test failed.\n";
      exit(EXIT_FAILURE);
    }
  
    // Read in context,
    std::unique_ptr<FHEcontext> context = buildContextFromBinary(inFile);
    readContextBinary(inFile, *context);  

    // Read in SecKey and PubKey.
    // Got to insert pubKey into seckey obj first.
    std::unique_ptr<FHESecKey> secKey(new FHESecKey(*context));
    FHEPubKey* pubKey = (FHEPubKey*) secKey.get();
  
    readPubKeyBinary(inFile, *pubKey);
    readSecKeyBinary(inFile, *secKey);
    inFile.close();
 
    // ASCII 
    cout << "\tWriting other endian file." << endl;
    writeContextBase(outFile, *context);
    outFile << *context << endl << endl;
    outFile << *pubKey << endl << endl;
    outFile << *secKey << endl << endl;
    outFile.close();

    cout << "Test successful.\n\n";

    // Compare byte-wise the two ASCII files

    cout << "Comparing the two ASCII files\n"; 
  
    long differ = compareFiles(otherEndianASCII, otherEndianFileOut); 

    if(differ != 0){
      cout << "\tFAIL - Files differ. Return Code: " << differ << endl;
      cout << "Test failed.\n";
      exit(EXIT_FAILURE);
    } else {
      cout << "\tSUCCESS - Files are identical.\n";
    }

    if(cleanup){
      cout << "Clean up. Deleting created files." << endl;
      if(unlink(otherEndianFileOut.c_str())) 
        cerr << "Delete of "<<otherEndianFileOut<<" failed."<<endl; 
    }
    cout << "Test successful.\n\n";
  }

  return 0;
}
