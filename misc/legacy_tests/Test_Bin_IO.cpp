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
#include <cassert>
#include <cstring>
#include <fstream>
#include <utility>
#include <unistd.h>

#include <NTL/ZZX.h>
#include <NTL/vector.h>

#include <helib/helib.h>
#include <helib/ArgMap.h>

#include <helib/debugging.h>

NTL_CLIENT
using namespace helib;

bool isLittleEndian()
{
    int i=1;
    return static_cast<bool>(*reinterpret_cast<char *>(&i));
}

void cleanupFiles(const char* file){
  if(unlink(file)) cerr<< "Delete of "<< file <<" failed."<<endl;
}

template<class... Files>
void cleanupFiles(const char * file, Files... files){
  cleanupFiles(file);
  cleanupFiles(files...);
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
  ArgMap amap;

  bool noPrint=true;
  long m=7;
  long r=1;
  long p=2;
  long c=2;
  long w=64;
  long L=300;
  long cleanup=1;
  string sampleFilePrefix;

  amap.arg("m", m, "order of cyclotomic polynomial");
  amap.arg("p", p, "plaintext base");
  amap.arg("r", r, "lifting");
  amap.arg("c", c, "number of columns in the key-switching matrices");
  amap.arg("L", L, "number of levels wanted");
  amap.arg("sample", sampleFilePrefix, "sample file prefix e.g. <prefix>_BE.txt");
  amap.arg("cleanup", cleanup, "cleanup files created");
  amap.arg("noPrint", noPrint, "suppress printouts");
  amap.parse(argc, argv);

  // FIXME: this is wrong!
  const char* asciiFile1 = "../misc/iotest_ascii1.txt";
  const char* asciiFile2 = "../misc/iotest_ascii2.txt";
  const char* binFile1 = "../misc/iotest_bin.bin";
  const char* otherEndianFileOut = "../misc/iotest_ascii3.txt";

  { // 1. Write ASCII and bin files.
    ofstream asciiFile(asciiFile1);
    ofstream binFile(binFile1, ios::binary);
    assert(asciiFile.is_open());

    std::unique_ptr<Context> context(new Context(m, p, r));
    buildModChain(*context, L, c);  // Set the modulus chain

    if (!noPrint) {
      cout << "Test to write out ASCII and Binary Files.\n";
      context->zMStar.printout(); // Printout context params
      cout << "\tSecurity Level: " << context->securityLevel() << endl;
    }
    std::unique_ptr<SecKey> secKey(new SecKey(*context));
    PubKey* pubKey = (PubKey*) secKey.get();
    secKey->GenSecKey(w);
    addSome1DMatrices(*secKey);
    addFrbMatrices(*secKey);

#ifdef HELIB_DEBUG
        dbgEa = context->ea;
        dbgKey = secKey.get();
#endif

    // ASCII
    if (!noPrint)
      cout << "\tWriting ASCII1 file " << asciiFile1 << endl;
    writeContextBase(asciiFile, *context);
    asciiFile << *context << endl << endl;
    asciiFile << *pubKey << endl << endl;
    asciiFile << *secKey << endl << endl;

    // Bin
    if (!noPrint)
      cout << "\tWriting Binary file " << binFile1<< endl;
    writeContextBaseBinary(binFile, *context);
    writeContextBinary(binFile, *context);
    writePubKeyBinary(binFile, *pubKey);
    writeSecKeyBinary(binFile, *secKey);

    asciiFile.close();
    binFile.close();
    cout << "GOOD\n";
  }
  { // 2. Read in bin files and write out ASCII.
    if (!noPrint)
      cout << "Test to read binary file and write it out as ASCII" << endl;

    ifstream inFile(binFile1, ios::binary);
    ofstream outFile(asciiFile2);

    // Read in context,
    std::unique_ptr<Context> context = buildContextFromBinary(inFile);
    readContextBinary(inFile, *context);

    // Read in SecKey and PubKey.
    std::unique_ptr<SecKey> secKey(new SecKey(*context));

#ifdef HELIB_DEBUG
        dbgEa = context->ea;
        dbgKey = secKey.get();
#endif

    PubKey* pubKey = (PubKey*) secKey.get();

    readPubKeyBinary(inFile, *pubKey);
    readSecKeyBinary(inFile, *secKey);

    // ASCII
    if (!noPrint)
      cout << "\tWriting ASCII2 file." << endl;
    writeContextBase(outFile, *context);
    outFile << *context << endl << endl;
    outFile << *pubKey << endl << endl;
    outFile << *secKey << endl << endl;

    inFile.close();
    outFile.close();

    cout << "GOOD\n";
  }
  { // 3. Compare byte-wise the two ASCII files
    if (!noPrint)
      cout << "Comparing the two ASCII files\n";

    long differ = compareFiles(asciiFile1, asciiFile2);

    if(differ != 0){
      cout << "BAD\n";
      if (!noPrint)
        cout << "\tFAIL - Files differ. Return Code: " << differ << endl;
      exit(EXIT_FAILURE);
    }
    cout << "GOOD\n";
  }
  { // 4. Read in binary and perform operation.
    if (!noPrint)
      cout << "Test reading in Binary files and performing an operation between two ctxts\n";

    ifstream inFile(binFile1, ios::binary);

    // Read in context,
    std::unique_ptr<Context> context = buildContextFromBinary(inFile);
    readContextBinary(inFile, *context);

    // Read in PubKey.
    std::unique_ptr<SecKey> secKey(new SecKey(*context));
    PubKey* pubKey = (PubKey*) secKey.get();

#ifdef HELIB_DEBUG
        dbgEa = context->ea;
        dbgKey = secKey.get();
#endif

        readPubKeyBinary(inFile, *pubKey);
    readSecKeyBinary(inFile, *secKey);
    inFile.close();

    // Get the ea
    const EncryptedArray& ea = *context->ea;

    // Setup some ptxts and ctxts.
    Ctxt c1(*pubKey), c2(*pubKey);
    PlaintextArray p1(ea),  p2(ea);

    random(ea, p1);
    random(ea, p2);

    ea.encrypt(c1, *pubKey, p1);
    ea.encrypt(c2, *pubKey, p2);

    // Operation multiply and add.
    mul(ea, p1, p2);
    c1.multiplyBy(c2);
    //c1 *= c2;

    // Decrypt and Compare.
    PlaintextArray pp1(ea);
    ea.decrypt(c1, *secKey, pp1);

    if(!equals(ea, p1, pp1)) {
      cout << "BAD\n";
      exit(EXIT_FAILURE);
    }
    cout << "GOOD\n";

    if(cleanup) {
      if (!noPrint)
        cout << "Clean up. Deleting created files." << endl;
      cleanupFiles(asciiFile1, asciiFile2, binFile1);
    }
  }
  { // 5. Read in binary from opposite little endian and print ASCII and compare
    bool littleEndian = isLittleEndian();

    string otherEndianFileIn
      = sampleFilePrefix + (littleEndian? "_BE.bin" : "_LE.bin");
    string otherEndianASCII
      = sampleFilePrefix + (littleEndian? "_BE.txt" : "_LE.txt");

    if (!noPrint)
      cout << "Test to read in" << (littleEndian? " BE ":" LE ")
           << "binary file and write it out as ASCII" << endl;

    if(sampleFilePrefix.empty()) {
      if (!noPrint)
        cout << "\tSample prefix not provided, test not done." << endl;
    } else {
      if (!noPrint)
        cout << "\tSample file used: " << otherEndianFileIn << endl;

      ifstream inFile(otherEndianFileIn, ios::binary);

      if(!inFile.is_open()) {
        cout << "BAD boo!\n";
        if (!noPrint)
          cout << "  file " << otherEndianFileIn
               << " could not be opened.\n";
        exit(EXIT_FAILURE);
      }
      ofstream outFile(otherEndianFileOut);

      // Read in context,
      std::unique_ptr<Context> context = buildContextFromBinary(inFile);
      readContextBinary(inFile, *context);

      // Read in SecKey and PubKey.
      std::unique_ptr<SecKey> secKey(new SecKey(*context));
      PubKey* pubKey = (PubKey*) secKey.get();

#ifdef HELIB_DEBUG
        dbgEa = context->ea;
        dbgKey = secKey.get();
#endif

        readPubKeyBinary(inFile, *pubKey);
      readSecKeyBinary(inFile, *secKey);
      inFile.close();

      // ASCII
      if (!noPrint)
        cout << "\tWriting other endian file." << endl;
      writeContextBase(outFile, *context);
      outFile << *context << endl << endl;
      outFile << *pubKey << endl << endl;
      outFile << *secKey << endl << endl;
      outFile.close();

      // Compare byte-wise the two ASCII files
      if (!noPrint)
        cout << "Comparing the two ASCII files\n";

      long differ = compareFiles(otherEndianASCII, otherEndianFileOut);

      if(differ != 0) {
        cout << "BAD\n";
        exit(EXIT_FAILURE);
      }
      cout << "GOOD\n";

      if(cleanup) {
        if (!noPrint)
          cout << "Clean up. Deleting created files." << endl;
        cleanupFiles(otherEndianFileOut);
      }
    }
  }
  return 0;
}
