/* Copyright (C) 2020-2021 IBM Corp.
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

// To do useful CKKS work, we need the ability to write
// our objects to disk or to transmit them.
// This is a sample program for education purposes only.
// It attempts to show the serialization/deserialization
// operations that can be performed on contexts, public
// keys, secret keys, plaintexts, and ciphertexts.

// NOTE The serialization used in this program reads and
// writes to JSON for demonstration. However, once you
// feel more confident you can also try the binary
// serialization APIs which are similar to these JSON
// ones.

// NOTE In this tutorial program we introduce a different
// API for creating plaintext objects. These objects can be
// (de)serialised.

#include <fstream>

#include <helib/helib.h>

int main(int argc, char* argv[])
{
  // CKKS context created with a builder.
  helib::Context context = helib::ContextBuilder<helib::CKKS>()
                               .m(128)
                               .precision(20)
                               .bits(30)
                               .c(3)
                               .build();

  // NOTE These chosen parameters are for demonstration only. They do not the
  // provide security level that might be required for real use/application
  // scenarios.

  std::ofstream outContextFile;
  outContextFile.open("context.json", std::ios::out);
  if (outContextFile.is_open()) {
    // Write the context to a file
    context.writeToJSON(outContextFile);
    // Close the ofstream
    outContextFile.close();
  } else {
    throw std::runtime_error("Could not open file 'context.json'.");
  }

  std::ifstream inContextFile;
  inContextFile.open("context.json");
  if (inContextFile.is_open()) {
    // Read in the context from the file
    helib::Context deserializedContext =
        helib::Context::readFromJSON(inContextFile);
    // Close the ifstream
    inContextFile.close();
  } else {
    throw std::runtime_error("Could not open file 'context.json'.");
  }
  // Remove the context file. Comment out the line below if you would like to
  // inspect the file.
  std::remove("context.json");

  // Create a secret key associated with the CKKS context
  helib::SecKey secretKey(context);

  // Generate the secret key
  secretKey.GenSecKey();

  // Compute key-switching matrices that we need
  helib::addSome1DMatrices(secretKey);

  std::ofstream outSecretKeyFile;
  outSecretKeyFile.open("sk.json", std::ios::out);
  if (outSecretKeyFile.is_open()) {
    // Write the secret key to a file
    secretKey.writeToJSON(outSecretKeyFile);
    // Close the ofstream
    outSecretKeyFile.close();
  } else {
    throw std::runtime_error("Could not open file 'sk.json'.");
  }

  std::ifstream inSecretKeyFile;
  inSecretKeyFile.open("sk.json");
  if (inSecretKeyFile.is_open()) {
    // Read in the secret key from the file
    helib::SecKey deserializedSecretKey =
        helib::SecKey::readFromJSON(inSecretKeyFile, context);
    // Note there are alternative methods for deserialization of SecKey objects.
    // After initialization
    // helib::SecKey deserializedSecretKey(context);
    // One can write
    // inSecretKeyFile >> deserializedSecretKey;
    // Or alternatively
    // deserializedSecretKey.readJSON(inSecretKeyFile);

    // Close the ifstream
    inSecretKeyFile.close();
  } else {
    throw std::runtime_error("Could not open file 'sk.json'.");
  }
  // Remove the context file. Comment out the line below if you would like to
  // inspect the file.
  std::remove("sk.json");

  // Create an alias public key part from the secret key
  const helib::PubKey& publicKey = secretKey;

  std::ofstream outPublicKeyFile;
  outPublicKeyFile.open("pk.json", std::ios::out);
  if (outPublicKeyFile.is_open()) {
    // Write the public key to a file
    publicKey.writeToJSON(outPublicKeyFile);
    // Close the ofstream
    outPublicKeyFile.close();
  } else {
    throw std::runtime_error("Could not open file 'pk.json'.");
  }

  std::ifstream inPublicKeyFile;
  inPublicKeyFile.open("pk.json");
  if (inPublicKeyFile.is_open()) {
    // Read in the public key from the file
    helib::PubKey deserializedPublicKey =
        helib::PubKey::readFromJSON(inPublicKeyFile, context);
    // Note there are alternative methods for deserialization of PubKey objects.
    // After initialization
    // helib::PubKey deserializedPublicKey(context);
    // One can write
    // inPublicKeyFile >> deserializedPublicKey;
    // Or alternatively
    // deserializedPublicKey.readJSON(inPublicKeyFile);

    // Close the ifstream
    inPublicKeyFile.close();
  } else {
    throw std::runtime_error("Could not open file 'pk.json'.");
  }
  // Remove the context file. Comment out the line below if you would like to
  // inspect the file.
  std::remove("pk.json");

  // Create a data object
  std::vector<long> data(context.getNSlots());

  // Generate some data
  std::iota(data.begin(), data.end(), 0);

  // Create a ptxt.
  helib::PtxtArray ptxt(context, data);

  std::ofstream outPtxtFile;
  outPtxtFile.open("ptxt.json", std::ios::out);
  if (outPtxtFile.is_open()) {
    // Write the ptxt to a file
    ptxt.writeToJSON(outPtxtFile);
    // Alternatively one can use
    // outPtxtFile << ptxt;

    // Close the ofstream
    outPtxtFile.close();
  } else {
    throw std::runtime_error("Could not open file 'ptxt.json'.");
  }

  std::ifstream inPtxtFile;
  inPtxtFile.open("ptxt.json");
  if (inPtxtFile.is_open()) {
    // Read in the ptxt from the file
    helib::PtxtArray deserializedPtxt =
        helib::PtxtArray::readFromJSON(inPtxtFile, context);
    // Note there are alternative methods for deserialization of Ptxt objects.
    // After initialization
    // helib::PtxtArray deserializedPtxt(publicKey);
    // One can write
    // inPtxtFile >> deserializedPtxt;
    // Or alternatively
    // deserializedPtxt.readJSON(inPtxtFile);

    // Close the ifstream
    inPtxtFile.close();
  } else {
    throw std::runtime_error("Could not open file 'ptxt.json'.");
  }
  // Remove the context file. Comment out the line below if you would like to
  // inspect the file.
  std::remove("ptxt.json");

  // Create a ctxt
  helib::Ctxt ctxt(publicKey);

  // Encrypt `data` into the ciphertext
  ptxt.encrypt(ctxt);

  std::ofstream outCtxtFile;
  outCtxtFile.open("ctxt.json", std::ios::out);
  if (outCtxtFile.is_open()) {
    // Write the ctxt to a file
    ctxt.writeToJSON(outCtxtFile);
    // Close the ofstream
    outCtxtFile.close();
  } else {
    throw std::runtime_error("Could not open file 'ctxt.json'.");
  }

  std::ifstream inCtxtFile;
  inCtxtFile.open("ctxt.json", std::ios::in);
  if (inCtxtFile.is_open()) {
    // Read in the ctxt from the file
    helib::Ctxt deserializedCtxt =
        helib::Ctxt::readFromJSON(inCtxtFile, publicKey);
    // Note there are alternative methods for deserialization of Ctxt objects.
    // After initialization
    // helib::Ctxt deserializedCtxt(publicKey);
    // One can write
    // inCtxtFile >> deserializedCtxt;
    // Or alternatively
    // deserializedCtxt.readJSON(inCtxtFile);

    // Close the fstream
    inCtxtFile.close();
  } else {
    throw std::runtime_error("Could not open file 'ctxt.json'.");
  }
  // Remove the context file. Comment out the line below if you would like to
  // inspect the file.
  std::remove("ctxt.json");

  return 0;
}
