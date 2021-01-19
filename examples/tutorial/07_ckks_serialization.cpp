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

#include <iostream>

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

  // Print context to stdout
  std::cout << "*** CKKS context:\n";
  // Below we pretty print the JSON. If you do not wish to pretty print an
  // alternative is to call `context.writeJSON(std::cout);`
  std::cout << context.writeToJSON().pretty() << std::endl;

  // Create a secret key associated with the CKKS context
  helib::SecKey secret_key(context);

  // Generate the secret key
  secret_key.GenSecKey();

  // Compute key-switching matrices that we need
  helib::addSome1DMatrices(secret_key);

  // Print secret key to stdout
  std::cout << "\n\n*** Secret Key:\n";
  std::cout << secret_key.writeToJSON().pretty() << std::endl;

  // Create an alias public key part from the secret key
  const helib::PubKey& public_key = secret_key;

  // Print the public key to stdout
  std::cout << "\n\n*** Public Key:\n";
  std::cout << public_key.writeToJSON().pretty() << std::endl;

  // Create a data object
  std::vector<long> data(context.getNSlots());

  // Generate some data
  std::iota(data.begin(), data.end(), 0);

  // Create a PtxtArray.
  helib::PtxtArray ptxt(context, data);

  // Print the ptxt to stdout
  std::cout << "\n\n*** Ptxt:\n";
  std::cout << ptxt.writeToJSON().pretty() << std::endl;

  // Create a ctxt
  helib::Ctxt ctxt(public_key);

  // Encrypt `data` into the ciphertext
  ptxt.encrypt(ctxt);

  // Print the ctxt to stdout
  std::cout << "\n\n*** Ctxt:\n";
  std::cout << ctxt.writeToJSON().pretty() << std::endl;

  return 0;
}
