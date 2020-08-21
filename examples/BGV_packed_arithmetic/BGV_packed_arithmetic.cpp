/* Copyright (C) 2019-2020 IBM Corp.
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

// This is a sample program for education purposes only.
// It attempts to show the various basic mathematical
// operations that can be performed on both ciphertexts
// and plaintexts.

#include <iostream>

#include <helib/helib.h>

int main(int argc, char* argv[])
{
  /*  Example of BGV scheme  */

  // Plaintext prime modulus
  unsigned long p = 4999;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = 32109;
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 500;
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 2;

  std::cout << "\n*********************************************************";
  std::cout << "\n*         Basic Mathematical Operations Example         *";
  std::cout << "\n*         =====================================         *";
  std::cout << "\n*                                                       *";
  std::cout << "\n* This is a sample program for education purposes only. *";
  std::cout << "\n* It attempts to show the various basic mathematical    *";
  std::cout << "\n* operations that can be performed on both ciphertexts  *";
  std::cout << "\n* and plaintexts.                                       *";
  std::cout << "\n*                                                       *";
  std::cout << "\n*********************************************************";
  std::cout << std::endl;

  std::cout << "Initialising context object..." << std::endl;
  // Initialize context
  // This object will hold information about the algebra created from the
  // previously set parameters
  helib::Context context(m, p, r);
  // Modify the context, adding primes to the modulus chain
  // This defines the ciphertext space
  std::cout << "Building modulus chain..." << std::endl;
  buildModChain(context, bits, c);

  // Print the context
  context.zMStar.printout();
  std::cout << std::endl;

  // Print the security level
  std::cout << "Security: " << context.securityLevel() << std::endl;

  // Secret key management
  std::cout << "Creating secret key..." << std::endl;
  // Create a secret key associated with the context
  helib::SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();
  std::cout << "Generating key-switching matrices..." << std::endl;
  // Compute key-switching matrices that we need
  helib::addSome1DMatrices(secret_key);

  // Public key management
  // Set the secret key (upcast: SecKey is a subclass of PubKey)
  const helib::PubKey& public_key = secret_key;

  // Get the EncryptedArray of the context
  const helib::EncryptedArray& ea = *(context.ea);

  // Get the number of slot (phi(m))
  long nslots = ea.size();
  std::cout << "Number of slots: " << nslots << std::endl;

  // Create a vector of long with nslots elements
  helib::Ptxt<helib::BGV> ptxt(context);
  // Set it with numbers 0..nslots - 1
  // ptxt = [0] [1] [2] ... [nslots-2] [nslots-1]
  for (int i = 0; i < ptxt.size(); ++i) {
    ptxt[i] = i;
  }

  // Print the plaintext
  std::cout << "Initial Plaintext: " << ptxt << std::endl;

  // Create a ciphertext object
  helib::Ctxt ctxt(public_key);
  // Encrypt the plaintext using the public_key
  public_key.Encrypt(ctxt, ptxt);

  /********** Operations **********/
  // Ciphertext and plaintext operations are performed
  // "entry-wise".

  // Square the ciphertext
  // [0] [1] [2] [3] [4] ... [nslots-1]
  // -> [0] [1] [4] [9] [16] ... [(nslots-1)*(nslots-1)]
  ctxt.multiplyBy(ctxt);
  // Plaintext version
  ptxt.multiplyBy(ptxt);

  // Divide the ciphertext by itself
  // To do this we must calculate the multiplicative inverse using Fermat's
  // Little Theorem.  We calculate a^{-1} = a^{p-2} mod p, where a is non-zero
  // and p is our plaintext prime.
  // First make a copy of the ctxt using copy constructor
  helib::Ctxt ctxt_divisor(ctxt);
  // Raise the copy to the exponent p-2
  // [0] [1] [4] ... [16] -> [0] [1] [1] ... [1]
  // Note: 0 is a special case because 0^n = 0 for any power n
  ctxt_divisor.power(p - 2);
  // a^{p-2}*a = a^{-1}*a = a / a = 1;
  ctxt.multiplyBy(ctxt_divisor);

  // Plaintext version
  helib::Ptxt<helib::BGV> ptxt_divisor(ptxt);
  ptxt_divisor.power(p - 2);
  ptxt.multiplyBy(ptxt_divisor);

  // Double it (using additions)
  // [0] [1] [1] ... [1] [1] -> [0] [2] [2] ... [2] [2]
  ctxt += ctxt;
  // Plaintext version
  ptxt += ptxt;

  // Subtract it from itself (result should be 0)
  // i.e. [0] [0] [0] [0] ... [0] [0]
  ctxt -= ctxt;
  // Plaintext version
  ptxt -= ptxt;

  // Create a plaintext for decryption
  helib::Ptxt<helib::BGV> plaintext_result(context);
  // Decrypt the modified ciphertext
  secret_key.Decrypt(plaintext_result, ctxt);

  // Print the decrypted plaintext
  // Should be [0] [0] [0] ... [0] [0]
  std::cout << "Decrypted Plaintext: " << plaintext_result << std::endl;
  // Print the plaintext version result, should be the same as the ctxt version
  std::cout << "Plaintext version: " << ptxt << std::endl;

  // We can also add constants
  // [0] [0] [0] ... [0] [0] -> [1] [1] [1] ... [1] [1]
  ctxt.addConstant(NTL::ZZX(1l));
  // Plaintext version
  ptxt.addConstant(NTL::ZZX(1l));

  // And multiply by constants
  // [1] [1] [1] ... [1] [1]
  // -> [1*1] [1*1] [1*1] ... [1*1] [1*1] = [1] [1] [1] ... [1] [1]
  ctxt *= NTL::ZZX(1l);
  // Plaintext version
  ptxt *= NTL::ZZX(1l);

  // We can also perform ciphertext-plaintext operations
  // ctxt = [1] [1] [1] ... [1] [1], ptxt = [1] [1] [1] ... [1] [1]
  // ctxt + ptxt = [2] [2] [2] ... [2] [2]
  // Note: the output of this is also a ciphertext
  ctxt += ptxt;

  // Decrypt the modified ciphertext into a new plaintext
  helib::Ptxt<helib::BGV> new_plaintext_result(context);
  secret_key.Decrypt(new_plaintext_result, ctxt);

  // Print the decrypted plaintext
  // Should be [2] [2] [2] ... [2] [2]
  std::cout << "Decrypted Plaintext: " << new_plaintext_result << std::endl;

  return 0;
}
