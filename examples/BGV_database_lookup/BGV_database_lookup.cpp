/* Copyright (C) 2020 IBM Corp.
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
#include <iostream>

#include <helib/helib.h>
#include <helib/EncryptedArray.h>

// Print the polynomial
void printPoly(NTL::ZZX& poly)
{
  for (int i = NTL::deg(poly); i >= 0; i--) {
    std::cout << poly[i] << "x^" << i;
    if (i > 0)
      std::cout << " + ";
    else
      std::cout << "\n";
  }
}

// This program is not a performant implementation but for teaching purposes
// only.

int main(int argc, char* argv[])
{
  /************ HElib boiler plate ************/

  // Note: The parameters have been chosen to provide fast running times as
  // opposed to a realistic security level.
  // Plaintext prime modulus
  unsigned long p = 131;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = 33;
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 1000;
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 2;

  std::cout << "Initialising context object..." << std::endl;
  // Intialise context
  helib::Context context(m, p, r);
  // Modify the context, adding primes to the modulus chain
  std::cout << "Building modulus chain..." << std::endl;
  helib::buildModChain(context, bits, c);

  // Print the context
  context.zMStar.printout();
  std::cout << std::endl;

  // Print the security level
  // Note: This will be negligible to improve performance time.
  std::cout << "Security: " << context.securityLevel() << std::endl;

  // Secret key management
  std::cout << "Creating secret key..." << std::endl;
  // Create a secret key associated with the context
  helib::SecKey secret_key = helib::SecKey(context);
  // Generate the secret key
  secret_key.GenSecKey();
  std::cout << "Generating key-switching matrices..." << std::endl;
  // Compute key-switching matrices that we need
  helib::addSome1DMatrices(secret_key);

  // Public key management
  // Set the secret key (upcast: FHESecKey is a subclass of FHEPubKey)
  const helib::PubKey& public_key = secret_key;

  // Get the EncryptedArray of the context
  const helib::EncryptedArray& ea = *(context.ea);

  // Get the number of slot (phi(m))
  long nslots = ea.size();
  std::cout << "\nNumber of slots: " << nslots << std::endl;

  /************ Create the database ************/

  std::vector<std::pair<std::string, std::string>> address_book = {
      {"Alice", "Amsterdam"},
      {"Bob", "London"},
      {"Charlie", "Liverpool"},
      {"Diana", "Washington"},
      {"Eve", "Paris"}};

  // Convert strings into numerical vectors
  std::vector<std::pair<helib::Ptxt<helib::BGV>, helib::Ptxt<helib::BGV>>>
      address_book_ptxt;
  for (const auto& name_addr_pair : address_book) {
    helib::Ptxt<helib::BGV> name(context);
    for (long i = 0; i < name_addr_pair.first.size(); ++i)
      name[i] = name_addr_pair.first[i];
    helib::Ptxt<helib::BGV> addr(context);
    for (long i = 0; i < name_addr_pair.second.size(); ++i)
      addr[i] = name_addr_pair.second[i];
    address_book_ptxt.emplace_back(std::move(name), std::move(addr));
  }

  // Encrypt the address book
  std::vector<std::pair<helib::Ctxt, helib::Ctxt>> encrypted_address_book;
  for (const auto& name_addr_pair : address_book_ptxt) {
    helib::Ctxt encrypted_name(public_key);
    helib::Ctxt encrypted_addr(public_key);
    public_key.Encrypt(encrypted_name, name_addr_pair.first);
    public_key.Encrypt(encrypted_addr, name_addr_pair.second);
    encrypted_address_book.emplace_back(encrypted_name, encrypted_addr);
  }

  /************ Create the query ************/

  // Read in query from the command line
  std::string query_string;
  std::cout << "\nPlease enter a name: ";
  std::cin >> query_string;
  std::cout << "Looking for the address of " << query_string << std::endl;

  // Convert query to a numerical vector
  helib::Ptxt<helib::BGV> query_ptxt(context);
  for (long i = 0; i < query_string.size(); ++i)
    query_ptxt[i] = query_string[i];

  // Encrypt the query
  helib::Ctxt query(public_key);
  public_key.Encrypt(query, query_ptxt);

  /************ Perform the database search ************/

  std::vector<helib::Ctxt> mask;
  mask.reserve(address_book.size());
  for (const auto& encrypted_pair : encrypted_address_book) {
    helib::Ctxt mask_entry = encrypted_pair.first; // Copy of database key
    mask_entry -= query;                           // Calculate the difference
    mask_entry.power(p - 1);                       // FLT
    mask_entry.negate();                           // Negate the ciphertext
    mask_entry.addConstant(NTL::ZZX(1));           // 1 - mask = 0 or 1
    // Create a vector of copies of the mask
    std::vector<helib::Ctxt> rotated_masks(ea.size(), mask_entry);
    for (int i = 1; i < rotated_masks.size(); i++)
      ea.rotate(rotated_masks[i], i);             // Rotate each of the masks
    totalProduct(mask_entry, rotated_masks);      // Multiply each of the masks
    mask_entry.multiplyBy(encrypted_pair.second); // multiply mask with values
    mask.push_back(mask_entry);
  }

  // Aggregate the results into a single ciphertext
  // Note: This code is for educational purposes and thus we try to refrain
  // from using the STL and do not use std::accumulate
  helib::Ctxt value = mask[0];
  for (int i = 1; i < mask.size(); i++)
    value += mask[i];

  /************ Decrypt and print result ************/

  helib::Ptxt<helib::BGV> plaintext_result(context);
  secret_key.Decrypt(plaintext_result, value);

  // Convert from ASCII to a string
  std::string string_result;
  for (long i = 0; i < plaintext_result.size(); ++i)
    string_result.push_back(static_cast<long>(plaintext_result[i]));

  std::cout << "\nQuery result: " << string_result << std::endl;

  return 0;
}
