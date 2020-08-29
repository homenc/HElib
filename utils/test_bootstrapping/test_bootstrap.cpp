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
#include <helib/helib.h>
#include <helib/debugging.h>

#include <helib/ArgMap.h>

#include <random>

#include <common.h>

helib::Ptxt<helib::BGV> generateRandomPtxt(const helib::Context& context)
{
  // Note: this function uses a single coefficient in each slot.  For this
  // reason there is no difference between THIN and THICK bootstrapping.
  // Consider generating random polynomials if THICK
  helib::Ptxt<helib::BGV> ptxt(context);
  std::mt19937 gen(231087);
  std::uniform_int_distribution<int> coinFlipDist(0, context.zMStar.getP() - 1);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    ptxt[i] = coinFlipDist(gen);
  }
  return ptxt;
}

int main(int argc, char** argv)
{
  std::string sk_file_name;

  bool thick = false;
  bool quiet = false;

  // clang-format off
  helib::ArgMap()
        .toggle()
          .arg("--thick", thick, "perform thick bootstrapping", nullptr)
          .arg("--quiet", quiet, "Suppress most of the output", nullptr)
        .positional()
        .required()
          .arg("<sk-file>", sk_file_name,
               "The file containing context and secret key.", nullptr)
        .parse(argc, argv);
  // clang-format on

  std::unique_ptr<helib::Context> contextp;
  std::unique_ptr<helib::SecKey> skp;

  try {
    std::tie(contextp, skp) = loadContextAndKey<helib::SecKey>(sk_file_name);
  } catch (const helib::RuntimeError& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  // Load Context and PubKey
  helib::Context& context(*contextp);
  helib::SecKey& secretKey(*skp);
  helib::PubKey& publicKey(secretKey);

  if (!context.isBootstrappable()) {
    std::cerr << "Context is not bootstrappable" << std::endl;
    return EXIT_FAILURE;
  }
  const helib::EncryptedArray& ea(*context.ea);
  helib::Ptxt<helib::BGV> ptxt =
      generateRandomPtxt(context); // Random in [0, p)
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);

  // Start doing bootstrapping
  if (!quiet) {
    std::cout << "Initial bootstrap" << std::endl;
  }
  if (thick) {
    secretKey.reCrypt(ctxt);
  } else {
    secretKey.thinReCrypt(ctxt);
  }
  if (!quiet) {
    helib::CheckCtxt(ctxt, "After recryption");
    std::cout << "ctxt.bitCapacity: " << ctxt.bitCapacity() << std::endl;
  }

  int round = 0;
  for (std::size_t round = 0; round < 2; ++round) {
    if (!quiet) {
      std::cout << "Round: " << round << std::endl;
      helib::CheckCtxt(ctxt, "Before multiplication");
    }
    // Multiply the ciphertext with itself n times
    // until number of bits falls below threshold
    HELIB_NTIMER_START(RoundTotal);
    while (ctxt.bitCapacity() >= ctxt.getContext().BPL() * 3) {
      long count = 0; // number of multiplications
      if (!quiet) {
        std::cout << "multiplication: " << count++ << std::endl;
        std::cout << "ctxt.bitCapacity: " << ctxt.bitCapacity() << std::endl;
      }
      HELIB_NTIMER_START(Multiplications);
      ctxt.square();
      HELIB_NTIMER_STOP(Multiplications);
      ptxt.square();
    }

    if (!quiet) {
      helib::CheckCtxt(ctxt, "Before recryption");
      std::cout << "ctxt.bitCapacity: " << ctxt.bitCapacity() << std::endl;
    }

    // Time the recryption step
    HELIB_NTIMER_START(Bootstrap);
    // Recrypt/Bootstrap the ctxt
    if (thick) {
      publicKey.reCrypt(ctxt);
    } else {
      publicKey.thinReCrypt(ctxt);
    }
    HELIB_NTIMER_STOP(Bootstrap);

    if (!quiet) {
      helib::CheckCtxt(ctxt, "After recryption");
      std::cout << "ctxt.bitCapacity " << ctxt.bitCapacity() << std::endl;
    }
    HELIB_NTIMER_STOP(RoundTotal);

    // Plaintext operation
    // Multiply with itself n times
  }
  helib::Ptxt<helib::BGV> decrypted(context);
  secretKey.Decrypt(decrypted, ctxt);

  if (ptxt != decrypted) {
    std::cerr << "Decryption error." << std::endl;
    return EXIT_FAILURE;
  }

  if (!quiet) {
    for (const auto& timerName :
         {"Setup", "Multiplications", "Bootstrap", "RoundTotal"}) {
      helib::printNamedTimer(std::cout << std::endl, timerName);
    }
  }
}
