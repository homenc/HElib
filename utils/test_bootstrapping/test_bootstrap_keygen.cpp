/* Copyright (C) 2020 IBM Corp.
 * Copyright (C) 2022 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
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
  std::uniform_int_distribution<int> coinFlipDist(0, context.getP() - 1);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    ptxt[i] = coinFlipDist(gen);
  }
  return ptxt;
}

int main(int argc, char** argv)
{
  std::string sk_file_name;
  std::string pk_file_name;

  bool thick = false;
  bool quiet = false;
  bool sk_only = false;

  // clang-format off
  helib::ArgMap()
        .toggle()
          .arg("--thick", thick, "perform thick bootstrapping", nullptr)
          .arg("--quiet", quiet, "Suppress most of the output", nullptr)
          .arg("-s", sk_only, "only the secret key polynomial is written to the secret key file.")
        .positional()
        .required()
          .arg("<sk-file>", sk_file_name,
               "The file containing context and secret key.", nullptr)
          .arg("<pk-file>", pk_file_name, "The file containing context and evaluation public key.", nullptr)
        .parse(argc, argv);
  // clang-format on

  std::unique_ptr<helib::Context> pk_contextp;
  std::unique_ptr<helib::PubKey> pkp;

  try {
    std::tie(pk_contextp, pkp) = loadContextAndKey<helib::PubKey>(pk_file_name);
  } catch (const helib::RuntimeError& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  // Load Context and PubKey
  helib::Context& context(*pk_contextp);
  helib::PubKey& publicKey(*pkp);

  if (!context.isBootstrappable()) {
    std::cerr << "Context is not bootstrappable" << std::endl;
    return EXIT_FAILURE;
  }
  const helib::EncryptedArray& ea(context.getEA());
  helib::Ptxt<helib::BGV> ptxt =
      generateRandomPtxt(context); // Random in [0, p)
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);

  // Start doing bootstrapping
  if (!quiet) {
    std::cout << "Initial bootstrap" << std::endl;
  }
  if (thick) {
    publicKey.reCrypt(ctxt);
  } else {
    publicKey.thinReCrypt(ctxt);
  }
  if (!quiet) {
    helib::CheckCtxt(ctxt, "After recryption");
    std::cout << "ctxt.bitCapacity: " << ctxt.bitCapacity() << std::endl;
  }

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

  std::unique_ptr<helib::Context> sk_contextp;
  std::unique_ptr<helib::SecKey> skp;

  try {
    std::tie(sk_contextp, skp) =
        loadContextAndKey<helib::SecKey>(sk_file_name, sk_only);
  } catch (const helib::RuntimeError& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  // Load secret key context and secret key
  helib::SecKey& secretKey(*skp);
  helib::Context& sk_context(*sk_contextp);

  // to avoid a context mismatch, we write ctxt to stream and then read it back
  // with the secret key
  std::stringstream ss;
  ctxt.writeTo(ss);
  helib::Ctxt ctxt_copy = helib::Ctxt::readFrom(ss, secretKey);
  helib::Ptxt<helib::BGV> decrypted(sk_context);
  secretKey.Decrypt(decrypted, ctxt_copy);

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
