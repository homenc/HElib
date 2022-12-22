/* Copyright (C) 2020 IBM Corp.
 * Copyright (C) 2022 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
 */

#include <fstream>
#include <ctime>

#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/debugging.h>

#include <NTL/BasicThreadPool.h>

#include "common.h"

struct CmdLineOpts
{
  std::string paramFileName;
  std::string outputPrefixPath;
  bool write_only_sk = 0; // Default to false for backward compatibility.
  std::string scheme = "BGV";
  std::string bootstrappable = "NONE"; // NONE | THIN | FAT
  bool noSKM = false;
  bool frobSKM = false;
  bool infoFile = false;
};

// Captures parameters of both BGV and CKKS
struct ParamsFileOpts
{
  long m = 0;
  long p = 0;
  long r = 0;
  long c = 0;
  long Qbits = 0;
  long scale = 4;
  long c_m = 100;
  NTL::Vec<long> mvec;
  NTL::Vec<long> gens;
  NTL::Vec<long> ords;
};

// Write context.printout to file.out
void printoutToStream(const helib::Context& context,
                      std::ostream& out,
                      bool noSKM,
                      bool frobSKM,
                      bool bootstrappable)
{
  if (!noSKM || bootstrappable)
    out << "Key switching matrices created.\n";
  if (frobSKM || bootstrappable)
    out << "Frobenius matrices created.\n";
  if (bootstrappable)
    out << "Recrypt data created.\n";

  // write the algebra info
  context.printout(out);
}

// sk is child of pk in HElib.
void writeKeyToFile(std::string& pathPrefix,
                    helib::Context& context,
                    helib::SecKey& secretKey,
                    bool pkNotSk,
                    bool write_only_sk = false)
{
  std::string path = pathPrefix + (pkNotSk ? ".pk" : ".sk");
  std::ofstream keysFile(path, std::ios::binary);
  if (!keysFile.is_open()) {
    std::runtime_error("Cannot write keys to file at '" + path);
  }

  // write the context
  context.writeTo(keysFile);

  // write the keys
  if (pkNotSk) {
    const helib::PubKey& pk = secretKey;
    pk.writeTo(keysFile);
  } else {
    secretKey.writeTo(keysFile, write_only_sk);
  }
}

int main(int argc, char* argv[])
{
  CmdLineOpts cmdLineOpts;

  // clang-format off
  helib::ArgMap()
        .toggle()
         .arg("--no-skm", cmdLineOpts.noSKM,
               "disable switch-key matrices.", nullptr)
         .arg("--frob-skm", cmdLineOpts.frobSKM,
               "generate Frobenius switch-key matrices.", nullptr)
         .arg("--info-file", cmdLineOpts.infoFile,
               "print algebra info to file.", nullptr)
         .arg("-s",cmdLineOpts.write_only_sk,"write only the secret key polynomial to the secret key file.")
        .separator(helib::ArgMap::Separator::WHITESPACE)
        .named()
          .arg("--scheme", cmdLineOpts.scheme,
               "choose scheme BGV | CKKS.")
          .arg("-o", cmdLineOpts.outputPrefixPath,
               "choose an output prefix path.", nullptr)
          .arg("--bootstrap", cmdLineOpts.bootstrappable,
               "choose boostrapping option NONE | THIN | THICK.")
        .required()
        .positional()
          .arg("<params-file>", cmdLineOpts.paramFileName,
               "the parameters file.", nullptr)
        .parse(argc, argv);
  // clang-format on

  ParamsFileOpts paramsOpts;

  try {
    // clang-format off
    helib::ArgMap()
          .arg("p", paramsOpts.p, "require p.", "")
          .arg("m", paramsOpts.m, "require m.", "")
          .arg("r", paramsOpts.r, "require r.", "")
          .arg("c", paramsOpts.c, "require c.", "")
          .arg("Qbits", paramsOpts.Qbits, "require Q bits.", "")
          .optional()
            .arg("scale", paramsOpts.scale, "require scale for CKKS")
            .arg("c_m", paramsOpts.c_m, "require c_m for bootstrapping.", "")
            .arg("mvec", paramsOpts.mvec, "require mvec for bootstrapping.", "")
            .arg("gens", paramsOpts.gens, "require gens for bootstrapping.", "")
            .arg("ords", paramsOpts.ords, "require ords for bootstrapping.", "")
          .parse(cmdLineOpts.paramFileName);
    // clang-format on
  } catch (const helib::RuntimeError& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  // Create the FHE context
  long p;
  if (cmdLineOpts.scheme.empty() || cmdLineOpts.scheme == "BGV") {
    if (paramsOpts.p < 2) {
      std::cerr << "BGV invalid plaintext modulus. "
                   "In BGV it must be a prime number greater than 1."
                << std::endl;
      return EXIT_FAILURE;
    }
    p = paramsOpts.p;
  } else if (cmdLineOpts.scheme == "CKKS") {
    if (paramsOpts.p != -1) {
      std::cerr << "CKKS invalid plaintext modulus. "
                   "In CKKS it must be set to -1."
                << std::endl;
      return EXIT_FAILURE;
    }
    p = -1;
    if (cmdLineOpts.bootstrappable != "NONE") {
      std::cerr << "CKKS does not currently support bootstrapping."
                << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    std::cerr << "Unrecognized scheme '" << cmdLineOpts.scheme << "'."
              << std::endl;
    return EXIT_FAILURE;
  }

  if (cmdLineOpts.noSKM && cmdLineOpts.frobSKM) {
    std::cerr << "Frobenius matrices reqires switch-key matrices to be "
                 "generated."
              << std::endl;
    return EXIT_FAILURE;
  }

  if (cmdLineOpts.bootstrappable != "NONE" &&
      cmdLineOpts.bootstrappable != "THIN" &&
      cmdLineOpts.bootstrappable != "THICK") {
    std::cerr << "Bad boostrap option: " << cmdLineOpts.bootstrappable
              << ".  Allowed options are NONE, THIN, THICK." << std::endl;
    return EXIT_FAILURE;
  }

  if (cmdLineOpts.bootstrappable != "NONE") {
    if (cmdLineOpts.noSKM) {
      std::cerr << "Cannot generate bootstrappable context without switch-key "
                   "and frobenius matrices."
                << std::endl;
      return EXIT_FAILURE;
    }
    if (paramsOpts.mvec.length() == 0) {
      std::cerr << "Missing mvec parameter for bootstrapping in "
                << cmdLineOpts.paramFileName << "." << std::endl;
      return EXIT_FAILURE;
    }
    if (paramsOpts.gens.length() == 0) {
      std::cerr << "Missing gens parameter for bootstrapping in "
                << cmdLineOpts.paramFileName << "." << std::endl;
      return EXIT_FAILURE;
    }
    if (paramsOpts.ords.length() == 0) {
      std::cerr << "Missing ords parameter for bootstrapping in "
                << cmdLineOpts.paramFileName << "." << std::endl;
      return EXIT_FAILURE;
    }
  }

  try {
    helib::Context* contextp;

    if (cmdLineOpts.scheme == "BGV") {
      helib::ContextBuilder<helib::BGV> cb;
      cb.m(paramsOpts.m)
          .p(p)
          .r(paramsOpts.r)
          .gens(helib::convert<std::vector<long>>(paramsOpts.gens))
          .ords(helib::convert<std::vector<long>>(paramsOpts.ords))
          .bits(paramsOpts.Qbits)
          .c(paramsOpts.c);

      if (cmdLineOpts.bootstrappable != "NONE") {
        cb.bootstrappable(true).mvec(paramsOpts.mvec);
        if (cmdLineOpts.bootstrappable == "THICK")
          cb.thickboot();
        else if (cmdLineOpts.bootstrappable == "THIN")
          cb.thinboot();
      }
      contextp = cb.buildPtr();
    } else if (cmdLineOpts.scheme == "CKKS") {
      contextp = helib::ContextBuilder<helib::CKKS>()
                     .m(paramsOpts.m)
                     .precision(paramsOpts.r)
                     .bits(paramsOpts.Qbits)
                     .scale(paramsOpts.scale)
                     .buildPtr();
    }

    // and a new secret/public key
    helib::SecKey secretKey(*contextp);
    secretKey.GenSecKey(); // A +-1/0 secret key

    // If not set by user, returns params file name with truncated UTC
    if (cmdLineOpts.outputPrefixPath.empty()) {
      cmdLineOpts.outputPrefixPath =
          stripExtension(cmdLineOpts.paramFileName) +
          std::to_string(std::time(nullptr) % 100000);
      std::cout << "File prefix: " << cmdLineOpts.outputPrefixPath << std::endl;
    }

    // Printout important info
    std::ostream* outp = &std::cout;
    std::ofstream fout;
    if (cmdLineOpts.infoFile) {
      // outputPrefixPath should be set further up main.
      std::string path = cmdLineOpts.outputPrefixPath + ".info";
      fout.open(path);
      if (!fout.is_open()) {
        throw std::runtime_error("Cannot write keys to file at '" + path +
                                 "'.");
      }
      outp = &fout;
    }

    printoutToStream(*contextp,
                     *outp,
                     cmdLineOpts.noSKM,
                     cmdLineOpts.frobSKM,
                     cmdLineOpts.bootstrappable != "NONE");

    NTL::SetNumThreads(2);
    // write the secret key to <outputPrefixPath>.sk and the encryption public
    // key to <outputPrefixPath>Enc.pk
    NTL_EXEC_INDEX(2, skOrPk)
    std::string path = cmdLineOpts.outputPrefixPath + (skOrPk ? "Enc" : "");
    writeKeyToFile(path,
                   *contextp,
                   secretKey,
                   skOrPk,
                   cmdLineOpts.write_only_sk);
    NTL_EXEC_INDEX_END
    // now compute the evaluation key material

    // compute key-switching matrices
    if (!cmdLineOpts.noSKM || cmdLineOpts.bootstrappable != "NONE") {
      helib::addSome1DMatrices(secretKey);
      if (cmdLineOpts.frobSKM || cmdLineOpts.bootstrappable != "NONE") {
        helib::addFrbMatrices(secretKey);
      }
    }

    if (cmdLineOpts.bootstrappable != "NONE") {
      secretKey.genRecryptData();
    }
    // and write to file <outputPrefixPath>Eval.pk
    std::string path = cmdLineOpts.outputPrefixPath + "Eval";
    writeKeyToFile(path, *contextp, secretKey, 1);

  } catch (const std::invalid_argument& e) {
    std::cerr << "Exit due to invalid argument thrown:\n"
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const helib::IOError& e) {
    std::cerr << "Exit due to IOError thrown:\n" << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::runtime_error& e) {
    std::cerr << "Exit due to runtime error thrown:\n" << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::logic_error& e) {
    std::cerr << "Exit due to logic error thrown:\n" << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "Exit due to unknown exception thrown:\n"
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
