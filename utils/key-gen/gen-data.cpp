/* Copyright (C) 2022 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
 */
 
// a file to generate random plaintexts for testing the key-gen pipeline.

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
  std::string ptxtCount = "1";
};

// Captures parameters for BGV
struct ParamsFileOpts
{
  long m = 0;
  long p = 0;
  long r = 0;
  long c = 0;
  long Qbits = 0;
};

int main(int argc, char* argv[])
{
  CmdLineOpts cmdLineOpts;

  // clang-format off
  helib::ArgMap()
        .toggle()
        .separator(helib::ArgMap::Separator::WHITESPACE)        
        .named()
          .arg("-o", cmdLineOpts.outputPrefixPath,
               "choose an output prefix path.", nullptr)
          .arg("-n", cmdLineOpts.ptxtCount, "the number of plaintexts to generate.")
        .required()
        .positional()
          .arg("<params-file>", cmdLineOpts.paramFileName,
               "the parameters file.", nullptr)          
        .parse(argc, argv);
  // clang-format on
  // If not set by user, returns params file name without extension and "params"
  if (cmdLineOpts.outputPrefixPath.empty()) {
    std::string prefix = stripExtension(cmdLineOpts.paramFileName);
    std::size_t params_pos = prefix.rfind("params");
    cmdLineOpts.outputPrefixPath = (params_pos == std::string::npos)
                                       ? prefix
                                       : prefix.substr(0, params_pos);
    std::cout << "File prefix: " << cmdLineOpts.outputPrefixPath << std::endl;
  }
  ParamsFileOpts paramsOpts;

  try {
    // clang-format off
    helib::ArgMap()
          .arg("p", paramsOpts.p, "require p.", "")
          .arg("m", paramsOpts.m, "require m.", "")
          .arg("r", paramsOpts.r, "require r.", "")
          .arg("c", paramsOpts.c, "require c.", "")
          .arg("Qbits", paramsOpts.Qbits, "require Q bits.", "")
          .parse(cmdLineOpts.paramFileName);
    // clang-format on
  } catch (const helib::RuntimeError& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }


  if (paramsOpts.p < 2) {
    std::cerr << "BGV invalid plaintext modulus. "
                 "In BGV it must be a prime number greater than 1."
              << std::endl;
    return EXIT_FAILURE;
  }

  try {
    const helib::Context* contextp = helib::ContextBuilder<helib::BGV>()
                                         .m(paramsOpts.m)
                                         .p(paramsOpts.p)
                                         .r(paramsOpts.r)
                                         .bits(paramsOpts.Qbits)
                                         .c(paramsOpts.c)
                                         .buildPtr();

    std::string ptxt_path = cmdLineOpts.outputPrefixPath + ".json";
    std::ofstream outPtxtFile(ptxt_path, std::ios::out);
    if (!outPtxtFile.is_open()) {
      throw std::runtime_error("Could not open file '" + ptxt_path + "'.");
    }
    // to call encrypt, need to have number of plaintexts at top of file
    outPtxtFile << stoi(cmdLineOpts.ptxtCount) << std::endl;
    for (int i = 0; i < stoi(cmdLineOpts.ptxtCount); i++) {
      helib::Ptxt<helib::BGV> ptxt(*contextp);
      ptxt.random();
      ptxt.writeToJSON(outPtxtFile);
      outPtxtFile << std::endl;
    }

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
