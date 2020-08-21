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
#include <fstream>

#include <helib/helib.h>
#include <helib/ArgMap.h>

#include <NTL/BasicThreadPool.h>

#include "Reader.h"
#include "common.h"

struct CmdLineOpts
{
  std::string skFilePath;
  std::string ctxtFilePath;
  std::string outFilePath;
  long batchSize = 0;
  long nthreads = 0; // Default is 0 for number of cpus.
};

void writeDimsHeader(std::ostream& os, std::pair<long, long>& dims)
{
  if (dims.second == 1) {
    os << dims.first << std::endl;
  } else {
    os << dims.first << " " << dims.second << std::endl;
  }
}

template <typename SCHEME>
void decryptFromTo(const CmdLineOpts& cmdLineOpts,
                   const helib::Context& context,
                   const helib::SecKey& sk)
{
  // Read in data file header
  std::ifstream dataFile(cmdLineOpts.ctxtFilePath);
  if (!dataFile.is_open()) {
    throw std::runtime_error("Could not open file '" +
                             cmdLineOpts.ctxtFilePath + "'.");
  }

  // Read in a batch
  std::vector<helib::Ptxt<SCHEME>> ptxts;
  std::vector<helib::Ctxt> ctxts;
  helib::Ctxt zero_ctxt(sk);
  helib::Ptxt<SCHEME> zero_ptxt(context);

  Reader<helib::Ctxt> reader(cmdLineOpts.ctxtFilePath, zero_ctxt);

  std::pair<long, long> dims = {reader.getTOC().getRows(),
                                reader.getTOC().getCols()};

  std::ofstream outFile;
  std::ostream* out;

  if (!cmdLineOpts.outFilePath.empty()) {
    outFile.open(cmdLineOpts.outFilePath);
    if (!outFile.is_open()) {
      throw std::runtime_error("Could not open file '" +
                               cmdLineOpts.outFilePath + "'.");
    }
    out = &outFile;
  } else {
    out = &std::cout;
  }

  writeDimsHeader(*out, dims);

  for (long remaining = dims.first * dims.second, readBatches = 0;
       remaining > 0;
       remaining -= cmdLineOpts.batchSize, ++readBatches) {

    // Read in a batch
    long bsz =
        (remaining > cmdLineOpts.batchSize) ? cmdLineOpts.batchSize : remaining;
    ptxts.resize(bsz, zero_ptxt);
    ctxts.resize(bsz, zero_ctxt);

    NTL_EXEC_RANGE(ctxts.size(), first, last)
    Reader<helib::Ctxt> threadReader(reader);
    for (long i = first; i < last; ++i) {
      if (dims.second == 1) {
        threadReader.readDatum(ctxts[i],
                               readBatches * cmdLineOpts.batchSize + i,
                               0);
      } else {
        ldiv_t qr = ldiv(readBatches * cmdLineOpts.batchSize + i, dims.second);
        threadReader.readDatum(ctxts[i], qr.quot, qr.rem);
      }
    }

    // Decrypt using NTL threads
    for (long i = first; i < last; ++i) {
      sk.Decrypt(ptxts[i], ctxts[i]);
    }
    NTL_EXEC_RANGE_END

    // Write out to stream
    for (const auto& ptxt : ptxts)
      *out << ptxt << std::endl;
  }
}

int main(int argc, char* argv[])
{
  CmdLineOpts cmdLineOpts;

  // clang-format off
  helib::ArgMap()
    .required()
    .positional()
      .arg("<sk>", cmdLineOpts.skFilePath,
           "secret key.", nullptr)
      .arg("<ctxt>", cmdLineOpts.ctxtFilePath,
           "path of the ciphertext to read in.", nullptr)
    .separator(helib::ArgMap::Separator::WHITESPACE)
    .named()
    .optional()
      .arg("-o", cmdLineOpts.outFilePath,
           "the output file name.", nullptr)
      .arg("-b", cmdLineOpts.batchSize,
           "batch size, how many ctxts in memory. If not set or 0 defaults to the number of threads used.")
      .arg("-n", cmdLineOpts.nthreads,
           "number of threads to use. If not set or 0 defaults to the number of concurrent threads supported.", "num. of cores")
    .parse(argc, argv);
  // clang-format on

  // Set NTL nthreads
  if (cmdLineOpts.nthreads == 0) {
    cmdLineOpts.nthreads = std::thread::hardware_concurrency();
    // hardware_concurrency may still return 0 if not supported on
    // implementation.
    if (cmdLineOpts.nthreads == 0) {
      std::cerr << "C++ `hardware_concurrency` call not available on this"
                   "platform.\nnthreads must be explicitly provided."
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (cmdLineOpts.nthreads > 0)
    NTL::SetNumThreads(cmdLineOpts.nthreads);
  else {
    std::cerr << "Number of threads must a be positive integer." << std::endl;
    return EXIT_FAILURE;
  }

  // Set default batch size.
  if (cmdLineOpts.batchSize == 0) {
    cmdLineOpts.batchSize = cmdLineOpts.nthreads;
  }

  // Check batch size above 1.
  if (cmdLineOpts.batchSize < 1) {
    std::cerr << "Batch size must be a positive integer." << std::endl;
    return EXIT_FAILURE;
  }

  // Warn if batch size is less than thread count
  if (cmdLineOpts.batchSize < cmdLineOpts.nthreads) {
    std::cerr << "WARNING: Not enough elements in the batch ("
              << cmdLineOpts.batchSize << ") to run with "
              << cmdLineOpts.nthreads << " threads." << std::endl;
    std::cerr << "To achieve better performance set batch size (-b) to be "
                 "equal to or more than the number of threads (-n). Ideally "
                 "the batch size should be a multiple of the number of threads."
              << std::endl;
  }

  // Load Context and SecKey
  std::unique_ptr<helib::Context> contextp;
  std::unique_ptr<helib::SecKey> skp;

  std::tie(contextp, skp) =
      loadContextAndKey<helib::SecKey>(cmdLineOpts.skFilePath);

  // Read in, decrypt, output.
  try {
    if (contextp->zMStar.getP() == -1) {
      decryptFromTo<helib::CKKS>(cmdLineOpts, *contextp, *skp);
    } else if (contextp->zMStar.getP() > 0) {
      decryptFromTo<helib::BGV>(cmdLineOpts, *contextp, *skp);
    } else {
      std::cerr << "Unrecognized scheme from context." << std::endl;
      return EXIT_FAILURE;
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
