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
#include <cstdio>
#include <cstdlib> // ldiv, system

#include <helib/helib.h>
#include <helib/ArgMap.h>

#include "Writer.h"
#include "common.h"

#include <NTL/BasicThreadPool.h>

struct CmdLineOpts
{
  std::string dataFilePath;
  std::string pkFilePath;
  std::string outFilePath;
  long batchSize = 0;
  long nthreads = 0; // Default is 0 for number of cpus.
  long offset = 0;
};

std::pair<long, long> parseDimsHeader(const std::string& s)
{
  std::stringstream iss(s);
  std::istream_iterator<long> issit(iss);
  std::vector<long> vl(issit, {});

  switch (vl.size()) {
  case 1:
    return {vl[0], 1};
  case 2:
    return {vl[0], vl[1]};
  default:
    std::ostringstream oss;
    oss << "Dimensions in header is wrong.\n";
    for (const auto& l : vl)
      oss << l << " ";
    throw std::runtime_error(oss.str());
  }
}

long estimateCtxtSize(const helib::Context& context, long offset)
{
  // Return in bytes.

  // We assume that the size of each element in the DCRT is BINIO_64BIT

  // sizeof(BINIO_EYE_CTXT_BEGIN) = 4;
  // BINIO_32BIT = 4
  // sizeof(long) = BINIO_64BIT = 8
  // xdouble = s * sizeof(long) = 2 * BINIO_64BIT = 16

  // We assume that primeSet after encryption is context.ctxtPrimes
  // We assume we have exactly 2 parts after encryption
  // We assume that the DCRT prime set is the same as the ctxt one

  long size = 0;

  size += 4;

  // Begin Ctxt metadata
  // 64 = header_size = ptxtSpace (long) + intFactor (long) + ptxtMag (xdouble)
  //                    + ratFactor (xdouble) + noiseBound (xdouble)
  size += 64;

  // primeSet.write(str);
  // size of set (long) + each prime (long)
  size += 8 + context.ctxtPrimes.card() * 8;

  // Begin Ctxt content size
  // write_raw_vector(str, parts);
  // Size of the parts vector (long)
  size += 8;

  long part_size = 0;
  // Begin CtxtPart size

  // skHandle.write(str);
  // powerOfS (long) + powerOfX (long) + secretKeyID (long)
  part_size += 24;

  // Begin DCRT size computation

  // this->DoubleCRT::write(str);
  // map.getIndexSet().write(str);
  // size of set (long) + each prime (long)
  part_size += 8 + context.ctxtPrimes.card() * 8;

  // DCRT data write as write_ntl_vec_long(str, map[i]);
  // For each prime in the ctxt modulus chain
  //    size of DCRT column (long) + size of each element (long) +
  //    size of all the slots (column in DCRT) (PhiM long elements)
  long dcrt_size =
      (8 + 8 * context.zMStar.getPhiM()) * context.ctxtPrimes.card();

  part_size += dcrt_size;

  // End DCRT size
  // End CtxtPart size

  size += 2 * part_size; // 2 * because we assumed 2 parts
  // End Ctxt content size

  // End eye-catcher
  size += 4;

  return size + offset;
}

template <typename SCHEME>
void encryptFromTo(const CmdLineOpts& cmdLineOpts,
                   const helib::Context& context,
                   const helib::PubKey& pk)
{

  // Read in data file header
  std::fstream dataFile(cmdLineOpts.dataFilePath);
  if (!dataFile.is_open()) {
    throw std::runtime_error("Could not open file '" +
                             cmdLineOpts.dataFilePath + "'.");
  }

  std::pair<long, long> dims = parseDimsHeader(readline(dataFile));

  // Here we 'batch'. Read in the batch size into memory.
  // Then, process with n threads. Repeat.
  // Write the header to file
  Writer<helib::Ctxt> writer(cmdLineOpts.outFilePath,
                             dims.first,
                             dims.second,
                             estimateCtxtSize(context, cmdLineOpts.offset));

  // Setting the stage
  std::vector<helib::Ptxt<SCHEME>> ptxts;
  std::vector<helib::Ctxt> ctxts;
  helib::Ctxt zero_ctxt(pk);
  helib::Ptxt<SCHEME> zero_ptxt(context);

  // Read in a batch
  for (long remaining = dims.first * dims.second, writtenBatches = 0;
       remaining > 0;
       remaining -= cmdLineOpts.batchSize, ++writtenBatches) {

    // Local batchSize
    long bsz =
        (remaining > cmdLineOpts.batchSize) ? cmdLineOpts.batchSize : remaining;
    ptxts.resize(bsz, zero_ptxt);
    ctxts.resize(bsz, zero_ctxt);

    std::vector<std::string> ptxt_strings(ptxts.size());
    for (std::size_t j = 0; j < ptxts.size(); j++) {
      std::getline(dataFile, ptxt_strings[j], '\n');
    }

    // Spread across n threads
    NTL_EXEC_RANGE(ctxts.size(), first, last)
    Writer<helib::Ctxt> threadWriter(writer);
    for (long i = first; i < last; ++i) {
      std::istringstream istr(ptxt_strings[i]);
      istr >> ptxts[i];
      pk.Encrypt(ctxts[i], ptxts[i]);
    }

    // Write to file
    for (long i = first; i < last; ++i) {
      if (dims.second == 1) {
        threadWriter.writeByLocation(ctxts[i],
                                     writtenBatches * cmdLineOpts.batchSize + i,
                                     0);
      } else {
        ldiv_t qr =
            ldiv(writtenBatches * cmdLineOpts.batchSize + i, dims.second);
        threadWriter.writeByLocation(ctxts[i], qr.quot, qr.rem);
      }
    }
    NTL_EXEC_RANGE_END
  }
}

int main(int argc, char* argv[])
{
  CmdLineOpts cmdLineOpts;

  // clang-format off
  helib::ArgMap()
    .required()
    .positional()
      .arg("<pk-file>", cmdLineOpts.pkFilePath,
           "the public key file.", nullptr)
      .arg("<input-file>", cmdLineOpts.dataFilePath,
           "the input file.", nullptr)
    .separator(helib::ArgMap::Separator::WHITESPACE)
    .named()
    .optional()
      .arg("-o", cmdLineOpts.outFilePath,
           "the output file name.", nullptr)
      .arg("-b", cmdLineOpts.batchSize,
           "batch size, how many ctxts in memory. If not set or 0 defaults to the number of threads used.")
      .arg("-n", cmdLineOpts.nthreads,
           "number of threads to use. If not set or 0 defaults to the number of concurrent threads supported.", "num. of cores")
      .arg("--offset", cmdLineOpts.offset,
           "byte packing offset in output file.")
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

  if (cmdLineOpts.nthreads > 0) {
    NTL::SetNumThreads(cmdLineOpts.nthreads);
  } else {
    std::cerr << "Number of threads must be a positive integer." << std::endl;
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

  // Set default outFilePath
  if (cmdLineOpts.outFilePath.empty()) {
    cmdLineOpts.outFilePath =
        stripExtension(cmdLineOpts.dataFilePath) + ".ctxt";
  }

  // Load Context and PubKey
  std::unique_ptr<helib::Context> contextp;
  std::unique_ptr<helib::PubKey> pkp;

  std::tie(contextp, pkp) =
      loadContextAndKey<helib::PubKey>(cmdLineOpts.pkFilePath);

  try {
    // Read in, encrypt, output.
    if (contextp->zMStar.getP() == -1) { // CKKS
      encryptFromTo<helib::CKKS>(cmdLineOpts, *contextp, *pkp);
    } else if (contextp->zMStar.getP() > 0) { // BGV
      encryptFromTo<helib::BGV>(cmdLineOpts, *contextp, *pkp);
    } else {
      std::cerr << "Unrecognized scheme from context." << std::endl;
      return EXIT_FAILURE;
    }
  } catch (const std::invalid_argument& e) {
    std::cerr << "Exit due to invalid argument thrown:\n"
              << e.what() << std::endl;
    std::remove(cmdLineOpts.outFilePath.c_str());
    return EXIT_FAILURE;
  } catch (const helib::IOError& e) {
    std::cerr << "Exit due to IOError thrown:\n" << e.what() << std::endl;
    std::remove(cmdLineOpts.outFilePath.c_str());
    return EXIT_FAILURE;
  } catch (const std::runtime_error& e) {
    std::cerr << "Exit due to runtime error thrown:\n" << e.what() << std::endl;
    std::remove(cmdLineOpts.outFilePath.c_str());
    return EXIT_FAILURE;
  } catch (const std::logic_error& e) {
    std::cerr << "Exit due to logic error thrown:\n" << e.what() << std::endl;
    std::remove(cmdLineOpts.outFilePath.c_str());
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "Exit due to unknown exception thrown:\n"
              << e.what() << std::endl;
    std::remove(cmdLineOpts.outFilePath.c_str());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
