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
#ifndef PSIIO_H
#define PSIIO_H

#include <iostream>
#include <fstream>

#include <helib/helib.h>
#include <helib/Matrix.h>

#include <common.h>
#include <Reader.h>
#include <Writer.h>

using sharedContext = std::shared_ptr<helib::Context>;

// Reads in encrypted database from file
inline helib::Database<helib::Ctxt> readDbFromFile(
                                          const std::string& databaseFilePath,
                                          const sharedContext& contextp,
                                          const helib::PubKey& pk)
{
  // Read in ctxt file header
  std::ifstream databaseFile(databaseFilePath);
  if (!databaseFile.is_open()) {
    throw std::runtime_error("Could not open file '" +
                             databaseFilePath + "'.");
  }

  helib::Ctxt zero_ctxt(pk);
  Reader<helib::Ctxt> reader(databaseFilePath, zero_ctxt);

  std::pair<long, long> dims = {reader.getTOC().getRows(),
                                reader.getTOC().getCols()};

  helib::Matrix<helib::Ctxt> data(zero_ctxt, dims.first, dims.second);

  // Read in ctxts
  if (dims.second == 1) { // Column vector
    for (long i = 0; i < dims.first * dims.second; ++i)
      reader.readDatum(data(i, 0), i, 0);
  } else {
    for (long i = 0; i < dims.first; ++i) {
      for (long j = 0; j < dims.second; ++j) {
        reader.readDatum(data(i, j), i, j);
      }
    }
  }

  return helib::Database<helib::Ctxt>(data, contextp);
}

// Reads in encrypted query from file
inline helib::Matrix<helib::Ctxt> readQueryFromFile(
                                             const std::string& queryFilePath,
                                             const helib::PubKey& pk)
{
  // Read in ctxt file header
  std::ifstream queryFile(queryFilePath);
  if (!queryFile.is_open()) {
    throw std::runtime_error("Could not open file '" +
                             queryFilePath + "'.");
  }

  helib::Ctxt zero_ctxt(pk);
  Reader<helib::Ctxt> reader(queryFilePath, zero_ctxt);

  std::pair<long, long> dims = {reader.getTOC().getRows(),
                                reader.getTOC().getCols()};

  if (dims.first != 1 && dims.second != 1) {
    std::cerr << "Query must be either a row or column vector." << std::endl;
    exit(EXIT_FAILURE);
  }

  helib::Matrix<helib::Ctxt> query(zero_ctxt, dims.first, dims.second);

  // Read in ctxts
  // FIXME: Always builds query as a row vector
  if (dims.first == 1) { // Read in row vector
    for (long i = 0; i < dims.first * dims.second; ++i)
      reader.readDatum(query(0, i), 0, i);
  } else if (dims.second == 1) { // Read in column vector
    for (long i = 0; i < dims.first * dims.second; ++i) {
      reader.readDatum(query(i, 0), i, 0);
    }
    query.transpose();
  } else {
    throw std::runtime_error("Trying to read in query that is not a vector. Dimensions " + std::to_string(dims.first) + " by " + std::to_string(dims.second) + ".");
  }

  return query;
}

// Writes out a matrix to file
inline void writeResultsToFile(const std::string& outFilePath,
                               const helib::Matrix<helib::Ctxt>& results,
                               const long offset = 0)
{
  Writer<helib::Ctxt> writer(outFilePath,
                             results.dims(0),
                             results.dims(1),
                             estimateCtxtSize(results(0,0).getContext(), offset));

  // Write the data
  NTL_EXEC_RANGE(results.dims(0) * results.dims(1), first, last)
  Writer<helib::Ctxt> threadWriter(writer);
  for (long i = first; i < last; ++i) {
    long row = i / results.dims(1);
    long col = i % results.dims(1);
    threadWriter.writeByLocation(results(row, col), row, col);
  }
  NTL_EXEC_RANGE_END
}

#endif
