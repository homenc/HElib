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

/* Copyright (C) 2022 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
 *
 * Extended existing functions found in helib/misc/psi/psiio/psiio.h to
 * work with Ptxt objects.
 */

#ifndef IO_H_
#define IO_H_

#include <helib/Matrix.h>
#include <helib/helib.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Reader.h"
#include "Writer.h"
#include "common.h"

using sharedContext = std::shared_ptr<helib::Context>;
using Ptxt = helib::Ptxt<helib::BGV>;

template <typename TXT>
helib::Database<TXT> readDbFromFile(const std::string& databaseFilePath,
                                    const sharedContext& contextp,
                                    const helib::PubKey& pk)
{
  // Read in TXT file header
  std::ifstream databaseFile(databaseFilePath);
  if (!databaseFile.is_open()) {
    throw std::runtime_error("Could not open file '" + databaseFilePath + "'.");
  }

  TXT zero_txt(pk);
  // This is only needed for TXT = Ctxt
  std::optional<Reader<TXT>> reader;
  long nrow, ncol;
  if constexpr (std::is_same_v<TXT, Ptxt>) {
    std::tie(nrow, ncol) = parseDimsHeader(readline(databaseFile));
  } else {
    reader.emplace(Reader<helib::Ctxt>(databaseFilePath, zero_txt));
    nrow = reader.value().getTOC().getRows();
    ncol = reader.value().getTOC().getCols();
  }

  helib::Matrix<TXT> data(zero_txt, nrow, ncol);

  if constexpr (std::is_same_v<TXT, Ptxt>) { // Ptxt query
    // Read in ptxts
    std::vector<std::string> ptxt_strings(nrow * ncol);
    for (auto& ptxt : ptxt_strings) {
      std::getline(databaseFile, ptxt, '\n');
    }
    // Populate Matrix
    for (long i = 0; i < nrow; ++i) {
      for (long j = 0; j < ncol; ++j) {
        std::istringstream iss(ptxt_strings[(i * ncol) + j]);
        iss >> data(i, j);
      }
    }
  } else { // Ctxt query
    NTL_EXEC_RANGE(nrow * ncol, first, last)
    Reader<TXT> threadReader(reader.value());
    for (long i = first; i < last; ++i) {
      long row = i / ncol;
      long col = i % ncol;
      threadReader.readDatum(data(row, col), row, col);
    }
    NTL_EXEC_RANGE_END
  }

  return helib::Database<TXT>(data, contextp);
}

template <typename TXT>
helib::Matrix<TXT> readQueryFromFile(const std::string& queryFilePath,
                                     const helib::PubKey& pk)
{
  // Read in TXT file header
  std::ifstream queryFile(queryFilePath);
  if (!queryFile.is_open()) {
    throw std::runtime_error("Could not open file '" + queryFilePath + "'.");
  }

  TXT zero_txt(pk);
  // This is only needed for TXT = Ctxt
  std::optional<Reader<TXT>> reader;
  long nrow, ncol;
  if constexpr (std::is_same_v<TXT, Ptxt>) { // Ptxt query
    std::tie(nrow, ncol) = parseDimsHeader(readline(queryFile));
  } else { // Ctxt query
    reader.emplace(Reader<helib::Ctxt>(queryFilePath, zero_txt));
    nrow = reader.value().getTOC().getRows();
    ncol = reader.value().getTOC().getCols();
  }

  if (nrow != 1 && ncol != 1) {
    throw std::runtime_error("Query must be either a row or column vector.\n");
  }

  helib::Matrix<TXT> query(zero_txt, nrow, ncol);

  // NOTE: Always builds query as a row vector
  if constexpr (std::is_same_v<TXT, Ptxt>) { // Ptxt query
    // Read in ptxts
    std::vector<std::string> ptxt_strings(nrow * ncol);
    for (auto& ptxt : ptxt_strings) {
      std::getline(queryFile, ptxt, '\n');
    }
    // Populate Matrix
    for (long i = 0; i < ptxt_strings.size(); ++i) {
      std::istringstream iss(ptxt_strings[i]);
      iss >> query(0, i);
    }
  } else { // Ctxt query
    // Read in ctxts
    NTL_EXEC_RANGE(nrow * ncol, first, last)
    Reader<TXT> threadReader(reader.value());
    for (long i = first; i < last; ++i) {
      long row = i / ncol;
      long col = i % ncol;
      threadReader.readDatum(query(row, col), row, col);
    }
    NTL_EXEC_RANGE_END
    if (ncol == 1) { // Transpose to make row vector
      query.transpose();
    }
  }

  return query;
}

// Writes out a matrix to file
template <typename TXT>
inline void writeResultsToFile(const std::string& outFilePath,
                               const helib::Matrix<TXT>& results,
                               const long offset = 0)
{
  if constexpr (std::is_same_v<TXT, Ptxt>) { // BGV Ptxt results
    // Open file
    std::ofstream outFile(outFilePath);
    if (!outFile.is_open()) {
      throw std::runtime_error("Could not open file '" + outFilePath + "'.");
    }

    long nrow = results.dims(0);
    long ncol = results.dims(1);

    // Write dims header
    if (ncol == 1) {
      outFile << nrow << std::endl;
    } else {
      outFile << nrow << " " << ncol << std::endl;
    }

    // Write the data
    // TODO(jlhcrawford): Is the result supposed to be row or col order?
    for (long i = 0; i < nrow; ++i) {
      for (long j = 0; j < ncol; ++j) {
        outFile << results(i, j) << std::endl;
      }
    }
  } else { // Ctxt results
    Writer<TXT> writer(outFilePath,
                       results.dims(0),
                       results.dims(1),
                       estimateCtxtSize(results(0, 0).getContext(), offset));

    // Write the data
    NTL_EXEC_RANGE(results.dims(0) * results.dims(1), first, last)
    Writer<TXT> threadWriter(writer);
    for (long i = first; i < last; ++i) {
      long row = i / results.dims(1);
      long col = i % results.dims(1);
      threadWriter.writeByLocation(results(row, col), row, col);
    }
    NTL_EXEC_RANGE_END
  }
}

#endif // IO_H_
