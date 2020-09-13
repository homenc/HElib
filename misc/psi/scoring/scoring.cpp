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
#include <helib/ArgMap.h>
#include <helib/Matrix.h>
#include <helib/partialMatch.h>
#include <helib/timing.h>

#include <psiio.h>

#if (defined(__unix__) || defined(__unix) || defined(unix))
  #include <sys/time.h>
  #include <sys/resource.h>
#endif

// Struct to hold command line arguments
struct CmdLineOpts
{
  std::string pkFilePath;
  std::string databaseFilePath;
  std::string queryFilePath;
  std::string outFilePath;
  bool isColumn = false;
  long nthreads = 1;
  long offset = 0;
};

int main(int argc, char *argv[])
{
  // PSI STUFF
  // 1. Read in context and pk file
  // 2. Python scripts to gen data (might use utils)
  // 3. Read in numbers (all numbers in a single row) and (all numbers in a single column)
  // 4. Create the query etc. Similar to TestPartial

  CmdLineOpts cmdLineOpts;

  // clang-format off
  helib::ArgMap()
    .required()
    .positional()
    .arg("<pkFile>", cmdLineOpts.pkFilePath, "Public Key file.", nullptr)
    .arg("<databaseFile>", cmdLineOpts.databaseFilePath,
      "Database file.", nullptr)
    .arg("<queryFile>", cmdLineOpts.queryFilePath, "Query file.", nullptr)
    .arg("<outFile>", cmdLineOpts.outFilePath, "Output file.", nullptr)
    .optional()
    .named()
    .arg("-n", cmdLineOpts.nthreads, "Number of threads.")
    .optional()
    .named()
    .arg("--offset", cmdLineOpts.offset, "Offset in bytes when writing to file.")
    .toggle()
    .arg("--column", cmdLineOpts.isColumn,
      "Flag to signify input is in column format.", nullptr)
    .parse(argc, argv);
  // clang-format on

  if (cmdLineOpts.nthreads < 1) {
    std::cerr << "Number of threads must be a postive integer. Setting n = 1." << std::endl;
    cmdLineOpts.nthreads = 1;
  }

  NTL::SetNumThreads(cmdLineOpts.nthreads);

  HELIB_NTIMER_START(readKey);
  // Load Context and PubKey
  std::shared_ptr<helib::Context> contextp;
  std::unique_ptr<helib::PubKey> pkp;
  std::tie(contextp, pkp) =
    loadContextAndKey<helib::PubKey>(cmdLineOpts.pkFilePath);
  HELIB_NTIMER_STOP(readKey);

  HELIB_NTIMER_START(readDatabase);
  // Read in database
  helib::Database<helib::Ctxt> database = readDbFromFile(cmdLineOpts.databaseFilePath, contextp, *pkp);
  HELIB_NTIMER_STOP(readDatabase);

  HELIB_NTIMER_START(readQuery);
  // Read in the query data
  helib::Matrix<helib::Ctxt> query = readQueryFromFile(cmdLineOpts.queryFilePath, *pkp);
  HELIB_NTIMER_STOP(readQuery);

  HELIB_NTIMER_START(buildQuery);
  // Build query info, Fs, mus, and taus
  // TODO: Read this in from the user rather than having a hard coded query
  // Fs = index_sets
  std::vector<std::vector<long>> Fs = {{0, 1, 2, 3}, {2, 3}, {0}, {1}};
  // mus = offsets
  std::vector<long> mus = {1, 5, 0, 1};
  // taus = weights
  std::vector<helib::Matrix<long>> taus = { {{0}, {7}, {0}, {1}},
                                            {{1}, {2}},
                                            {{1}},
                                            {{1}} };

  // Build the actual query
  helib::Query_t weighted_query(Fs, mus, taus, false);
  HELIB_NTIMER_STOP(buildQuery);

  HELIB_NTIMER_START(scoring);
  // Calculate scores
  // FIXME: Query currently must always be a row vector.
  auto clean = [](auto& x){x.cleanUp();};
  auto scores = database.getScore(weighted_query, query).apply(clean);
  HELIB_NTIMER_STOP(scoring);

  HELIB_NTIMER_START(writeResults);
  // Write results to file
  writeResultsToFile(cmdLineOpts.outFilePath, scores);
  HELIB_NTIMER_STOP(writeResults);

  std::ofstream timers("times.log");
  if (timers.is_open()) {
    helib::printAllTimers(timers);
  }

#if (defined(__unix__) || defined(__unix) || defined(unix))
  std::ofstream usage("usage.log");
  if (usage.is_open()) {
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    usage << "\n  rusage.ru_utime="<<rusage.ru_utime.tv_sec << '\n';
    usage << "  rusage.ru_stime="<<rusage.ru_utime.tv_sec << '\n';
    usage << "  rusage.ru_maxrss="<<rusage.ru_maxrss << '\n';
    usage << "  rusage.ru_minflt="<<rusage.ru_minflt << '\n';
    usage << "  rusage.ru_majflt="<<rusage.ru_majflt << '\n';
    usage << "  rusage.ru_inblock="<<rusage.ru_inblock << '\n';
    usage << "  rusage.ru_oublock="<<rusage.ru_majflt << '\n';
    usage << "  rusage.ru_nvcsw="<<rusage.ru_nvcsw << '\n';
    usage << "  rusage.ru_nivcsw="<<rusage.ru_minflt << std::endl;
  }
#endif

  return 0;
}
