/* Copyright (C) 2020-2021 IBM Corp.
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

#include <helib/log.h>
#include "test_common.h"
#include "gtest/gtest.h"

namespace {

// Find matching line in file. Match substr. Currently, just scans file.
// Fine as for test log files will be relatively small.
static bool findMatchingLine(std::istream& is, const std::string& str)
{
  std::string line;
  while (std::getline(is, line)) {
    if (line.find(str) != std::string::npos) {
      return true;
    }
  }
  return false;
}

// Similar to findMatchingLine above however it looks for s1 first and if found
// continues to look for s2. so s1 comes before s2.
static bool findMatchingLines(std::istream& is,
                              const std::string& s1,
                              const std::string& s2)
{
  if (findMatchingLine(is, s1)) {
    return findMatchingLine(is, s2);
  }

  return false;
}

// First we run the tests that do not change the default `ostream`.
// This is to test the default behaviour.

TEST(TestLogging, testWarningToDefaultLogfile)
{
  const char* filepath = "helib.log";

  const std::string warningMsg = "Warning message 900!";
  helib::Warning(warningMsg);
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);
  EXPECT_TRUE(findMatchingLine(filestream, warningMsg));

  // Should not really delete the default log file
}

// Now we run the tests below that explicitly set the stream.

TEST(TestLogging, testWarningToSetLogfile)
{
  const char* filepath = "other.log";

  helib::helog.setLogToFile(filepath);

  const std::string warningMsg = "Warning message 700!";
  helib::Warning(warningMsg);
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);
  EXPECT_TRUE(findMatchingLine(filestream, warningMsg));
  ASSERT_FALSE(std::remove(filepath));
}

TEST(TestLogging, testWarningToSetLogfileAppend)
{
  const char* filepath = "other_append.log";

  // Create file first
  std::ofstream fout(filepath);

  ASSERT_TRUE(fout);
  const std::string firstLine = "This line was first.";
  fout << firstLine << std::endl;
  fout.close();

  // Now set helib to it.
  helib::helog.setLogToFile(filepath);

  const std::string warningMsg = "Warning message 600!";
  helib::Warning(warningMsg);
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);
  ASSERT_TRUE(findMatchingLines(filestream, firstLine, warningMsg));
  ASSERT_FALSE(std::remove(filepath));
}

TEST(TestLogging, testWarningToSetLogfileOverwrite)
{
  const char* filepath = "other_overwrite.log";

  // Create file first
  std::ofstream fout(filepath);

  ASSERT_TRUE(fout);
  const std::string firstLine = "This line was first.";
  fout << firstLine << std::endl;
  fout.close();

  // Now set helib to it.
  helib::helog.setLogToFile(filepath, true);

  const std::string warningMsg = "Warning message 500!";
  helib::Warning(warningMsg);
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);
  ASSERT_FALSE(findMatchingLine(filestream, firstLine));
  filestream.clear();  // clear stream flags.
  filestream.seekg(0); // file reset.
  ASSERT_TRUE(findMatchingLine(filestream, warningMsg));
  ASSERT_FALSE(std::remove(filepath));
}

TEST(TestLogging, testWarningToStderr)
{
  helib::helog.setLogToStderr();

  // Redirect the read buffer of std cerr to stringstream
  std::stringstream capturedStderr;
  // This saves the original buffer to restore later
  std::streambuf* sbuf = std::cerr.rdbuf();
  std::cerr.rdbuf(capturedStderr.rdbuf());

  const std::string warningMsg = "Warning message 400!";
  helib::Warning(warningMsg);

  // Reset std cerr buffer
  std::cerr.rdbuf(sbuf);

  EXPECT_TRUE(capturedStderr.str().find(warningMsg) != std::string::npos);
}

} // namespace
