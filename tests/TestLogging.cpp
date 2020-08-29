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

#include <helib/log.h>
#include "test_common.h"
#include "gtest/gtest.h"

namespace {

static std::string penultimateLine(std::istream& is)
{
  std::string line;
  std::string penultimateLine;
  while (!is.eof()) {
    penultimateLine = line;
    std::getline(is, line);
  }
  return penultimateLine;
}

// First we run the tests that do not change the default `ostream`.
// This is to test the default behaviour.

TEST(TestLogging, testWarningToDefaultLogfile)
{
  const char* filepath = "helib.log";

  helib::Warning("Warning message!");
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);

  std::string warningMsg;

  // The default file will append so need to check last line.
  // Technically the penultimate as the last line is blank due to a final '\n'.
  warningMsg = penultimateLine(filestream);

  EXPECT_EQ("WARNING: Warning message!", warningMsg);

  // Should not really delete the default log file
}

// Now we run the tests below that explicitly set the stream.

TEST(TestLogging, testWarningToSetLogfile)
{
  const char* filepath = "other.log";

  helib::helog.setLogToFile(filepath);

  helib::Warning("Warning message!");
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);

  std::string warningMsg;
  std::getline(filestream, warningMsg);

  EXPECT_EQ("WARNING: Warning message!", warningMsg);

  ASSERT_FALSE(std::remove(filepath));
}

TEST(TestLogging, testWarningToSetLogfileAppend)
{
  const char* filepath = "other_append.log";

  // Create file first
  std::ofstream fout(filepath);

  ASSERT_TRUE(fout);

  fout << "This line was first.\n";
  fout.close();

  // Now set helib to it.
  helib::helog.setLogToFile(filepath);

  helib::Warning("Warning message!");
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);

  std::string warningMsg;
  std::getline(filestream, warningMsg);

  EXPECT_EQ("This line was first.", warningMsg);

  std::getline(filestream, warningMsg);

  EXPECT_EQ("WARNING: Warning message!", warningMsg);

  ASSERT_FALSE(std::remove(filepath));
}

TEST(TestLogging, testWarningToSetLogfileOverwrite)
{
  const char* filepath = "other_overwrite.log";

  // Create file first
  std::ofstream fout(filepath);

  ASSERT_TRUE(fout);

  fout << "This line was first.\n";
  fout.close();

  // Now set helib to it.
  helib::helog.setLogToFile(filepath, true);

  helib::Warning("Warning message!");
  std::ifstream filestream(filepath);

  ASSERT_TRUE(filestream);

  std::string warningMsg;
  std::getline(filestream, warningMsg);

  EXPECT_NE("This line was first.", warningMsg);
  EXPECT_EQ("WARNING: Warning message!", warningMsg);

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

  helib::Warning("Warning message!");

  // Reset std cerr buffer
  std::cerr.rdbuf(sbuf);

  EXPECT_EQ("WARNING: Warning message!\n", capturedStderr.str());
}

} // namespace
