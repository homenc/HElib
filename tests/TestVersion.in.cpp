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

#include <helib/version.h>
#include <tuple>
#include <string>
#include <regex>
#include <sstream>
#include <fstream>
#include "test_common.h"
#include "gtest/gtest.h"

namespace {

// Util function to grab version from cmake file.
std::tuple<int, int, int, std::string> readVersionFromCmakeFile()
{
  // File will be in src/, we are in tests/
  const std::string cmakeFilePath = "@HELIB_PROJECT_ROOT_DIR@/VERSION";

  // Slurp the file. It isn't a large file.
  std::string fileStr;

  {
    std::ifstream cmakeFile(cmakeFilePath);
    if (!cmakeFile.is_open())
      throw std::runtime_error("Could not find '" + cmakeFilePath + "'.");

    std::ostringstream oss;
    oss << cmakeFile.rdbuf();
    fileStr = oss.str();
  }

  // Find the version.
  // e.g.
  // x.y.z

  std::regex re_version(R"((\d+)\.(\d+)\.(\d))");
  std::smatch match;
  std::regex_search(fileStr, match, re_version);
  if (match.size() != 4) {
    std::ostringstream oss;
    oss << "Expected 4 matches, got " << match.size() << ".";
    throw std::runtime_error(oss.str());
  }

  return std::tuple<int, int, int, std::string>(std::stoi(match[1]),
                                                std::stoi(match[2]),
                                                std::stoi(match[3]),
                                                match[0].str());
}

TEST(TestVersion, versionMatchesThatFoundInCMakelists)
{
  int major;
  int minor;
  int patch;
  std::string verStr;

  std::tie(major, minor, patch, verStr) = readVersionFromCmakeFile();

  // EXPECT_EQ does not like these class static members.
  int cur_maj = helib::version::major;
  int cur_min = helib::version::minor;
  int cur_pat = helib::version::patch;

  EXPECT_EQ(cur_maj, major);
  EXPECT_EQ(cur_min, minor);
  EXPECT_EQ(cur_pat, patch);
  EXPECT_STREQ(helib::version::asString, verStr.c_str());
}

TEST(TestVersion, versionGreaterThanOrEqualTo)
{
  EXPECT_FALSE(helib::version::greaterEquals(3));
  EXPECT_FALSE(helib::version::greaterEquals(3, 1));
  EXPECT_FALSE(helib::version::greaterEquals(3, 1, 4));
  // Shouldn't work for negative numbers
  EXPECT_FALSE(helib::version::greaterEquals(1, -1, 0));
  EXPECT_FALSE(helib::version::greaterEquals(0, -1, 0));
  // Older version
  EXPECT_TRUE(helib::version::greaterEquals(1, 0, 0));
  // Current version
  // clang-format off
  EXPECT_TRUE(helib::version::greaterEquals(@PROJECT_VERSION_MAJOR@, @PROJECT_VERSION_MINOR@, @PROJECT_VERSION_PATCH@));
  // clang-format on
}

TEST(TestVersion, versionLibString)
{
  const char* libString = helib::version::libString();
  EXPECT_STREQ(libString, "v@PROJECT_VERSION@");
}

} // namespace
