/* Copyright (C) 2019-2020 IBM Corp.
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

#include <cstring> // strcpy
#include <set>
#include <array>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <chrono>
#include <random>
#include "test_common.h"
#include "gtest/gtest.h"
#include <helib/assertions.h>
#include <NTL/vector.h>

namespace {

// Util func to convert a string to vector of type T.
template <typename T = std::string>
std::vector<T> lineToVector(const std::string& line)
{
  std::istringstream iss(line);
  std::istream_iterator<T> isiter(iss);
  return std::vector<T>(isiter, {});
}

class TestArgMapCmdLine : public ::testing::Test
{

protected:
  int argc = 0;
  char** argv = nullptr;

public:
  // Helper func to generate argc and argv
  virtual void mockCmdLineArgs(const std::string& s)
  {

    std::vector<std::string> words;
    // pre-process the single quotes.
    if (std::count(s.begin(), s.end(), '\'')) {
      std::istringstream iss(s);
      std::string line;
      for (bool append_only = false; std::getline(iss, line, '\'');
           append_only = !append_only) {
        if (append_only) {
          words.push_back(line);
        } else {
          std::vector<std::string> tmp_hold = lineToVector<std::string>(line);
          words.insert(words.end(), tmp_hold.begin(), tmp_hold.end());
        }
      }
    } else {
      words = lineToVector<std::string>(s);
    }

    this->argc = words.size();
    this->argv = new char*[argc + 1];

    // Copy strings to argv
    for (int i = 0; i < this->argc; i++) {
      this->argv[i] = new char[words[i].length() + 1];
      strcpy(this->argv[i], words[i].c_str());
    }

    this->argv[argc] = nullptr;
  };

  // Need to delete argv nicely.
  ~TestArgMapCmdLine() override
  {
    for (int i = 0; i < this->argc; i++) {
      delete[] argv[i];
    }
    delete[] argv;
  };
};

class TestArgMapSampleFile : public ::testing::Test
{

private:
  bool del_file_flag = false;

protected:
  std::string filepath = "sample.params";

public:
  bool createSampleTestFile(const std::string& s,
                            const std::string& sample_file = "")
  {

    this->del_file_flag = true;
    if (!sample_file.empty())
      this->filepath = sample_file;

    std::ofstream file(this->filepath);
    if (!file.is_open())
      return false;

    file << s;

    return true;
  }

  void doNotDeleteFile() { del_file_flag = false; }

  ~TestArgMapSampleFile() override
  {
    if (!this->filepath.empty() && del_file_flag) {
      // Delete the tmp test file
      if (remove(filepath.c_str()) != 0) {
        std::cerr << "[***WARN***] Error removing file: " << this->filepath
                  << std::endl;
      }
    }
  }
};

// For death tests naming convention
using DeathTestArgMapCmdLine = TestArgMapCmdLine;

TEST_F(DeathTestArgMapCmdLine, documentationShownIfDefaultHelpSelectedCmdLine)
{
  mockCmdLineArgs("./prog --help");

  helib::ArgMap amap1;

  EXPECT_EXIT(amap1.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Usage");

  // And again got '-h'
  mockCmdLineArgs("./prog -h");

  helib::ArgMap amap2;

  EXPECT_EXIT(amap2.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Usage.*");
}

TEST_F(DeathTestArgMapCmdLine, documentationCheckCmdLine)
{
  mockCmdLineArgs("./prog -h");

  struct Opts
  {
    int i = 0;
    bool a = false;
    bool b = true;
    std::string s;
    std::string file;
    std::vector<std::string> v;
  } opts;

  helib::ArgMap amap;
  amap.required()
      .positional()
      .arg("string", opts.s, "a string to set.")
      .optional()
      .arg("file", opts.file, "a string to set.")
      .note("A small note.")
      .named()
      .arg("i", opts.i, "an int to place.")
      .toggle()
      .arg("-a", opts.a, "a bool to set true.")
      .toggle(false)
      .arg("-b", opts.b, "b bool to set false.")
      .dots(opts.v, "file");

  EXPECT_EXIT(
      amap.parse(argc, argv),
      ::testing::ExitedWithCode(EXIT_FAILURE),
      R"(^Usage: \./prog \[i=<arg>\] \[-a\] \[-b\] string \[file\] \[file \.\.\.\].)");
}

TEST_F(DeathTestArgMapCmdLine, documentationShownIfCustomHelpSelectedCmdLine)
{
  mockCmdLineArgs("./prog --foo");

  struct Opts
  {
    int arg1 = 5;
    float arg2 = 34.5;
    std::string arg3 = "Hello World";
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      // This now should not show default
      .arg("chris", opts.arg3, "message string", "")
      .helpArgs({"--bar", "--foo"});

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Usage.*");
}

TEST_F(DeathTestArgMapCmdLine,
       documentationShownIfSingleCustomHelpSelectedCmdLine)
{
  mockCmdLineArgs("./prog --bar");

  struct Opts
  {
    int arg1 = 5;
    float arg2 = 34.5;
    std::string arg3 = "Hello World";
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      // This now should not show default
      .arg("chris", opts.arg3, "message string", "")
      .helpArgs("--bar");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Usage.*");
}

TEST_F(TestArgMapSampleFile, documentationShownIfHelpSelectedFromFile)
{
  std::ostringstream oss;
  oss << "-h\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1 = 5;
    float arg2 = 34.5;
    std::string arg3 = "Hello World";
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      // This now should not show default
      .arg("chris", opts.arg3, "message string", "");

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST_F(TestArgMapSampleFile, documentationShownIfCustomHelpSelectedFromFile)
{
  std::ostringstream oss;
  oss << "--foo";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1 = 5;
    float arg2 = 34.5;
    std::string arg3 = "Hello World";
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      // This now should not show default
      .arg("chris", opts.arg3, "message string", "")
      .helpArgs({"--bar", "--foo"});

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST_F(TestArgMapSampleFile,
       documentationShownIfSingleCustomHelpSelectedFromFile)
{
  std::ostringstream oss;
  oss << "--bar";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1 = 5;
    float arg2 = 34.5;
    std::string arg3 = "Hello World";
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      // This now should not show default
      .arg("chris", opts.arg3, "message string", "")
      .helpArgs("--bar");

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST_F(TestArgMapSampleFile, emptyLinesFromFile)
{
  std::ostringstream oss;
  oss << "alice=3\n\n\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    int arg2;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      .parse(filepath);
}

TEST_F(DeathTestArgMapCmdLine, illFormedCmdLine)
{
  mockCmdLineArgs("./prog alic");

  struct Opts
  {
    int arg1;
    float arg2;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Unrecognised argument 'alic'\nUsage.*");
}

TEST_F(DeathTestArgMapCmdLine, danglingSeparatorCmdLine)
{
  mockCmdLineArgs("./prog alice=");

  struct Opts
  {
    int arg1;
    float arg2;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "");

  EXPECT_EXIT(
      amap.parse(argc, argv),
      ::testing::ExitedWithCode(EXIT_FAILURE),
      "^Dangling value for named argument 'alice' after separator.\nUsage.*");
}

TEST_F(TestArgMapSampleFile, danglingSeparatorFromFile)
{
  std::ostringstream oss;
  oss << "./prog alice=";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "");

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST_F(TestArgMapSampleFile, illFormedSampleFile)
{
  std::ostringstream oss;
  oss << "alice=5\n"
      << "bo\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "");

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST_F(DeathTestArgMapCmdLine, nullptrAndEmptyStringsForNoDefaultsCmdLine)
{
  mockCmdLineArgs("./prog -h");

  struct Opts
  {
    int arg1;
    float arg2;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string1", "")
      .arg("bob", opts.arg2, "message string2", nullptr);

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Usage.*string1\n.*string2\n");
}

TEST_F(TestArgMapCmdLine, namedArgsCmdLine)
{
  mockCmdLineArgs("./prog alice=1 bob=2.2 chris=NotIn");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "")
      .arg("chris", opts.arg3, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "NotIn");
}

TEST_F(TestArgMapCmdLine, argumentValuesWithSpacesCmdLine)
{
  mockCmdLineArgs(
      "./prog --alice=1 --chris='Not in here' --bob=2.2  'file name.txt'");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    std::string arg4;
  } opts;

  helib::ArgMap()
      .arg("--alice", opts.arg1, "message string", "")
      .arg("--bob", opts.arg2, "message string", "")
      .arg("--chris", opts.arg3, "message string", "")
      .positional()
      .arg("dave", opts.arg4, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "Not in here");
  EXPECT_EQ(opts.arg4, "file name.txt");
}

TEST_F(TestArgMapSampleFile, argumentValuesWithSpacesFromFile)
{
  std::ostringstream oss;
  oss << "--alice=1\n"
      << "--chris=Not in here\n"
      << "--bob=2.2\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    std::string arg4;
  } opts;

  helib::ArgMap()
      .arg("--alice", opts.arg1, "message string", "")
      .arg("--bob", opts.arg2, "message string", "")
      .arg("--chris", opts.arg3, "message string", "")
      .parse(filepath);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "Not in here");
}

TEST_F(TestArgMapCmdLine, argumentValuesWithSpacesAndWhitespaceSeparatorCmdLine)
{
  mockCmdLineArgs(
      "./prog --alice   1 --chris 'Not in here' --bob  2.2  'file name.txt'");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    std::string arg4;
  } opts;

  helib::ArgMap()
      .separator(helib::ArgMap::Separator::WHITESPACE)
      .arg("--alice", opts.arg1, "message string", "")
      .arg("--bob", opts.arg2, "message string", "")
      .arg("--chris", opts.arg3, "message string", "")
      .positional()
      .arg("dave", opts.arg4, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "Not in here");
  EXPECT_EQ(opts.arg4, "file name.txt");
}

TEST_F(DeathTestArgMapCmdLine, settingSameNameTwiceCmdLine)
{
  mockCmdLineArgs("./prog bob=2 bob=1");

  struct Opts
  {
    int arg1;
  } opts;

  helib::ArgMap amap;
  amap.arg("bob", opts.arg1, "message string");
  amap.parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
}

TEST_F(TestArgMapSampleFile, settingSameNameTwiceFromFile)
{
  std::ostringstream oss;
  oss << "bob=2\n"
      << "bob=1\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
  } opts;

  helib::ArgMap amap;
  amap.arg("bob", opts.arg1, "message string");

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST(TestArgMap, settingSameVariableTwice)
{
  struct Opts
  {
    int arg1;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string");

  EXPECT_THROW(amap.arg("bob", opts.arg1, "message string"), helib::LogicError);
}

TEST_F(TestArgMapCmdLine, spacedArgsCmdLine)
{
  mockCmdLineArgs("./prog alice= 1 bob = 2.2 chris =NotIn");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "")
      .arg("chris", opts.arg3, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "NotIn");
}

TEST_F(TestArgMapSampleFile, spacedArgsFromFile)
{
  std::ostringstream oss;
  oss << "alice= 1\n"
      << "bob = 2.2\n"
      << "chris =NotIn\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "")
      .arg("chris", opts.arg3, "message string", "")
      .parse(filepath);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "NotIn");
}

TEST_F(DeathTestArgMapCmdLine, unrecognisedArgsCmdLine)
{
  mockCmdLineArgs("./prog lice=1 bob=2.2 chris=NotIn");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "")
      .arg("chris", opts.arg3, "message string", "");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Unrecognised argument 'lice'\nUsage.*");
}

TEST_F(TestArgMapSampleFile, unrecognisedArgsFromFile)
{
  std::ostringstream oss;
  oss << "lice=1\n"
      << "bob=2.2\n"
      << "chris=NotIn\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string", "")
      .arg("bob", opts.arg2, "message string", "")
      .arg("chris", opts.arg3, "message string", "");

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST_F(TestArgMapCmdLine, changingSeparatorCmdLine)
{
  mockCmdLineArgs("./prog alice:1 bob:7.5 chris:Hi");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .separator(helib::ArgMap::Separator::COLON)
      .arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      .arg("chris", opts.arg3, "message string")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 7.5);
  EXPECT_EQ(opts.arg3, "Hi");
}

TEST_F(TestArgMapCmdLine, changingToWhitespaceSeparatorCmdLine)
{
  mockCmdLineArgs("./prog -alice  1 -bob   7.5 -chris Hi");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .separator(helib::ArgMap::Separator::WHITESPACE)
      .arg("-alice", opts.arg1, "message string")
      .arg("-bob", opts.arg2, "message string")
      .arg("-chris", opts.arg3, "message string")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 7.5);
  EXPECT_EQ(opts.arg3, "Hi");
}

TEST_F(DeathTestArgMapCmdLine, wrongSeparatorNonWhitespaceCaseCmdLine)
{
  mockCmdLineArgs("./prog alice:1 bob:7.5 chris:Hi");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap amap;
  amap.separator(helib::ArgMap::Separator::EQUALS)
      .arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      .arg("chris", opts.arg3, "message string");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Unrecognised argument 'alice:1'\nUsage.*");
}

TEST_F(DeathTestArgMapCmdLine, wrongSeparatorWhitespaceCaseCmdLine)
{
  mockCmdLineArgs("./prog alice=1 bob=7.5 chris=Hi");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap amap;
  amap.separator(helib::ArgMap::Separator::WHITESPACE)
      .arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      .arg("chris", opts.arg3, "message string");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              "^Unrecognised argument 'alice=1'\nUsage.*");
}

TEST_F(TestArgMapSampleFile, changingSeparatorFromFile)
{
  std::ostringstream oss;
  oss << "alice:1\n"
      << "bob:7.5\n"
      << "chris:Hi\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .separator(helib::ArgMap::Separator::COLON)
      .arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      .arg("chris", opts.arg3, "message string")
      .parse(filepath);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 7.5);
  EXPECT_EQ(opts.arg3, "Hi");
}

TEST_F(TestArgMapCmdLine, requiredNamedArgumentGivenCmdLine)
{
  mockCmdLineArgs("./prog alice=1 bob=7.5");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string")
      .required()
      .arg("bob", opts.arg2, "message string")
      .optional()
      .arg("chris", opts.arg3, "message string")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 7.5);
  EXPECT_EQ(opts.arg3, "");
}

TEST_F(TestArgMapSampleFile, requiredNamedArgumentGivenFromFile)
{
  std::ostringstream oss;
  oss << "alice=1\n"
      << "bob=7.5\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string")
      .required()
      .arg("bob", opts.arg2, "message string")
      .optional()
      .arg("chris", opts.arg3, "message string")
      .parse(filepath);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 7.5);
  EXPECT_EQ(opts.arg3, "");
}

TEST_F(DeathTestArgMapCmdLine, requiredPositionalNotGivenArgsCmdLine)
{
  mockCmdLineArgs("./prog");

  struct Opts
  {
    std::string dave;
  } opts;

  helib::ArgMap amap;
  amap.required().positional().arg("dave", opts.dave, "message string", "");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              R"(^Required argument\(s\) not given:.*)");
}

TEST_F(DeathTestArgMapCmdLine, tooManyPositionalArgsCmdLine)
{
  mockCmdLineArgs("./prog dave1 notEve bob5");

  struct Opts
  {
    std::string dave;
    std::string eve;
  } opts;

  helib::ArgMap amap;
  amap.required()
      .positional()
      .arg("dave", opts.dave, "message string", "")
      .arg("eve", opts.eve, "message string", "");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              ("^Unrecognised argument 'bob5'\n"
               "There could be too many positional arguments\n"
               "Usage.*"));
}

TEST_F(DeathTestArgMapCmdLine, requiredNamedArgumentNotGivenCmdLine)
{
  mockCmdLineArgs("./prog alice=1");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string")
      .required()
      .arg("bob", opts.arg2, "message string")
      .optional()
      .arg("chris", opts.arg3, "message string");

  EXPECT_EXIT(amap.parse(argc, argv),
              ::testing::ExitedWithCode(EXIT_FAILURE),
              R"(^Required argument\(s\) not given:.*)");
}

TEST_F(TestArgMapSampleFile, requiredArgumentNotGivenFromFile)
{
  std::ostringstream oss;
  oss << "alice=1\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string")
      .required()
      .arg("bob", opts.arg2, "message string")
      .optional()
      .arg("chris", opts.arg3, "message string");

  EXPECT_THROW(amap.parse(filepath), helib::RuntimeError);
}

TEST_F(TestArgMapCmdLine, readInAVectorCmdLine)
{
  mockCmdLineArgs("./prog alice='[1 2]'");

  struct Opts
  {
    NTL::Vec<long> arg1;
  } opts;

  helib::ArgMap()
      .required()
      .arg("alice", opts.arg1, "message string", "")
      .parse(argc, argv);

  NTL::Vec<long> test_v;
  std::stringstream ss("[1 2]");
  ss >> test_v;

  EXPECT_EQ(opts.arg1, test_v);
}

TEST_F(TestArgMapSampleFile, readInAVectorFromFile)
{
  std::ostringstream oss;
  oss << "alice=[1 2]\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    NTL::Vec<long> arg1;
  } opts;

  helib::ArgMap()
      .required()
      .arg("alice", opts.arg1, "message string", "")
      .parse(filepath);

  NTL::Vec<long> test_v;
  std::stringstream ss("[1 2]");
  ss >> test_v;

  EXPECT_EQ(opts.arg1, test_v);
}

TEST_F(TestArgMapSampleFile, namedArgsFromFile)
{
  std::ostringstream oss;
  oss << "alice = 1\n"
      << "bob=7.5\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      .arg("chris", opts.arg3, "message string")
      .parse(filepath);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 7.5);
  EXPECT_EQ(opts.arg3, "");
}

TEST_F(TestArgMapSampleFile, handlingCommentsFromFile)
{
  std::ostringstream oss;
  oss << "# An initial comment line.\n"
      << "alice = 1\n"
      << "# A comment line. bob=1.0\n"
      << "bob=7.5\n"
      << " # A space and a comment line.\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
  } opts;

  helib::ArgMap()
      .arg("alice", opts.arg1, "message string")
      .arg("bob", opts.arg2, "message string")
      .arg("chris", opts.arg3, "message string")
      .parse(filepath);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 7.5);
  EXPECT_EQ(opts.arg3, "");
}

TEST_F(TestArgMapSampleFile, fileDoesNotExist)
{
  struct Opts
  {
    int arg1;
  } opts;

  // File must not be made.

  helib::ArgMap amap;
  amap.arg("alice", opts.arg1, "message string");

  EXPECT_THROW(amap.parse("aaaa.bbb"), helib::RuntimeError);
}

TEST_F(TestArgMapSampleFile, toggleArgsFromFile)
{
  std::ostringstream oss;
  oss << "-f\n"
      << "-alice=1\n"
      << "-bob=2.2\n"
      << "-t\n"
      << "-chris=NotIn\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    bool toggle_t = true;
    bool toggle_f = false;
  } opts;

  helib::ArgMap amap;
  amap.arg("-alice", opts.arg1, "message string", "")
      .arg("-bob", opts.arg2, "message string", "")
      .toggle()
      .arg("-t", opts.toggle_t, "message string", "")
      .toggle(false)
      .arg("-f", opts.toggle_f, "message string", "")
      .named()
      .required()
      .arg("-chris", opts.arg3, "message string", "");

  EXPECT_THROW(amap.parse(filepath), helib::LogicError);
}

TEST_F(TestArgMapCmdLine, toggleArgsCmdLine)
{
  mockCmdLineArgs("./prog -f -alice=1 -bob=2.2 -t -chris=NotIn");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    bool toggle_t = true;
    bool toggle_f = false;
  } opts;

  helib::ArgMap()
      .arg("-alice", opts.arg1, "message string", "")
      .arg("-bob", opts.arg2, "message string", "")
      .toggle()
      .arg("-t", opts.toggle_t, "message string", "")
      .toggle(false)
      .arg("-f", opts.toggle_f, "message string", "")
      .named()
      .required()
      .arg("-chris", opts.arg3, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "NotIn");
  EXPECT_TRUE(opts.toggle_t);
  EXPECT_FALSE(opts.toggle_f);
}

TEST_F(TestArgMapCmdLine, requiredToggleArgsAreOptionalCmdLine)
{
  mockCmdLineArgs("./prog");

  struct Opts
  {
    bool toggle_t = true;
    bool toggle_f = false;
  } opts;

  helib::ArgMap amap1;
  amap1.required().toggle();

  // For TOGGLE_TRUE
  EXPECT_THROW(amap1.arg("-t", opts.toggle_t), helib::LogicError);

  helib::ArgMap amap2;
  amap2.required().toggle(false);

  // For TOGGLE_FALSE
  EXPECT_THROW(amap2.arg("-f", opts.toggle_f), helib::LogicError);
}

TEST_F(TestArgMapCmdLine, positionalArgsCmdLine)
{
  mockCmdLineArgs("./prog -alice=1 -bob=2.2 -t dave1 -chris=NotIn eve2");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    bool toggle_t = true;
    std::string dave;
    std::string eve;
  } opts;

  helib::ArgMap()
      .arg("-alice", opts.arg1, "message string", "")
      .arg("-bob", opts.arg2, "message string", "")
      .toggle()
      .arg("-t", opts.toggle_t, "message string", "")
      .named()
      .required()
      .arg("-chris", opts.arg3, "message string", "")
      .positional()
      .arg("dave", opts.dave, "message string", "")
      .optional()
      .arg("eve", opts.eve, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "NotIn");
  EXPECT_TRUE(opts.toggle_t);
  EXPECT_EQ(opts.dave, "dave1");
  EXPECT_EQ(opts.eve, "eve2");
}

TEST_F(TestArgMapSampleFile, positionalArgsFromFile)
{
  std::ostringstream oss;
  oss << "-alice=1\n"
      << "-bob=2.2\n"
      << "-t\n"
      << "dave1\n"
      << "-chris=NotIn\n"
      << "eve2\n";

  // Cannot perform test without file.
  ASSERT_TRUE(createSampleTestFile(oss.str()));

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    bool toggle_t = true;
    std::string dave;
    std::string eve;
  } opts;

  helib::ArgMap amap;
  amap.arg("-alice", opts.arg1, "message string", "")
      .arg("-bob", opts.arg2, "message string", "")
      .toggle()
      .arg("-t", opts.toggle_t, "message string", "")
      .named()
      .required()
      .arg("-chris", opts.arg3, "message string", "")
      .positional()
      .arg("dave", opts.dave, "message string", "")
      .optional()
      .arg("eve", opts.eve, "message string", "");

  EXPECT_THROW(amap.parse(filepath), helib::LogicError);
}

TEST_F(TestArgMapCmdLine, requiredThenOptionalPositionalArgsCmdLine)
{
  mockCmdLineArgs("./prog -alice=1 -bob=2.2 -t dave1 -chris=NotIn ");

  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    bool toggle_t = true;
    std::string dave;
    std::string eve;
  } opts;

  helib::ArgMap()
      .arg("-alice", opts.arg1, "message string", "")
      .arg("-bob", opts.arg2, "message string", "")
      .toggle()
      .arg("-t", opts.toggle_t, "message string", "")
      .named()
      .required()
      .arg("-chris", opts.arg3, "message string", "")
      .positional()
      .arg("dave", opts.dave, "message string", "")
      .optional()
      .arg("eve", opts.eve, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.arg1, 1);
  EXPECT_FLOAT_EQ(opts.arg2, 2.2);
  EXPECT_EQ(opts.arg3, "NotIn");
  EXPECT_TRUE(opts.toggle_t);
  EXPECT_EQ(opts.dave, "dave1");
  EXPECT_EQ(opts.eve, "");
}

TEST_F(TestArgMapCmdLine, positionalArgsPassedWithNameOfArgCmdLine)
{
  mockCmdLineArgs("./prog dave notEve");

  struct Opts
  {
    std::string dave;
    std::string eve;
  } opts;

  helib::ArgMap()
      .required()
      .positional()
      .arg("dave", opts.dave, "message string", "")
      .arg("eve", opts.eve, "message string", "")
      .parse(argc, argv);

  EXPECT_EQ(opts.dave, "dave");
  EXPECT_EQ(opts.eve, "notEve");
}

// Create from 1 to 10 inclusive random order list for testing variable
// positional args.
struct VariablePositionalArgs
{

  static const int upto = 10;
  unsigned seed;
  int random = 0;
  std::array<std::string, upto> default_word_list{{"apple.txt",
                                                   "banana.txt",
                                                   "cat.dat",
                                                   "dog.dat",
                                                   "eagle.md",
                                                   "fox.sh",
                                                   "gorilla.dat",
                                                   "herring",
                                                   "imp.cpp",
                                                   "jaguar.car"}};

  std::vector<std::string> word_list;

  VariablePositionalArgs()
  {
    // between 1 and 10 inclusive is good enough.
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::uniform_int_distribution<int> distribution(1, upto);
    std::default_random_engine generator(seed);
    random = distribution(generator);

    std::shuffle(default_word_list.begin(),
                 default_word_list.end(),
                 std::default_random_engine(seed));
  }

  std::string genStr(std::string first,
                     std::multiset<std::string> priority_mixin = {})
  {

    // Random number of positional dots args form default_word_list
    std::ostringstream string_list;
    string_list << first;

    std::size_t max_size = random + priority_mixin.size();
    word_list.resize(max_size);

    // In case this isn't our first call.
    if (word_list.size() > max_size) {
      word_list.resize(max_size);
    }

    // random [1, upto]
    std::copy(default_word_list.begin(),
              default_word_list.begin() + random,
              word_list.begin());

    if (priority_mixin.size() > 0) {
      std::copy(priority_mixin.begin(),
                priority_mixin.end(),
                word_list.begin() + random);
      // Need to mix them.
      std::shuffle(word_list.begin(),
                   word_list.end(),
                   std::default_random_engine(seed));
    }

    for (const auto& word : word_list)
      string_list << ' ' << word;

    return string_list.str();
  }

  int getRandom() const { return random; }

  const std::vector<std::string>& getWordList() const { return word_list; }
};

TEST_F(TestArgMapCmdLine, onlyVariableNumbersOfOptionalPositional)
{

  VariablePositionalArgs var_pos;

  mockCmdLineArgs(var_pos.genStr("./prog"));

  struct Opts
  {
    std::vector<std::string> dots;
  } opts;

  helib::ArgMap().dots(opts.dots, "file").parse(argc, argv);

  // Has it been populated first?
  ASSERT_TRUE(opts.dots.size() > 0);
  // Now check the entries.
  for (size_t i = 0; i < opts.dots.size(); ++i)
    EXPECT_EQ(opts.dots[i], var_pos.getWordList()[i])
        << "i = " << i << std::endl;
}

TEST_F(TestArgMapCmdLine,
       variableNumbersOfOptionalPositionalWithMixedInOtherArgTypes)
{

  VariablePositionalArgs var_pos;

  std::multiset<std::string> otherArgs = {"-t", "-t", "k = v", "-f"};

  std::string cmdline = var_pos.genStr("./prog", otherArgs);

  mockCmdLineArgs(cmdline);

  struct Opts
  {
    bool t = false;
    bool f = true;
    std::string value;
    std::string bob;
    std::vector<std::string> dots;
  } opts;

  helib::ArgMap()
      .optional()
      .toggle()
      .arg("-t", opts.t, "")
      .toggle(false)
      .arg("-f", opts.f, "")
      .required()
      .named()
      .arg("k", opts.value, "")
      .positional()
      .arg("bob", opts.bob, "")
      .dots(opts.dots, "file")
      .parse(argc, argv);

  std::vector<std::string> word_list;
  std::remove_copy_if(
      begin(var_pos.getWordList()),
      end(var_pos.getWordList()),
      std::back_inserter(word_list),
      [&otherArgs](const std::string& s) { return otherArgs.count(s) > 0; });

  ASSERT_TRUE(word_list.size() > 0);
  ASSERT_EQ(opts.dots.size(), word_list.size() - 1);
  // Now check the entries.
  EXPECT_TRUE(opts.t);
  EXPECT_FALSE(opts.f);
  EXPECT_EQ(opts.value, "v");
  EXPECT_EQ(opts.bob, word_list.front());
  // and the dots.

  // Start form 1 as zeroth is the normal positional arg.
  for (size_t i = 1; i < opts.dots.size(); ++i)
    EXPECT_EQ(opts.dots[i], word_list[i + 1]) << "i = " << i << std::endl;
}

TEST(TestArgMap, secondTimeOptionalPositionalArgs)
{
  struct Opts
  {
    int arg1;
    float arg2;
    std::string arg3;
    bool toggle_t = true;
    std::string dave;
    std::string eve;
    std::string freya;
  } opts;

  helib::ArgMap amap;
  amap.arg("-alice", opts.arg1, "message string", "")
      .arg("-bob", opts.arg2, "message string", "")
      .toggle()
      .arg("-t", opts.toggle_t, "message string", "")
      .required()
      .named()
      .arg("-chris", opts.arg3, "message string", "")
      .positional()
      .arg("dave", opts.dave, "message string", "")
      .optional()
      .arg("eve", opts.eve, "message string", "")
      .required();

  EXPECT_THROW(amap.arg("freya", opts.freya, "message string", ""),
               helib::LogicError);
}

TEST(TestArgMap, emptyNameOfArg)
{
  int i;

  helib::ArgMap amap;

  EXPECT_THROW(amap.arg("", i, ""), helib::LogicError);
}

TEST(TestArgMap, whitespaceInNameOfArg)
{
  int i;

  helib::ArgMap amap1;
  // At beginning
  EXPECT_THROW(amap1.arg(" alice", i, ""), helib::LogicError);

  helib::ArgMap amap2;
  // In middle
  EXPECT_THROW(amap2.arg("ali ce", i, ""), helib::LogicError);

  helib::ArgMap amap3;
  // At end
  EXPECT_THROW(amap3.arg("alice ", i, ""), helib::LogicError);
}

} // namespace
