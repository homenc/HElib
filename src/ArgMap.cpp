/* Copyright (C) 2012-2019 IBM Corp.
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
#include <regex>
#include <fstream>
#include "ArgMap.h"

namespace helib {

ArgMap& ArgMap::note(const std::string& s) {
  docStream << "\t\t   " << s << "\n";
  return *this;
}

void ArgMap::usage(const std::string& msg) const {
  if (!msg.empty())
    std::cerr << msg << std::endl;
  std::cerr << "Usage: " << this->progname << " [ name=value ]...\n";
  std::cerr << doc();
  exit(EXIT_FAILURE);
}

std::string ArgMap::doc() const { return docStream.str(); }

ArgMap& ArgMap::kvSep(char c) {
  this->kv_separator = c;
  return *this;
}

ArgMap& ArgMap::optional() {
  this->required_flag = false;
  return *this;
}

ArgMap& ArgMap::required() {
  this->required_flag = true;
  return *this;
}

// static void printMatchResults(
//  int size,
//  std::string token,
//  std::string key,
//  std::string value)
//{
//    std::cerr << "Match size: " << size << '\n'
//              << "Token (Full Match): " << token << '\n'
//              << "Key: "   << key << '\n'
//              << "Value: " << value << std::endl;
//}

void ArgMap::simpleRegexParse(const std::string& line,
                              bool duplicates,
                              std::function<void(const std::string&)> stop) {
  if (stop == nullptr) {
    stop = std::bind(&ArgMap::usage, this, std::placeholders::_1);
  }

  // Accepts blanks between the word and the kv separator.
  // TODO Values can be grouped by '[]'. May require a different approach.
  const std::string pattern =
      R"((\S+)\s*)" + std::string(1, kv_separator) + R"(\s*(\[.*?\]|\S+)|\S+)";
  std::regex re(pattern);

  // regex iterator
  auto words = std::sregex_iterator(line.begin(), line.end(), re);
  auto words_end = std::sregex_iterator();

  // iterate through matches
  for (; words != words_end; ++words) {

    const std::string token(words->str(0));
    const std::string key(words->str(1));
    const std::string value(words->str(2));

    // Useful for debug this under some print debug
    //    printMatchResults(words->size(), token, key, value);

    // Check if not called before
    if (!duplicates && this->previous_call_set.count(key)) {
      stop("Attempting to set same variable '" + key + "' twice.");
    }

    // Special case built-in '-h' flag for help
    if (token == "-h") {
      stop("");
    }

    // Select ArgProcessor
    std::shared_ptr<ArgProcessor> ap = map[key];

    // Is it a registered variable?
    if (!ap) {
      stop("Unrecognised argument " + ((key.empty()) ? token : key));
    }

    // Process value (parse and set)
    if (!ap->process(value))
      stop("");

    // Remove from required_set (if it is there)
    this->required_set.erase(key);

    // Previously called.
    this->previous_call_set.insert(key);
  }
}

ArgMap& ArgMap::parse(int argc, char** argv) {
  this->progname = std::string(argv[0]);

  if (argc > 1) {
    // put cmdline back together (without progname)
    std::string line;
    for (long i = 1; i < argc; i++) {
      line += argv[i];
      line += " ";
    }

    simpleRegexParse(line);
  }

  // Have the required args been provided.
  if (!this->required_set.empty()) {
    std::ostringstream oss;
    oss << "Required argument(s) not given:\n";
    for (auto& e : required_set)
      oss << "\t" << e << '\n';
    usage(oss.str()); // exits
  }

  return *this;
}

ArgMap& ArgMap::parse(const std::string& filepath) {
  std::ifstream file(filepath);
  this->progname = filepath;

  if (!file.is_open()) {
    throw helib::RuntimeError("Could not open file " + filepath);
  }

  std::string single_line;
  std::string line;
  std::regex re_comment_lines(R"((^\s*\#)|(^\s+$))");
  while (getline(file, line)) {
    if (std::regex_search(line, re_comment_lines)) {
      continue; // ignore comment lines and empties.
    }
    // create single line
    single_line += line;
    single_line += " ";
  }

  file.close();

  simpleRegexParse(single_line, false, [&filepath](const std::string& msg) {
    throw helib::RuntimeError("Could not parse params file: " + filepath);
  });

  // Have the required args been provided.
  if (!this->required_set.empty()) {
    std::ostringstream oss;
    oss << "Required argument(s) not given:\n";
    for (auto& e : required_set)
      oss << "\t" << e << '\n';
    throw helib::RuntimeError(oss.str());
  }

  return *this;
}

}
