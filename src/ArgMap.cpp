/* Copyright (C) 2012-2020 IBM Corp.
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
#include <algorithm>
#include <regex>
#include <fstream>
#include <cctype>
#include <helib/ArgMap.h>

namespace helib {

// Three functions strip whitespaces before and after strings.
static void lstrip(std::string& s)
{
  auto it = std::find_if(
      s.begin(), s.end(), [](unsigned char c) { return !std::isspace(c); });

  s.erase(s.begin(), it);
}

static void rstrip(std::string& s)
{
  auto it = std::find_if(
      s.rbegin(), s.rend(), [](unsigned char c) { return !std::isspace(c); });

  s.erase(it.base(), s.end());
}

static void strip(std::string& s)
{
  lstrip(s);
  rstrip(s);
}

ArgMap& ArgMap::note(const std::string& s)
{
  docStream << "\t\t" << s << "\n";
  return *this;
}

static const std::string str_if_cond(bool cond, const char sep, const char* txt)
{
  std::string ext;
  if (cond)
    ext.append(1, sep).append(txt);

  return ext;
}

void ArgMap::usage(const std::string& msg) const
{
  if (!msg.empty())
    std::cerr << msg << std::endl;

  std::cerr << "Usage: " << this->progname;
  for (const auto& n : this->optional_set) {
    bool named = (this->map.at(n)->getArgType() == ArgType::NAMED);
    std::cerr << " [" << n << str_if_cond(named, this->kv_separator, "<v>")
              << "]";
  }

  for (const auto& n : this->required_set) {
    bool named = (this->map.at(n)->getArgType() == ArgType::NAMED);
    std::cerr << " " << n << str_if_cond(named, this->kv_separator, "<v>");
  }

  if (this->dots_enabled)
    std::cerr << " [" << this->dots_name << " ...]"
              << "\n";
  else
    std::cerr << "\n";

  std::cerr << doc() << std::endl;

  exit(EXIT_FAILURE);
}

std::string ArgMap::doc() const { return docStream.str(); }

ArgMap& ArgMap::helpArgs(const std::initializer_list<std::string> s)
{
  this->help_tokens = s;
  return *this;
}

ArgMap& ArgMap::helpArgs(const std::string s)
{
  this->help_tokens = {s};
  return *this;
}

ArgMap& ArgMap::separator(Separator s)
{

  switch (s) {
  case Separator::EQUALS:
    this->kv_separator = '=';
    break;
  case Separator::COLON:
    this->kv_separator = ':';
    break;
  case Separator::WHITESPACE:
    this->kv_separator = ' ';
    break;
  default:
    // Use of class enums means it should never reach here.
    throw helib::LogicError("Unrecognised option for kv seperator.");
  }

  return *this;
}

ArgMap& ArgMap::optional()
{
  this->required_mode = false;
  return *this;
}

ArgMap& ArgMap::required()
{
  this->required_mode = true;
  return *this;
}

ArgMap& ArgMap::toggle(bool t)
{
  this->named_args_only = false;
  this->arg_type = t ? ArgType::TOGGLE_TRUE : ArgType::TOGGLE_FALSE;
  return *this;
}

ArgMap& ArgMap::named()
{
  this->arg_type = ArgType::NAMED;
  return *this;
}

ArgMap& ArgMap::positional()
{
  this->named_args_only = false;
  this->arg_type = ArgType::POSITIONAL;
  return *this;
}

ArgMap& ArgMap::diagnostics(std::ostream& ostrm)
{
  this->diagnostics_strm = &ostrm;
  return *this;
}

void ArgMap::printDiagnostics(const std::forward_list<std::string>& args) const
{
  if (this->diagnostics_strm != nullptr) {
    // argv as seen by ArgMap
    *this->diagnostics_strm << "Args pre-parse:\n";
    for (const auto& e : args) {
      *this->diagnostics_strm << e << std::endl;
    }
    // required set
    *this->diagnostics_strm << "Required args set:\n";
    for (const auto& e : required_set) {
      *this->diagnostics_strm << e << std::endl;
    }
    // optional set
    *this->diagnostics_strm << "Optional args set:\n";
    for (const auto& e : optional_set) {
      *this->diagnostics_strm << e << std::endl;
    }
  }
}

// Correct the list from argv by splitting on the seperator.
static void splitOnSeparator(std::forward_list<std::string>& args_lst, char sep)
{
  if (sep == ' ')
    return;

  for (auto it = args_lst.begin(); it != args_lst.end(); ++it) {
    if (it->size() != 1) {
      std::size_t pos = it->find(sep);
      if (pos != std::string::npos) {
        if (pos == 0) {
          std::string sub = it->substr(1, std::string::npos);
          *it = sep;
          args_lst.insert_after(it, sub);
        } else {
          std::string sub = it->substr(pos);
          *it = it->substr(0, pos);
          args_lst.insert_after(it, sub);
        }
      }
    }
  }
}

void ArgMap::simpleParse(const std::forward_list<std::string>& args,
                         bool duplicates,
                         std::function<void(const std::string&)> stop)
{
  if (stop == nullptr) {
    stop = std::bind(&ArgMap::usage, this, std::placeholders::_1);
  }

  auto pos_args_it = this->positional_args_list.begin();
  for (auto it = args.begin(); it != args.end(); ++it) {

    const std::string token = *it;

    if (this->help_tokens.count(token))
      stop("");

    // Check if not called before
    if (!duplicates && this->previous_call_set.count(token))
      stop("Attempting to set same variable '" + token + "' twice.");

    // Select ArgProcessor
    auto map_it = this->map.find(token);
    std::shared_ptr<ArgProcessor> ap =
        (map_it == this->map.end()) ? nullptr : map_it->second;

    if (ap && ap->getArgType() != ArgType::POSITIONAL) {

      switch (ap->getArgType()) {
      case ArgType::NAMED:
        // Process value (parse and set)
        if ((++it) == args.end())
          stop("Dangling value for named argument '" + token + "'.");

        if (this->kv_separator == ' ') {
          if (!ap->process(*it))
            stop("Whitespace separator issue. Value:'" + *it + "'");
        } else {
          if ((++it) == args.end())
            stop("Dangling value for named argument '" + token +
                 "' after separator.");
          if (!ap->process(*it))
            stop("Not a valid value '" + *it + "'.");
        }
        break;
      case ArgType::TOGGLE_TRUE:
        if (!ap->process("1"))
          stop("");
        break;
      case ArgType::TOGGLE_FALSE:
        if (!ap->process("0"))
          stop("");
        break;
      default:
        // Should never get here.
        throw helib::LogicError("Unrecognised ArgType.");
        break;
      }

      // Remove from required_set (if it is there)
      this->required_set.erase(token);

      // Previously called.
      this->previous_call_set.insert(token);
    } else if (pos_args_it != this->positional_args_list.end()) {
      // POSITIONAL args are treated differently as it is technically
      // never a recognised token.
      std::shared_ptr<ArgProcessor> pos_ap = map[*pos_args_it];
      if (!pos_ap->process(*it))
        throw helib::LogicError(
            "Positional name does not match a ArgMap name.");
      // Remove from required_set (if it is there)
      this->required_set.erase(*pos_args_it);
      ++pos_args_it;
    } else if (this->dots_enabled) {
      this->dots_ap->process(token);
    } else {
      std::string msg = "Unrecognised argument \'" + token + "\'";
      if (!this->positional_args_list.empty())
        msg += "\nThere could be too many positional arguments";
      stop(msg);
    }
  }
}

ArgMap& ArgMap::parse(int argc, char** argv)
{

  this->progname = std::string(argv[0]);

  std::forward_list<std::string> args(argv + 1, argv + argc);

  splitOnSeparator(args, this->kv_separator);

  // Take any leading and trailing whitespace away.
  std::for_each(args.begin(), args.end(), strip);

  printDiagnostics(args);

  simpleParse(args);

  // Have the required args been provided - if not exit
  if (!this->required_set.empty()) {
    std::ostringstream oss;
    oss << "Required argument(s) not given:\n";
    for (const auto& e : this->required_set)
      oss << "\t" << e << '\n';
    usage(oss.str()); // exits
  }

  return *this;
}

ArgMap& ArgMap::parse(const std::string& filepath)
{

  if (this->kv_separator == ' ') { // Not from files.
    throw helib::LogicError("Whitespace separator not possible from files.");
  }

  if (!this->named_args_only) {
    throw helib::LogicError("Toggle and Positional arguments not possible from "
                            "files. Only named arguments.");
  }

  std::ifstream file(filepath);
  this->progname = filepath;

  if (!file.is_open()) {
    throw helib::RuntimeError("Could not open file " + filepath);
  }

  std::forward_list<std::string> args;
  auto it = args.before_begin();
  std::regex re_comment_lines(R"((^\s*\#)|(^\s+$))");
  std::string line;
  while (getline(file, line)) {
    if (std::regex_search(line, re_comment_lines)) {
      continue; // ignore comment lines and empties.
    }
    it = args.insert_after(it, line);
  }

  splitOnSeparator(args, this->kv_separator);

  // Take any leading and trailing whitespace away.
  std::for_each(args.begin(), args.end(), strip);

  printDiagnostics(args);

  simpleParse(args, false, [&filepath](const std::string& msg) {
    throw helib::RuntimeError("Could not parse params file: " + filepath +
                              ". " + msg);
  });

  // Have the required args been provided - if not throw
  if (!this->required_set.empty()) {
    std::ostringstream oss;
    oss << "Required argument(s) not given:\n";
    for (const auto& e : this->required_set)
      oss << "\t" << e << '\n';
    throw helib::RuntimeError(oss.str());
  }

  return *this;
}

} // namespace helib
