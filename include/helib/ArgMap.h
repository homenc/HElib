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

#ifndef HELIB_ARGMAP_H
#define HELIB_ARGMAP_H

#include <iostream>
#include <iomanip>
#include <forward_list>
#include <initializer_list>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <string>
#include <sstream>
#include <fstream>
#include <cctype>
#include <memory>
#include <vector>
#include <tuple>
#include <type_traits>
#include <regex>
#include <helib/assertions.h>

/**
 * @file ArgMap.h
 * @brief Easier arg parsing.
 **/

namespace helib {

/**
 * @brief Basic class for arg parsing.
 * Example use:
 * @code
 *   // Variables to be set by command line.
 *   long p = 2;                                 // default values.
 *   long m = 19;
 *   bool t = false;
 *   bool f = true;
 *   std::string k = "Hello World";
 *
 *   ArgMap()                                    // (*) marks default.
 *     .required()                               // set args to required.
 *     .positional()                             //
 *       .arg("p", p, "doc for p")               //
 *       .arg("m", m, "doc for m", "undefined")  // special default info.
 *     .optional()                               // swap to optional args (*).
 *     .named()                                  // named args (*) e.g.k=v.
 *     .separator(ArgMap::Separator::WHITESPACE) // change separator to
 *       .arg("-k", k, "doc for k", "")          // whitespace ('=' is (*)).
 *       .note("an extra note")                  // no default value info.
 *     .toggle()                                 // add extra doc/note.
 *        .arg("-t", t, "doc for t", "")         // toggle flag sets bool true.
 *     .toggle(false)                            // toggle flag sets bool false.
 *        .arg("-f", f, "doc for f", "")         //
 *     .helpArgs({"--myhelp"})                   // changes default help flags
 *     .parse(argc, argv);                       // (*) is {"-h", "--help"}.
 *                                               // parses and overwrites values
 * @endcode
 **/
class ArgMap
{
private:
  enum class ArgType
  {
    NAMED,
    TOGGLE_TRUE,
    TOGGLE_FALSE,
    POSITIONAL,
    DOTS
  };

  // requires latching logic.
  class PositionalArgsList
  {
  private:
    std::vector<std::string> positional_args;
    bool optional_flag = false;

  public:
    void insert(std::string name, bool optional)
    {
      if (optional) {
        this->optional_flag = true;
        positional_args.push_back(name);
      } else if (!optional && !optional_flag) {
        positional_args.push_back(name);
      } else {
        throw LogicError("Attempting to have argument '" + name +
                         "' required after optional positional args given.");
      }
    }

    std::vector<std::string>::iterator begin()
    {
      return this->positional_args.begin();
    }

    std::vector<std::string>::iterator end()
    {
      return this->positional_args.end();
    }

    bool empty() { return this->positional_args.empty(); }
  }; // end of PositionalArgsList

  struct ArgProcessor;

  template <typename C>
  class ArgProcessorContainer;

  /* ArgProcessorValue: templated subclasses */
  template <typename T>
  class ArgProcessorValue;

private:
  ArgType arg_type = ArgType::NAMED;
  char kv_separator = '=';
  // Modes and other flags.
  bool required_mode = false;
  bool dots_enabled = false;
  bool named_args_only = true;

  std::string progname;
  std::string dots_name;

  // Track addresses to stop assigning same variable to more than one .arg(...)
  std::set<void*> addresses_used;

  // Track what has been called previously whilst parsing.
  std::set<std::string> previous_call_set;

  // Store the args.
  std::unordered_map<std::string, std::shared_ptr<ArgProcessor>> map;

  // Docs held in vector until called by methods such as doc and usage
  // Triple (arg name (+ value),  docString, whether required)
  std::vector<std::tuple<std::string, std::string, bool>> docVec;

  PositionalArgsList positional_args_list;

  std::unique_ptr<ArgProcessor> dots_ap;

  // Set for tracking required.
  std::set<std::string> required_set;

  // Set for tracking optional.
  std::set<std::string> optional_set;

  std::ostream* diagnostics_strm = nullptr;

  std::set<std::string> help_tokens = {"-h", "--help"};

  // Private for diagnostics
  void printDiagnostics(const std::forward_list<std::string>& args) const;

  //! @brief Arg parsing function
  /**
   * Arg parsing function (private)
   * @param line Space-separated argument line to parse
   * @param duplicates If true does not fail in case of duplicated arguments
   * @param stop Callback function called in case of failure. (Default is Usage)
   */
  void simpleParse(const std::forward_list<std::string>& args,
                   bool duplicates = true,
                   std::function<void(const std::string&)> stop = {});

public:
  enum class Separator
  {
    COLON,
    EQUALS,
    WHITESPACE
  };

  /**
   * @brief Add a new argument description
   * Adds a new argument description with value of type T.
   * Throws helib::RuntimeError if the arg key is duplicated or if the storing
   * variable is used more than once
   * @tparam T The type of the argument
   * @param name The argument name (key)
   * @param value a variable where the argument will be stored. Also used as
   * default value
   * @return A reference to the modified ArgMap object
   */
  template <typename T>
  ArgMap& arg(const std::string& name, T& value);

  /**
   * @brief Add a new argument with docs
   * Adds a new argument description with value of type T and docs.
   * Throws helib::RuntimeError if the arg key is duplicated or if the storing
   * variable is used more than once
   * @tparam T The type of the argument
   * @param name The argument name (key)
   * @param value a variable where the argument will be stored. Also used as
   * default value
   * @param doc1 Description of the argument used when displaying usage
   * @return A reference to the modified ArgMap object
   */
  template <typename T>
  ArgMap& arg(const std::string& name, T& value, const std::string& doc);

  /**
   * @brief Add a new argument with docs and default description
   * Adds a new argument description with value of type T, with docs and
   * default description. Throws helib::RuntimeError if the arg key is
   * duplicated or if the storing variable is used more than once
   * @tparam T The type of the argument
   * @param name The argument name (key)
   * @param value a variable where the argument will be stored. Also used as
   * default value
   * @param doc1 Description of the argument used when displaying usage
   * @param info The default value description (ignored if nullptr or "")
   * @return A reference to the modified ArgMap object
   */
  template <typename T>
  ArgMap& arg(const std::string& name,
              T& value,
              const std::string& doc,
              const char* info);

  /**
   * @brief Adds variable number of positional arg types after defined arg types
   * are exhausted. These are treated as optional.
   * @param container holds the variable positional args. It must have a
   * push_back method for insertion
   * @return A reference to the ArgMap object
   */
  template <typename C>
  ArgMap& dots(C& container, const char* name);

  /**
   * @brief Parse the argv array
   * Parse the argv array
   * If it fails or -h is an argument it prints the usage and exits the program
   * @param argc number of entries in argv
   * @param argv array containing the arguments
   * @return A reference to the ArgMap object
   */
  ArgMap& parse(int argc, char** argv);

  /**
   * @brief Parse the configuration/parameters file
   * Parsing a configuration file only functions with named arguments
   * Parse the config file
   * Throws RuntimeError on failure
   * @param filepath the config file path
   * @return A reference to the ArgMap object
   */
  ArgMap& parse(const std::string& filepath);

  /**
   * @brief Swaps to optional arg mode (default)
   * Swaps to optional arg mode. Following arguments will be considered optional
   * @return A reference to the ArgMap object
   */
  ArgMap& optional();

  /**
   * @brief Swaps to required arg mode
   * Swaps to required arg mode. Following arguments will be considered required
   * @return A reference to the ArgMap object
   */
  ArgMap& required();

  /**
   * @brief Swaps to toggle arg type
   * Swaps to required arg mode. Following arguments will be considered of
   * toggle type
   * @return A reference to the ArgMap object
   */
  ArgMap& toggle(bool t = true);

  /**
   * @brief Swaps to named arg type (default)
   * Swaps to required arg mode. Following arguments will be considered of named
   * type
   * @return A reference to the ArgMap object
   */
  ArgMap& named();

  /**
   * @brief Swaps to positional arg type
   * Swaps to required arg mode. Following arguments will be considered of
   * positional type
   * @return A reference to the ArgMap object
   */
  ArgMap& positional();

  /**
   * @brief Provide custom help toggle args. (defaults are "-h", "--help")
   * Overwrite default help toggle args to custom ones for parsing.
   * @return A reference to the ArgMap object
   */
  ArgMap& helpArgs(const std::initializer_list<std::string> s);
  ArgMap& helpArgs(const std::string s);

  /**
   * @brief Turns on diagnostics printout when parsing
   * Swaps to required arg mode. Following arguments will be considered of
   * positional type
   * @return A reference to the ArgMap object
   */
  ArgMap& diagnostics(std::ostream& ostrm = std::cout);

  /**
   * @brief Sets the key-value separator
   * Sets the named args key-value pair separator character
   * @param s the separator enum must be set either to COLON or EQUALS(default).
   * @return A reference to the ArgMap object
   */
  ArgMap& separator(Separator s);

  /**
   * @brief Adds a note to usage
   * Adds a note to the arg usage description
   * @param s The note string
   * @return A reference to the ArgMap object
   */
  ArgMap& note(const std::string& s);

  /**
   * @brief Print usage and exit
   * Prints the usage and exits the program
   * @param msg An additional message to print before showing usage
   */
  void usage(const std::string& msg = "") const;

  /**
   * @brief Return arg docs
   * Returns the argument documentation as a string
   * @return the argument documentation string
   */
  std::string doc() const;
}; // End of class ArgMap

// Three functions strip whitespaces before and after strings.
static inline void lstrip(std::string& s)
{
  auto it = std::find_if(s.begin(), s.end(), [](unsigned char c) {
    return !std::isspace(c);
  });

  s.erase(s.begin(), it);
}

static inline void rstrip(std::string& s)
{
  auto it = std::find_if(s.rbegin(), s.rend(), [](unsigned char c) {
    return !std::isspace(c);
  });

  s.erase(it.base(), s.end());
}

static inline void strip(std::string& s)
{
  lstrip(s);
  rstrip(s);
}

// Correct the list from argv by splitting on the separator.
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

/* ArgProcessor: A virtual base class that acts as the interface to hold
 * args in a map of different types.
 */
struct ArgMap::ArgProcessor
{
  virtual ~ArgProcessor() = default;
  virtual ArgType getArgType() = 0;
  virtual bool process(const std::string& s) = 0;
}; // end of ArgProcessor

template <typename C>
class ArgMap::ArgProcessorContainer : public ArgProcessor
{
private:
  C* container;
  ArgType arg_type = ArgType::DOTS;

  using T = typename C::value_type;

  // For strings. Avoids a stream breaking on whitespace.
  template <typename U = T,
            typename S,
            typename std::enable_if_t<std::is_same<U, S>::value, int> = 0>
  bool do_process(const S& s)
  {
    container->push_back(s);
    return true;
  }

  // For non-string types, at the mercy of stringstream.
  template <typename U = T,
            typename S,
            typename std::enable_if_t<!std::is_same<U, S>::value, int> = 0>
  bool do_process(const S& s)
  {
    std::istringstream iss(s);
    U tmp_value;
    bool rt = (iss >> tmp_value);
    container->push_back(tmp_value);
    return rt;
  }

public:
  ArgType getArgType() override { return arg_type; }

  bool process(const std::string& s) override { return this->do_process(s); }

  explicit ArgProcessorContainer(C* c) : container(c) {}

}; // end of ArgProcessorContainer

template <typename T>
class ArgMap::ArgProcessorValue : public ArgProcessor
{
private:
  T* value;
  ArgType arg_type;

  // For strings. Avoids a stream breaking on whitespace.
  template <typename U = T,
            typename S,
            typename std::enable_if_t<std::is_same<U, S>::value, int> = 0>
  bool do_process(const S& s)
  {
    *value = s;
    return true;
  }

  // For non-string types, at the mercy of stringstream.
  template <typename U = T,
            typename S,
            typename std::enable_if_t<!std::is_same<U, S>::value, int> = 0>
  bool do_process(const S& s)
  {
    std::istringstream iss(s);
    return bool(iss >> *value);
  }

public:
  ArgType getArgType() override { return arg_type; }

  bool process(const std::string& s) override { return this->do_process(s); }

  explicit ArgProcessorValue(T* v, ArgType at) : value(v), arg_type(at) {}
}; // end of ArgProcessorValue

template <typename T>
ArgMap& ArgMap::arg(const std::string& name, T& value)
{
  // trying to add empty or whitespace name?
  assertTrue<LogicError>(
      !name.empty() &&
          std::none_of(name.begin(),
                       name.end(),
                       [](unsigned char c) { return std::isspace(c); }),
      "Attempting to register an empty string or string with whitespace");

  // has this name already been added?
  assertTrue<LogicError>(map.count(name) == 0,
                         "Key already in arg map (key: " + name + ")");

  // have we seen this addr before?
  assertEq<LogicError>(addresses_used.count(&value),
                       0ul,
                       "Attempting to register variable twice");

  addresses_used.insert(&value);

  map[name] = std::make_shared<ArgProcessorValue<T>>(&value, this->arg_type);

  if (this->arg_type == ArgType::POSITIONAL) {
    this->positional_args_list.insert(name, !this->required_mode);
  }

  if (this->required_mode) {
    // The user should make toggles correct in the code with optional
    if (this->arg_type == ArgType::TOGGLE_TRUE ||
        this->arg_type == ArgType::TOGGLE_FALSE)
      throw LogicError("Toggle argument types cannot be required.");
    this->required_set.insert(name);
  } else {
    // It is optional
    this->optional_set.insert(name);
  }

  return *this;
}

template <typename T>
ArgMap& ArgMap::arg(const std::string& name, T& value, const std::string& doc)
{
  arg(name, value);
  std::stringstream ss;
  ss << doc << " [ default=" << value << " ]";
  docVec.push_back({name, ss.str(), required_mode});

  return *this;
}

template <typename T>
ArgMap& ArgMap::arg(const std::string& name,
                    T& value,
                    const std::string& doc,
                    const char* info)
{
  arg(name, value);

  docVec.push_back({name, doc, required_mode});
  if (info != nullptr && info[0] != '\0')
    std::get<1>(docVec.back()).append(" [ default=").append(info).append(" ]");

  return *this;
}

template <typename C>
ArgMap& ArgMap::dots(C& container, const char* name)
{
  if (this->dots_enabled)
    throw LogicError(".dots() can only be called once.");

  this->dots_enabled = true;
  this->dots_name = name;

  // Have it out of the map as it may be called many times and has no
  // name/token.
  this->dots_ap = std::make_unique<ArgProcessorContainer<C>>(&container);

  return *this;
}

inline ArgMap& ArgMap::note(const std::string& s)
{
  std::get<1>(docVec.back()).append("\t\t").append(s);
  return *this;
}

inline void ArgMap::usage(const std::string& msg) const
{
  if (!msg.empty())
    std::cerr << msg << '\n';

  decltype(docVec) docVecCopy = docVec;
  std::stable_partition(
      docVecCopy.begin(),
      docVecCopy.end(),
      [this](const auto& item) {
        const auto& it = map.find(std::get<0>(item));
        if (it == map.end())
          throw LogicError("Not found in map '" + std::get<0>(item) + "'.");
        else
          return it->second->getArgType() != ArgType::POSITIONAL;
      });

  std::cerr << "Usage: " << this->progname;
  for (const auto& doc : docVecCopy) {
    std::string name_ext(std::get<0>(doc));
    auto it = map.find(name_ext);
    if (it == map.end())
      throw LogicError("Not found in map '" + name_ext + "'.");
    if (it->second->getArgType() == ArgType::NAMED)
      name_ext.append(1, this->kv_separator).append("<arg>");
    if (std::get<2>(doc)) {
      std::cerr << " " << name_ext;
    } else {
      std::cerr << " [" << name_ext << "]";
    }
  }

  if (this->dots_enabled)
    std::cerr << " [" << this->dots_name << " ...]\n";
  else
    std::cerr << "\n";

  std::cerr << doc() << std::endl;

  exit(EXIT_FAILURE);
}

inline ArgMap& ArgMap::helpArgs(const std::initializer_list<std::string> s)
{
  this->help_tokens = s;
  return *this;
}

inline ArgMap& ArgMap::helpArgs(const std::string s)
{
  this->help_tokens = {s};
  return *this;
}

inline ArgMap& ArgMap::separator(Separator s)
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
    throw LogicError("Unrecognised option for kv separator.");
  }

  return *this;
}

inline ArgMap& ArgMap::optional()
{
  this->required_mode = false;
  return *this;
}

inline ArgMap& ArgMap::required()
{
  this->required_mode = true;
  return *this;
}

inline ArgMap& ArgMap::toggle(bool t)
{
  this->named_args_only = false;
  this->arg_type = t ? ArgType::TOGGLE_TRUE : ArgType::TOGGLE_FALSE;
  return *this;
}

inline ArgMap& ArgMap::named()
{
  this->arg_type = ArgType::NAMED;
  return *this;
}

inline ArgMap& ArgMap::positional()
{
  this->named_args_only = false;
  this->arg_type = ArgType::POSITIONAL;
  return *this;
}

inline ArgMap& ArgMap::diagnostics(std::ostream& ostrm)
{
  this->diagnostics_strm = &ostrm;
  return *this;
}

inline void ArgMap::printDiagnostics(
    const std::forward_list<std::string>& args) const
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

inline std::string ArgMap::doc() const
{
  std::stringstream ss;
  auto maxSzElem =
      std::max_element(docVec.begin(),
                       docVec.end(),
                       [](const auto& x, const auto& y) {
                         return std::get<0>(x).size() < std::get<0>(y).size();
                       });

  for (const auto& p : docVec) {
    ss << "  " << std::left << std::setw(std::get<0>(*maxSzElem).length() + 1)
       << std::get<0>(p) << std::setw(0) << std::get<1>(p) << '\n';
  }

  return ss.str();
}

inline void ArgMap::simpleParse(const std::forward_list<std::string>& args,
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
        throw LogicError("Unrecognised ArgType.");
        break;
      }

      // Remove from required_set (if it is there)
      this->required_set.erase(token);

      // Previously called.
      this->previous_call_set.insert(token);
    } else if (pos_args_it != this->positional_args_list.end()) {
      // POSITIONAL args are treated differently as it is technically
      // never a recognised token.
      std::shared_ptr<ArgProcessor> pos_ap = map.find(*pos_args_it)->second;
      if (!pos_ap->process(*it))
        throw LogicError("Positional name does not match a ArgMap name.");
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

inline ArgMap& ArgMap::parse(int argc, char** argv)
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

inline ArgMap& ArgMap::parse(const std::string& filepath)
{

  if (this->kv_separator == ' ') { // Not from files.
    throw LogicError("Whitespace separator not possible from files.");
  }

  if (!this->named_args_only) {
    throw LogicError("Toggle and Positional arguments not possible from "
                     "files. Only named arguments.");
  }

  std::ifstream file(filepath);
  this->progname = filepath;

  if (!file.is_open()) {
    throw RuntimeError("Could not open file '" + filepath + "'.");
  }

  std::forward_list<std::string> args;
  auto it = args.before_begin();
  std::regex re_comment_lines(R"((^\s*\#)|(^\s+$))");
  std::string line;
  while (getline(file, line)) {
    if (line.empty() || std::regex_search(line, re_comment_lines)) {
      continue; // ignore comment lines and empties.
    }
    it = args.insert_after(it, line);
  }

  splitOnSeparator(args, this->kv_separator);

  // Take any leading and trailing whitespace away.
  std::for_each(args.begin(), args.end(), strip);

  printDiagnostics(args);

  simpleParse(args, false, [&filepath](const std::string& msg) {
    throw RuntimeError("Could not parse params file: '" + filepath + "'. " +
                       msg);
  });

  // Have the required args been provided - if not throw
  if (!this->required_set.empty()) {
    std::ostringstream oss;
    oss << "Required argument(s) not given:\n";
    for (const auto& e : this->required_set)
      oss << "\t" << e << '\n';
    throw RuntimeError(oss.str());
  }

  return *this;
}

} // namespace helib

#endif // ifndef HELIB_ARGMAP_H
