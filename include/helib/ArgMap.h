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
#include <forward_list>
#include <initializer_list>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <string>
#include <sstream>
#include <memory>
#include <vector>
#include <type_traits>
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
 *     .separator(ArgMap::Separator::WHITESPACE) // change seperator to
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
    POSITIONAL
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
        throw helib::LogicError(
            "Attempting to have argument '" + name +
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

  /* ArgProcessor: A virtual base class that acts as the interface to hold
   * args in a map of different types.
   */
  struct ArgProcessor
  {
    virtual ~ArgProcessor() = default;
    virtual ArgType getArgType() = 0;
    virtual bool process(const std::string& s) = 0;
  }; // end of ArgProcessor

  /* ArgProcessorDerived: templated subclasses */
  template <class T>
  class ArgProcessorDerived : public ArgProcessor
  {

  private:
    T* value;
    ArgType arg_type;

    // For strings. Avoids a stream breaking on whitespace.
    template <typename U = T,
              typename S,
              typename std::enable_if<std::is_same<U, S>::value, int>::type = 0>
    bool do_process(const S& s)
    {
      *value = s;
      return true;
    }

    // For non-string types, at the mercy of stringstream.
    template <
        typename U = T,
        typename S,
        typename std::enable_if<!std::is_same<U, S>::value, int>::type = 0>
    bool do_process(const S& s)
    {
      std::stringstream ss(s);
      return bool(ss >> *value);
    }

  public:
    ArgType getArgType() override { return arg_type; }

    bool process(const std::string& s) override { return this->do_process(s); }

    explicit ArgProcessorDerived(T* v, ArgType at) : value(v), arg_type(at) {}
  }; // end of ArgProcessorDerived

  char kv_separator = '=';
  std::string progname;

  // Track addresses to stop assigning same variable to more than one .arg(...)
  std::unordered_set<void*> addresses_used;

  // Track what has been called previously whilst parsing.
  std::unordered_set<std::string> previous_call_set;

  // Store the args.
  std::unordered_map<std::string, std::shared_ptr<ArgProcessor>> map;

  // Docs are held in a stream.
  std::stringstream docStream;

  PositionalArgsList positional_args_list;

  std::unordered_set<std::string> help_tokens = {"-h", "--help"};

  // Modes and other flags.
  bool required_mode = false;
  bool named_args_only = true;
  ArgType arg_type = ArgType::NAMED;

  std::ostream* diagnostics_strm = nullptr;

  // Set for tracking required.
  std::unordered_set<std::string> required_set;

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
  template <class T>
  ArgMap& arg(const char* name, T& value)
  {
    // has this name already been added?
    helib::assertTrue<helib::RuntimeError>(
        map[name] == nullptr,
        "Key already in arg map (key: " + std::string(name) + ")");

    // have we seen this addr before?
    helib::assertEq<helib::RuntimeError>(
        addresses_used.count(&value),
        0ul,
        "Attempting to register variable twice");

    addresses_used.insert(&value);

    map[name] = std::shared_ptr<ArgProcessor>(
        new ArgProcessorDerived<T>(&value, this->arg_type));

    if (this->arg_type == ArgType::POSITIONAL) {
      this->positional_args_list.insert(name, !this->required_mode);
    }

    if (this->required_mode) {
      this->required_set.insert(name);
    }

    return *this;
  }

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
  template <class T>
  ArgMap& arg(const char* name, T& value, const char* doc1)
  {
    arg(name, value);
    docStream << "\t" << name << " \t" << doc1 << "  [ default=" << value
              << " ]"
              << "\n";

    return *this;
  }

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
  template <class T>
  ArgMap& arg(const char* name, T& value, const char* doc1, const char* info)
  {
    arg(name, value);
    docStream << "\t" << name << " \t" << doc1;
    if (info != nullptr && info[0] != '\0')
      docStream << "  [ default=" << info << " ]"
                << "\n";
    else
      docStream << "\n";

    return *this;
  }

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
};

} // namespace helib

#endif // ifndef HELIB_ARGMAP_H
