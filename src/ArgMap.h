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

#ifndef HELIB_ARGMAP_H
#define HELIB_ARGMAP_H

#include <unordered_map>
#include <functional>
#include <string>
#include <sstream>
#include <memory>
#include <set>
#include "assertions.h"

//! @brief Easier arg parsing
/**
 * Example use:
 *
 *   // Variables to be set by command line.
 *   long p = 2;                               // default values.
 *   float f = 5.5;
 *   std::string k = "Hello World";
 *
 *   ArgMap()
 *     .required()                             // set args to required
 *     .arg("p", p, "doc for p")               //
 *     .arg("m", m, "doc for m", "undefined")  // special default info
 *     .optional()                             // back to optional args
 *     .arg("k", k, "doc for k", "")           // no default info
 *     .note("an extra note")                  // add extra doc/note
 *     .parse(argc, argv);       // parses and overrides initial values
 *
 **/

namespace helib {

class ArgMap {
private:
  /* ArgProcessor: virtual base class */

  class ArgProcessor {

  public:
    virtual ~ArgProcessor() = default;
    virtual bool process(const std::string& s) = 0;
  }; // end of ArgProcessor

  /* ArgProcessorDerived: templated subclasses */

  template <class T>
  class ArgProcessorDerived : public ArgProcessor {

  private:
    T* value;

  public:
    bool process(const std::string& s) override {
      std::stringstream ss(s);
      return bool(ss >> *value);
    }

    explicit ArgProcessorDerived(T* _value) : value(_value) {}
  }; // end of ArgProcessorDerived

  char kv_separator = '=';
  std::string progname;

  std::set<void*> addresses_used;
  std::set<std::string> previous_call_set;

  // map to store the args.
  std::unordered_map<std::string, std::shared_ptr<ArgProcessor>> map;
  std::stringstream docStream;
  bool required_flag = false;
  std::set<std::string> required_set;

  //! @brief Arg parsing function
  /**
   * Regex-based low level arg parsing function (private)
   * @param line Space-separated argument line to parse
   * @param duplicates If true does not fail in case of duplicated arguments
   * @param stop Callback function called in case of failure. (Default is Usage)
   */
  void simpleRegexParse(const std::string& line,
                        bool duplicates = true,
                        std::function<void(const std::string&)> stop = {});

public:
  //! @brief Add a new argument description
  /**
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
  ArgMap& arg(const char* name, T& value) {
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

    map[name] =
        std::shared_ptr<ArgProcessor>(new ArgProcessorDerived<T>(&value));
    if (required_flag) {
      required_set.insert(name);
    }

    return *this;
  }

  //! @brief Add a new argument with docs
  /**
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
  ArgMap& arg(const char* name, T& value, const char* doc1) {
    arg(name, value);
    docStream << "\t" << name << " \t" << doc1 << "  [ default=" << value
              << " ]"
              << "\n";

    return *this;
  }

  //! @brief Add a new argument with docs and default description
  /**
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
  ArgMap& arg(const char* name, T& value, const char* doc1, const char* info) {
    arg(name, value);
    docStream << "\t" << name << " \t" << doc1;
    if (info != nullptr && info[0] != '\0')
      docStream << "  [ default=" << info << " ]"
                << "\n";
    else
      docStream << "\n";

    return *this;
  }

  //! @brief Parse the argv array
  /**
   * Parse the argv array
   * If it fails or -h is an argument it prints the usage and exits the program
   * @param argc number of entries in argv
   * @param argv array containing the arguments
   * @return A reference to the ArgMap object
   */
  ArgMap& parse(int argc, char** argv);

  //! @brief Parse the config file
  /**
   * Parse the config file
   * Throws RuntimeError on failure
   * @param filepath the config file path
   * @return A reference to the ArgMap object
   */
  ArgMap& parse(const std::string& filepath);

  //! @brief Swaps to optional arg mode
  /**
   * Swaps to optional arg mode. Following arguments will be considered optional
   * @return A reference to the ArgMap object
   */
  ArgMap& optional();

  //! @brief Swaps to required arg mode
  /**
   * Swaps to required arg mode. Following arguments will be considered required
   * @return A reference to the ArgMap object
   */
  ArgMap& required();

  //! @brief Sets the key-value separator
  /**
   * Sets the key-value separator character
   * @param c the separator character. It must be a non-whitespace character.
   * @return A reference to the ArgMap object
   */
  ArgMap& kvSep(char c);

  //! @brief Adds a note to usage
  /**
   * Adds a note to the arg usage description
   * @param s The note string
   * @return A reference to the ArgMap object
   */
  ArgMap& note(const std::string& s);

  //! @brief Print usage and exit
  /**
   * Prints the usage and exits the program
   * @param msg An additional message to print before showing usage
   */
  void usage(const std::string& msg = "") const;

  //! @brief Return arg docs
  /**
   * Returns the argument documentation as a string
   * @return the argument documentation string
   */
  std::string doc() const;
};

}

#endif // ifndef HELIB_ARGMAP_H
