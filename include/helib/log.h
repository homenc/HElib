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

#ifndef HELIB_LOG_H
#define HELIB_LOG_H

#include <ostream>
#include <sstream>
#include <string>

/**
 * @file log.h
 * @brief The location of where the internal logging functions and classes are
 * defined.
 * @note HElib internal only. Not a designed as a general logging utility.
 * Subject to complete re-implementation.
 **/

namespace helib {

/**
 * @brief Logger class that handles warning printouts.
 **/
class Logger
{
private:
  // A stream pointer to be set for the logs to go to.
  std::ostream* logStream_p;

public:
  /**
   * @brief Set the logger object to write to `stderr`.
   **/
  void setLogToStderr();

  /**
   * @brief Set the logger object to write to specified file.
   * @param filepath The name of the file to write to.
   * @param overwrite Flag to tell the logger to overwrite the file.
   * @note Appends to file by default.
   **/
  void setLogToFile(const std::string& filepath, bool overwrite = false);

  /**
   * @brief Default constructor creates a logger object that does not point to
   * any target/destination.
   **/
  Logger() = default;

  /**
   * @brief Copy constructor, creates a copy of a logger object.
   **/
  Logger(const Logger& other) = default;

  /**
   * @brief Move constructor, can be used with `std::move` but does the same as
   * the copy constructor.
   **/
  Logger(Logger&& other) = default;

  /**
   * @brief Copy assignment operator, copies a logger object.
   **/
  Logger& operator=(Logger& other) = default;

  /**
   * @brief Move assignment operator, does the same as the copy assignment
   * operator.
   **/
  Logger& operator=(Logger&& other) = default;

  /**
   * @brief Destructor that closes and deletes the log stream object if
   * required i.e. if the log stream is a file.
   **/
  ~Logger();

  friend inline void Warning(const char* msg);
};

/**
 * @brief Internal global logger.
 **/
extern Logger helog;

/**
 * @brief Function for logging a warning message.
 * @param msg The warning message.
 **/
inline void Warning(const char* msg)
{
  *helog.logStream_p << "WARNING: " << msg << std::endl;
}

/**
 * @brief Function for logging a warning message.
 * @param msg The warning message.
 **/
inline void Warning(const std::string& msg) { Warning(msg.c_str()); }

// Errors should raise exceptions through throw or assertion functions
// in helib/assertions.h

// TODO Info Debug
// inline void Info(const std::string& msg)
// inline void Debug(const std::string& msg)

} // namespace helib

#endif // ifndef HELIB_LOG_H
