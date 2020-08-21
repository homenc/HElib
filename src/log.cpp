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

#include <iostream>
#include <sstream>
#include <fstream>
#include <helib/assertions.h>
#include <helib/log.h>
#include <helib/apiAttributes.h>

namespace helib {

static bool checkDeletable(std::ostream* os)
{
  if (os != nullptr && os != &std::cerr)
    return true;
  return false;
}

Logger helog = []() -> Logger {
  Logger defaultLog;
  helog.setLogToFile("helib.log");
  return defaultLog;
}();

Logger::~Logger()
{
  if (checkDeletable(logStream_p))
    delete logStream_p;
}

void Logger::setLogToStderr()
{
  if (checkDeletable(logStream_p))
    delete logStream_p;

  logStream_p = &std::cerr;
}

void Logger::setLogToFile(const std::string& filepath, bool overwrite)
{
  if (checkDeletable(logStream_p))
    delete logStream_p;

  if (overwrite)
    logStream_p = new std::ofstream(filepath);
  else
    logStream_p = new std::ofstream(filepath, std::ios::app);

  assertNotNull<helib::IOError>(logStream_p,
                                "Could not open file '" + filepath +
                                    "' when setting the logStream.");
}

} // namespace helib
