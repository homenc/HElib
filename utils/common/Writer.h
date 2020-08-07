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

#ifndef WRITER_H
#define WRITER_H

#include <string>
#include <fstream>
#include <memory>
#include <exception>

#include "TOC.h"

template <typename D>
class Writer
{

private:
  std::shared_ptr<TOC> toc;
  std::fstream writeStream;
  const std::string filepath;

  void allocFile(long sizeInBytes) const
  {
    if (sizeInBytes < 0)
      throw std::invalid_argument("Cannot allocate negative sized file. "
                                  "Size requested is " +
                                  std::to_string(sizeInBytes));

    std::ostringstream oss;
#if __linux__
    oss << "fallocate -l " << sizeInBytes << " " << filepath;
#elif __APPLE__
    oss << "mkfile -n " << sizeInBytes << " " << filepath;
#else
    throw std::logic_error("Unsupported system (not __APPLE__ or __linux__).");
#endif

    // return for system is implementation dependent
    std::system(oss.str().c_str());

    std::ifstream test(filepath, std::ios::ate);
    long actualFileSize = test.tellg();
    if (actualFileSize != sizeInBytes) {
      std::ostringstream err_msg;
      err_msg << "Could not allocate file '" << filepath
              << "'.\nRequested size " << sizeInBytes << ", actual size "
              << actualFileSize;
      throw std::runtime_error(err_msg.str());
    }
  }

public:
  // A fresh writer is responsible for setting the TOC.
  Writer(const std::string& fpath,
         uint64_t rows,
         uint64_t cols,
         long recordSizeInBytes) :
      // File needs alloc in body before stream is opened.
      toc(std::make_shared<TOC>(rows, cols)),
      filepath(fpath)
  {
    // Estimate for ctxt sizes.
    const long tocSize = toc->memorySize();
    const long fileSize = (cols * rows * recordSizeInBytes) + tocSize;
    allocFile(fileSize);

    // Set the TOC
    for (uint64_t j = 0; j < cols; ++j)
      for (uint64_t i = 0; i < rows; ++i)
        toc->setIdx(i, j, (i + j * rows) * recordSizeInBytes + tocSize);

    // Write it to file
    writeStream.open(fpath, std::ios::in | std::ios::out | std::ios::binary);
    if (!writeStream.is_open())
      throw std::runtime_error("Could not open '" + fpath +
                               "' for writing out TOC.");
    writeStream.seekp(0);
    toc->write(writeStream);
  }

  // A copied writer is not responsible for setting the TOC.
  Writer(const Writer& other) :
      toc(other.toc),
      writeStream(other.filepath,
                  std::ios::in | std::ios::out | std::ios::binary),
      filepath(other.filepath)
  {
    if (!writeStream.is_open())
      throw std::runtime_error("Could not open file '" + other.filepath +
                               "' for copied Writer.");
  }

  void writeByLocation(const D& data, uint64_t row, uint64_t col)
  {
    writeStream.seekp(toc->getIdx(row, col));
    data.write(writeStream);
  }

  TOC& getTOC() { return *toc; }
};

#endif // WRITER_H
