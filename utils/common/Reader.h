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

#ifndef READER_H
#define READER_H

#include <vector>
#include <string>
#include <fstream>
#include <exception>

#include "TOC.h"

template <typename D>
class Reader
{

private:
  const std::string filepath;
  std::ifstream readStream;
  D& scratch;
  std::shared_ptr<TOC> toc;

public:
  Reader(const std::string& fname, D& init) :
      filepath(fname),
      readStream(filepath, std::ios::binary),
      scratch(init),
      toc(std::make_shared<TOC>())
  {
    if (!readStream.is_open())
      throw std::runtime_error("Could not open '" + filepath + "'.");
    toc->read(readStream);
  }

  Reader(const Reader& rdr) :
      filepath(rdr.filepath),
      readStream(filepath, std::ios::binary),
      scratch(rdr.scratch),
      toc(rdr.toc)
  {
    if (!readStream.is_open())
      throw std::runtime_error("Could not open '" + rdr.filepath + "'.");
  }

  void readDatum(D& dest, int i, int j)
  {
    if (readStream.eof())
      readStream.clear();

    readStream.seekg(toc->getIdx(i, j));
    dest.read(readStream);
  }

  std::unique_ptr<D> readDatum(int i, int j)
  {
    if (readStream.eof())
      readStream.clear();

    std::unique_ptr<D> ptr = std::make_unique<D>(scratch);
    readStream.seekg(toc->getIdx(i, j));
    ptr->read(readStream);

    return std::move(ptr);
  }

  std::unique_ptr<std::vector<std::vector<D>>> readAll()
  {
    if (readStream.eof())
      readStream.clear();

    auto m_ptr = std::make_unique<std::vector<std::vector<D>>>(
        toc->getRows(),
        std::vector<D>(toc->getCols(), scratch));

    for (int i = 0; i < toc->getRows(); i++) {
      for (int j = 0; j < toc->getCols(); j++) {
        readStream.seekg(toc->getIdx(i, j));
        (*m_ptr)[i][j].read(readStream);
      }
    }

    return std::move(m_ptr);
  }

  std::unique_ptr<std::vector<D>> readRow(int i)
  {

    if (readStream.eof())
      readStream.clear();

    auto v_ptr = std::make_unique<std::vector<D>>(toc->getCols(), scratch);
    for (int n = 0; n < toc->getCols(); n++) {
      readStream.seekg(toc->getIdx(i, n));
      (*v_ptr)[n].read(readStream);
    }

    return std::move(v_ptr);
  }

  std::unique_ptr<std::vector<D>> readCol(int j)
  {
    if (readStream.eof())
      readStream.clear();

    auto v_ptr = std::make_unique<std::vector<D>>(toc->getRows(), scratch);
    for (int n = 0; n < toc->getRows(); n++) {
      readStream.seekg(toc->getIdx(n, j));
      (*v_ptr)[n].read(readStream);
    }

    return std::move(v_ptr);
  }

  const TOC& getTOC() const { return *toc; }
};

template <typename D>
static std::vector<Reader<D>> createReaders(long count,
                                            const std::string& dataFilePath,
                                            D& dummy)
{
  std::vector<Reader<D>> readers;
  readers.reserve(count);
  for (long i = 0; i < count; ++i) {
    readers.emplace_back(dataFilePath.c_str(), dummy);
  }
  return readers;
}

#endif // READER_H
