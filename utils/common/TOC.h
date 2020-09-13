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

#ifndef TOC_H
#define TOC_H

#include <iostream>
#include <fstream>

// Table of Contents or Index or Catalog
class TOC
{

private:
  uint64_t rows;
  uint64_t cols;
  std::unique_ptr<uint64_t[]> idx;

public:
  TOC() = default;

  TOC(uint64_t r, uint64_t c) :
      rows(r), cols(c), idx(std::make_unique<uint64_t[]>(r * c))
  {}

  void print(std::ostream& out = std::cout) const
  {
    out << "Rows, Cols: " << rows << ", " << cols << std::endl;
    for (uint64_t i = 0; i < rows; ++i) {
      for (uint64_t j = 0; j < cols; ++j) {
        out << getIdx(i, j) << " ";
      }
      out << std::endl;
    }
  }

  uint64_t getRows() const { return this->rows; }

  uint64_t getCols() const { return this->cols; }

  uint64_t getIdx(int i, int j) const { return this->idx[i * cols + j]; }

  void setIdx(int i, int j, uint64_t value) { this->idx[i * cols + j] = value; }

  void write(std::ostream& s)
  {
    s.write(reinterpret_cast<char*>(&rows), sizeof(uint64_t));
    s.write(reinterpret_cast<char*>(&cols), sizeof(uint64_t));
    s.write(reinterpret_cast<char*>(idx.get()), sizeof(uint64_t) * rows * cols);
  }

  void read(std::istream& s)
  {
    s.read(reinterpret_cast<char*>(&rows), sizeof(uint64_t));
    s.read(reinterpret_cast<char*>(&cols), sizeof(uint64_t));
    // Deallocate existing, allocate new.
    idx.reset(new uint64_t[rows * cols]);
    s.read(reinterpret_cast<char*>(idx.get()), sizeof(uint64_t) * rows * cols);
  }

  long memorySize() { return sizeof(uint64_t) * (2 + rows * cols); }
};

#endif // TOC_H
