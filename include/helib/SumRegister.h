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

#ifndef HELIB_SUMREGISTER_H
#define HELIB_SUMREGISTER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

namespace helib {

/**
 * @class SumRegister
 * @brief Class to do a binary tree summation as results appear to keep memory
 * usage to a minimum.
 * @tparam T The type of object to be summed.
 **/
template <typename T>
class SumRegister
{
private:
  std::vector<std::unique_ptr<T>> intermediateResults;

  unsigned int maxNumOfInputs;
  unsigned int remainingInputs;
  bool flushRequiredFlag = false;
  bool resultFlag = false;
  size_t depth = 0;

public:
  /**
   * @brief Constructor
   * @param _maxNumOfInputs The maximum number expected to be summed together.
   **/
  SumRegister(unsigned int _maxNumOfInputs) :
      maxNumOfInputs(_maxNumOfInputs), remainingInputs(_maxNumOfInputs)
  {
    // log2 returns -inf if arg is zero!
    if (_maxNumOfInputs != 0)
      this->depth = std::ceil(std::log2(_maxNumOfInputs));

    // Set if NOT power of 2
    if ((_maxNumOfInputs & (_maxNumOfInputs - 1)) || (_maxNumOfInputs == 0)) {
      this->flushRequiredFlag = true;
    }

    // Default size will be 1.
    this->intermediateResults = std::vector<std::unique_ptr<T>>(depth + 1);
  }

  /**
   * @brief Add to the sum another object of type `T`.
   * @param t The object of type `T` to be added to the sum.
   * @note This method is destructive on the input argument.
   **/
  void add(std::unique_ptr<T>& t)
  {
    if (this->remainingInputs == 0) {
      return;
    }
    this->remainingInputs--;

    if (this->intermediateResults.at(0) != nullptr) {
      *this->intermediateResults.at(0) += *t;
      if (depth > 0) {
        for (size_t i = 1, j = 0; i < this->intermediateResults.size();
             i++, j++) {
          if (intermediateResults[i] != nullptr) {
            *this->intermediateResults[i] += *this->intermediateResults[j];
            this->intermediateResults[j].reset();
          } else {
            this->intermediateResults[i] =
                std::move(this->intermediateResults[j]);
            break;
          }
        }
      }

      if (intermediateResults.at(depth) != nullptr)
        resultFlag = true;

    } else {
      intermediateResults.at(0) = std::move(t);
    }

    if ((this->remainingInputs == 0) && this->flushRequiredFlag) {
      this->flush();
    }
  }

  /**
   * @brief Get the result of the summation.
   * @return The result of the summation as a `std::unique_ptr<T>`.
   * @note This is destructive on the result, thus should only be called once.
   **/
  std::unique_ptr<T> getResult()
  {
    return std::move(intermediateResults.at(depth));
  }

  /**
   * @brief Check result exists.
   * @return `true` if there is a result `false` otherwise.
   **/
  bool hasResult() const { return resultFlag; }

  /**
   * @brief Get depth of summation binary tree.
   * @return The size of the depth of the summation binary tree.
   **/
  size_t getDepth() const { return depth; }

  /**
   * @brief Flush the binary tree to force producing a result on current tree.
   **/
  void flush()
  {
    // Flushing with a result should do nothing.
    if (this->hasResult()) {
      return;
    }

    for (size_t i = 0, j = 1; i < depth; i++, j++) {
      if (this->intermediateResults[i] == nullptr)
        continue;

      if (this->intermediateResults[j] == nullptr) {
        this->intermediateResults[j] = std::move(this->intermediateResults[i]);
      } else {
        *this->intermediateResults[j] += *this->intermediateResults[i];
        this->intermediateResults[i].reset();
      }
    }

    if (intermediateResults.at(depth) != nullptr)
      resultFlag = true;
  }

  /**
   * @brief Print the information in the binary tree.
   **/
  void print()
  {
    std::cout << "Current values\n";
    for (size_t i = 0; i < this->intermediateResults.size(); i++) {
      std::cout << "[" << i << "]: "
                << (this->intermediateResults[i] != nullptr
                        ? *this->intermediateResults[i]
                        : 0)
                << " (" << this->intermediateResults[i] << ")" << '\n';
    }
  }

  /**
   * @brief Get the intermediate results.
   * @return The intermediate results in a `std::vector`.
   **/
  std::vector<std::unique_ptr<T>>& getIntermediates()
  {
    return intermediateResults;
  };
};

} // namespace helib
#endif // HELIB_SUMREGISTER_H
