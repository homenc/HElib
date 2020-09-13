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
#ifndef HELIB_MATRIX_H
#define HELIB_MATRIX_H

#include <iostream>
#include <utility> // swap
#include <vector>
#include <numeric>
#include <array>
#include <string>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <initializer_list>

#include <NTL/BasicThreadPool.h>

#include "assertions.h"
#include "helib.h"
#include "zeroValue.h"

/*
 * The code in this file should be considered very alpha.
 * It is in flux, so not advised for public use, yet.
 * These APIs differ from matmul; here we think of vectors of Ctxts
 * not a Ctxt encoding a vector.
 */

namespace helib {

template <std::size_t N>
struct TensorSlice
{
  std::array<std::size_t, N> lengths; // the dims.
  std::array<std::size_t, N> strides;
  std::vector<long> start = {0};
  std::size_t size = 0;

  // Constructor used when taking a 'view'
  template <typename Iter1, typename Iter2>
  TensorSlice(const Iter1& firstLength,
              const Iter1& lastLength,
              const Iter2& firstStride,
              const Iter2& lastStride,
              const std::vector<long>& st) :
      start(st)
  {
    std::copy(firstLength, lastLength, this->lengths.begin());
    std::copy(firstStride, lastStride, this->strides.begin());
    this->size = (std::accumulate(lengths.begin(),
                                  lengths.end(),
                                  1,
                                  std::multiplies<std::size_t>()));
  }

  // Constructor used when taking a 'view'
  template <typename Iter1, typename Iter2>
  TensorSlice(const Iter1& firstLength,
              const Iter1& lastLength,
              const Iter2& firstStride,
              const Iter2& lastStride,
              unsigned long pos) :
      start(1, pos)
  {
    std::copy(firstLength, lastLength, this->lengths.begin());
    std::copy(firstStride, lastStride, this->strides.begin());
    this->size = (std::accumulate(lengths.begin(),
                                  lengths.end(),
                                  1,
                                  std::multiplies<std::size_t>()));
  }

  // Constructor for setting up a Tensor.
  template <typename... Dims>
  TensorSlice(Dims... dims) :
      lengths{std::size_t(dims)...},
      size(std::accumulate(lengths.begin(),
                           lengths.end(),
                           1,
                           std::multiplies<std::size_t>()))
  {
    // TODO add helib::assert or static checkhere w.r.t. N
    this->strides.back() = 1;
    if (N > 1) {
      std::copy(lengths.begin() + 1, lengths.end(), strides.begin());
    }
  }

  template <typename... Dims>
  std::size_t operator()(Dims... dims) const
  {
    static_assert(sizeof...(Dims) == N, "Wrong number of indices given.");

    std::array<std::size_t, N> args{std::size_t(dims)...};
    // Check indices are within bounds.
    for (long i = 0; i < long(N); ++i) {
      if (args[i] >= this->lengths[i]) {
        throw helib::OutOfRangeError(
            "Index given: " + std::to_string(args[i]) +
            ". Max value is: " + std::to_string(this->lengths[i]));
      }
    }

    if (this->start.size() == 1) {
      return std::inner_product(args.begin(),
                                args.end(),
                                this->strides.begin(),
                                this->start.at(0));
    } else {
      // FIXME - this will only work with Matrices
      return std::inner_product(args.begin(),
                                args.end(),
                                this->strides.begin(),
                                this->start.at(args.at(1)) * strides.back());
    }
  }

  std::size_t order() const { return N; }

  bool operator==(const TensorSlice& rhs) const
  {
    if (this == &rhs)
      return true;
    else if (this->size == rhs.size && this->start == rhs.start &&
             this->lengths == rhs.lengths && this->strides == rhs.strides)
      return true;
    else
      return false;
  }

  bool operator!=(const TensorSlice& rhs) const { return !(*this == rhs); }
};

// N is the dimension of the tensor
// At the moment we care about N = 1 or 2 only.

template <typename T, std::size_t N>
class Tensor
{

private:
  TensorSlice<N> subscripts;
  std::shared_ptr<std::vector<T>> elements_ptr;
  bool full_view = true;

public:
  // Construct with all entries equal to obj, often useful if T does not have a
  // default constructor or if default-constructed Ts are invalid.
  // This has to be disabled if dims is convertible to size_t for now,
  // otherwise a call to e.g.  Tensor(obj, 5) will use the below version
  // instead of this.  TODO: think of a way around this.
  template <typename U = T,
            typename... Dims,
            typename std::enable_if_t<
                !std::is_convertible<U, std::size_t>::value>* = nullptr>
  Tensor(const T& obj, Dims... dims) :
      subscripts{dims...},
      elements_ptr(std::make_shared<std::vector<T>>(subscripts.size, obj))
  {}

  template <typename... Dims>
  Tensor(Dims... dims) :
      subscripts{std::size_t(dims)...},
      elements_ptr(std::make_shared<std::vector<T>>(subscripts.size))
  {}

  // TODO generalise this
  Tensor(std::initializer_list<std::vector<T>> lst) :
      subscripts{lst.size(), lst.begin()->size()},
      elements_ptr(
          std::make_shared<std::vector<T>>(lst.size() * lst.begin()->size()))
  {
    int column_length = lst.begin()->size();
    int cnt = 0;
    for (const auto& v : lst) {
      if (column_length != long(v.size()))
        throw helib::LogicError(
            "Column dimensions do not match on initializer list.");

      std::copy(v.begin(),
                v.end(),
                this->elements_ptr->begin() + (column_length * cnt++));
    }
  }

  Tensor(const TensorSlice<N>& ts,
         const std::shared_ptr<std::vector<T>>& elems) :
      subscripts(ts), elements_ptr(elems), full_view(false)
  {}

  Tensor(const Tensor& other) = default;
  Tensor(Tensor&& other) = default;

  Tensor& operator=(const Tensor& rhs) = default;
  Tensor& operator=(Tensor&& rhs) = default;

  ~Tensor() = default;

  std::size_t order() const { return N; }

  template <typename... Args>
  T& operator()(Args... args)
  {
    return this->elements_ptr->at(subscripts(args...));
  }

  template <typename... Args>
  const T& operator()(Args... args) const
  {
    return this->elements_ptr->at(subscripts(args...));
  }

  std::size_t size() const { return this->subscripts.size; }

  std::size_t dims(int i) const { return this->subscripts.lengths.at(i); }

  bool fullView() const { return this->full_view; }

  bool operator==(const Tensor& rhs) const
  {
    if (this == &rhs) {
      return true;
    } else if (this->subscripts != rhs.subscripts) {
      return false;
    } else if (this->fullView() && rhs.fullView() &&
               *elements_ptr == *(rhs.elements_ptr)) {
      return true;
    } else {
      for (size_t i = 0; i < dims(0); ++i)
        for (size_t j = 0; j < dims(1); ++j)
          if (this->operator()(i, j) != rhs(i, j)) {
            return false;
          }
      return true;
    }
  }

  bool operator!=(const Tensor& rhs) const { return !(*this == rhs); }

  Tensor<T, N - 1> row(std::size_t i) const
  {
    TensorSlice<N - 1> ts(this->subscripts.lengths.begin() + 1,
                          this->subscripts.lengths.end(),
                          this->subscripts.strides.begin() + 1,
                          this->subscripts.strides.end(),
                          /*row X strides*/ i * this->subscripts.strides.at(0));
    return Tensor<T, N - 1>(ts, this->elements_ptr);
  }

  Tensor<T, N - 1> column(std::size_t j) const
  {
    // TODO - can only exist with N > 1
    TensorSlice<N - 1> ts(this->subscripts.lengths.begin(),
                          this->subscripts.lengths.end() - 1,
                          this->subscripts.strides.begin(),
                          this->subscripts.strides.end() - 1,
                          /*col*/ j);
    return Tensor<T, N - 1>(ts, this->elements_ptr);
  }

  // FIXME - returns several columns, same order as this. Only Matrices.
  Tensor<T, N> columns(const std::vector<long>& js) const
  {
    // Check the js are valid
    for (const auto& j : js)
      assertInRange<LogicError>(
          j,
          0l,
          static_cast<long>(this->dims(1)),
          "Index for column does not exist. Given index " + std::to_string(j) +
              ". Expected index in " + "range [0, " +
              std::to_string(this->dims(1)) + ").");

    // FIXME: lengths - works only for matrices at mo.
    std::vector<std::size_t> lengths = {this->dims(0), js.size()};

    std::vector<long> offsets(js);
    for (std::size_t i = 0; i < offsets.size(); ++i) {
      offsets[i] -= i;
    }

    TensorSlice<N> ts(lengths.begin(),
                      lengths.end(),
                      this->subscripts.strides.begin(),
                      this->subscripts.strides.end(),
                      offsets);

    return Tensor<T, N>(ts, this->elements_ptr);
  }

  template <typename T2>
  Tensor<T, N>& entrywiseOperation(const Tensor<T2, N>& rhs,
                                   std::function<T&(T&, const T2&)> operation)
  {
    // rhs is not of the same type thus need the dims.
    std::array<std::size_t, N> rhs_subscripts;
    for (std::size_t i = 0; i < N; ++i) {
      rhs_subscripts[i] = rhs.dims(i);
    }

    // Sanity Check: Are they the same dimensions?
    if (!std::equal(this->subscripts.lengths.begin(),
                    this->subscripts.lengths.end(),
                    rhs_subscripts.begin())) {
      throw helib::LogicError("Matrix dimensions do not match.");
    }

    // TODO For now, we do not allow views to operate on views to the same data.
    if (static_cast<const void*>(&this->data()) ==
        static_cast<const void*>(&rhs.data())) {
      throw helib::LogicError("Views point to same underlying data.");
    }

    // Optimisation if they have full view of underlying memmory.
    if (this->full_view && rhs.fullView()) {
      const std::vector<T2>& rhs_v = rhs.data();
      for (std::size_t i = 0; i < this->elements_ptr->size(); ++i) {
        operation((*this->elements_ptr)[i], rhs_v[i]);
      }
    } else {
      // TODO - again will only work for Matrices.
      for (std::size_t j = 0; j < this->dims(1); ++j) {
        for (std::size_t i = 0; i < this->dims(0); ++i) {
          operation(this->operator()(i, j), rhs(i, j));
        }
      }
    }

    return *this;
  }

  template <typename T2>
  Tensor<T, N>& operator+=(const Tensor<T2, N>& rhs)
  {
    return entrywiseOperation<T2>(
        rhs,
        [](auto& lhs, const auto& rhs) -> decltype(auto) {
          return lhs += rhs;
        });
  }

  template <typename T2>
  Tensor<T, N>& operator-=(const Tensor<T2, N>& rhs)
  {
    return entrywiseOperation<T2>(
        rhs,
        [](auto& lhs, const auto& rhs) -> decltype(auto) {
          return lhs -= rhs;
        });
  }

  template <typename T2>
  Tensor<T, N>& hadamard(const Tensor<T2, N>& rhs)
  {
    return entrywiseOperation<T2>(
        rhs,
        [](auto& lhs, const auto& rhs) -> decltype(auto) {
          return lhs *= rhs;
        });
  }

  Tensor<T, N>& apply(std::function<void(T& x)> fn)
  {
    // Optimisation if they have full view of underlying memory.
    if (this->full_view) {
      NTL_EXEC_RANGE(long(this->elements_ptr->size()), first, last)
      for (long i = first; i < last; ++i)
        fn((*elements_ptr)[i]);
      NTL_EXEC_RANGE_END

    } else {
      // TODO - again will only work for Matrices.
      NTL_EXEC_RANGE(this->dims(1), first, last)
      for (long j = first; j < last; ++j)
        for (std::size_t i = 0; i < this->dims(0); ++i)
          fn(this->operator()(i, j));
      NTL_EXEC_RANGE_END
    }
    return *this;
  }

  // Matrix special
  Tensor<T, 2>& transpose()
  {
    if (fullView()) {
      // In this case we can perform a transpose copy-free using swaps
      // First express the desired permutation as a list of numbers
      std::vector<int> permutation(size());
      std::iota(permutation.begin(), permutation.end(), 0);
      for (int& num : permutation)
        num = (num % subscripts.lengths[1]) * subscripts.lengths[0] +
              num / subscripts.lengths[1];
      // Now write the permutation as a list of disjoint cycles
      std::vector<std::vector<int>> cycles;
      std::vector<bool> seen(size(), false);
      int num_processed = 0;
      int current_pos = 0;
      while (num_processed < long(size())) {
        // Build the cycle starting at current_pos
        std::vector<int> cycle = {current_pos};
        seen[current_pos] = true;
        while (permutation.at(cycle.back()) != cycle.front()) {
          seen[permutation.at(cycle.back())] = true;
          cycle.push_back(permutation.at(cycle.back()));
        }
        num_processed += cycle.size();
        cycles.push_back(std::move(cycle));
        // Move current_pos to the next non-seen value
        while (current_pos < long(size()) && seen[current_pos])
          ++current_pos;
      }
      // Now turn this product of disjoint cycles into a list of pairs
      std::vector<std::pair<int, int>> swaps;
      for (const auto& cycle : cycles)
        if (cycle.size() >= 2)
          for (int i = cycle.size() - 1; i > 0; --i)
            swaps.emplace_back(cycle[i], cycle[i - 1]);
      // Now do the swaps
      for (const auto& swap : swaps)
        std::swap(elements_ptr->at(swap.first), elements_ptr->at(swap.second));
    } else {
      auto new_elements =
          std::make_shared<std::vector<T>>(size(), data().front());
      int count = size();
      int i = 0;
      int j = 0;
      int new_i = 0;
      int new_j = 0;
      while (count--) {
        // Put the new element in the correct place
        new_elements->at(new_i * subscripts.lengths[0] + new_j) = operator()(i,
                                                                             j);

        // Increment everything.
        // (i,j) wraps around as per the current configuration.
        j = (j + 1) % subscripts.lengths[1];
        if (j == 0)
          i = (i + 1) % subscripts.lengths[0];

        // (new_i,new_j) increments and wraps around as per the new, transposed
        // configuration.
        new_i = (new_i + 1) % subscripts.lengths[1];
        if (new_i == 0)
          new_j = (new_j + 1) % subscripts.lengths[0];
      }
      elements_ptr = new_elements;
      full_view = true;
    }
    std::reverse(subscripts.lengths.begin(), subscripts.lengths.end());
    subscripts.strides.front() = subscripts.lengths.back();
    subscripts.start = {0};
    return *this;
  }

  const std::vector<T>& data() const { return *this->elements_ptr; }
};

// Matrix special - Different types
template <typename T,
          typename T2,
          typename std::enable_if_t<
              std::is_convertible<T, std::size_t>::value>* = nullptr>
inline Tensor<T, 2> operator*(const Tensor<T, 2>& M1, const Tensor<T2, 2>& M2)
{
  HELIB_NTIMER_START(MatrixMultiplicationConv);
  // rows from M1 x cols from M2
  Tensor<T, 2> R(M1.dims(0), M2.dims(1));

  // Sanity check: M1 cols should match M2 rows
  if (M1.dims(1) != M2.dims(0)) {
    throw helib::LogicError("Matrix inner dimensions do not match.");
  }

  // TODO add NTL thread pool.
  // TODO swap += for some binary sum.
  NTL_EXEC_RANGE(M1.dims(0), first, last)
  for (long i = first; i < last; ++i)
    // for (std::size_t i = 0; i < M1.dims(0); ++i)
    for (std::size_t j = 0; j < M2.dims(1); ++j)
      for (std::size_t k = 0; k < M2.dims(0); ++k) {
        T R_tmp = M1(i, k);
        R_tmp *= M2(k, j);
        R(i, j) += R_tmp;
      }
  NTL_EXEC_RANGE_END

  HELIB_NTIMER_STOP(MatrixMultiplicationConv);
  return R;
}

// Matrix special - Different types
template <typename T,
          typename T2,
          typename std::enable_if_t<
              !std::is_convertible<T, std::size_t>::value>* = nullptr>
inline Tensor<T, 2> operator*(const Tensor<T, 2>& M1, const Tensor<T2, 2>& M2)
{
  HELIB_NTIMER_START(MatrixMultiplicationNotConv);
  // rows from M1 x cols from M2
  Tensor<T, 2> R(helib::zeroValue(M1(0, 0)), M1.dims(0), M2.dims(1));

  // Sanity check: M1 cols should match M2 rows
  if (M1.dims(1) != M2.dims(0)) {
    throw helib::LogicError(
        "The number of columns in left matrix (" + std::to_string(M1.dims(1)) +
        ") do not match the number of rows of the right matrix (" +
        std::to_string(M2.dims(0)) + ").");
  }

  // TODO swap += for some binary sum.
  NTL_EXEC_RANGE(M1.dims(0), first, last)
  for (long i = first; i < last; ++i)
    // For reference before NTL threads: for (std::size_t i = 0; i < M1.dims(0);
    // ++i)
    for (std::size_t j = 0; j < M2.dims(1); ++j)
      for (std::size_t k = 0; k < M2.dims(0); ++k) {
        T R_tmp = M1(i, k);
        R_tmp *= M2(k, j);
        R(i, j) += R_tmp;
      }
  NTL_EXEC_RANGE_END

  HELIB_NTIMER_STOP(MatrixMultiplicationNotConv);
  return R;
}

// These are the alias we will actually use.
template <typename T>
using Vector = Tensor<T, 1>;

template <typename T>
using Matrix = Tensor<T, 2>;

// TODO Columns Tuple
template <typename T>
class MatrixView
{

private:
  std::vector<Matrix<T>> columns;

public:
  MatrixView(const std::initializer_list<Matrix<T>> lst) : columns(lst) {}
};

template <typename T>
void printMatrix(const Matrix<T>& M, std::ostream& out = std::cout)
{
  for (std::size_t i = 0; i < M.dims(0); ++i) {
    for (std::size_t j = 0; j < M.dims(1); ++j)
      out << M(i, j) << " ";
    out << "\n";
  }
}

} // namespace helib

#endif // ifndef HELIB_MATRIX_H
