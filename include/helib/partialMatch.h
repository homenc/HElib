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
/* Copyright (C) 2022 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
 *
 * Extended HElib to support the Not operator in queries.
 * Contributions include
 *
 * Added:
 *   Database:
 *                 template <typename TXT2>
              auto contains(const std::string& query_string, const Matrix<TXT2>&
 query_data) const
 */
#ifndef HELIB_PARTIALMATCH_H
#define HELIB_PARTIALMATCH_H

#include <sstream>

#include <helib/Matrix.h>
#include <helib/PolyMod.h>
#include <helib/query.h>

// This code is in flux and should be considered very alpha.
// Not recommended for public use.

namespace helib {

/**
 * @brief Given a query set and a database, calculates a mask of {0,1} where 1
 * signifies a matching element and 0 otherwise.
 * @tparam TXT type of the query set. Must be a `Ptxt` or `Ctxt`.
 * @param ea The encrypted array object holding information about the scheme.
 * @param query The query set to mask against the database. Must be a row
 * vector of the same dimension as the second dimension of the database matrix.
 * @param database The matrix holding the plaintext database.
 * @return The calculated mask. Is the same size as the database.
 * @note This is an overloaded function for when the database is not encrypted.
 **/
template <typename TXT>
inline Matrix<TXT> calculateMasks(const EncryptedArray& ea,
                                  Matrix<TXT> query,
                                  const Matrix<Ptxt<BGV>>& database)
{
  if (query.dims(0) != 1)
    throw InvalidArgument("Query must be a row vector");
  if (query.dims(1) != database.dims(1))
    throw InvalidArgument(
        "Database and query must have same number of columns");
  // TODO: case where query.dims(0) != database.dims(0)

  // Replicate the query once per row of the database
  // TODO: Some such replication will be needed once blocks/bands exist
  std::vector<long> columns(database.dims(0), 0l);
  Matrix<TXT>& mask = query;
  mask.inPlaceTranspose();
  mask = mask.columns(columns);
  mask.inPlaceTranspose();

  (mask -= database)
      .apply([&](auto& entry) { mapTo01(ea, entry); })
      .apply([](auto& entry) { entry.negate(); })
      .apply([](auto& entry) { entry.addConstant(NTL::ZZX(1l)); });

  return mask;
}

/**
 * @brief Given a query set and a database, calculates a mask of {0,1} where 1
 * signifies a matching element and 0 otherwise.
 * @tparam TXT type of the query set. Must be a `Ptxt` or `Ctxt`.
 * @param ea The encrypted array object holding information about the scheme.
 * @param query The query set to mask against the database. Must be a row
 * vector of the same dimension as the second dimension of the database matrix.
 * @param database The matrix holding the encrypted database.
 * @return The calculated mask. Is the same size as the database.
 * @note This is an overloaded function for when the database is encrypted.
 **/
template <typename TXT>
Matrix<Ctxt> calculateMasks(const EncryptedArray& ea,
                            Matrix<TXT> query,
                            const Matrix<Ctxt>& database)
{
  if (query.dims(0) != 1)
    throw InvalidArgument("Query must be a row vector");
  if (query.dims(1) != database.dims(1))
    throw InvalidArgument(
        "Database and query must have same number of columns");
  // TODO: case where query.dims(0) != database.dims(0)

  // Replicate the query once per row of the database
  // TODO: Some such replication will be needed once blocks/bands exist
  std::vector<long> columns(database.dims(0), 0l);
  Matrix<TXT>& mask = query;
  mask.inPlaceTranspose();
  mask = mask.columns(columns);
  mask.inPlaceTranspose();

  // FIXME: Avoid deep copy
  // Ptxt Query
  if constexpr (std::is_same_v<TXT, Ptxt<BGV>>) {
    auto tmp = database.deepCopy();
    (tmp -= mask)
        .apply([&](auto& entry) { mapTo01(ea, entry); })
        .apply([](auto& entry) { entry.negate(); })
        .apply([](auto& entry) { entry.addConstant(NTL::ZZX(1l)); });

    return tmp;
  } else { // Ctxt Query
    (mask -= database)
        .apply([&](auto& entry) { mapTo01(ea, entry); })
        .apply([](auto& entry) { entry.negate(); })
        .apply([](auto& entry) { entry.addConstant(NTL::ZZX(1l)); });

    return mask;
  }
}

/**
 * @brief Given a mask and information about the query to be performed,
 * calculates a score for each matching element signified by the mask.
 * @tparam TXT type of the mask matrix. Must be a `Ptxt` or `Ctxt`.
 * @param index_sets The set of indices signifying which columns of the mask
 * to query.
 * @param offsets The constant term to be added to the final score of each
 * queried column.
 * @param weights The weighted importance assigned to each queried column.
 * @param mask The mask with which to calculate the score from.
 * @return A single `Ctxt` or `Ptxt` containing the total score for each
 * queried column.
 **/
template <typename TXT>
inline Matrix<TXT> calculateScores(
    const std::vector<std::vector<long>> index_sets,
    const std::vector<long>& offsets,
    const std::vector<Matrix<long>>& weights,
    const Matrix<TXT>& mask)
{
  assertEq<InvalidArgument>(index_sets.size(),
                            offsets.size(),
                            "index_sets and offsets must have matching size");
  assertEq<InvalidArgument>(index_sets.size(),
                            weights.size(),
                            "index_sets and weights must have matching size");
  auto ones(mask(0, 0));
  ones.clear();
  ones.addConstant(NTL::ZZX(1L));
  Matrix<TXT> result(ones, mask.dims(0), 1l);
  for (std::size_t i = 0; i < index_sets.size(); ++i) {
    const auto& index_set = index_sets.at(i);
    const auto& weight_set = weights.at(i);
    long offset = offsets.at(i);

    assertEq<InvalidArgument>(
        weight_set.dims(0),
        index_set.size(),
        "found mismatch between index set size and weight set size");

    assertEq<InvalidArgument>(weight_set.dims(1),
                              1lu,
                              "all weight sets must be column vectors");

    Matrix<TXT> submatrix = mask.columns(index_set);
    Matrix<TXT> factor(submatrix * weight_set);
    // factor should in fact be a 1*1 matrix
    factor.apply([&](auto& entry) { entry.addConstant(NTL::ZZX(offset)); });
    result.template entrywiseOperation<TXT>(
        factor,
        [](auto& lhs, const auto& rhs) -> decltype(auto) {
          lhs.multiplyBy(rhs);
          return lhs;
        });
  }
  return result;
}

/**
 * @brief Given a value, encode the value across the coefficients of a
 * polynomial.
 * @param input The value of which to encode.
 * @param context The context object holding information on how to encode the
 * value.
 * @return A polynomial representing the encoded value.
 **/
inline PolyMod partialMatchEncode(uint32_t input, const Context& context)
{
  const long p = context.getP();
  std::vector<long> coeffs(context.getOrdP());
  // TODO - shouldn't keep checking input.
  for (long i = 0; i < long(coeffs.size()) && input != 0; ++i) {
    coeffs[i] = input % p;
    input /= p;
  }
  return PolyMod(coeffs, context.getSlotRing());
}

/**
 * @class Database
 * @tparam TXT The database is templated on `TXT` which can either be a `Ctxt`
 * or a `Ptxt<BGV>`
 * @brief An object representing a database which is a `HElib::Matrix<TXT>`.
 **/
template <typename TXT>
class Database
{
public:
  // FIXME: Generally, should Database own the Matrix uniquely?
  // Should we force good practice and ask that Context always be shared_ptr?

  // FIXME: Should probably move Matrix or make it unique_ptr or both?
  /**
   * @brief Constructor.
   * @param M The `Matrix<TXT>` containing the data of the database.
   * @param c A shared pointer to the context used to create the data.
   **/
  Database(const Matrix<TXT>& M, std::shared_ptr<const Context> c) :
      data(M), context(c)
  {}

  // FIXME: Should this option really exist?
  /**
   * @brief Constructor.
   * @param M The `Matrix<TXT>` containing the data of the database.
   * @param c The context object used to create the data.
   * @note This version accepts a `Context` that this object is not responsible
   * for i.e. if it is on the stack. The programmer is responsible in this case
   * for scope.
   **/
  Database(const Matrix<TXT>& M, const Context& c) :
      data(M),
      context(std::shared_ptr<const helib::Context>(&c, [](auto UNUSED p) {}))
  {}

  /**
   * @brief Overloaded function for performing a database lookup given a query
   *  expression and query data.
   * @tparam TXT2 The type of the query data, can be either a `Ctxt` or
   * `Ptxt<BGV>`.
   * @param lookup_query The lookup query expression to perform.
   * @param query_data The lookup query data to compare with the database.
   * @return A `Matrix<TXT>` containing 1s and 0s in slots where there was a
   * match or no match respectively. The return type `TXT` is `Ctxt` if either
   * query of database is encrypted, and `Ptxt<BGV>` otherwise.
   **/
  template <typename TXT2>
  auto contains(const QueryType& lookup_query,
                const Matrix<TXT2>& query_data) const;

  /**
   * @brief Overloaded function for performing a database lookup given a query
   *  string and query data.
   * @tparam TXT2 The type of the query data, can be either a `Ctxt` or
   * `Ptxt<BGV>`.
   * @param lookup_query A query string in Reverse Polish Notation with only
   *`ANDs` and `NOTs`.
   * @param query_data The lookup query data to compare with the database.
   * @return A `Matrix<TXT>` containing 1s and 0s in slots where there was a
   * match or no match respectively. The return type `TXT` is `Ctxt` if either
   * query of database is encrypted, and `Ptxt<BGV>` otherwise.
   **/
  template <typename TXT2>
  auto contains(const std::string& query_string,
                const Matrix<TXT2>& query_data) const;

  /**
   * @brief Function for performing a weighted partial match given a query
   * expression and query data.
   * @tparam TXT2 The type of the query data, can be either a `Ctxt` or
   * `Ptxt<BGV>`.
   * @param weighted_query The weighted lookup query expression to perform.
   * @param query_data The query data to compare with the database.
   * @return A `Matrix<TXT>` containing a score on weighted matches. The return
   * type `TXT` is `Ctxt` if either query of database is encrypted, and
   * `Ptxt<BGV>` otherwise.
   **/
  template <typename TXT2>
  auto getScore(const QueryType& weighted_query,
                const Matrix<TXT2>& query_data) const;

  // TODO - correct name?
  /**
   * @brief Returns number of columns in the database.
   * @return The number of columns in the database.
   **/
  long columns() const { return data.dims(1); }

  Matrix<TXT>& getData();

private:
  Matrix<TXT> data;
  std::shared_ptr<const Context> context;
};

template <typename TXT>
template <typename TXT2>
inline auto Database<TXT>::contains(const std::string& query_string,
                                    const Matrix<TXT2>& query_data) const
{
  // resolve type names: if both are ptxt, return ptxt, otherwise return ctxt
  using RTXT = typename std::conditional<(std::is_same_v<TXT, Ptxt<BGV>>)&&(
                                             std::is_same_v<TXT2, Ptxt<BGV>>),
                                         Ptxt<BGV>,
                                         Ctxt>::type;

  auto mask = calculateMasks(context->getEA(), query_data, this->data);

  std::istringstream input(query_string);

  std::stack<Matrix<RTXT>> txtStack;
  std::string symbol;

  while (input >> symbol) {
    if (symbol == "!") {
      auto& top = txtStack.top();
      top.apply([](auto& entry) { entry.negate(); });
      top.apply([](auto& entry) { entry.addConstant(NTL::ZZX(1L)); });
    } else if (symbol == "&&") {
      auto rhs = txtStack.top();
      txtStack.pop();

      auto& lhs = txtStack.top();
      lhs.template entrywiseOperation<RTXT>(
          rhs,
          [](auto& l, const auto& r) -> decltype(auto) {
            l.multiplyBy(r);
            return l;
          });
    } else {
      // If symbol is ||, throw a specific error
      assertFalse(symbol == "||",
                  "Cannot evaluate contains() on a string which contains Or "
                  "operator: first call removeOr()");
      // Otherwise, verify symbol is a number
      assertTrue(isNumber(symbol), "String is not a number: '" + symbol + "'");
      // And verify the number is in [0,...,columns - 1]
      assertTrue((std::stol(symbol) >= 0) &&
                     ((std::stol(symbol) < this->columns())),
                 "column is out of range");
      // push a copy of mask[symbol]. To take a deep copy, first make a matrix
      // of the right type and dimension with zeroes in every entry.
      auto ones(mask(0, 0));
      ones.clear();
      ones.addConstant(NTL::ZZX(0L));
      Matrix<RTXT> col(ones, mask.dims(0), 1l);
      // now add the entries we want to the RHS
      col += mask.getColumn(std::stol(symbol));
      txtStack.push(col);
    }
  }
  assertEq<LogicError>(1UL,
                       txtStack.size(),
                       "Size of stack after evaluation should be 1");
  return std::move(txtStack.top());
}

template <typename TXT>
template <typename TXT2>
inline auto Database<TXT>::contains(const QueryType& lookup_query,
                                    const Matrix<TXT2>& query_data) const
{
  auto result = getScore<TXT2>(lookup_query, query_data);

  if (lookup_query.containsOR) {
    // FLT on the scores
    result.apply([&](auto& txt) {
      txt.power(context->getAlMod().getPPowR() - 1);
      return txt;
    });
  }

  return result;
}

template <typename TXT>
template <typename TXT2>
inline auto Database<TXT>::getScore(const QueryType& weighted_query,
                                    const Matrix<TXT2>& query_data) const
{
  auto mask = calculateMasks(context->getEA(), query_data, this->data);

  auto result = calculateScores(weighted_query.Fs,
                                weighted_query.mus,
                                weighted_query.taus,
                                mask);

  return result;
}

template <typename TXT>
inline Matrix<TXT>& Database<TXT>::getData()
{
  return data;
}

} // namespace helib

#endif
