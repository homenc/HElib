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

#ifndef HELIB_PARTIALMATCH_H
#define HELIB_PARTIALMATCH_H

#include <sstream>
#include <stack>

#include <helib/Matrix.h>
#include <helib/PolyMod.h>

// This code is in flux and should be considered bery alpha.
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
  mask.transpose();
  mask = mask.columns(columns);
  mask.transpose();

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
  mask.transpose();
  mask = mask.columns(columns);
  mask.transpose();

  (mask -= database)
      .apply([&](auto& entry) { mapTo01(ea, entry); })
      .apply([](auto& entry) { entry.negate(); })
      .apply([](auto& entry) { entry.addConstant(NTL::ZZX(1l)); });

  return mask;
}

/**
 * @brief Given a mask and information about the query to be performed,
 * calculates a score for each matching element signified by the mask.
 * @tparam TXT type of the mask matrix. Must be a `Ptxt` or `Ctxt`.
 * @param index_sets The set of indicies signifying which columns of the mask
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
  const long p = context.zMStar.getP();
  std::vector<long> coeffs(context.zMStar.getOrdP());
  // TODO - shouldn't keep checking input.
  for (long i = 0; i < long(coeffs.size()) && input != 0; ++i) {
    coeffs[i] = input % p;
    input /= p;
  }
  return PolyMod(coeffs, context.slotRing);
}

struct Expr;
class ColNumber;

/**
 * @brief An alias for a shared pointer to an `Expr` object.
 **/
using QueryExpr = std::shared_ptr<Expr>;

/**
 * @brief Utility function for creating a shared pointer to a specified column
 * in a query.
 * @param cl The index of the column to be used in the query.
 * @return Shared pointer to the class `ColNumber`.
 **/
inline std::shared_ptr<ColNumber> makeQueryExpr(int cl)
{
  return std::make_shared<ColNumber>(cl);
}

/**
 * @struct Expr
 * @brief Base structure for logical expressions.
 * @note This is pure virtual.
 **/
struct Expr
{
  virtual std::string eval() const = 0;
  virtual ~Expr() = default;
};

/**
 * @class ColNumber
 * @brief An object representing a column of a database as an expression which
 * inherits from `Expr`.
 **/
class ColNumber : public Expr
{
public:
  /**
   * @brief Function for returning the column number of the object.
   * @return A string representation of the column number.
   **/
  std::string eval() const override { return std::to_string(column); }

  /**
   * @brief Constructor.
   * @param c The column number.
   **/
  ColNumber(int c) : column(c) {}

private:
  int column;
};

/**
 * @class And
 * @brief An object representing the logical AND expression which inherits from
 * `Expr`.
 **/
class And : public Expr
{
public:
  /**
   * @brief Function for returning the logical AND expression in reverse polish
   * notation where the AND operation is represented by `&&` and each operand
   * is a column number.
   * @return A string representing the AND expression in reverse polish
   * notation.
   **/
  std::string eval() const override
  {
    return lhs->eval() + " " + rhs->eval() + " &&";
  }

  /**
   * @brief Constructor.
   * @param l The left operand of the expression.
   * @param r The right operand of the expression.
   **/
  And(const QueryExpr& l, const QueryExpr& r) : lhs(l), rhs(r) {}

private:
  QueryExpr lhs;
  QueryExpr rhs;
};

/**
 * @class Or
 * @brief An object representing the logical OR expression which inherits from
 * `Expr`.
 **/
class Or : public Expr
{
public:
  /**
   * @brief Function for returning the logical OR expression in reverse polish
   * notation where the OR operation is represented by `||` and each operand
   * is a column number.
   * @return A string representing the OR expression in reverse polish
   * notation.
   **/
  std::string eval() const override
  {
    return lhs->eval() + " " + rhs->eval() + " ||";
  }

  /**
   * @brief Constructor.
   * @param l The left operand of the expression.
   * @param r The right operand of the expression.
   **/
  Or(const QueryExpr& l, const QueryExpr& r) : lhs(l), rhs(r) {}

private:
  QueryExpr lhs;
  QueryExpr rhs;
};

/**
 * @brief Overloaded operator for creating a shared pointer to an AND
 * expression.
 * @param lhs Left operand of the AND expression.
 * @param rhs Right operand of the AND expression.
 * @return Shared pointer to the class `And`.
 **/
inline std::shared_ptr<And> operator&&(const QueryExpr& lhs,
                                       const QueryExpr& rhs)
{
  return std::make_shared<And>(lhs, rhs);
}

/**
 * @brief Overloaded operator for creating a shared pointer to an OR
 * expression.
 * @param lhs Left operand of the OR expression.
 * @param rhs Right operand of the OR expression.
 * @return Shared pointer to the class `Or`.
 **/
inline std::shared_ptr<Or> operator||(const QueryExpr& lhs,
                                      const QueryExpr& rhs)
{
  return std::make_shared<Or>(lhs, rhs);
}

/**
 * @struct Query_t
 * @brief Structure containing all information required for an HE query.
 **/
struct Query_t
{
  /**
   * @brief `std::vector` of index sets. These index sets specify the indexes
   * of the columns in each column subset.
   **/
  std::vector<std::vector<long>> Fs;

  /**
   * @brief `std::vector` of offsets. Each offset is a constant value. There
   * should be a single offset for each index set.
   **/
  std::vector<long> mus;

  /**
   * @brief `std::vector` of a set of weights. Each weight set corresponds to a
   * single index set where each individual weight corresponds to the index of
   * the index set.
   **/
  std::vector<Matrix<long>> taus;

  /**
   * @brief Flag indicating if the query contains a logical OR operation. This
   * is used for optimization purposes.
   **/
  bool containsOR = false;

  /**
   * @brief Constructor.
   * @param index_sets The set of column subsets.
   * @param offsets The set of offset constants.
   * @param weights The set of weight sets.
   * @param isThereAnOR Boolean value indicating if the query contains an OR
   * operation.
   **/
  Query_t(const std::vector<std::vector<long>>& index_sets,
          const std::vector<long>& offsets,
          const std::vector<Matrix<long>>& weights,
          const bool isThereAnOR) :
      Fs(index_sets), mus(offsets), taus(weights), containsOR(isThereAnOR)
  {}

  /**
   * @brief Constructor.
   * @param index_sets The set of column subsets.
   * @param offsets The set of offset constants.
   * @param weights The set of weight sets.
   * @param isThereAnOR Boolean value indicating if the query contains an OR
   * operation.
   **/
  Query_t(std::vector<std::vector<long>>&& index_sets,
          std::vector<long>&& offsets,
          std::vector<Matrix<long>>&& weights,
          bool isThereAnOR) :
      Fs(index_sets), mus(offsets), taus(weights), containsOR(isThereAnOR)
  {}
};

/**
 * @class QueryBuilder
 * @brief An object used to construct a `Query_t` object from a logical
 * expression.
 **/
class QueryBuilder
{

  // 'outer' vec are the and groups and 'inner' are the or groups
  using vecvec = std::vector<std::vector<long>>;

public:
  /**
   * @brief Constructor.
   * @param expr The logical expression to build.
   * @note The expression is evaluated to reverse polish notation.
   **/
  QueryBuilder(const QueryExpr& expr) : query_str(expr->eval()) {}

  /**
   * @brief Function for building the `Query_t` object from the expression.
   * @param columns The total number of columns in the column set.
   * @return The resultant `Query_t` object containing information relating to
   * the query.
   **/
  Query_t build(long columns) const
  {

    // Convert the query to "type 1" by expanding out necessary ORs
    vecvec expr = expandOr(query_str);
    bool containsOR = false;

    vecvec Fs(expr.size());
    {
      std::vector<long> v(columns);
      std::iota(v.begin(), v.end(), 0);
      std::fill(Fs.begin(), Fs.end(), v);
    }
    std::vector<long> mus(expr.size(), 1);
    std::vector<Matrix<long>> taus;
    taus.reserve(expr.size());

    // Create the taus
    for (long i = 0; i < long(expr.size()); ++i) { // Each tau
      mus[i] = 0;                                  // Set mu to zero.
      Matrix<long> M(columns, 1);                  // Create temp tau matrix
      containsOR = (expr[i].size() > 1) ? true : false;
      for (long j = 0; j < long(expr[i].size()); ++j) // Each column index
        M(expr[i][j], 0) = 1;                         // Mark those columns as 1
      taus.push_back(std::move(M));
    }

    return Query_t(std::move(Fs), std::move(mus), std::move(taus), containsOR);
  }

private:
  std::string query_str;

  void printStack(std::stack<vecvec> stack)
  {
    while (!stack.empty()) {
      printVecVec(stack.top());
      stack.pop();
    }
  }

  void printVecVec(const vecvec& vv)
  {
    for (const auto& v : vv) {
      std::cout << "[ ";
      for (const auto& e : v) {
        std::cout << e << " ";
      }
      std::cout << "]";
    }
    std::cout << "\n";
  }

  bool isNumber(const std::string& s) const
  {
    // Positive only
    return std::all_of(s.begin(), s.end(), ::isdigit);
  }

  vecvec expandOr(const std::string& s) const
  {
    std::stack<vecvec> convertStack;

    std::istringstream input{s};
    std::ostringstream output{};

    std::string symbol;

    while (input >> symbol) {
      if (!symbol.compare("&&")) {
        // Squash the top into penultimate.
        auto op = convertStack.top();
        convertStack.pop();
        auto& top = convertStack.top();
        top.insert(top.end(), op.begin(), op.end());
      } else if (!symbol.compare("||")) {
        // Cartesian-esque product
        auto op1 = convertStack.top();
        convertStack.pop();
        auto op2 = convertStack.top();
        convertStack.pop();

        vecvec prod;
        prod.reserve(op1.size() * op2.size());
        for (const auto& i : op1)
          for (const auto& j : op2) {
            auto x = i;
            x.insert(x.end(), j.begin(), j.end());
            prod.push_back(std::move(x));
          }

        convertStack.push(std::move(prod));
      } else {
        // Assume it is a number. But sanity check anyway.
        assertTrue(isNumber(symbol),
                   "String is not a number: '" + symbol + "'");
        convertStack.emplace(vecvec(1, {std::stol(symbol)}));
      }
    }

    // Now read answer off stack (should be size == 1).
    assertEq<LogicError>(1UL,
                         convertStack.size(),
                         "Size of stack after expandOr should be 1");

    return std::move(convertStack.top());
  }
};

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

  // FIXME: Combination of TXT = ctxt and TXT2 = ptxt does not work
  /**
   * @brief Function for performing a database lookup given a query expression
   * and query data.
   * @tparam TXT2 The type of the query data, can be either a `Ctxt` or
   * `Ptxt<BGV>`.
   * @param lookup_query The lookup query expression to perform.
   * @param query_data The lookup query data to compare with the database.
   * @return A `Matrix<TXT2>` containing 1s and 0s in slots where there was a
   * match or no match respectively.
   **/
  template <typename TXT2>
  Matrix<TXT2> contains(const Query_t& lookup_query,
                        const Matrix<TXT2>& query_data) const;

  // FIXME: Combination of TXT = ctxt and TXT2 = ptxt does not work
  /**
   * @brief Function for performing a weighted partial match given a query
   * expression and query data.
   * @tparam TXT2 The type of the query data, can be either a `Ctxt` or
   * `Ptxt<BGV>`.
   * @param weighted_query The weighted lookup query expression to perform.
   * @param query_data The query data to compare with the database.
   * @return A `Matrix<TXT2>` containing a score on weighted matches.
   **/
  template <typename TXT2>
  Matrix<TXT2> getScore(const Query_t& weighted_query,
                        const Matrix<TXT2>& query_data) const;

  // TODO - correct name?
  /**
   * @brief Returns number of columns in the database.
   * @return The number of columns in the database.
   **/
  long columns() { return data.dims(1); }

private:
  Matrix<TXT> data;
  std::shared_ptr<const Context> context;
};

template <typename TXT>
template <typename TXT2>
inline Matrix<TXT2> Database<TXT>::contains(
    const Query_t& lookup_query,
    const Matrix<TXT2>& query_data) const
{
  auto result = getScore<TXT2>(lookup_query, query_data);

  if (lookup_query.containsOR) {
    // FLT on the scores
    result.apply([&](auto& txt) {
      txt.power(context->alMod.getPPowR() - 1);
      return txt;
    });
  }

  return result;
}

template <typename TXT>
template <typename TXT2>
inline Matrix<TXT2> Database<TXT>::getScore(
    const Query_t& weighted_query,
    const Matrix<TXT2>& query_data) const
{
  auto mask = calculateMasks(*(context->ea), query_data, this->data);

  auto result = calculateScores(weighted_query.Fs,
                                weighted_query.mus,
                                weighted_query.taus,
                                mask);

  return result;
}

} // namespace helib

#endif
