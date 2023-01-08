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
 * Extended HElib to support runtime queries and Table objects.
 * Contributions include
 *
 * Added:
 *   class Table "and its full contents"
 *   inline QueryType pseudoParser(const std::string& s)
 *   inline QueryType pseudoParserFromFile(const std::string& filename)
 *   inline void addSpacesAroundChars(std::stringstream& ss, const std::string&
 * s)
 *
 * Modified:
 *   Functions added to QueryBuilder class
 *     constructor QueryBuilder(const std::string& query)
 *     inline std::string convertToPostFix(const std::string& s)
 *     std::string getQueryString() const
 *
 *   Query_t -> QueryType
 *
 *   inline bool isNumber(const std::string& s) "now a free function"
 *
 *   QueryExpr is now a class instead of an alias
 *
 * Extended HElib to support the Not operator in queries.
 * Contributions include
 *
 * Added:
 *   class Not "and its full contents"
 *   inline QueryExpr operator!(const QueryExpr& p)
 *   QueryBuilder:
 *       void removeOr()
 *       void tidy(vecvec& expr) const
 *       vecvec negate(const vecvec& clause) const
 *       void tidyClause(std::vector<long>& clause) const
 *
 * Modified:
 *     QueryBuilder:
 *        Moved code out of QueryType build(long columns) const into a utility
 * function QueryType buildWeights(const vecvec& expr, const long columns) const
 *        modified the code in buildWeights()
 *        vecvec expandOr(const std::string& s) const
 */

#ifndef HELIB_QUERY_H
#define HELIB_QUERY_H

#include <helib/Matrix.h>

#include <stack>
#include <regex>
#include <algorithm> // std::count
#include <unordered_map>
#include <unordered_set>

/**
 * @file query.h
 * @brief Implementation of classes and functions related to HE queries.
 * @note This code is to be considered experimental and APIs can be expected to
 * change without prior warning.
 **/

namespace helib {

// Forward declarations
struct Expr;
class ColNumber;

/**
 * @brief A class wrapping a shared pointer to an `Expr` object.
 **/
class QueryExpr
{
public:
  std::shared_ptr<Expr> exp;
  /**
   * @brief Constructor.
   * @param p The `std::shared_ptr` to wrap into a `QueryExpr`.
   **/
  QueryExpr(std::shared_ptr<Expr> p) : exp(p) {}
};

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
  ColNumber(long c) : column(c) {}

private:
  long column;
};

/**
 * @brief Utility function for creating a `QueryExpr` for a specified column
 * in a query.
 * @param cl The index of the column to be used in the query.
 * @return `QueryExpr` with the specified `ColNumber` as the `Expr`.
 **/
inline QueryExpr makeQueryExpr(long cl)
{
  return QueryExpr(std::make_shared<ColNumber>(ColNumber(cl)));
}

/**
 * @class Not
 * @brief An object representing the logical NOT expression which inherits from
 * `Expr`.
 **/
class Not : public Expr
{
public:
  /**
   * @brief Function for returning the logical `NOT` expression in reverse
   * polish notation where the `NOT` operation is represented by `!` and each
   * operand is a column number.
   * @return A string representing the `NOT` expression in reverse polish
   * notation.
   **/
  std::string eval() const override { return p.exp->eval() + " !"; }
  /**
   * @brief Constructor.
   * @param exp The operand of the expression.
   **/
  Not(const QueryExpr& exp) : p(exp) {}

private:
  QueryExpr p;
};

/**
 * @class And
 * @brief An object representing the logical `AND` expression which inherits
 * from `Expr`.
 **/
class And : public Expr
{
public:
  /**
   * @brief Function for returning the logical `AND` expression in reverse
   * polish notation where the `AND` operation is represented by `&&` and each
   * operand is a column number.
   * @return A string representing the `AND` expression in reverse polish
   * notation.
   **/
  std::string eval() const override
  {
    return lhs.exp->eval() + " " + rhs.exp->eval() + " &&";
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
    return lhs.exp->eval() + " " + rhs.exp->eval() + " ||";
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
 * @brief Overloaded operator for creating a `QueryExpr` to a `NOT`
 * expression.
 * @param p Operand of the `NOT` expression.
 * @return `QueryExpr` with the `NOT` as the `Expr`.
 **/
inline QueryExpr operator!(const QueryExpr& p)
{
  return QueryExpr(std::make_shared<Not>(p));
}

/**
 * @brief Overloaded operator for creating a `QueryExpr` to an `AND`
 * expression.
 * @param lhs Left operand of the `AND` expression.
 * @param rhs Right operand of the `AND` expression.
 * @return `QueryExpr` with the `AND` as the `Expr`.
 **/
inline QueryExpr operator&&(const QueryExpr& lhs, const QueryExpr& rhs)
{
  return QueryExpr(std::make_shared<And>(lhs.exp, rhs.exp));
}

/**
 * @brief Overloaded operator for creating a `QueryExpr` of an `OR`
 * expression.
 * @param lhs Left operand of the `OR` expression.
 * @param rhs Right operand of the `OR` expression.
 * @return `QueryExpr` with the OR as the `Expr`.
 **/
inline QueryExpr operator||(const QueryExpr& lhs, const QueryExpr& rhs)
{
  return QueryExpr(std::make_shared<Or>(lhs.exp, rhs.exp));
}

/**
 * @brief Check if string is a number. i.e. All characters are digits.
 * @param s The input string to be checked.
 * @return bool `True` if all chars are digits, `False` otherwise.
 */
inline bool isNumber(const std::string& s)
{
  // Positive only
  return std::all_of(s.begin(), s.end(), ::isdigit);
}

/**
 * @brief Utility function that adds spaces around the chars '(', ')', ',' given
 * a string.
 * @param ss The output `std::stringstream` which will contain the space
 * separated items.
 * @param s The input string that contains characters to add spaces around.
 * @note This will add a trailing space if the string ends with one of the
 * specified characters.
 */
inline void addSpacesAroundChars(std::stringstream& ss, const std::string& s)
{
  std::unordered_set<char> set = {')', '(', ','};
  for (const auto& c : s) {
    if (set.count(c) == 1) {
      ss << " " << c << " ";
      continue;
    }
    ss << c;
  }
}

/**
 * @struct QueryType
 * @brief Structure containing all information required for an HE query.
 **/
struct QueryType
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
  QueryType(const std::vector<std::vector<long>>& index_sets,
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
  QueryType(std::vector<std::vector<long>>&& index_sets,
            std::vector<long>&& offsets,
            std::vector<Matrix<long>>&& weights,
            bool isThereAnOR) :
      Fs(index_sets), mus(offsets), taus(weights), containsOR(isThereAnOR)
  {}
};

/**
 * @class QueryBuilder
 * @brief An object used to construct a `QueryType` object from a logical
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
  QueryBuilder(const QueryExpr& expr) : query_str(expr.exp->eval()) {}

  /**
   * @brief Constructor.
   * @param expr Query string to build, in infix notation.
   * @note The expression is evaluated to reverse polish notation.
   **/
  QueryBuilder(const std::string& query) : query_str(convertToPostFix(query)) {}

  /**
   * @brief Function for building the `QueryType` object from the expression.
   * @param columns The total number of columns in the column set.
   * @return The resultant `QueryType` object containing information relating to
   * the query.
   **/
  QueryType build(long columns) const
  {
    // First convert the query to "type 1" by expanding out necessary ORs
    // to allow for zero ordering, (i+1) corresponds to i, negatives
    // correspond to not
    vecvec expr = expandOr(query_str);
    // then eliminate duplicates from inner clauses and any empty clauses
    vecvec tidyExpr = tidy(expr);
    // lastly form weights
    return this->buildWeights(tidyExpr, columns);
  }

  /**
   * @brief Function for replacing a `QueryBuilder` object with a logically
   * equivalent 'QueryBuilder` that uses only `AND`s and `NOT`s.
   */
  void removeOr()
  {
    std::stack<std::string> convertStack;

    std::istringstream input{query_str};
    std::string symbol;
    while (input >> symbol) {
      if (symbol == "&&") {
        auto rhs = convertStack.top();
        convertStack.pop();
        auto& lhs = convertStack.top();
        lhs += " " + rhs + " &&";
      } else if (symbol == "||") {
        // A B || is logically equivalent to A ! B ! && !
        auto rhs = convertStack.top();
        convertStack.pop();
        auto& lhs = convertStack.top();
        lhs += " ! " + rhs + " ! && !";
      } else if (symbol == "!") {
        convertStack.top() += " !";
      } else {
        // Operand: should be a number
        assertTrue(isNumber(symbol),
                   "String is not a number: '" + symbol + "'");
        convertStack.push(symbol);
      }
    }
    assertEq<LogicError>(1UL,
                         convertStack.size(),
                         "Size of stack after removeOr should be 1");
    query_str = std::move(convertStack.top());
  }

  /**
   * @brief Getter function that returns the `query_str`.
   * @return The string representation of the query.
   **/
  std::string getQueryString() const { return query_str; }

private:
  std::string query_str;

  /**
   * @brief Function that converts a space separated infix query with ANDs and
   * ORs to reverse polish notation.
   * @param s Input string in infix notation.
   * @return A new string in polish notation.
   **/
  inline std::string convertToPostFix(const std::string& s)
  {
    // checks for well formed query:
    // no consecutive arguments
    std::regex consecArgs(R"(\b(?!(AND|OR)\b)\w+\s+\b(?!(AND|OR)\b)\w+)");
    assertFalse(regex_search(s.begin(), s.end(), consecArgs),
                "invalid query string, found consecutive arguments");
    // No consecutive operators separated by spaces
    std::regex consecOps(R"(\b(AND|OR)+?\s+?(AND|OR)+?)");
    assertFalse(regex_search(s.begin(), s.end(), consecOps),
                "invalid query string, found consecutive operators");
    // no empty brackets, no operators directly before right bracket,
    // no operators directly after left bracket
    std::regex bracketsAndOps(R"((\(\s*(AND|OR|\))\s+)|(\b(AND|OR)\s*\)))");
    assertFalse(regex_search(s.begin(), s.end(), bracketsAndOps),
                "invalid query, binary operator next to bracket of wrong type");
    std::regex startOrEndOp(R"(^\s*(AND|OR)|(AND|OR)\s*$)");
    assertFalse(regex_search(s.begin(), s.end(), startOrEndOp),
                "invalid query, query starts or ends with an operator");
    // number of left brackets = number of right brackets
    assertEq<LogicError>(std::count(s.begin(), s.end(), '('),
                         std::count(s.begin(), s.end(), ')'),
                         "invalid query, unbalanced brackets");
    std::string ret = "";
    std::stack<char> st;

    // Add spaces in the input stream before and after brackets
    std::stringstream input;
    addSpacesAroundChars(input, s);
    input << " )";
    st.push('(');

    std::string symbol;
    while (input >> symbol) {
      if (symbol == "(") {
        st.push('(');
      } else if (symbol == "AND" || symbol == "OR") {
        if (st.top() == '&' || (symbol == "OR" && st.top() == '|')) {
          // cases where precedence has to be handled
          ret += (st.top() == '&') ? " &&" : " ||";
          st.pop(); // pop old operator
        }
        st.push((symbol == "AND") ? '&' : '|'); // add operator to stack
      } else if (symbol == ")") {               // empty stack until hit a "("
        while (st.top() != '(') {
          // add either " &&" or " ||"
          // Note: For some reason this doesn't work with +=
          ret = ret + " " + st.top() + st.top();
          st.pop();
        }
        st.pop(); // pop the "("
      } else {    // operands
        if (ret.length() != 0)
          ret += " ";
        ret += symbol;
      }
    }

    assertEq<LogicError>(
        0UL,
        st.size(),
        "error -- after conversion stack should be empty: query badly formed");

    return ret;
  }

  void printStack(std::stack<vecvec> stack) const
  {
    while (!stack.empty()) {
      printVecVec(stack.top());
      stack.pop();
    }
  }

  void printVecVec(const vecvec& vv) const
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

  /**
   * @brief Take a string in Reverse Polish Notation (RPN) and produce a
   * vector of vectors representing a logically equivalent `AND` of `OR`s.
   * @param s String in RPN, operators given by `&&`, `||` or `!`.
   * @return A `vecvec` representing a logically equivalent conjunction of
   * disjunctions of either columns or their negations. `(i + 1)` corresponds to
   * column `i`, and negatives correspond to negations of columns.
   */
  vecvec expandOr(const std::string& s) const
  {
    std::stack<vecvec> convertStack;

    std::istringstream input{s};

    std::string symbol;

    while (input >> symbol) {
      if (symbol == "&&") {
        // Squash the top into penultimate.
        auto op = convertStack.top();
        convertStack.pop();
        auto& top = convertStack.top();
        top.insert(top.end(), op.begin(), op.end());
      } else if (symbol == "||") {
        // Cartesian-esque product
        auto op1 = convertStack.top();
        convertStack.pop();
        auto op2 = convertStack.top();
        convertStack.pop();
        vecvec prod;
        prod.reserve(op1.size() * op2.size());
        for (const auto& i : op1)
          for (const auto& j : op2) {
            auto x = i; // Make the set i U j for every i, j
            x.insert(x.end(), j.begin(), j.end());
            prod.push_back(std::move(x));
          }
        convertStack.push(std::move(prod));
      } else if (symbol == "!") {
        vecvec top = convertStack.top();
        vecvec clause = negate(top); // negate top of stack
        convertStack.pop();          // pop
        convertStack.push(clause);   // push negation
      } else {
        // Operand: should be a number
        assertTrue(isNumber(symbol),
                   "String is not a number: '" + symbol + "'");
        convertStack.emplace(
            vecvec(1,
                   {std::stol(symbol) + 1})); // using 1 ordering temporarily:
        // To allow negatives to correspond to nots of columns even when there
        // is a zero column, column i is temporarily written as (i+1): In
        // particular, in this line, pushing back {{i + 1}} rather than {{i}}.
      }
    }

    // Now read answer off stack (should be size == 1).
    assertEq<LogicError>(1UL,
                         convertStack.size(),
                         "Size of stack after expandOr should be 1");

    return std::move(convertStack.top());
  }

  /**
   * @brief Given a `vecvec` where outer clauses correspond to `AND`s and inner
   * clauses correspond to `OR`s, deletes duplicated columns from inner clauses
   * (i.e., `i OR i`), deletes a column and its negation when both appear in an
   * inner clause (i.e., `i OR NOT i`), and lastly deletes empty clauses.
   * @param expr A `vecvec` corresponding to an `AND` of `OR`s.
   * @return A `vecvec` corresponding to a logically equivalent `AND` of `OR`s.
   */
  vecvec tidy(const vecvec& expr) const
  {
    vecvec x;
    for (const auto& y : expr) {
      auto z = tidyClause(y);
      if (z.size() != 0)
        x.emplace_back(z);
    }
    return x;
  }

  /**
   * @brief Build the weights for a query given as an `AND` of `OR`s.
   * @param expr Inner groups correspond to `OR`s, and we take the `AND` across
   * all inner groups. `(i + 1)` corresponds to column `i`, and negatives to
   * negations.
   * @param columns The number of columns in the database.
   * @return `QueryType` with the correct weights.
   */
  QueryType buildWeights(const vecvec& expr, const long columns) const
  {
    assertTrue(expr.size() > 0, "Cannot build weights for an empty query.");
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
    for (size_t i = 0; i < expr.size(); ++i) { // Each tau
      assertTrue(expr[i].size() > 0,
                 "empty clause found at index " + std::to_string(i) +
                     ", function called without Tidy().");
      mus[i] = 0;                 // Set mu to zero.
      Matrix<long> M(columns, 1); // Create temp tau matrix
      containsOR = (expr[i].size() > 1) ? true : containsOR;
      for (size_t j = 0; j < expr[i].size(); ++j) // Each column index
      {
        // column should not already have been seen in this clause
        assertTrue(M(std::abs(expr[i][j]) - 1, 0) == 0,
                   "Column " + std::to_string(std::abs(expr[i][j]) - 1) +
                       " repeated in clause " + std::to_string(i) +
                       ", function called without Tidy().");
        // add one to the offset for each NOT
        if (expr[i][j] < 0)
          mus[i] += 1;
        // mark NOT columns as -1, columns as 1
        M(std::abs(expr[i][j]) - 1, 0) = ((expr[i][j] >= 0) ? 1 : -1);
      }
      taus.push_back(std::move(M));
    }
    return QueryType(std::move(Fs),
                     std::move(mus),
                     std::move(taus),
                     containsOR);
  }

  /**
   * @brief Given a `vecvec` interepreted as an `AND` of `OR`s, return a
   * `vecvec` corresponding to the negation. Called by `expandOr()`.
   *
   * @param clause A `vecvec` of column numbers and negations of column numbers,
   * with `(i+1)` corresponding to column `i`.
   * @return A `vecvec` of column numbers and negations of column
   * numbers, with `(i+1)` corresponding to column `i`.
   */
  vecvec negate(const vecvec& clause) const
  {
    // We build the negation clause by clause
    vecvec notclause = {{}};
    for (const auto& i : clause) {
      // i = {a,b,c,...}. Negate this, to give {{-a},{-b},{-c},...}
      vecvec nextclause;
      for (const auto& j : i) {
        nextclause.push_back({-1 * j});
      }
      vecvec notclausetemp;
      notclausetemp.reserve(nextclause.size() * notclause.size());
      // then Cartesian style product
      for (const auto& j : notclause) {
        for (const auto& k : nextclause) {
          std::vector<long> x = j;
          x.insert(x.end(), k.begin(), k.end());
          notclausetemp.push_back(x);
        }
      }
      notclause = notclausetemp;
    }
    return notclause;
  }

  /**
   * @brief Tidies an inner clause by deleting duplicate columns and
   * removing columns for which both the column and its negation appear. Called
   * by `Tidy()`.
   * @param clause A vector of column labels and negative column labels,
   * interpreted as an `OR` clause, where `(i+1)` corresponds to column `i`, and
   * negatives correspond to negations.
   * @return A `std::vector<long>` of unique column labels.
   */
  std::vector<long> tidyClause(const std::vector<long>& clause) const
  {
    std::unordered_set<long> vars;
    std::vector<long> newclause;

    // remove duplicate columns
    std::copy_if(clause.begin(),
                 clause.end(),
                 std::back_inserter(newclause),
                 [&vars](long i) {
                   if (vars.count(i) == 0) {
                     vars.insert(i);
                     return true;
                   }
                   return false;
                 });

    // remove columns for which the negation appears in the clause
    std::remove_if(newclause.begin(), newclause.end(), [&vars](long i) {
      return (vars.count(-1 * i) > 0);
    });
    return newclause;
  }
};

/**
 * @class Table
 * @brief A class for converting cleartext queries with column labels to PSI
 * queries with column numbers.
 * Given a string describing a Table of the form
 * "TABLE (col1 (dim1), col2 (dim2) ...)"
 * a Table object will be created. This object will keep track of the column
 * names along with their respective indices and dimensions.
 */
class Table
{
public:
  /**
   * @brief Construct a new `Table` object from a table string.
   * @param table The input string listing the names and dimensions of the
   * columns of the table.
   */
  Table(const std::string& table) { buildTableFromStr(table); }

  /*************************
   * This section contains future possible public APIs for the Table class.
   * This is considered experimental.
   **************************/
  //  /**
  //   * @brief Append a new column to the end of the `Table`.
  //   * @param col The column label.
  //   * @param dim The ciphertext dimension of the new column. Default = 1.
  //   */
  //  void append(const std::string& col, long dim = 1L)
  //  {
  //    col_names.push_back(col);
  //    col_dims.push_back(dim);
  //    col_count += dim;
  //  }
  //
  //  /**
  //   * @brief Insert a new column at index.
  //   * @param index The `Table` index to insert the new column, using zero
  //   * ordering.
  //   * @param col The column label.
  //   * @param dim The ciphertext dimension of the new column. Default = 1.
  //   */
  //  void insert(long index, const std::string& col, long dim = 1L)
  //  {
  //    // sanitize/check index?
  //    col_names.insert(col_names.begin() + index, col);
  //    col_dims.insert(col_dims.begin() + index, dim);
  //    col_count += dim;
  //  }
  //
  //  /**
  //   * @brief Resize the width of a column.
  //   * @param index The column to resize.
  //   * @param newdim The new column width.
  //   */
  //  void resize(long index, long newdim)
  //  {
  //    assertTrue<LogicError>(newdim > 0,
  //                           "error, column width must be a positive
  //                           integer");
  //    col_count += (newdim - col_dims[index]);
  //    col_dims[index] = newdim;
  //  }
  //
  //  /**
  //   * @brief Erase a column at specified index.
  //   * @param index The `Table` index to erase.
  //   */
  //  void erase(long index)
  //  {
  //    // sanitize/check index?
  //    col_count -= col_dims[index];
  //    col_names.erase(col_names.begin() + index);
  //    col_dims.erase(col_dims.begin() + index);
  //  }
  /*************************************
   * End of experimental section.
   ************************************/

  /**
   * @brief Given a Table query string, convert this into a query expression
   * string by replacing the column names with indices.
   * i.e. col1 AND col2 OR col3 -> 0 AND 1 OR 2
   * @param s The input query string with column names.
   * @return A string with column numbers that can be provided to a
   * `QueryBuilder`.
   */
  std::string buildQueryString(const std::string& s) const
  {
    // checks for valid query
    // no consecutive operands
    std::regex consecArgs(R"(\b(?!(AND|OR)\b)\w+\s+\b(?!(AND|OR)\b)\w+)");
    assertFalse(regex_search(s.begin(), s.end(), consecArgs),
                "invalid query, found consecutive arguments");
    // No consecutive operators
    std::regex consecOps(R"(\b(AND|OR)+?\s+?(AND|OR)+?)");
    assertFalse(regex_search(s.begin(), s.end(), consecOps),
                "invalid query, found consecutive operators");
    // no empty brackets, no operators directly before right bracket,
    // no operators directly after left bracket
    std::regex bracketsAndOps(R"((\(\s*(AND|OR|\))\s+)|(\b(AND|OR)\s*\)))");
    assertFalse(regex_search(s.begin(), s.end(), bracketsAndOps),
                "invalid query, binary operator next to bracket of wrong type");
    std::regex startOrEndOp(R"(^\s*(AND|OR)|(AND|OR)\s*$)");
    assertFalse(regex_search(s.begin(), s.end(), startOrEndOp),
                "invalid query, query starts or ends with an operator");
    // number of left brackets = number of right brackets
    assertEq<LogicError>(std::count(s.begin(), s.end(), '('),
                         std::count(s.begin(), s.end(), ')'),
                         "invalid query, unbalanced brackets");
    std::stringstream ss;
    addSpacesAroundChars(ss, s);
    std::ostringstream out;
    std::string symbol;
    while (ss >> symbol) {
      if (operators.count(symbol) == 1) { // Operator found
        out << symbol << " ";
      } else if (symbol == "(" || symbol == ")") {
        out << symbol << " ";
      } else { // Must be an operand. Replace name with its index
        assertEq(columns.count(symbol),
                 1UL,
                 "Column name \'" + symbol + "\' not found.");
        out << columns.at(symbol).first << " ";
        // Append any additional columns if part of a composite column
        for (long i = 1; i < columns.at(symbol).second; ++i) {
          out << "AND " << columns.at(symbol).first + i << " ";
        }
      }
    }
    return out.str();
  }

  /**
   * @brief Return the number of column entries/names registered in the `Table`.
   * @return The number of column entries.
   */
  long getColCount() const { return columns.size(); }

  /**
   * @brief Return the total number of columns in the `Table`.
   * @return The sum of the dimensions of each column.
   */
  long size() const { return sz; }

  /**
   * @brief Print out the contents of the `Table` to a stream of the format
   * Key:[column name] Value:[index, dim].
   */
  void printTable(std::ostream& os = std::cout) const
  {
    // copy the map into a vector
    std::vector<std::pair<std::string, std::pair<long, long>>> column_list(
        columns.begin(),
        columns.end());
    // sort on index
    std::sort(column_list.begin(),
              column_list.end(),
              [](const std::pair<std::string, std::pair<long, long>>& a,
                 const std::pair<std::string, std::pair<long, long>>& b) {
                return a.second.first < b.second.first;
              });
    for (const auto& [key, value] : column_list) {
      os << "Key:[" << key << "] Value:[" << value.first << ", " << value.second
         << "]\n";
    }
  }

private:
  // Reserved key-words. Check potential column names against these words
  static const inline std::unordered_set<std::string> operators = {"AND",
                                                                   "OR",
                                                                   "NOT"};
  // key = column name, value = (index, dim)
  std::unordered_map<std::string, std::pair<long, long>> columns;
  long sz = 0; // Total size of the table (includes multi-dims)

  // Utility function that build the Table object from a string description
  void buildTableFromStr(const std::string& tableStr)
  {
    // Look for strings of the form TABLE (...) or TABLE(...) and capture
    // everything in the outer brackets in a capture group.
    std::regex grabTable(R"(^\s*TABLE\s*\((.+)\)\s*$)");
    std::smatch match;
    // If match was found, match[0] will be string that matched and match[1]
    // will be the capture group
    std::regex_search(tableStr, match, grabTable);
    // Check if string is has the TABLE label and a pair of brackets
    assertFalse<InvalidArgument>(match.empty(),
                                 "Invalid table string detected. Must be of the"
                                 " form \'TABLE(...)\', received \'" +
                                     tableStr + "\'");

    // This regex looks for strings with the pattern "name (dim)" where
    // number is optional and it captures "name" and "(dim)"
    std::regex grabNameAndDim(R"(^\s*([a-zA-Z_]+)\s*(.*))");
    long index = 0;
    // Adding a trailing comma to the string will trigger the assert in the
    // while loop if original string contained a trailing comma
    std::istringstream columnStr(match.str(1) + ",");
    std::string col;

    while (std::getline(columnStr, col, ',')) {
      // Check column entry exists i.e. something before a comma
      assertFalse<InvalidArgument>(col.empty(),
                                   "Empty column entry detected. Check table"
                                   " string.");
      std::smatch m; // match
      std::regex_search(col, m, grabNameAndDim);
      assertFalse<InvalidArgument>(m.empty(),
                                   "Invalid string detected. Columns must be a "
                                   "comma separated list of the form \'name "
                                   "(dimension)\' or \'name\'. Received \'" +
                                       col + "\'");
      isValidName(m.str(1)); // Check for valid column name

      // m[1] will be name and m[2] will be dim if it exists otherwise it will
      // be an empty string.
      // Set dimension to 1 if m[2] is empty string, otherwise set to match
      long dim = (m.str(2).empty()) ? 1 : isValidDimension(m.str(2));
      assertTrue<InvalidArgument>(dim > 0,
                                  "Invalid dimension value. Dimension must be "
                                  "an integer greater than 0.");
      columns.insert({m.str(1), std::make_pair(index, dim)});
      index += dim;
    }
    sz = index; // Final index will be total size of the table
  }

  // Given a string of alphabetical characters checks if it is a valid column
  // name
  void isValidName(const std::string& s)
  {
    // Check for empty string
    assertFalse<InvalidArgument>(s.empty(),
                                 "Invalid column name: valid names must be "
                                 "underscore separated alphabetical words. "
                                 "Received \'" +
                                     s + "\'");
    // Check the value is not one of the reserved keywords
    assertNeq<InvalidArgument>(operators.count(s),
                               1UL,
                               "Invalid column name: the following words are "
                               "reserved {AND, OR, NOT}.");
    // Check for duplicate column names
    assertEq<InvalidArgument>(columns.count(s),
                              0UL,
                              "Duplicate column name detected.");
  }

  // Given a string of the form (dim) checks if dim is valid and returns it
  long isValidDimension(const std::string& s)
  {
    // Strip trailing whitespace
    auto pos = s.find_last_not_of(" \t");
    auto dim = (pos == std::string::npos) ? s : s.substr(0, pos + 1);
    // Check the string is surrounded by brackets
    assertTrue<InvalidArgument>(dim[0] == '(' && dim[dim.size() - 1] == ')',
                                "Invalid dimension. All dimensions must be "
                                "enclosed in round brackets. Received \'" +
                                    dim + "\'");
    dim = dim.substr(1, dim.size() - 2); // Trim the brackets
    // Check for empty '()'
    assertFalse<InvalidArgument>(dim.empty(),
                                 "Invalid dimension. Brackets must enclose a "
                                 "value. Received \'" +
                                     s + "\'");
    // Check it is a number
    assertTrue<InvalidArgument>(isNumber(dim),
                                "Dimension can only be an integer greater "
                                "than 0. Received \'" +
                                    dim + "\'");
    return std::stoi(dim);
  }
};

/**
 * @brief Given a string, it searches for a Table descriptor and a query
 * expression that will be used to build a `QueryType`.
 * @param s Input string to parse.
 * @return A `QueryType` representing the query against the table described in
 * the input string.
 */
inline QueryType pseudoParser(const std::string& s)
{
  // Ignore comment lines and empty lines
  std::regex ignoreLines(R"((^\s*\#)|(^\s+$))");
  std::regex findTable(R"(^\s*TABLE\s*\(.+\)\s*$)");
  std::string line, tableStr, queryStr;
  std::istringstream is(s);

  // Break up the lines into tableStr and queryStr
  while (std::getline(is, line)) {
    std::smatch match;
    if (line.empty() || std::regex_search(line, ignoreLines))
      continue; // Ignore empty lines or comments
    else if (std::regex_search(line, match, findTable)) {
      assertTrue<InvalidArgument>(tableStr.empty(),
                                  "Multiple table strings found.");
      tableStr = match.str(0);
    } else { // Must be a query string
      assertTrue<InvalidArgument>(queryStr.empty(),
                                  "Multiple query strings found.");
      queryStr = line;
    }
  }
  assertFalse<InvalidArgument>(tableStr.empty(), "Table string not found.");
  assertFalse<InvalidArgument>(queryStr.empty(), "Query string not found.");

  Table t(tableStr);
  auto queryExprStr = t.buildQueryString(queryStr);
  return QueryBuilder(queryExprStr).build(/*columns=*/t.size());
}

/**
 * @brief Given a filename, read the file in and build a `QueryType` object.
 * @param filename Input file to read.
 * @return A `QueryType` representing the query against the table described in
 * the input file.
 */
inline QueryType pseudoParserFromFile(const std::string& filename)
{
  std::ifstream file(filename);
  std::ostringstream os;

  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      os << line << "\n";
    }
    file.close();
  } else {
    throw std::runtime_error("Could not open file " + filename);
  }
  return pseudoParser(os.str());
}

} // namespace helib

#endif // HELIB_QUERY_H
