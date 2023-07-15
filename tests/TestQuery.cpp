// Copyright (C) 2022 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#include <helib/helib.h>
#include <helib/query.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

TEST(TestQuery, QueryBuilderFromStringInfixToPolishConversion)
{
  std::array<helib::QueryBuilder, 14> input_exprs{
      helib::QueryBuilder("a AND b"),
      helib::QueryBuilder("a OR b"),
      helib::QueryBuilder("( b OR c )"),
      helib::QueryBuilder("a AND ( b OR c )"),
      helib::QueryBuilder("( a AND b ) OR c"),
      helib::QueryBuilder("( a AND b ) OR ( c AND d )"),
      helib::QueryBuilder("( a AND b ) OR ( ( c AND d ) AND ( e OR f ) )"),
      /* test spaces, order of precedence, triple char operands */
      helib::QueryBuilder("a AND b AND c"),
      helib::QueryBuilder("( a OR b AND c )"),
      helib::QueryBuilder("a AND b OR c"),
      helib::QueryBuilder("aaa OR bbc OR ac"),
      helib::QueryBuilder("( a AND b ) OR ( c AND d ) AND ( e OR f )"),
      /* test brackets */
      helib::QueryBuilder("(a AND b) OR (c AND d)"),
      helib::QueryBuilder("(a AND (b AND (c AND d)))")};
  std::array<std::string, 14> expected{"a b &&",
                                       "a b ||",
                                       "b c ||",
                                       "a b c || &&",
                                       "a b && c ||",
                                       "a b && c d && ||",
                                       "a b && c d && e f || && ||",
                                       "a b && c &&",
                                       "a b c && ||",
                                       "a b && c ||",
                                       "aaa bbc || ac ||",
                                       "a b && c d && e f || && ||",
                                       "a b && c d && ||",
                                       "a b c d && && &&"};

  for (size_t i = 0; i < input_exprs.size(); ++i) {
    EXPECT_EQ(input_exprs[i].getQueryString(), expected[i]) << "*** i = " << i;
  }
}

TEST(TestQuery, ExpectedThrowsInfixToPostFixNotation)
{
  // throws because left bracket immediately before AND/OR
  EXPECT_THROW(helib::QueryBuilder("(AND a OR b"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder("( OR a OR b"), helib::LogicError);
  // throws because right bracket immediately after AND/OR
  EXPECT_THROW(helib::QueryBuilder("( a AND a OR )"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder("(a AND a AND)"), helib::LogicError);
  // throws because OR next to AND
  EXPECT_THROW(helib::QueryBuilder("c OR AND b"), helib::LogicError);
  // throws because left bracket next to right bracket
  EXPECT_THROW(helib::QueryBuilder("() c AND b"), helib::LogicError);
  // throws because consecutive operands
  EXPECT_THROW(helib::QueryBuilder("a b c"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder("date of birth AND c"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder("aAND c OR b"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder("c aAND OR b"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder("c ORb"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder(" ORb c"), helib::LogicError);
  // throws because unequal number of left and right brackets
  EXPECT_THROW(helib::QueryBuilder("(a AND b))"), helib::LogicError);
  EXPECT_THROW(helib::QueryBuilder("((b)"), helib::LogicError);
  // throws because expression starts with AND
  EXPECT_THROW(helib::QueryBuilder("AND a AND b"), helib::LogicError);
  // throws because expression ends with OR
  EXPECT_THROW(helib::QueryBuilder("y OR "), helib::LogicError);
}

TEST(TestQuery, TableFromString1DColumns)
{
  std::vector<helib::Table> tables;
  std::vector<std::string> expected_tablestrings;
  tables.emplace_back(helib::Table("TABLE(name,age,weight)"));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 1]\nKey:[weight] Value:[2, 1]\n");
  // a table with a space between TABLE and left bracket
  tables.emplace_back(helib::Table("TABLE (name,age,weight)"));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 1]\nKey:[weight] Value:[2, 1]\n");

  // spaces and tabs before and after
  tables.emplace_back(helib::Table("  TABLE (name,age,weight) "));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 1]\nKey:[weight] Value:[2, 1]\n");
  tables.emplace_back(helib::Table("  TABLE (name,age,weight)     "));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 1]\nKey:[weight] Value:[2, 1]\n");
  tables.emplace_back(helib::Table("        TABLE (name,age,weight)  "));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 1]\nKey:[weight] Value:[2, 1]\n");

  // a table with keywords as substrings of column names
  tables.emplace_back(helib::Table("TABLE(FLOOR, LANDS, NOTE, OTHER, and)"));
  expected_tablestrings.push_back(
      "Key:[FLOOR] Value:[0, 1]\nKey:[LANDS] "
      "Value:[1, 1]\nKey:[NOTE] Value:[2, 1]\nKey:[OTHER] Value:[3, "
      "1]\nKey:[and] Value:[4, 1]\n");

  for (std::size_t i = 0; i < tables.size(); i++) {
    std::ostringstream tablestring;
    tables[i].printTable(tablestring);
    EXPECT_EQ(tablestring.str(), expected_tablestrings[i]) << "*** i=" << i;
  }
}

TEST(TestQuery, TableFromStringMultiDimensionalColumns)
{
  std::vector<helib::Table> tables;
  std::vector<std::string> expected_tablestrings;

  tables.emplace_back(helib::Table("TABLE(name (2) ,age (3),weight (2))"));
  expected_tablestrings.push_back("Key:[name] Value:[0, 2]\nKey:[age] "
                                  "Value:[2, 3]\nKey:[weight] Value:[5, 2]\n");

  tables.emplace_back(helib::Table("TABLE(name,age (4),weight(1))"));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 4]\nKey:[weight] Value:[5, 1]\n");

  tables.emplace_back(helib::Table("  TABLE   (name ,age (4),weight(1))   "));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 4]\nKey:[weight] Value:[5, 1]\n");

  tables.emplace_back(helib::Table("TABLE(  name ,age   (4) ,   weight(1)  )"));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 4]\nKey:[weight] Value:[5, 1]\n");

  tables.emplace_back(helib::Table("TABLE(name,age(4),weight(1))"));
  expected_tablestrings.push_back("Key:[name] Value:[0, 1]\nKey:[age] "
                                  "Value:[1, 4]\nKey:[weight] Value:[5, 1]\n");

  tables.emplace_back(helib::Table("TABLE(first_name,age(4),weight_in_kg(1))"));
  expected_tablestrings.push_back(
      "Key:[first_name] Value:[0, 1]\nKey:[age] "
      "Value:[1, 4]\nKey:[weight_in_kg] Value:[5, 1]\n");

  for (std::size_t i = 0; i < tables.size(); i++) {
    std::ostringstream tablestring;
    tables[i].printTable(tablestring);
    EXPECT_EQ(tablestring.str(), expected_tablestrings[i]) << "*** i=" << i;
  }
}

TEST(TestQuery, ExpectedThrowsTableFromString)
{
  // Throw for empty string
  EXPECT_THROW(helib::Table(""), helib::InvalidArgument);
  // Throw for string does not start with TABLE(
  EXPECT_THROW(helib::Table("Tab(a,b,c)"), helib::InvalidArgument);
  // Throw for string does not end with )
  EXPECT_THROW(helib::Table("TABLE(a,b,c"), helib::InvalidArgument);
  // Throw for empty table
  EXPECT_THROW(helib::Table("TABLE()"), helib::InvalidArgument);
  // Throw for comma syntax
  EXPECT_THROW(helib::Table("TABLE(a,b,)"), helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE(,a,b)"), helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE(a,,b)"), helib::InvalidArgument);
  // Throw for numbers in column names
  EXPECT_THROW(helib::Table("TABLE(Col_1,col_two)"), helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE(Col_one, col2)"), helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE(firstCol, 2ndCol)"), helib::InvalidArgument);
  // Throw for repeated column
  EXPECT_THROW(helib::Table("TABLE(a,b,a)"), helib::InvalidArgument);
  // Throw for multiword column
  EXPECT_THROW(helib::Table("TABLE(a b, c, d)"), helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE(a (1) b (2), c, d)"),
               helib::InvalidArgument);
  // Throw for brackets in column name
  EXPECT_THROW(helib::Table("TABLE((a), c, d)"), helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE((a,b,c))"), helib::InvalidArgument);
  // Throw for column name with AND, OR, and NOT
  EXPECT_THROW(helib::Table("TABLE(AND (1), OR"), helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE(NOT (4)"), helib::InvalidArgument);
  // Throw for non-positive dimensions
  EXPECT_THROW(helib::Table("TABLE (name (-1), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name, age(0) , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name, age, weight (0.5))"),
               helib::InvalidArgument);
  // Throw for dimensions that are not numbers, empty dimensions, and multiple
  // dimensions
  EXPECT_THROW(helib::Table("TABLE (name (dim(7)), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name (), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name (2), age(dim) , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name (2)(7), age, weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name (!), age , weight)"),
               helib::InvalidArgument);
  // Throw for dimensions with incomplete/too many brackets and spaces inside
  // the brackets
  EXPECT_THROW(helib::Table("TABLE (name (1 ), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name ( 2 ), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name 3, age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name (4, age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name 5), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name ((6), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name ((7)), age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table("TABLE (name (1 2), age , weight)"),
               helib::InvalidArgument);
  // Throw for text before Table, noncomment text after Table, two tables in one
  // line
  EXPECT_THROW(helib::Table("this is a table TABLE (name, age , weight)"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::Table(" TABLE (name, age , weight) this is a table"),
               helib::InvalidArgument);
  EXPECT_THROW(
      helib::Table(" TABLE (name, age , weight) TABLE (name, age , weight)"),
      helib::InvalidArgument);
}

TEST(TestQuery, BuildQueryStringMultidimensionalColumns)
{
  helib::Table t0("   TABLE   (a (1) ,b ,  c (2),d)  ");
  std::array<std::string, 7> input_exprs{
      t0.buildQueryString("a AND b"),
      t0.buildQueryString("a OR b"),
      t0.buildQueryString("( b OR c )"),
      t0.buildQueryString("a OR b AND c"),
      t0.buildQueryString("( a OR b ) AND c"),
      t0.buildQueryString("( a AND b ) OR ( c AND d )"),
      t0.buildQueryString("(a AND b AND c AND d)")};
  std::array<std::string, 7> expected{"0 AND 1 ",
                                      "0 OR 1 ",
                                      "( 1 OR 2 AND 3 ) ",
                                      "0 OR 1 AND 2 AND 3 ",
                                      "( 0 OR 1 ) AND 2 AND 3 ",
                                      "( 0 AND 1 ) OR ( 2 AND 3 AND 4 ) ",
                                      "( 0 AND 1 AND 2 AND 3 AND 4 ) "};
  for (size_t i = 0; i < input_exprs.size(); ++i) {
    EXPECT_EQ(input_exprs[i], expected[i]) << "*** i = " << i;
  }
}

TEST(TestQuery, ExpectedThrowsBuildQueryString)
{
  helib::Table t("TABLE(a, b, c)");
  // Throw for consecutive operands
  EXPECT_THROW(t.buildQueryString("a b AND"), helib::LogicError);
  // Throws for consecutive operators and operators followed by right bracket,
  // operators preceded by left bracket
  EXPECT_THROW(t.buildQueryString("a AND AND b"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("(a AND) OR (b AND c)"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("(OR a)"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("() c AND b"), helib::LogicError);
  // Throws for query starting or ending with operator
  EXPECT_THROW(t.buildQueryString("a AND"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("OR a AND b"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("a OR "), helib::LogicError);
  // Throws for unbalanced brackets
  EXPECT_THROW(t.buildQueryString("(a AND b))"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("( (a OR c )"), helib::LogicError);
  // Throws for unrecognised operators (also handles no spaces around AND or OR)
  // Throws for unrecognized operands
  EXPECT_THROW(t.buildQueryString("a AND z"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("a ANDb"), helib::LogicError);
  EXPECT_THROW(t.buildQueryString("c ORb"), helib::LogicError);
}

TEST(TestQuery, BuildQueryBuilderFromTable1DColumns)
{
  helib::Table t("TABLE(name, age, weight, height)");
  std::array<helib::QueryBuilder, 7> input_exprs{
      helib::QueryBuilder(t.buildQueryString("name AND age")),
      helib::QueryBuilder(t.buildQueryString("name OR age")),
      helib::QueryBuilder(t.buildQueryString("( age OR weight )")),
      helib::QueryBuilder(t.buildQueryString("name OR age AND weight")),
      helib::QueryBuilder(t.buildQueryString("( name OR age ) AND weight")),
      helib::QueryBuilder(
          t.buildQueryString("( name AND age ) OR ( weight AND height )")),
      helib::QueryBuilder(
          t.buildQueryString("(name AND age AND weight AND height)"))};
  std::array<std::string, 7> expected{"0 1 &&",
                                      "0 1 ||",
                                      "1 2 ||",
                                      "0 1 2 && ||",
                                      "0 1 || 2 &&",
                                      "0 1 && 2 3 && ||",
                                      "0 1 && 2 && 3 &&"};
  for (size_t i = 0; i < input_exprs.size(); ++i) {
    EXPECT_EQ(input_exprs[i].getQueryString(), expected[i]) << "*** i= " << i;
  }
}

TEST(TestQuery, BuildQueryBuilderFromTableMultidimensionalColumns)
{
  // Multi-width columns
  helib::Table t1("TABLE(country (2), capital, population (1), currency (2))");

  std::array<helib::QueryBuilder, 7> input_exprs{
      helib::QueryBuilder(t1.buildQueryString("country AND capital")),
      helib::QueryBuilder(t1.buildQueryString("country OR capital")),
      helib::QueryBuilder(t1.buildQueryString("( capital OR population )")),
      helib::QueryBuilder(
          t1.buildQueryString("country OR capital AND population")),
      helib::QueryBuilder(
          t1.buildQueryString("( country OR capital ) AND population")),
      helib::QueryBuilder(t1.buildQueryString(
          "( country AND capital ) OR ( population AND currency )")),
      helib::QueryBuilder(t1.buildQueryString(
          "(country AND capital AND population AND currency)"))};

  std::array<std::string, 7> expected{"0 1 && 2 &&",
                                      "0 1 && 2 ||",
                                      "2 3 ||",
                                      "0 1 && 2 3 && ||",
                                      "0 1 && 2 || 3 &&",
                                      "0 1 && 2 && 3 4 && 5 && ||",
                                      "0 1 && 2 && 3 && 4 && 5 &&"};
  for (size_t i = 0; i < input_exprs.size(); ++i) {
    EXPECT_EQ(input_exprs[i].getQueryString(), expected[i]) << "*** i= " << i;
  }
}

TEST(TestQuery, QueryExprFromPseudoParser)
{
  long cases = 2;
  std::vector<std::string> tableStrings;
  std::vector<helib::QueryType> queries;
  std::vector<std::vector<std::vector<long>>> expected_Fs_vector;
  std::vector<std::vector<helib::Matrix<long>>> expected_taus_vector;
  std::vector<std::vector<long>> expected_mus_vector;
  std::vector<long> column_sizes = {3, 6};
  tableStrings.reserve(cases);
  queries.reserve(cases);
  expected_Fs_vector.reserve(cases);
  expected_taus_vector.reserve(cases);
  expected_mus_vector.reserve(cases);
  expected_mus_vector.reserve(cases);
  std::string s("TABLE(name,age,weight)\n"
                "name AND age OR weight\n"
                "\n"
                "# This is a comment line!\n"
                "\n");
  queries.push_back(helib::pseudoParser(s));
  expected_Fs_vector.push_back({{0, 1, 2}, {0, 1, 2}});
  expected_taus_vector.push_back({{{1}, {0}, {1}}, {{0}, {1}, {1}}});
  expected_mus_vector.push_back({0, 0});

  std::string s1("\n"
                 "TABLE(name(2),age (3),weight)\n"
                 "name AND age OR weight\n"
                 "# This is a comment line!\n");
  queries.push_back(helib::pseudoParser(s1));
  expected_Fs_vector.push_back({{0, 1, 2, 3, 4, 5},
                                {0, 1, 2, 3, 4, 5},
                                {0, 1, 2, 3, 4, 5},
                                {0, 1, 2, 3, 4, 5},
                                {0, 1, 2, 3, 4, 5}});
  expected_taus_vector.push_back({{{1}, {0}, {0}, {0}, {0}, {1}},
                                  {{0}, {1}, {0}, {0}, {0}, {1}},
                                  {{0}, {0}, {1}, {0}, {0}, {1}},
                                  {{0}, {0}, {0}, {1}, {0}, {1}},
                                  {{0}, {0}, {0}, {0}, {1}, {1}}});
  expected_mus_vector.push_back({0, 0, 0, 0, 0});

  for (long j = 0; j < cases; j++) {
    EXPECT_EQ(expected_Fs_vector[j].size(), queries[j].Fs.size())
        << "*** j = " << j;
    EXPECT_EQ(expected_mus_vector[j].size(), queries[j].mus.size())
        << "*** j = " << j;
    EXPECT_EQ(expected_taus_vector[j].size(), queries[j].taus.size())
        << "*** j = " << j;
    EXPECT_EQ(column_sizes[j], queries[j].taus[0].size()) << "*** j = " << j;
    for (size_t i = 0; i < expected_Fs_vector[j].size(); ++i)
      EXPECT_EQ(expected_Fs_vector[j][i], queries[j].Fs[i])
          << "*** j = " << j << ", i = " << i;
    for (size_t i = 0; i < expected_mus_vector[j].size(); ++i)
      EXPECT_EQ(expected_mus_vector[j][i], queries[j].mus[i])
          << "*** j = " << j << ", i= " << i;
    for (size_t i = 0; i < expected_taus_vector[j].size(); ++i)
      EXPECT_TRUE(expected_taus_vector[j][i] == queries[j].taus[i])
          << "*** j = " << j << ", i= " << i;
  }
}

TEST(TestQuery, ExpectedThrowsPseudoParser)
{
  // no query string provided
  EXPECT_THROW(helib::pseudoParser("TABLE(a, b, c)"), helib::InvalidArgument);
  EXPECT_THROW(helib::pseudoParser("TABLE(a, b, c)\n"
                                   "#this is a comment line"),
               helib::InvalidArgument);
  // no table string provided
  EXPECT_THROW(helib::pseudoParser("a AND b AND c"), helib::InvalidArgument);
  EXPECT_THROW(helib::pseudoParser("a AND b AND c\n"
                                   "#this is comment 1\n"
                                   "#this is comment 2"),
               helib::InvalidArgument);
  // multiple table strings provided
  EXPECT_THROW(helib::pseudoParser("TABLE(a, b, c)\n"
                                   "TABLE(x,y,z)\n"
                                   "x AND y OR z"),
               helib::InvalidArgument);
  // multiple query strings provided
  EXPECT_THROW(helib::pseudoParser("TABLE(a, b, c)\n"
                                   "a AND b AND c\n"
                                   "a AND b OR c\n"),
               helib::InvalidArgument);
  // two table entries on one line
  EXPECT_THROW(helib::pseudoParser("TABLE(a, b, c) TABLE(x, y, z)\n"
                                   "a AND c OR b\n"),
               helib::InvalidArgument);
  // words before or after table
  EXPECT_THROW(helib::pseudoParser("this is a table    TABLE(a, b, c) \n"
                                   "a AND b OR c\n"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::pseudoParser(" TABLE(a, b, c) this is a table\n"
                                   "a AND b OR c\n"),
               helib::InvalidArgument);
  // commented out table or query string
  EXPECT_THROW(helib::pseudoParser("# TABLE(a, b, c)\n"
                                   "# this is a comment\n"
                                   "a AND b OR c\n"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::pseudoParser("# this is a comment\n"
                                   " TABLE(a, b, c)\n"
                                   "# a AND b OR c\n"),
               helib::InvalidArgument);
  // Empty string or only comments/empty lines
  EXPECT_THROW(helib::pseudoParser("   \n"), helib::InvalidArgument);
  EXPECT_THROW(helib::pseudoParser("\n"
                                   " \n"),
               helib::InvalidArgument);
  EXPECT_THROW(helib::pseudoParser("# this is a comment\n"
                                   " \n"),
               helib::InvalidArgument);
}

TEST(TestQuery, QueryExprGeneratesPostFix)
{
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  const helib::QueryExpr& weight = helib::makeQueryExpr(3);

  helib::QueryExpr res = ((name && age) || height);
  EXPECT_EQ("0 1 && 2 ||", res.exp->eval());

  res = name || age || height;
  EXPECT_EQ("0 1 || 2 ||", res.exp->eval());

  res = (name || age) && height;
  EXPECT_EQ("0 1 || 2 &&", res.exp->eval());

  res = (name || age) && (name || height);
  EXPECT_EQ("0 1 || 0 2 || &&", res.exp->eval());

  res = (name && age) || (age && height);
  EXPECT_EQ("0 1 && 1 2 && ||", res.exp->eval());

  res = name && (age || (height && weight));
  EXPECT_EQ("0 1 2 3 && || &&", res.exp->eval());
}

TEST(TestQuery, QueryExprGeneratesPostFixWithNot)
{
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  const helib::QueryExpr& weight = helib::makeQueryExpr(3);

  helib::QueryExpr res = name || !name;
  EXPECT_EQ("0 0 ! ||", res.exp->eval());

  res = !(name || age || height);
  EXPECT_EQ("0 1 || 2 || !", res.exp->eval());

  res = !(name || age) && height;
  EXPECT_EQ("0 1 || ! 2 &&", res.exp->eval());

  res = !!name;
  EXPECT_EQ("0 ! !", res.exp->eval());

  res = age || !!weight;
  EXPECT_EQ("1 3 ! ! ||", res.exp->eval());

  res = !(name && age && height);
  EXPECT_EQ("0 1 && 2 && !", res.exp->eval());

  res = !((name || age) && (height || weight));
  EXPECT_EQ("0 1 || 2 3 || && !", res.exp->eval());
}

TEST(TestQuery, containsOrFlagInBuild)
{
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  long columns = 4;
  helib::QueryBuilder res1((name || age) && height);
  helib::QueryType query = res1.build(columns);
  EXPECT_EQ(query.containsOR, true);

  helib::QueryBuilder res2(height && (name || age));
  query = res2.build(columns);
  EXPECT_EQ(query.containsOR, true);

  helib::QueryBuilder res3(height && name && age);
  query = res3.build(columns);
  EXPECT_EQ(query.containsOR, false);
}

TEST(TestQuery, QueryBuilderGeneratesMusAndTaus)
{
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  const helib::QueryExpr& weight = helib::makeQueryExpr(3);

  // ((0 && 1) || 2) = (0 || 2) && (1 || 2)
  helib::QueryBuilder qbExpand0((name && age) || height);
  // (0 || 1 || 2)
  helib::QueryBuilder qbExpand1(name || age || height);
  // (0 || 1) && (2)
  helib::QueryBuilder qbExpand2((name || age) && height);
  //((0 || 1) && (0 || 2))
  helib::QueryBuilder qbExpand3((name || age) && (name || height));
  // ((0 && 1) || (1 && 2)) = (0 || 1) && 1 && (2 || 0) && (2 || 1)
  helib::QueryBuilder qbExpand4((name && age) || (age && height));
  // 0 && (1 || (2 && 3)) = 0 && (1 || 2) && (1 && 3)
  helib::QueryBuilder qbExpand5(name && (age || (height && weight)));
  std::array<helib::QueryBuilder, 6> qbs =
      {qbExpand0, qbExpand1, qbExpand2, qbExpand3, qbExpand4, qbExpand5};
  long columns = 5;
  long cases = 6;
  std::array<std::vector<std::vector<long>>, 6> expected_Fs_vector;
  std::vector<long> F(columns);
  std::iota(F.begin(), F.end(), 0);
  std::vector<std::vector<long>> Fs(4);
  std::fill(Fs.begin(), Fs.end(), F);
  std::fill(expected_Fs_vector.begin(), expected_Fs_vector.end(), Fs);
  expected_Fs_vector[0].resize(2); // query 0 has 2 conjunctions
  expected_Fs_vector[1].resize(1); // query 1 has 1 conjunction
  expected_Fs_vector[2].resize(2); // query 2 has 2 conjunctions
  expected_Fs_vector[3].resize(2); // query 3 has 2 conjunctions
  expected_Fs_vector[4].resize(4); // query 4 has 4 conjunctions
  expected_Fs_vector[5].resize(3); // query 5 has 3 conjunctions

  std::array<std::vector<helib::Matrix<long>>, 6> expected_taus_vector;
  expected_taus_vector[0] = {
      {{1}, {0}, {1}, {0}, {0}},  // Either 0-th or 2nd column
      {{0}, {1}, {1}, {0}, {0}}}; // Either 1st or 2nd column

  expected_taus_vector[1] = {
      {{1}, {1}, {1}, {0}, {0}}}; // Either 0-th or 1st or 2nd column

  expected_taus_vector[2] = {
      {{1}, {1}, {0}, {0}, {0}},  // Either 0-th or 1st column
      {{0}, {0}, {1}, {0}, {0}}}; // Only 2nd column

  expected_taus_vector[3] = {
      {{1}, {1}, {0}, {0}, {0}},  // Either 0-th or 1st column
      {{1}, {0}, {1}, {0}, {0}}}; // Either 0-th or 2nd column

  expected_taus_vector[4] = {
      {{1}, {1}, {0}, {0}, {0}},  // Either 0-th or 1st column
      {{0}, {1}, {0}, {0}, {0}},  // Only 1st column
      {{1}, {0}, {1}, {0}, {0}},  // Either 2nd or 0-th column
      {{0}, {1}, {1}, {0}, {0}}}; // Either 2nd or 1st column

  expected_taus_vector[5] = {
      {{1}, {0}, {0}, {0}, {0}},  // Only 0-th column
      {{0}, {1}, {1}, {0}, {0}},  // Either 1st or 2nd column
      {{0}, {1}, {0}, {1}, {0}}}; // Either 1st or 3rd column
  std::array<std::vector<long>, 6> expected_mus_vector{
      {{0, 0}, {0}, {0, 0}, {0, 0}, {0, 0, 0, 0}, {0, 0, 0}}};
  EXPECT_EQ(expected_Fs_vector.size(), cases);
  EXPECT_EQ(expected_taus_vector.size(), cases);
  EXPECT_EQ(expected_mus_vector.size(), cases);

  for (int j = 0; j < cases; j++) {
    helib::QueryType query = qbs[j].build(columns);
    EXPECT_EQ(expected_Fs_vector[j].size(), query.Fs.size()) << "*** j = " << j;
    EXPECT_EQ(expected_mus_vector[j].size(), query.mus.size())
        << "*** j = " << j;
    EXPECT_EQ(expected_taus_vector[j].size(), query.taus.size())
        << "*** j = " << j;
    EXPECT_EQ(columns, query.taus[0].size()) << "*** j = " << j;
    for (size_t i = 0; i < expected_Fs_vector[j].size(); ++i)
      EXPECT_EQ(expected_Fs_vector[j][i], query.Fs[i])
          << "*** j = " << j << ", i = " << i;
    for (size_t i = 0; i < expected_mus_vector[j].size(); ++i)
      EXPECT_EQ(expected_mus_vector[j][i], query.mus[i])
          << "*** j = " << j << ", i= " << i;
    for (size_t i = 0; i < expected_taus_vector[j].size(); ++i)
      EXPECT_TRUE(expected_taus_vector[j][i] == query.taus[i])
          << "*** j = " << j << ", i= " << i;
  }
}

TEST(TestQuery, QueryBuilderGeneratesMusAndTausWithNot)
{
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  const helib::QueryExpr& weight = helib::makeQueryExpr(3);

  // ! (0 || 1 || 2) -> {{-3},{-2},{-1}}
  helib::QueryBuilder qbNotOfOr1(!(name || age || height));
  // !(0 || 1) && 2 -> {{-2},{-1},{3}}
  helib::QueryBuilder qbNotOfOr2(!(name || age) && height);
  // ! ! 0 -> {{1}}
  helib::QueryBuilder qbDoubleNot1(!!name);
  // 0 || ! ! 3 ->{{4,2}}
  helib::QueryBuilder qbDoubleNot2(age || !!weight);
  // ! (0 && 1 && 2) -> {{-1,-2,-3}}
  helib::QueryBuilder qbNotOfAnd1(!(name && age && height));
  // !((0 || 1) && (2 || 3)) -> {{-2,-4},{-2,-3},{-1,-4},{-1,-3}}
  helib::QueryBuilder qbNotOfAnd2(!((name || age) && (height || weight)));
  std::array<helib::QueryBuilder, 6> qbs = {qbNotOfOr1,
                                            qbNotOfOr2,
                                            qbDoubleNot1,
                                            qbDoubleNot2,
                                            qbNotOfAnd1,
                                            qbNotOfAnd2};
  long columns = 5;
  long cases = 6;
  long max_clauses = 4;
  std::array<std::vector<std::vector<long>>, 6> expected_Fs_vector;
  std::vector<long> F(columns);
  std::iota(F.begin(), F.end(), 0);
  std::vector<std::vector<long>> Fs(max_clauses);
  std::fill(Fs.begin(), Fs.end(), F);
  std::fill(expected_Fs_vector.begin(), expected_Fs_vector.end(), Fs);
  expected_Fs_vector[0].resize(3); // query 0 has 3 conjunctions
  expected_Fs_vector[1].resize(3); // query 1 has 3 conjunctions
  expected_Fs_vector[2].resize(1); // query 2 has 1 conjunction
  expected_Fs_vector[3].resize(1); // query 3 has 1 conjunction
  expected_Fs_vector[4].resize(1); // query 4 has 1 conjunction
  expected_Fs_vector[5].resize(4); // query 5 has 4 conjunctions

  std::array<std::vector<helib::Matrix<long>>, 6> expected_taus_vector;
  expected_taus_vector[0] = {{{0}, {0}, {-1}, {0}, {0}},  // Not 2nd
                             {{0}, {-1}, {0}, {0}, {0}},  // Not 1st
                             {{-1}, {0}, {0}, {0}, {0}}}; // Not 0th

  expected_taus_vector[1] = {{{0}, {-1}, {0}, {0}, {0}}, // Not 1st
                             {{-1}, {0}, {0}, {0}, {0}}, // Not 0th
                             {{0}, {0}, {1}, {0}, {0}}}; // 2nd

  expected_taus_vector[2] = {{{1}, {0}, {0}, {0}, {0}}}, // 0th

      expected_taus_vector[3] = {{{0}, {1}, {0}, {1}, {0}}}; // 3rd or 1st

  expected_taus_vector[4] = {
      {{-1}, {-1}, {-1}, {0}, {0}}}; // Not 0th or not 1st or not 2nd

  expected_taus_vector[5] = {{{0}, {-1}, {0}, {-1}, {0}},  // Not 1st, not 3rd
                             {{0}, {-1}, {-1}, {0}, {0}},  // Not 1st, not 2nd
                             {{-1}, {0}, {0}, {-1}, {0}},  // Not 0th, not 3rd
                             {{-1}, {0}, {-1}, {0}, {0}}}; // Not 0th, not 2nd
  std::array<std::vector<long>, 6> expected_mus_vector{
      {{1, 1, 1}, {1, 1, 0}, {0}, {0}, {3}, {2, 2, 2, 2}}};
  EXPECT_EQ(expected_Fs_vector.size(), cases);
  EXPECT_EQ(expected_taus_vector.size(), cases);
  EXPECT_EQ(expected_mus_vector.size(), cases);

  for (int j = 0; j < cases; j++) {
    helib::QueryType query = qbs[j].build(columns);
    EXPECT_EQ(expected_Fs_vector[j].size(), query.Fs.size()) << "*** j = " << j;
    EXPECT_EQ(expected_mus_vector[j].size(), query.mus.size())
        << "*** j = " << j;
    EXPECT_EQ(expected_taus_vector[j].size(), query.taus.size())
        << "*** j = " << j;
    EXPECT_EQ(columns, query.taus[0].size()) << "*** j = " << j;
    for (size_t i = 0; i < expected_Fs_vector[j].size(); ++i)
      EXPECT_EQ(expected_Fs_vector[j][i], query.Fs[i])
          << "*** j = " << j << ", i = " << i;
    for (size_t i = 0; i < expected_mus_vector[j].size(); ++i)
      EXPECT_EQ(expected_mus_vector[j][i], query.mus[i])
          << "*** j = " << j << ", i= " << i;
    for (size_t i = 0; i < expected_taus_vector[j].size(); ++i)
      EXPECT_TRUE(expected_taus_vector[j][i] == query.taus[i])
          << "*** j = " << j << ", i= " << i;
  }
}

TEST(TestQuery, RemoveOrExpectedString)
{
  const helib::QueryExpr& a = helib::makeQueryExpr(0);
  const helib::QueryExpr& b = helib::makeQueryExpr(1);
  const helib::QueryExpr& c = helib::makeQueryExpr(2);
  const helib::QueryExpr& d = helib::makeQueryExpr(3);

  std::array<helib::QueryBuilder, 13> qbs = {
      helib::QueryBuilder((a && b) || c),
      helib::QueryBuilder(a || b || c),
      helib::QueryBuilder((a || b) && c),
      helib::QueryBuilder((a || b) && (c || d)),
      helib::QueryBuilder((a && b) || (c && d)),
      helib::QueryBuilder(a && (b || (c && d))),
      helib::QueryBuilder(a || !a),
      helib::QueryBuilder(!(a || b || c)),
      helib::QueryBuilder(!(a || b) && c),
      helib::QueryBuilder(!!a),
      helib::QueryBuilder(b || !!a),
      helib::QueryBuilder(!(a && b && c)),
      helib::QueryBuilder(!((a || b) && (c || d)))};
  std::array<std::string, 13> expected_strings = {
      "0 1 && ! 2 ! && !",
      "0 ! 1 ! && ! ! 2 ! && !",
      "0 ! 1 ! && ! 2 &&",
      "0 ! 1 ! && ! 2 ! 3 ! && ! &&",
      "0 1 && ! 2 3 && ! && !",
      "0 1 ! 2 3 && ! && ! &&",
      "0 ! 0 ! ! && !",
      "0 ! 1 ! && ! ! 2 ! && ! !",
      "0 ! 1 ! && ! ! 2 &&",
      "0 ! !",
      "1 ! 0 ! ! ! && !",
      "0 1 && 2 && !",
      "0 ! 1 ! && ! 2 ! 3 ! && ! && !"};

  for (std::size_t i = 0; i < qbs.size(); i++) {
    qbs[i].removeOr();
    EXPECT_EQ(expected_strings[i], qbs[i].getQueryString()) << "*** i = " << i;
  }
}
} // namespace
