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

#include <algorithm>
#include <functional>
#include <helib/Matrix.h>
#include <helib/helib.h>

#include "test_common.h"
#include "gtest/gtest.h"

// Helper classes
struct NoDefaultConstructor
{
  int data = 0;
  NoDefaultConstructor() = delete;
  NoDefaultConstructor(int i) : data(i) {}
};

struct DefaultConstructedThrows
{
  DefaultConstructedThrows() : defaulted(true) {}
  DefaultConstructedThrows(int i) : data(i), defaulted(false) {}
  int data;
  bool defaulted;
};

DefaultConstructedThrows operator+(const DefaultConstructedThrows& a,
                                   const DefaultConstructedThrows& b)
{
  if (a.defaulted || b.defaulted)
    throw helib::RuntimeError("Thrown due to default construction");
  return DefaultConstructedThrows(a.data + b.data);
}

DefaultConstructedThrows operator*(const DefaultConstructedThrows& a,
                                   const DefaultConstructedThrows& b)
{
  if (a.defaulted || b.defaulted)
    throw helib::RuntimeError("Thrown due to default construction");
  return DefaultConstructedThrows(a.data * b.data);
}

DefaultConstructedThrows operator-(const DefaultConstructedThrows& a,
                                   const DefaultConstructedThrows& b)
{
  if (a.defaulted || b.defaulted)
    throw helib::RuntimeError("Thrown due to default construction");
  return DefaultConstructedThrows(a.data - b.data);
}

DefaultConstructedThrows& operator+=(DefaultConstructedThrows& a,
                                     const DefaultConstructedThrows& b)
{
  if (a.defaulted || b.defaulted)
    throw helib::RuntimeError("Thrown due to default construction");
  a.data += b.data;
  return a;
}

DefaultConstructedThrows& operator-=(DefaultConstructedThrows& a,
                                     const DefaultConstructedThrows& b)
{
  if (a.defaulted || b.defaulted)
    throw helib::RuntimeError("Thrown due to default construction");
  a.data -= b.data;
  return a;
}

DefaultConstructedThrows& operator*=(DefaultConstructedThrows& a,
                                     const DefaultConstructedThrows& b)
{
  if (a.defaulted || b.defaulted)
    throw helib::RuntimeError("Thrown due to default construction");
  a.data *= b.data;
  return a;
}

bool operator==(const DefaultConstructedThrows& a,
                const DefaultConstructedThrows& b)
{
  return (a.defaulted && b.defaulted) ||
         (!a.defaulted && !b.defaulted && a.data == b.data);
}

struct TypeB
{
  int data = 0;
};

// TypeA has methods that take TypeB
struct TypeA
{
  int data = 0;

  TypeA() = default;
  TypeA(int i) : data(i) {}

  TypeA& operator+=(const TypeA& other)
  {
    this->data += other.data;
    return *this;
  }

  TypeA& operator-=(const TypeA& other)
  {
    this->data -= other.data;
    return *this;
  }

  TypeA& operator+=(const TypeB& b)
  {
    this->data += b.data;
    return *this;
  }

  TypeA& operator-=(const TypeB& b)
  {
    this->data -= b.data;
    return *this;
  }

  TypeA& operator*=(const TypeB& b)
  {
    this->data *= b.data;
    return *this;
  }
};

template <typename T>
void printVector(std::vector<T> v)
{
  for (const auto& e : v) {
    std::cout << e << ' ';
  }
  std::cout << std::endl;
}

template <typename T>
void printMatrix(helib::Matrix<T> M)
{
  for (long i = 0; i < M.dims(0); ++i) {
    for (long j = 0; j < M.dims(1); ++j) {
      std::cout << M(i, j) << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

struct BGVParameters
{
  BGVParameters(unsigned m, unsigned p, unsigned r) : m(m), p(p), r(r){};

  const unsigned m;
  const unsigned p;
  const unsigned r;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << "}";
  }
};

class TestMatrixWithCtxt : public ::testing::TestWithParam<BGVParameters>
{

protected:
  const unsigned long m;
  const unsigned long p;
  const unsigned long r;

  helib::Context context;
  helib::SecKey sk;
  const helib::PubKey& pk;
  const helib::EncryptedArray& ea;

  TestMatrixWithCtxt() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      context(m, p, r),
      sk((helib::buildModChain(context, /*bits*/ 300, /*c*/ 2), context)),
      pk((sk.GenSecKey(), sk)),
      ea(*context.ea)
  {}
};

namespace {

TEST(TestMatrix, ConstructMatrix)
{
  helib::Matrix<int> M(2, 3);

  EXPECT_EQ(M.size(), 2 * 3);
  EXPECT_EQ(M.dims(0), 2);
  EXPECT_EQ(M.dims(1), 3);
  EXPECT_EQ(M.order(), 2);
}

TEST(TestMatrix, PopulateMatrix)
{
  helib::Matrix<int> M(4, 2);

  M(0, 0) = 1;
  M(0, 1) = 2;
  M(1, 0) = 3;
  M(1, 1) = 4;
  M(2, 0) = 5;
  M(2, 1) = 6;
  M(3, 0) = 7;
  M(3, 1) = 8;

  EXPECT_EQ(M(0, 0), 1);
  EXPECT_EQ(M(0, 1), 2);
  EXPECT_EQ(M(1, 0), 3);
  EXPECT_EQ(M(1, 1), 4);
  EXPECT_EQ(M(2, 0), 5);
  EXPECT_EQ(M(2, 1), 6);
  EXPECT_EQ(M(3, 0), 7);
  EXPECT_EQ(M(3, 1), 8);
}

TEST(TestMatrix, PopulateMatrixWithInitializerList)
{
  helib::Matrix<int> M = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};

  EXPECT_EQ(M(0, 0), 1);
  EXPECT_EQ(M(0, 1), 2);
  EXPECT_EQ(M(1, 0), 3);
  EXPECT_EQ(M(1, 1), 4);
  EXPECT_EQ(M(2, 0), 5);
  EXPECT_EQ(M(2, 1), 6);
  EXPECT_EQ(M(3, 0), 7);
  EXPECT_EQ(M(3, 1), 8);

  EXPECT_EQ(M.size(), 4 * 2);
  EXPECT_EQ(M.dims(0), 4);
  EXPECT_EQ(M.dims(1), 2);
  EXPECT_EQ(M.order(), 2);
}

TEST(TestMatrix, PopulateMatrixWithInitializerListColumnMismatch)
{
  EXPECT_THROW((helib::Matrix<int>{{1, 2}, {3, 4}, {5, 6, 6}, {7, 8}}),
               helib::LogicError);
}

TEST(TestMatrix, PopulateMatrixWithNonDefaultConstructor)
{
  helib::Matrix<NoDefaultConstructor> M(NoDefaultConstructor(0), 4, 2);

  M(0, 0) = NoDefaultConstructor(1);
  M(0, 1) = NoDefaultConstructor(2);
  M(1, 0) = NoDefaultConstructor(3);
  M(1, 1) = NoDefaultConstructor(4);
  M(2, 0) = NoDefaultConstructor(5);
  M(2, 1) = NoDefaultConstructor(6);
  M(3, 0) = NoDefaultConstructor(7);
  M(3, 1) = NoDefaultConstructor(8);

  EXPECT_EQ(M(0, 0).data, 1);
  EXPECT_EQ(M(0, 1).data, 2);
  EXPECT_EQ(M(1, 0).data, 3);
  EXPECT_EQ(M(1, 1).data, 4);
  EXPECT_EQ(M(2, 0).data, 5);
  EXPECT_EQ(M(2, 1).data, 6);
  EXPECT_EQ(M(3, 0).data, 7);
  EXPECT_EQ(M(3, 1).data, 8);
}

TEST(TestMatrix, EqualsOpertorForMatrices)
{
  helib::Matrix<int> M1 = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
  helib::Matrix<int> M2 = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
  helib::Matrix<int> M3 = {{0, 1}, {2, 3}, {4, 5}, {6, 7}};
  helib::Matrix<int> M4 = {{0, 1}, {2, 3}, {4, 5}};
  helib::Matrix<int> M5 = M1;
  // FIXME: Transpose of M5 = M1
  // M5.transpose();

  EXPECT_TRUE(M1 == M1);
  EXPECT_TRUE(M1 == M2);
  EXPECT_FALSE(M1 == M3);
  EXPECT_FALSE(M1 == M4);
  // EXPECT_FALSE(M1 == M5);
}

TEST(TestMatrix, MoveMatrix)
{
  helib::Matrix<int> M1 = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
  helib::Matrix<int> M2 = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
  helib::Matrix<int> M3 = std::move(M2);

  EXPECT_TRUE(M1 == M3);
  EXPECT_FALSE(M1 == M2);
}

TEST(TestMatrix, NotEqualsOpertorForMatrices)
{
  helib::Matrix<int> M1 = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
  helib::Matrix<int> M2 = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
  helib::Matrix<int> M3 = {{0, 1}, {2, 3}, {4, 5}, {6, 7}};
  helib::Matrix<int> M4 = {{0, 1}, {2, 3}, {4, 5}};
  helib::Matrix<int> M5 = M1;
  // FIXME: Transpose of M5 = M1
  // M5.transpose();

  EXPECT_FALSE(M1 != M1);
  EXPECT_FALSE(M1 != M2);
  EXPECT_TRUE(M1 != M3);
  EXPECT_TRUE(M1 != M4);
  // EXPECT_TRUE(M1 != M5);
}

TEST(TestMatrix, AccessingMatrixElementOutOfBound)
{
  helib::Matrix<int> M(2, 2);

  EXPECT_EQ(M.size(), 2 * 2);
  EXPECT_EQ(M.dims(0), 2);
  EXPECT_EQ(M.dims(1), 2);
  EXPECT_EQ(M.order(), 2);

  EXPECT_THROW(M(0, 2), helib::OutOfRangeError);
  EXPECT_THROW(M(2, 0), helib::OutOfRangeError);
  EXPECT_THROW(M(2, 2), helib::OutOfRangeError);
}

TEST(TestMatrix, GetMatrixRow)
{
  helib::Matrix<int> M(2, 3);
  M(0, 0) = 1;
  M(0, 1) = 2;
  M(0, 2) = 3;
  M(1, 0) = 4;
  M(1, 1) = 5;
  M(1, 2) = 6;

  auto subM1 = M.row(0);
  auto subM2 = M.row(1);

  EXPECT_EQ(subM1.size(), 3);
  EXPECT_EQ(subM1.dims(0), 3);
  EXPECT_EQ(subM1.order(), 1);

  EXPECT_EQ(subM1(0), 1);
  EXPECT_EQ(subM1(1), 2);
  EXPECT_EQ(subM1(2), 3);

  EXPECT_EQ(subM2.size(), 3);
  EXPECT_EQ(subM2.dims(0), 3);
  EXPECT_EQ(subM2.order(), 1);

  EXPECT_EQ(subM2(0), 4);
  EXPECT_EQ(subM2(1), 5);
  EXPECT_EQ(subM2(2), 6);
}

TEST(TestMatrix, GetMatrixColumn)
{
  helib::Matrix<int> M(3, 2);
  M(0, 0) = 1;
  M(0, 1) = 2;
  M(1, 0) = 3;
  M(1, 1) = 4;
  M(2, 0) = 5;
  M(2, 1) = 6;

  helib::Vector<int> subM1 = M.column(0);
  helib::Vector<int> subM2 = M.column(1);

  EXPECT_EQ(subM1.size(), 3);
  EXPECT_EQ(subM1.dims(0), 3);
  EXPECT_EQ(subM1.order(), 1);

  EXPECT_EQ(subM1(0), 1);
  EXPECT_EQ(subM1(1), 3);
  EXPECT_EQ(subM1(2), 5);

  EXPECT_EQ(subM2.size(), 3);
  EXPECT_EQ(subM2.dims(0), 3);
  EXPECT_EQ(subM2.order(), 1);

  EXPECT_EQ(subM2(0), 2);
  EXPECT_EQ(subM2(1), 4);
  EXPECT_EQ(subM2(2), 6);
}

TEST(TestMatrix, AddMatrices)
{
  helib::Matrix<int> M1(2, 3);
  M1(0, 0) = 1;
  M1(0, 1) = 2;
  M1(0, 2) = 3;
  M1(1, 0) = 4;
  M1(1, 1) = 5;
  M1(1, 2) = 6;

  helib::Matrix<int> M2(2, 3);
  M2(0, 0) = 1;
  M2(0, 1) = -3;
  M2(0, 2) = 7;
  M2(1, 0) = 1;
  M2(1, 1) = 7;
  M2(1, 2) = 5;

  M1 += M2;

  EXPECT_EQ(M1(0, 0), 2);
  EXPECT_EQ(M1(0, 1), -1);
  EXPECT_EQ(M1(0, 2), 10);
  EXPECT_EQ(M1(1, 0), 5);
  EXPECT_EQ(M1(1, 1), 12);
  EXPECT_EQ(M1(1, 2), 11);
}

TEST(TestMatrix, SubtractMatrices)
{
  helib::Matrix<int> M1(2, 3);
  M1(0, 0) = 1;
  M1(0, 1) = 2;
  M1(0, 2) = 3;
  M1(1, 0) = 4;
  M1(1, 1) = 5;
  M1(1, 2) = 6;

  helib::Matrix<int> M2(2, 3);
  M2(0, 0) = 1;
  M2(0, 1) = -3;
  M2(0, 2) = 7;
  M2(1, 0) = 1;
  M2(1, 1) = 7;
  M2(1, 2) = 5;

  M1 -= M2;

  EXPECT_EQ(M1(0, 0), 0);
  EXPECT_EQ(M1(0, 1), 5);
  EXPECT_EQ(M1(0, 2), -4);
  EXPECT_EQ(M1(1, 0), 3);
  EXPECT_EQ(M1(1, 1), -2);
  EXPECT_EQ(M1(1, 2), 1);
}

TEST(TestMatrix, hadamardWorks)
{
  helib::Matrix<int> M1(2, 3);
  M1(0, 0) = 1;
  M1(0, 1) = 2;
  M1(0, 2) = 3;
  M1(1, 0) = 4;
  M1(1, 1) = 5;
  M1(1, 2) = 6;

  helib::Matrix<int> M2(2, 3);
  M2(0, 0) = 1;
  M2(0, 1) = -3;
  M2(0, 2) = 7;
  M2(1, 0) = 1;
  M2(1, 1) = 7;
  M2(1, 2) = 5;

  M1.hadamard(M2);

  EXPECT_EQ(M1(0, 0), 1);
  EXPECT_EQ(M1(0, 1), -6);
  EXPECT_EQ(M1(0, 2), 21);
  EXPECT_EQ(M1(1, 0), 4);
  EXPECT_EQ(M1(1, 1), 35);
  EXPECT_EQ(M1(1, 2), 30);
}

TEST(TestMatrix, entrywiseOperationWorks)
{
  helib::Matrix<int> M1(2, 3);
  M1(0, 0) = 1;
  M1(0, 1) = 2;
  M1(0, 2) = 3;
  M1(1, 0) = 4;
  M1(1, 1) = 5;
  M1(1, 2) = 6;

  helib::Matrix<int> M2(2, 3);
  M2(0, 0) = 1;
  M2(0, 1) = 3;
  M2(0, 2) = 7;
  M2(1, 0) = 1;
  M2(1, 1) = 7;
  M2(1, 2) = 5;

  std::function<int&(int&, const int&)> f =
      [](auto& left, const auto& right) -> decltype(auto) {
    return left %= right;
  };

  M1.entrywiseOperation(M2, f);

  EXPECT_EQ(M1(0, 0), 0);
  EXPECT_EQ(M1(0, 1), 2);
  EXPECT_EQ(M1(0, 2), 3);
  EXPECT_EQ(M1(1, 0), 0);
  EXPECT_EQ(M1(1, 1), 5);
  EXPECT_EQ(M1(1, 2), 1);
}

TEST(TestMatrix, ExceptionWhenMatrixDimensionsDoNotMatchDuringSubtraction)
{
  helib::Matrix<int> M1(2, 3);
  helib::Matrix<int> M2(2, 2);

  EXPECT_THROW(M1 -= M2, helib::LogicError);
}

TEST(TestMatrix, ApplyTransformToElementsWithoutFullView)
{
  std::cout << "Available Threads: " << NTL::AvailableThreads() << std::endl;

  helib::Matrix<int> M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  auto view = M.columns({0, 2});
  view.apply([](int& x) { x += 5; });

  EXPECT_EQ(view(0, 0), 6);
  EXPECT_EQ(view(0, 1), 8);
  EXPECT_EQ(view(1, 0), 9);
  EXPECT_EQ(view(1, 1), 11);
  EXPECT_EQ(view(2, 0), 12);
  EXPECT_EQ(view(2, 1), 14);
}

TEST(TestMatrix, ApplyTransformToElementsWithFullView)
{
  std::cout << "Available Threads: " << NTL::AvailableThreads() << std::endl;

  helib::Matrix<int> M(2, 3);
  M(0, 0) = 1;
  M(0, 1) = 2;
  M(0, 2) = 3;
  M(1, 0) = 4;
  M(1, 1) = 5;
  M(1, 2) = 6;

  M.apply([](int& x) { x += 5; });

  EXPECT_EQ(M(0, 0), 6);
  EXPECT_EQ(M(0, 1), 7);
  EXPECT_EQ(M(0, 2), 8);
  EXPECT_EQ(M(1, 0), 9);
  EXPECT_EQ(M(1, 1), 10);
  EXPECT_EQ(M(1, 2), 11);
}

TEST(TestMatrix, TransposeMatrix)
{

  helib::Matrix<int> M(2, 3);
  M(0, 0) = 1;
  M(0, 1) = 2;
  M(0, 2) = 3;
  M(1, 0) = 4;
  M(1, 1) = 5;
  M(1, 2) = 6;

  M.transpose();

  EXPECT_EQ(M(0, 0), 1);
  EXPECT_EQ(M(1, 0), 2);
  EXPECT_EQ(M(2, 0), 3);
  EXPECT_EQ(M(0, 1), 4);
  EXPECT_EQ(M(1, 1), 5);
  EXPECT_EQ(M(2, 1), 6);
}

TEST(TestMatrix, MultiplyTwoMatrices)
{
  helib::Matrix<int> M1(2, 3);
  M1(0, 0) = 1;
  M1(0, 1) = 2;
  M1(0, 2) = 3;
  M1(1, 0) = 4;
  M1(1, 1) = 5;
  M1(1, 2) = 6;

  helib::Matrix<int> M2(2, 3);
  M2(0, 0) = 1;
  M2(0, 1) = -3;
  M2(0, 2) = 7;
  M2(1, 0) = 1;
  M2(1, 1) = 7;
  M2(1, 2) = 5;

  helib::Matrix<int> M3 = M1 * M2.transpose();

  // Check dimensions (expect 2x2 Matrix)
  EXPECT_EQ(M3.size(), 2 * 2);
  EXPECT_EQ(M3.dims(0), 2);
  EXPECT_EQ(M3.dims(1), 2);
  EXPECT_EQ(M3.order(), 2);
  // Check result
  EXPECT_EQ(M3(0, 0), 16);
  EXPECT_EQ(M3(0, 1), 30);
  EXPECT_EQ(M3(1, 0), 31);
  EXPECT_EQ(M3(1, 1), 69);
}

TEST(TestMatrix, MultiplyTwoMatricesWhoseTypeThrowsIfDefaultConstructed)
{
  helib::Matrix<DefaultConstructedThrows> M1 = {{5, 2}, {-1, 3}, {-4, 0}};
  helib::Matrix<DefaultConstructedThrows> M2 = {{-1, 8}, {0, -2}};
  helib::Matrix<DefaultConstructedThrows> M3 = M1 * M2;
  helib::Matrix<DefaultConstructedThrows> expected = {{-5, 36},
                                                      {1, -14},
                                                      {4, -32}};

  EXPECT_EQ(M3.size(), expected.size());
  EXPECT_EQ(M3.dims(0), expected.dims(0));
  EXPECT_EQ(M3.dims(1), expected.dims(1));
  EXPECT_EQ(M3.order(), expected.order());

  EXPECT_EQ(M3(0, 0), expected(0, 0));
  EXPECT_EQ(M3(0, 1), expected(0, 1));
  EXPECT_EQ(M3(1, 0), expected(1, 0));
  EXPECT_EQ(M3(1, 1), expected(1, 1));
  EXPECT_EQ(M3(2, 0), expected(2, 0));
  EXPECT_EQ(M3(2, 1), expected(2, 1));
}

TEST(TestMatrix, ExceptionWhenMatrixDimensionsDoNotMatchDuringMultiplication)
{
  helib::Matrix<int> M1(2, 3);
  helib::Matrix<int> M2(2, 2);

  EXPECT_THROW(M1 * M2, helib::LogicError);
}

TEST(TestMatrix, ValidSubtractTwoMatricesOfDifferentTypes)
{
  // TypeA can be subtracted by TypeB
  helib::Matrix<TypeA> A(2, 3);
  A(0, 0).data = 1;
  A(0, 1).data = 2;
  A(0, 2).data = 3;
  A(1, 0).data = 4;
  A(1, 1).data = 5;
  A(1, 2).data = 6;

  helib::Matrix<TypeB> B(2, 3);
  B(0, 0).data = 6;
  B(0, 1).data = 5;
  B(0, 2).data = 4;
  B(1, 0).data = 3;
  B(1, 1).data = 2;
  B(1, 2).data = 1;

  A -= B;

  EXPECT_EQ(A(0, 0).data, -5);
  EXPECT_EQ(A(0, 1).data, -3);
  EXPECT_EQ(A(0, 2).data, -1);
  EXPECT_EQ(A(1, 0).data, 1);
  EXPECT_EQ(A(1, 1).data, 3);
  EXPECT_EQ(A(1, 2).data, 5);
}

TEST(TestMatrix, ValidMultiplyTwoMatricesOfDifferentTypes)
{
  // TypeA can be multiplied by TypeB
  helib::Matrix<TypeA> A(2, 3);
  A(0, 0).data = 1;
  A(0, 1).data = 2;
  A(0, 2).data = 3;
  A(1, 0).data = 4;
  A(1, 1).data = 5;
  A(1, 2).data = 6;

  helib::Matrix<TypeB> B(2, 3);
  B(0, 0).data = 1;
  B(0, 1).data = -3;
  B(0, 2).data = 7;
  B(1, 0).data = 1;
  B(1, 1).data = 7;
  B(1, 2).data = 5;

  helib::Matrix<TypeA> C = A * B.transpose();

  // Check dimensions (expect 2x2 Matrix)
  EXPECT_EQ(C.size(), 2 * 2);
  EXPECT_EQ(C.dims(0), 2);
  EXPECT_EQ(C.dims(1), 2);
  EXPECT_EQ(C.order(), 2);
  // Check result
  EXPECT_EQ(C(0, 0).data, 16);
  EXPECT_EQ(C(0, 1).data, 30);
  EXPECT_EQ(C(1, 0).data, 31);
  EXPECT_EQ(C(1, 1).data, 69);
}

TEST(TestMatrix, ConstructMatrixView)
{
  helib::Matrix<int> M = {{9, 8, 7, -9, -8, -7, 10},
                          {6, 5, 4, -6, -5, -4, 20},
                          {3, 2, 1, -3, -2, -1, 30}};

  helib::Matrix<int> view = M.columns({1, 2, 3, 6, 0, 2});

  // Check dimensions
  EXPECT_EQ(view.size(), 3 * 6);
  EXPECT_EQ(view.dims(0), 3);
  EXPECT_EQ(view.dims(1), 6);
  EXPECT_EQ(view.order(), 2);
  EXPECT_FALSE(view.fullView());

  // Check result
  // row 1
  EXPECT_EQ(view(0, 0), 8);  // 1
  EXPECT_EQ(view(0, 1), 7);  // 2
  EXPECT_EQ(view(0, 2), -9); // 3
  EXPECT_EQ(view(0, 3), 10); // 6
  EXPECT_EQ(view(0, 4), 9);  // 0
  EXPECT_EQ(view(0, 5), 7);  // 2
  // row 2
  EXPECT_EQ(view(1, 0), 5);  // 1
  EXPECT_EQ(view(1, 1), 4);  // 2
  EXPECT_EQ(view(1, 2), -6); // 3
  EXPECT_EQ(view(1, 3), 20); // 6
  EXPECT_EQ(view(1, 4), 6);  // 0
  EXPECT_EQ(view(1, 5), 4);  // 2
  // row 3
  EXPECT_EQ(view(2, 0), 2);  // 1
  EXPECT_EQ(view(2, 1), 1);  // 2
  EXPECT_EQ(view(2, 2), -3); // 3
  EXPECT_EQ(view(2, 3), 30); // 6
  EXPECT_EQ(view(2, 4), 3);  // 0
  EXPECT_EQ(view(2, 5), 1);  // 2
}

TEST(TestMatrix, TransposeMatrixView)
{
  helib::Matrix<int> M = {{10, 0, 3, 20, 6, 9, 30},
                          {11, 1, 4, 21, 7, 10, 31},
                          {12, 2, 5, 22, 8, 11, 32}};

  helib::Matrix<int> view = M.columns({1, 2, 4, 5});

  view.transpose();

  EXPECT_EQ(view.size(), 4 * 3);
  EXPECT_EQ(view.dims(0), 4);
  EXPECT_EQ(view.dims(1), 3);
  EXPECT_EQ(view.order(), 2);
  // When transpose is called, view becomes a fullView and copies in elements
  EXPECT_TRUE(view.fullView());

  // row1
  EXPECT_EQ(view(0, 0), 0);
  EXPECT_EQ(view(0, 1), 1);
  EXPECT_EQ(view(0, 2), 2);
  // row2
  EXPECT_EQ(view(1, 0), 3);
  EXPECT_EQ(view(1, 1), 4);
  EXPECT_EQ(view(1, 2), 5);
  // row3
  EXPECT_EQ(view(2, 0), 6);
  EXPECT_EQ(view(2, 1), 7);
  EXPECT_EQ(view(2, 2), 8);
  // row4
  EXPECT_EQ(view(3, 0), 9);
  EXPECT_EQ(view(3, 1), 10);
  EXPECT_EQ(view(3, 2), 11);
}

TEST(TestMatrix, SubtractMatrixViewFromMatrix)
{
  helib::Matrix<int> S = {{9, 8, 7, -9, -8, -7, 10},
                          {6, 5, 4, -6, -5, -4, 20},
                          {3, 2, 1, -3, -2, -1, 30}};

  helib::Matrix<int> view = S.columns({2, 2, 5, 6});

  helib::Matrix<int> M = {{10, 20, 30, 40}, {5, 6, 7, 8}, {50, 60, 70, 80}};

  M -= view;

  // row1
  EXPECT_EQ(M(0, 0), 3);
  EXPECT_EQ(M(0, 1), 13);
  EXPECT_EQ(M(0, 2), 37);
  EXPECT_EQ(M(0, 3), 30);
  // row2
  EXPECT_EQ(M(1, 0), 1);
  EXPECT_EQ(M(1, 1), 2);
  EXPECT_EQ(M(1, 2), 11);
  EXPECT_EQ(M(1, 3), -12);
  // row3
  EXPECT_EQ(M(2, 0), 49);
  EXPECT_EQ(M(2, 1), 59);
  EXPECT_EQ(M(2, 2), 71);
  EXPECT_EQ(M(2, 3), 50);
}

TEST(TestMatrix, MatrixSubtractionOfTwoViewsFromSameData)
{
  helib::Matrix<int> M(3, 7);

  helib::Matrix<int> view1 = M.columns({2, 2, 5, 6});
  helib::Matrix<int> view2 = M.columns({1, 2, 3, 5});

  // TODO Decide on policy here.
  // But for now, not allowed.
  EXPECT_THROW(view1 -= view2, helib::LogicError);
}

TEST(TestMatrix, MatrixMultiplicationOfTwoViewsFromSameData)
{
  helib::Matrix<int> M = {{9, 8, 7, -9, -8, -7, 10},
                          {6, 5, 4, -6, -5, -4, 20},
                          {3, 2, 1, -3, -2, -1, 30}};

  helib::Matrix<int> view1 = M.columns({2, 2, 5, 6});
  helib::Matrix<int> view2 = M.columns({1, 2, 3, 5});

  helib::Matrix<int> R = view1 * view2.transpose();

  // Check dimensions
  EXPECT_EQ(view1.size(), 3 * 4);
  EXPECT_EQ(view1.dims(0), 3);
  EXPECT_EQ(view1.dims(1), 4);
  EXPECT_EQ(view1.order(), 2);
  EXPECT_FALSE(view1.fullView());

  EXPECT_EQ(view2.size(), 4 * 3);
  EXPECT_EQ(view2.dims(0), 4);
  EXPECT_EQ(view2.dims(1), 3);
  EXPECT_EQ(view2.order(), 2);
  // transposed => data copied
  EXPECT_TRUE(view2.fullView());

  // Check dimensions
  EXPECT_EQ(R.size(), 3 * 3);
  EXPECT_EQ(R.dims(0), 3);
  EXPECT_EQ(R.dims(1), 3);
  EXPECT_EQ(R.order(), 2);
  EXPECT_TRUE(R.fullView());

  // row1
  EXPECT_EQ(R(0, 0), 98);
  EXPECT_EQ(R(0, 1), 65);
  EXPECT_EQ(R(0, 2), 32);
  // row2
  EXPECT_EQ(R(1, 0), -44);
  EXPECT_EQ(R(1, 1), -20);
  EXPECT_EQ(R(1, 2), 4);
  // row3
  EXPECT_EQ(R(2, 0), -186);
  EXPECT_EQ(R(2, 1), -105);
  EXPECT_EQ(R(2, 2), -24);
}

TEST(TestMatrix, TransposeIssue)
{
  helib::Matrix<int> M = {{1, 2, 3, 4, 5}};
  M.transpose();
  helib::Matrix<int> view = M.columns({0, 0, 0});
  // Make sure that the 'transposed' property of a matrix is preserved even
  // after taking a view
  EXPECT_EQ(view.size(), 5 * 3);
  EXPECT_EQ(view.dims(0), 5);
  EXPECT_EQ(view.dims(1), 3);
  EXPECT_EQ(view.order(), 2);
  EXPECT_FALSE(view.fullView());
  for (size_t i = 0; i < 5; ++i)
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(view(i, j), i + 1);
    }
}

TEST(TestMatrix, ColumnsWorksWithUnmodifiedPositions)
{
  helib::Matrix<int> M = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};
  helib::Matrix<int> view = M.columns({0, 1, 2});
  EXPECT_EQ(view.size(), 9);
  EXPECT_EQ(view.dims(0), 3);
  EXPECT_EQ(view.dims(1), 3);
  EXPECT_FALSE(view.fullView());
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(view(i, j), M(i, j));
    }
}

TEST(TestMatrix, AccessingMatrixViewElementOutOfBound)
{
  helib::Matrix<int> S(2, 7);
  helib::Matrix<int> M = S.columns({5, 6});

  EXPECT_EQ(M.size(), 2 * 2);
  EXPECT_EQ(M.dims(0), 2);
  EXPECT_EQ(M.dims(1), 2);
  EXPECT_EQ(M.order(), 2);

  EXPECT_THROW(M(0, 2), helib::OutOfRangeError);
  EXPECT_THROW(M(2, 0), helib::OutOfRangeError);
  EXPECT_THROW(M(2, 2), helib::OutOfRangeError);
}

TEST_P(TestMatrixWithCtxt, ConstructMatrixWithCtxt)
{
  helib::Matrix<helib::Ctxt> M(helib::Ctxt(pk), 2, 3);

  EXPECT_EQ(M.size(), 2 * 3);
  EXPECT_EQ(M.dims(0), 2);
  EXPECT_EQ(M.dims(1), 3);
  EXPECT_EQ(M.order(), 2);
}

TEST_P(TestMatrixWithCtxt, PopulateMatrixWithCtxt)
{
  helib::Matrix<helib::Ctxt> M(helib::Ctxt(pk), 2, 2);
  helib::Matrix<std::vector<long>> P(2, 2);

  // some values
  std::vector<long> d1 = {0, 1};
  std::vector<long> d2 = {2, 3};
  std::vector<long> d3 = {4, 5};
  std::vector<long> d4 = {6, 7};

  // encrypt some values
  ea.encrypt(M(0, 0), pk, d1);
  ea.encrypt(M(0, 1), pk, d2);
  ea.encrypt(M(1, 0), pk, d3);
  ea.encrypt(M(1, 1), pk, d4);

  // decrypt some values
  ea.decrypt(M(0, 0), sk, P(0, 0));
  ea.decrypt(M(0, 1), sk, P(0, 1));
  ea.decrypt(M(1, 0), sk, P(1, 0));
  ea.decrypt(M(1, 1), sk, P(1, 1));

  // check basics
  EXPECT_EQ(M.size(), 2 * 2);
  EXPECT_EQ(M.dims(0), 2);
  EXPECT_EQ(M.dims(1), 2);
  EXPECT_EQ(M.order(), 2);
  EXPECT_TRUE(M.fullView());

  EXPECT_EQ(P.size(), 2 * 2);
  EXPECT_EQ(P.dims(0), 2);
  EXPECT_EQ(P.dims(1), 2);
  EXPECT_EQ(P.order(), 2);
  EXPECT_TRUE(P.fullView());

  // check decrypted values
  EXPECT_TRUE(std::equal(P(0, 0).begin(), P(0, 0).end(), d1.begin()));
  EXPECT_TRUE(std::equal(P(0, 1).begin(), P(0, 1).end(), d2.begin()));
  EXPECT_TRUE(std::equal(P(1, 0).begin(), P(1, 0).end(), d3.begin()));
  EXPECT_TRUE(std::equal(P(1, 1).begin(), P(1, 1).end(), d4.begin()));
}

TEST_P(TestMatrixWithCtxt, TransposeCtxtMatrix)
{
  helib::Matrix<helib::Ctxt> M(helib::Ctxt(pk), 2, 3);
  helib::Matrix<std::vector<long>> P(3, 2);

  // some values
  std::vector<long> d1 = {0, 1};
  std::vector<long> d2 = {2, 3};
  std::vector<long> d3 = {4, 5};
  std::vector<long> d4 = {6, 7};
  std::vector<long> d5 = {8, 9};
  std::vector<long> d6 = {10, 11};

  // encrypt some values
  ea.encrypt(M(0, 0), pk, d1);
  ea.encrypt(M(0, 1), pk, d2);
  ea.encrypt(M(0, 2), pk, d3);
  ea.encrypt(M(1, 0), pk, d4);
  ea.encrypt(M(1, 1), pk, d5);
  ea.encrypt(M(1, 2), pk, d6);

  M.transpose();

  // decrypt some values
  ea.decrypt(M(0, 0), sk, P(0, 0));
  ea.decrypt(M(1, 0), sk, P(1, 0));
  ea.decrypt(M(2, 0), sk, P(2, 0));
  ea.decrypt(M(0, 1), sk, P(0, 1));
  ea.decrypt(M(1, 1), sk, P(1, 1));
  ea.decrypt(M(2, 1), sk, P(2, 1));

  // check decrypted values
  EXPECT_TRUE(std::equal(P(0, 0).begin(), P(0, 0).end(), d1.begin()));
  EXPECT_TRUE(std::equal(P(1, 0).begin(), P(1, 0).end(), d2.begin()));
  EXPECT_TRUE(std::equal(P(2, 0).begin(), P(2, 0).end(), d3.begin()));
  EXPECT_TRUE(std::equal(P(0, 1).begin(), P(0, 1).end(), d4.begin()));
  EXPECT_TRUE(std::equal(P(1, 1).begin(), P(1, 1).end(), d5.begin()));
  EXPECT_TRUE(std::equal(P(2, 1).begin(), P(2, 1).end(), d6.begin()));
}

TEST_P(TestMatrixWithCtxt, MultiplyTwoCtxtMatrices)
{
  helib::Ctxt blank(pk);

  helib::Matrix<helib::Ctxt> M1(blank, 2, 3);
  helib::Matrix<helib::Ctxt> M2(blank, 2, 3);

  // some values
  std::vector<long> d1 = {0, 1};
  std::vector<long> d2 = {2, 3};
  std::vector<long> d3 = {4, 5};
  std::vector<long> d4 = {6, 7};
  std::vector<long> d5 = {8, 9};
  std::vector<long> d6 = {10, 11};

  // some answers
  std::vector<long> a1 = {20, 35};
  std::vector<long> a2 = {56, 89};
  std::vector<long> a3 = {200, 251};

  auto mod_p = [this](long& x) { x %= this->p; };
  std::for_each(a1.begin(), a1.end(), mod_p);
  std::for_each(a2.begin(), a2.end(), mod_p);
  std::for_each(a3.begin(), a3.end(), mod_p);

  // encrypt some values
  ea.encrypt(M1(0, 0), pk, d1);
  ea.encrypt(M1(0, 1), pk, d2);
  ea.encrypt(M1(0, 2), pk, d3);
  ea.encrypt(M1(1, 0), pk, d4);
  ea.encrypt(M1(1, 1), pk, d5);
  ea.encrypt(M1(1, 2), pk, d6);

  // other matrix
  ea.encrypt(M2(0, 0), pk, d1);
  ea.encrypt(M2(0, 1), pk, d2);
  ea.encrypt(M2(0, 2), pk, d3);
  ea.encrypt(M2(1, 0), pk, d4);
  ea.encrypt(M2(1, 1), pk, d5);
  ea.encrypt(M2(1, 2), pk, d6);

  // multiply
  helib::Matrix<helib::Ctxt> M3 = M1 * M2.transpose();
  helib::Matrix<std::vector<long>> P(2, 2);

  // decrypt some values
  ea.decrypt(M3(0, 0), sk, P(0, 0));
  ea.decrypt(M3(0, 1), sk, P(0, 1));
  ea.decrypt(M3(1, 0), sk, P(1, 0));
  ea.decrypt(M3(1, 1), sk, P(1, 1));

  // Check dimensions (expect 2x2 Matrix)
  EXPECT_EQ(M3.size(), 2 * 2);
  EXPECT_EQ(M3.dims(0), 2);
  EXPECT_EQ(M3.dims(1), 2);
  EXPECT_EQ(M3.order(), 2);

  // Check result
  EXPECT_TRUE(std::equal(P(0, 0).begin(), P(0, 0).end(), a1.begin()));
  EXPECT_TRUE(std::equal(P(0, 1).begin(), P(0, 1).end(), a2.begin()));
  EXPECT_TRUE(std::equal(P(1, 0).begin(), P(1, 0).end(), a2.begin()));
  EXPECT_TRUE(std::equal(P(1, 1).begin(), P(1, 1).end(), a3.begin()));
}

TEST_P(TestMatrixWithCtxt, TestBasicCtxt)
{
  helib::Ctxt blank(pk);
  ea.encrypt(blank, pk, std::vector<long>{0, 0});
}

INSTANTIATE_TEST_SUITE_P(various_parameters,
                         TestMatrixWithCtxt,
                         ::testing::Values(BGVParameters(5, 19, 1)));

} // namespace
