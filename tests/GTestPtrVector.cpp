/* Copyright (C) 2012-2020 IBM Corp.
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
#include <iostream>
#include <NTL/tools.h>
#include <helib/NumbTh.h>
#include <helib/PtrVector.h>
#include <helib/PtrMatrix.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

// A class with no default constructor
class MyClass
{
  int myInt;
  MyClass(){}; // private
public:
  MyClass(int i) : myInt(i) {}
  int get() const { return myInt; }
  void set(int i) { myInt = i; }
};

class GTestPtrVector : public ::testing::Test
{
protected:
  GTestPtrVector() : vLength(6), zero(0){};

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }

  const int vLength;
  MyClass zero;
};

typedef helib::PtrVector<MyClass> MyPtrVec;
typedef helib::PtrVector_VecT<MyClass> MyPtrVec_Vec;
typedef helib::PtrVector_VecPt<MyClass> MyPtrVec_VecPt;
typedef helib::PtrVector_vectorT<MyClass> MyPtrVec_vector;
typedef helib::PtrVector_vectorPt<MyClass> MyPtrVec_vectorPt;

typedef helib::PtrVector_slice<MyClass> MyPtrVec_slice; // A slice of MyPtrVec

typedef helib::PtrMatrix<MyClass> MyPtrMat;
typedef helib::PtrMatrix_Vec<MyClass> MyPtrMat_Vec;
typedef helib::PtrMatrix_ptVec<MyClass> MyPtrMat_ptVec;
typedef helib::PtrMatrix_vector<MyClass> MyPtrMat_vector;
typedef helib::PtrMatrix_ptvector<MyClass> MyPtrMat_ptvector;

// compare a "generic" vectors to pointers to vector to objects
template <typename T2>
::testing::AssertionResult pointersEqual(const MyPtrVec& a, T2& b)
{
  if (helib::lsize(a) != helib::lsize(b)) {
    return ::testing::AssertionFailure()
           << "sizes do not match (" << helib::lsize(a) << " vs "
           << helib::lsize(b) << ")";
  }
  for (long i = 0; i < helib::lsize(b); i++) {
    if (a[i] != &b[i]) {
      return ::testing::AssertionFailure()
             << "difference found in the " << i << "th position: " << a[i]
             << " vs " << &b[i];
    }
  }
  return ::testing::AssertionSuccess();
}

void test1(MyClass array[], int length, const MyPtrVec& ptrs)
{
  if (helib_test::verbose) {
    std::cout << "test1 " << std::flush;
  }
  for (int i = 0; i < length; i++)
    ASSERT_EQ(ptrs[i], &(array[i]));
  ASSERT_EQ(ptrs.numNonNull(), length);
  ASSERT_EQ(ptrs.size(), length);
  const MyClass* pt = ptrs.ptr2nonNull();
  if (length > 0) {
    ASSERT_NE(pt, nullptr) << "but length > 0";
  } else if (length <= 0) {
    ASSERT_EQ(pt, nullptr) << "however length <= 0";
  }
}

void test2(MyClass* array[], int length, const MyPtrVec& ptrs)
{
  if (helib_test::verbose) {
    std::cout << "test2 " << std::flush;
  }
  for (int i = 0; i < length; i++)
    ASSERT_EQ(ptrs[i], array[i]);
  ASSERT_EQ(ptrs.size(), length);

  ASSERT_EQ(ptrs.numNonNull(), std::min(4, length));
  ASSERT_NE(ptrs.ptr2nonNull(), nullptr);
}

void printPtrVector(const MyPtrVec& ptrs)
{
  if (helib_test::verbose) {
    for (int i = 0; i < ptrs.size(); i++) {
      std::cout << ((i == 0) ? '[' : ',');
      MyClass* pt = ptrs[i];
      if (pt == nullptr)
        std::cout << "null";
      else
        std::cout << pt->get();
    }
    std::cout << ']';
  }
}
void test3(MyPtrVec& ptrs)
{
  if (helib_test::verbose) {
    std::cout << "\nBefore resize: ";
    printPtrVector(ptrs);

    int length = ptrs.size();
    ptrs.resize(length + 1);
    ptrs[length]->set(length + 1);

    std::cout << "\n After resize: ";
    printPtrVector(ptrs);
    std::cout << std::endl;
  }
}

template <typename T>
void test4(const MyPtrMat& mat, const T& array)
{
  if (helib_test::verbose) {
    std::cout << "test4 " << std::flush;
  }
  ASSERT_EQ(mat.size(), helib::lsize(array));
  for (int i = 0; i < helib::lsize(array); i++)
    ASSERT_TRUE(pointersEqual(mat[i], array[i]));
}

template <typename T>
void test5(const MyPtrMat& mat, const T& array)
{
  if (helib_test::verbose) {
    std::cout << "test5 " << std::flush;
  }
  ASSERT_EQ(mat.size(), helib::lsize(array));
  for (int i = 0; i < helib::lsize(array); i++)
    ASSERT_TRUE(pointersEqual(mat[i], *array[i]));
}

TEST_F(GTestPtrVector, pointerVectorsRemainConsistent)
{
  MyClass zero(0);
  std::vector<MyClass> v1(vLength, zero);
  NTL::Vec<MyClass> v2(NTL::INIT_SIZE, vLength, zero);

  std::vector<MyClass*> v3(vLength, nullptr);
  for (int i = 1; i < vLength - 1; i++)
    v3[i] = &(v1[i]);

  NTL::Vec<MyClass*> v4(NTL::INIT_SIZE, vLength, nullptr);
  for (int i = 1; i < vLength - 1; i++)
    v4[i] = &(v2[i]);

  MyPtrVec_vector vv1(v1);
  MyPtrVec_VecPt vv4(v4);

  ASSERT_NO_FATAL_FAILURE(test1(&v1[0], 6, vv1));
  ASSERT_NO_FATAL_FAILURE(test1(&v2[0], 6, MyPtrVec_Vec(v2)));

  MyPtrVec_slice vs1(vv1, 1);
  ASSERT_NO_FATAL_FAILURE(test1(&v1[1], 5, vs1));
  MyPtrVec_slice vss1(vs1, 1, 3);
  ASSERT_NO_FATAL_FAILURE(test1(&v1[2], 3, vss1));

  ASSERT_NO_FATAL_FAILURE(test2(&v3[0], 6, MyPtrVec_vectorPt(v3)));
  ASSERT_NO_FATAL_FAILURE(test2(&v4[0], 6, vv4));

  MyPtrVec_slice vs4(vv4, 1);
  ASSERT_NO_FATAL_FAILURE(test2(&v4[1], 5, vs4));
  MyPtrVec_slice vss4(vs4, 1, 3);
  ASSERT_NO_FATAL_FAILURE(test2(&v4[2], 3, vss4));

  ASSERT_NO_FATAL_FAILURE(test3(vv1));

  std::vector<std::vector<MyClass>> mat1(6);
  std::vector<std::vector<MyClass>*> mat2(6);
  NTL::Vec<NTL::Vec<MyClass>> mat3(NTL::INIT_SIZE, 6);
  NTL::Vec<NTL::Vec<MyClass>*> mat4(NTL::INIT_SIZE, 6);
  for (long i = 0; i < 6; i++) {
    mat1[i].resize(4, MyClass(i));
    mat2[5 - i] = &mat1[i];

    mat3[i].SetLength(3, MyClass(i + 10));
    mat4[5 - i] = &mat3[i];
  }

  ASSERT_NO_FATAL_FAILURE(test4(MyPtrMat_vector(mat1), mat1));
  ASSERT_NO_FATAL_FAILURE(test4(MyPtrMat_Vec(mat3), mat3));

  ASSERT_NO_FATAL_FAILURE(test5(MyPtrMat_ptvector(mat2), mat2));
  ASSERT_NO_FATAL_FAILURE(test5(MyPtrMat_ptVec(mat4), mat4));
  if (helib_test::verbose) {
    // Tidy up newline-free previous output
    std::cout << std::endl;
  }
}

} // namespace
