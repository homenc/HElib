/* Copyright (C) 2012-2017 IBM Corp.
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
#include <cassert>
#include <iostream>
#include <NTL/tools.h>
#include "NumbTh.h"
#include "PtrVector.h"
#include "PtrMatrix.h"


// A class with no default constructor
class MyClass {
  int myInt;
  MyClass(){}; // private
public:
  MyClass(int i): myInt(i){}
  int get() const { return myInt; }
  void set(int i) { myInt=i; }
};

typedef PtrVector<MyClass> MyPtrVec;
typedef PtrVector_VecT<MyClass> MyPtrVec_Vec;
typedef PtrVector_VecPt<MyClass> MyPtrVec_VecPt;
typedef PtrVector_vectorT<MyClass> MyPtrVec_vector;
typedef PtrVector_vectorPt<MyClass> MyPtrVec_vectorPt;

typedef PtrVector_slice<MyClass> MyPtrVec_slice; // A slice of MyPtrVec

typedef PtrMatrix<MyClass> MyPtrMat;
typedef PtrMatrix_Vec<MyClass> MyPtrMat_Vec;
typedef PtrMatrix_ptVec<MyClass> MyPtrMat_ptVec;
typedef PtrMatrix_vector<MyClass> MyPtrMat_vector;
typedef PtrMatrix_ptvector<MyClass> MyPtrMat_ptvector;

// compare a "generic" vectors to pointers to vector to objects
template<typename T2>
bool compPointers(const MyPtrVec& a, T2& b)
{
  if (lsize(a)!=lsize(b)) return false;
  for (long i=0; i<lsize(b); i++)
    if (a[i]!=&b[i]) return false;
  return true;
}

void test1(MyClass array[], int length, const MyPtrVec& ptrs)
{
  std::cout <<"test1 "<<std::flush;
  for (int i=0; i<length; i++) assert(ptrs[i]==&(array[i]));
  assert (ptrs.numNonNull()==length);
  assert (ptrs.size()==length);
  const MyClass* pt=ptrs.ptr2nonNull();
  assert ((length>0 && pt!=nullptr) || (length<=0 && pt==nullptr));
}

void test2(MyClass* array[], int length, const MyPtrVec& ptrs)
{
  std::cout <<"test2 "<<std::flush;
  for (int i=0; i<length; i++) assert(ptrs[i]==array[i]);
  assert (ptrs.size()==length);
  
  assert (ptrs.numNonNull()==std::min(4,length));
  assert (ptrs.ptr2nonNull() != nullptr);
}

void printPtrVector(const MyPtrVec& ptrs)
{
  for (int i=0; i<ptrs.size(); i++) {
    std::cout << ((i==0)? '[' : ',');
    MyClass* pt = ptrs[i];
    if (pt==nullptr) std::cout << "null";
    else             std::cout << pt->get();
  }
  std::cout << ']';
}
void test3(MyPtrVec& ptrs)
{
  std::cout << "\nBefore resize: ";
  printPtrVector(ptrs);

  int length = ptrs.size();
  ptrs.resize(length+1);
  ptrs[length]->set(length+1);
  
  std::cout << "\n After resize: ";
  printPtrVector(ptrs);
  std::cout << std::endl;
}

template<typename T>
void test4(const MyPtrMat& mat, const T& array)
{
  std::cout <<"test4 "<<std::flush;
  assert(mat.size()==lsize(array));
  for (int i=0; i<lsize(array); i++)
    assert(compPointers(mat[i], array[i]));
}

template<typename T>
void test5(const MyPtrMat& mat, const T& array)
{
  std::cout <<"test5 "<<std::flush;
  assert(mat.size()==lsize(array));
  for (int i=0; i<lsize(array); i++)
    assert(compPointers(mat[i], *array[i]));
}


int main()
{
  const int vLength=6;

  MyClass zero(0);
  std::vector<MyClass> v1(vLength, zero);
  NTL::Vec<MyClass> v2(NTL::INIT_SIZE, vLength, zero);

  std::vector<MyClass*> v3(vLength, nullptr);
  for (int i=1; i<vLength-1; i++) v3[i] = &(v1[i]);

  NTL::Vec<MyClass*> v4(NTL::INIT_SIZE, vLength, nullptr);
  for (int i=1; i<vLength-1; i++) v4[i] = &(v2[i]);

  MyPtrVec_vector vv1(v1);
  MyPtrVec_VecPt  vv4(v4);

  test1(&v1[0], 6, vv1);
  test1(&v2[0], 6, MyPtrVec_Vec(v2));

  MyPtrVec_slice vs1(vv1,1);
  test1(&v1[1], 5, vs1);
  MyPtrVec_slice vss1(vs1,1,3);
  test1(&v1[2], 3, vss1);

  test2(&v3[0], 6, MyPtrVec_vectorPt(v3));
  test2(&v4[0], 6, vv4);

  MyPtrVec_slice vs4(vv4,1);
  test2(&v4[1], 5, vs4);
  MyPtrVec_slice vss4(vs4,1,3);
  test2(&v4[2], 3, vss4);

  test3(vv1);

  std::vector< std::vector<MyClass> > mat1(6);
  std::vector< std::vector<MyClass>* > mat2(6);
  NTL::Vec< NTL::Vec<MyClass> > mat3(NTL::INIT_SIZE, 6);
  NTL::Vec< NTL::Vec<MyClass>* > mat4(NTL::INIT_SIZE, 6);
  for (long i=0; i<6; i++) {
    mat1[i].resize(4, MyClass(i));
    mat2[5-i] = &mat1[i];

    mat3[i].SetLength(3, MyClass(i+10));
    mat4[5-i] = &mat3[i];
  }

  test4(MyPtrMat_vector(mat1), mat1);
  test4(MyPtrMat_Vec(mat3), mat3);

  test5(MyPtrMat_ptvector(mat2), mat2);
  test5(MyPtrMat_ptVec(mat4), mat4);

  std::cout <<"All done\n";
}
