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

#include <helib/ClonedPtr.h>
#include "test_common.h"
#include "gtest/gtest.h"

namespace {

struct SimpleClonable
{
  int x;
  int y;

  SimpleClonable(int x_, int y_) : x(x_), y(y_) {}

  bool operator==(const SimpleClonable& other) const
  {
    return x == other.x && y == other.y;
  }

  SimpleClonable* clone() const { return new SimpleClonable(x, y); }
};

struct WithCloneMethod
{

  int someInt;
  int* somePtr1;
  std::unique_ptr<int> somePtr2;
  std::shared_ptr<int> somePtr3;

  WithCloneMethod() = delete;
  ~WithCloneMethod()
  {
    if (somePtr1 != nullptr)
      delete somePtr1;
  }
  WithCloneMethod(const WithCloneMethod& other) = delete;
  WithCloneMethod(WithCloneMethod&& other) = delete;
  WithCloneMethod& operator=(const WithCloneMethod& other) = delete;
  WithCloneMethod& operator=(WithCloneMethod&& other) = delete;

  WithCloneMethod(int i, int r, int u, int s) :
      someInt(i),
      somePtr1(new int(r)),
      somePtr2(std::make_unique<int>(u)),
      somePtr3(std::make_shared<int>(s))
  {}

  WithCloneMethod* clone() const
  {
    return new WithCloneMethod(someInt, *somePtr1, *somePtr2, *somePtr3);
  }
};

struct Copyable
{
  int someInt;
  int* rawPtr;
  std::shared_ptr<int> shrPtr;

  Copyable() = delete;
  Copyable(int i, int r, int s) :
      someInt(i), rawPtr(new int(r)), shrPtr(std::make_shared<int>(s))
  {}
  Copyable(const Copyable&) = default;
  Copyable(Copyable&&) = delete;
  Copyable& operator=(const Copyable&) = delete;
  Copyable& operator=(Copyable&&) = delete;
  ~Copyable() = default; // Avoid double free.
};

TEST(TestClonedPtr, deepCopyObject)
{
  helib::ClonedPtr<WithCloneMethod> to_clone(new WithCloneMethod(1, 2, 3, 4));
  helib::ClonedPtr<WithCloneMethod> cloned(to_clone);

  // Check pointers
  EXPECT_NE(&to_clone, &cloned);
  EXPECT_NE(to_clone->somePtr1, cloned->somePtr1);
  EXPECT_NE(to_clone->somePtr2.get(), cloned->somePtr2.get());
  EXPECT_NE(to_clone->somePtr3.get(), cloned->somePtr3.get());

  // Check values same
  EXPECT_EQ(to_clone->someInt, cloned->someInt);
  EXPECT_EQ(*to_clone->somePtr1, *cloned->somePtr1);
  EXPECT_EQ(*to_clone->somePtr2, *cloned->somePtr2);
  EXPECT_EQ(*to_clone->somePtr3, *cloned->somePtr3);
}

TEST(TestClonedPtr, shallowCopyObject)
{
  helib::CopiedPtr<Copyable> to_copy(new Copyable(1, 2, 3));
  helib::CopiedPtr<Copyable> copied(to_copy);

  // Check pointers
  EXPECT_NE(&to_copy, &copied);
  EXPECT_EQ(to_copy->rawPtr, copied->rawPtr);
  EXPECT_EQ(to_copy->shrPtr.get(), copied->shrPtr.get());

  // Check values same
  EXPECT_EQ(to_copy->someInt, copied->someInt);
  EXPECT_EQ(*to_copy->rawPtr, *copied->rawPtr);
  EXPECT_EQ(*to_copy->shrPtr, *copied->shrPtr);
}

TEST(TestClonedPtr, moveConstructor)
{
  helib::ClonedPtr<WithCloneMethod> clone1(new WithCloneMethod(1, 2, 3, 4));
  auto clone2(std::move(clone1));

  // Check values same
  EXPECT_EQ(clone2->someInt, 1);
  EXPECT_EQ(*clone2->somePtr1, 2);
  EXPECT_EQ(*clone2->somePtr2, 3);
  EXPECT_EQ(*clone2->somePtr3, 4);

  // Check original has been nulled
  EXPECT_EQ(clone1.get(), nullptr);
}

TEST(TestClonedPtr, moveAssignment)
{
  helib::ClonedPtr<WithCloneMethod> clone1(new WithCloneMethod(1, 2, 3, 4));
  helib::ClonedPtr<WithCloneMethod> clone2;
  clone2 = std::move(clone1);

  // Check values same
  EXPECT_EQ(clone2->someInt, 1);
  EXPECT_EQ(*clone2->somePtr1, 2);
  EXPECT_EQ(*clone2->somePtr2, 3);
  EXPECT_EQ(*clone2->somePtr3, 4);

  // Check original has been nulled
  EXPECT_EQ(clone1.get(), nullptr);
}

TEST(TestClonedPtr, relationalOperators)
{
  helib::ClonedPtr<WithCloneMethod> clone1(new WithCloneMethod(1, 2, 3, 4));
  auto clone2 = clone1;

  // Check values same
  EXPECT_EQ(clone1 == clone2, clone1.get() == clone2.get());
  EXPECT_EQ(clone1 != clone2, clone1.get() != clone2.get());
  EXPECT_EQ(clone1 < clone2, clone1.get() < clone2.get());
  EXPECT_EQ(clone1 > clone2, clone1.get() > clone2.get());
  EXPECT_EQ(clone1 <= clone2, clone1.get() <= clone2.get());
  EXPECT_EQ(clone1 >= clone2, clone1.get() >= clone2.get());
}

TEST(TestClonedPtr, swapClonedPtrsByStdSwap)
{
  helib::ClonedPtr<WithCloneMethod> clone1(new WithCloneMethod(1, 2, 3, 4));
  helib::ClonedPtr<WithCloneMethod> clone2(new WithCloneMethod(5, 6, 7, 8));

  std::swap(clone1, clone2);

  // Check values same
  EXPECT_EQ(clone1->someInt, 5);
  EXPECT_EQ(*clone1->somePtr1, 6);
  EXPECT_EQ(*clone1->somePtr2, 7);
  EXPECT_EQ(*clone1->somePtr3, 8);

  // Check values same
  EXPECT_EQ(clone2->someInt, 1);
  EXPECT_EQ(*clone2->somePtr1, 2);
  EXPECT_EQ(*clone2->somePtr2, 3);
  EXPECT_EQ(*clone2->somePtr3, 4);
}

TEST(TestClonedPtr, swapClonedPtrsByMethod)
{
  helib::ClonedPtr<WithCloneMethod> clone1(new WithCloneMethod(1, 2, 3, 4));
  helib::ClonedPtr<WithCloneMethod> clone2(new WithCloneMethod(5, 6, 7, 8));

  clone1.swap(clone2);

  // Check values same
  EXPECT_EQ(clone1->someInt, 5);
  EXPECT_EQ(*clone1->somePtr1, 6);
  EXPECT_EQ(*clone1->somePtr2, 7);
  EXPECT_EQ(*clone1->somePtr3, 8);

  // Check values same
  EXPECT_EQ(clone2->someInt, 1);
  EXPECT_EQ(*clone2->somePtr1, 2);
  EXPECT_EQ(*clone2->somePtr2, 3);
  EXPECT_EQ(*clone2->somePtr3, 4);
}

TEST(TestClonedPtr, makeClonedPtr)
{
  auto clone = helib::makeClonedPtr<WithCloneMethod>(1, 2, 3, 4);

  // Check values same
  EXPECT_EQ(clone->someInt, 1);
  EXPECT_EQ(*clone->somePtr1, 2);
  EXPECT_EQ(*clone->somePtr2, 3);
  EXPECT_EQ(*clone->somePtr3, 4);
}

TEST(TestClonedPtr, resetClonedPtr)
{
  auto clone = helib::makeClonedPtr<WithCloneMethod>(1, 2, 3, 4);

  // Reset to new values.
  clone.reset(new WithCloneMethod(7, 8, 9, 0));

  // Check values same
  EXPECT_EQ(clone->someInt, 7);
  EXPECT_EQ(*clone->somePtr1, 8);
  EXPECT_EQ(*clone->somePtr2, 9);
  EXPECT_EQ(*clone->somePtr3, 0);

  // Reset to null.
  clone.reset();

  EXPECT_EQ(clone.get(), nullptr);
}

TEST(TestClonedPtr, releaseClonedPtr)
{
  auto clone = helib::makeClonedPtr<WithCloneMethod>(1, 2, 3, 4);

  // Release pointer.
  auto ptr = clone.release();

  // Check values same.
  EXPECT_EQ(ptr->someInt, 1);
  EXPECT_EQ(*ptr->somePtr1, 2);
  EXPECT_EQ(*ptr->somePtr2, 3);
  EXPECT_EQ(*ptr->somePtr3, 4);

  // Clone should be null.
  EXPECT_EQ(clone.get(), nullptr);
}

TEST(TestClonedPtr, cloneOtherPtrTypesRawPtr)
{
  SimpleClonable* a = new SimpleClonable(1, 2);
  SimpleClonable* a_clone = helib::clonePtr(a);

  // Point to same values not same object.
  EXPECT_EQ(*a, *a_clone);
  EXPECT_NE(&a, &a_clone);
  // Objects at different addresses.
  EXPECT_NE(a, a_clone);
}

TEST(TestClonedPtr, cloneOtherPtrTypesUniquePtr)
{
  std::unique_ptr<SimpleClonable> a = std::make_unique<SimpleClonable>(3, 4);
  std::unique_ptr<SimpleClonable> a_clone = helib::clonePtr(a);

  // a and a_clone the same values.
  EXPECT_EQ(*a, *a_clone);
  // ptr different addresses.
  EXPECT_NE(&a, &a_clone);
  // obj different addresses.
  EXPECT_NE(a.get(), a_clone.get());
}

TEST(TestClonedPtr, cloneOtherPtrTypesSharedPtr)
{
  std::shared_ptr<SimpleClonable> a = std::make_shared<SimpleClonable>(3, 4);
  std::shared_ptr<SimpleClonable> a_clone = helib::clonePtr(a);

  // a and a_clone the same values.
  EXPECT_EQ(*a, *a_clone);
  // ptr different addresses.
  EXPECT_NE(&a, &a_clone);
  // obj different addresses.
  EXPECT_NE(a.get(), a_clone.get());
}

TEST(TestClonedPtr, copyConstructorPromotingNullptrSpeacialCase)
{
  helib::ClonedPtr<int> a = nullptr;
  EXPECT_EQ(a, nullptr);
}

} // namespace
