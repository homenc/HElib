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
#ifndef HELIB_CLONEDPTR_H
#define HELIB_CLONEDPTR_H

#include <utility> // std::swap
/**
 * @file ClonedPtr.h
 * @brief Implementation of smart pointer `ClonedPtr` with "cloning" semantics.
 * The API should be familiar to those that use the standard library smart
 * pointers (unique_ptr and shared_ptr).
 *
 * Two "cloning" policies are provided; one for "deep cloning" (default) and the
 * other "shallow cloning". For an object to be managed by `ClonedPtr`, a class
 * needs a "clone" method for deep cloning that creates a raw pointer to a copy
 * of the object or a copy constructor for shallow cloning.  An alias
 * `CopiedPtr` is provided for convenience which is a `ClonedPtr` with shallow
 * clone policy.
 **/

namespace helib {

/**
 * @brief Clone any type of pointer.
 * @param ptr Pointer to be cloned. That is the object that it points to is
 *`cloned`.
 * @return A new pointer object that points to a freshly cloned object.
 * @note Requires that the object being cloned has a `clone` method.
 **/
template <typename PTR>
PTR clonePtr(PTR& ptr)
{
  static_assert(std::is_pointer<decltype(ptr->clone())>::value,
                "`clone` method must return raw pointer.");

  return PTR(ptr->clone());
}

/**
 * @struct DeepCopy
 * @brief Deep copy an object by using its clone method.
 * @tparam T The type of object to be cloned.
 **/
template <typename T>
struct DeepCopy
{
  /**
   * @brief Apply the cloning policy.
   * @param ptr The pointer to apply clone policy.
   * @return The pointer to apply clone policy.
   **/
  static T* apply(const T* ptr)
  {
    static_assert(std::is_pointer<decltype(ptr->clone())>::value,
                  "`clone` method must return raw pointer.");

    return ptr->clone();
  }
};

/**
 * @struct ShallowCopy
 * @brief Shallow copy an object by using its copy constructor.
 * @tparam T The type of object to be cloned.
 **/
template <typename T>
struct ShallowCopy
{
  /**
   * @brief Apply the cloning policy.
   * @param ptr The pointer to apply clone policy.
   * @return The pointer to apply clone policy.
   **/
  static T* apply(const T* ptr)
  {
    static_assert(std::is_copy_constructible<T>::value,
                  "Must be copy constructable.");

    return new T(*ptr);
  }
};

/**
 * @class ClonedPtr
 * @brief A smart pointer that `clones` the object it holds when it is copied.
 * @tparam T The type of object managed by the pointer.
 * @tparam Cloner A policy of how to clone the object pointed to, either `deep`
 * or `shallow`. The default is `deep`.
 **/
template <typename T, typename Cloner = DeepCopy<T>>
class ClonedPtr
{
public:
  /**
   * @brief The type of the managed object.
   **/
  using value_type = T;

  /**
   * @brief Explicit constructor to create new `ClonedPtr` object.
   * @param p The raw pointer to the object that `ClonedPtr` will manage.
   **/
  explicit ClonedPtr(T* p = nullptr) : ptr(p) {}

  /**
   * @brief Special constructor to create new `ClonedPtr` object promoting a
   * nullptr.
   **/
  ClonedPtr(std::nullptr_t) : ptr(nullptr) {}

  /**
   * @brief The destructor deletes the managed object.
   **/
  ~ClonedPtr()
  {
    if (ptr != nullptr)
      delete ptr;
  }

  /**
   * @brief Copy construct by cloning the managed object based on the
   * cloning policy.
   * @param other The `ClonedPtr` object to be cloned.
   **/
  ClonedPtr(const ClonedPtr& other) : ptr(copy(other.ptr)) {}

  /**
   * @brief Move construct by taking ownership of the managed object by the
   * `ClonedPtr` object being `moved` and setting its pointer to `nullptr`.
   * @param other The `ClonedPtr` object to be moved.
   **/
  ClonedPtr(ClonedPtr&& other) noexcept : ptr{std::exchange(other.ptr, nullptr)}
  {}

  /**
   * @brief Copy assign by cloning the managed object based on the
   * cloning policy.
   * @param other The `ClonedPtr` object to be cloned.
   * @return Reference to `*this`.
   **/
  ClonedPtr& operator=(const ClonedPtr& other)
  {
    if (this != &other) {
      delete ptr;
      ptr = copy(other.ptr);
    }
    return *this;
  }

  /**
   * @brief Move construct by taking ownership of the managed object by the
   * `ClonedPtr` object being `moved` and setting its pointer to `nullptr`.
   * @param other The `ClonedPtr` object to be moved.
   * @return Reference to `*this`.
   **/
  ClonedPtr& operator=(ClonedPtr&& other) noexcept
  {
    if (this != &other) {
      delete ptr;
      ptr = std::exchange(other.ptr, nullptr);
    }
    return *this;
  }

  /**
   * @brief Reset method deletes the object that it currently managed and
   * manages the object given by a raw pointer.
   * @param p The raw pointer to replace the raw pointer currently managed by
   * the `ClonedPtr` object. The default is `nullptr`.
   **/
  void reset(T* p = nullptr)
  {
    if (ptr != p) {
      if (ptr != nullptr)
        delete ptr;
      ptr = p;
    }
  }

  /**
   * @brief Release the ownership of the currently managed object by
   * returning the pointer and setting the held pointer to nullptr.
   * @return the pointer currently held by the `ClonedPtr` object.
   **/
  T* release() noexcept
  {
    T* ptr_to_release = ptr;
    ptr = nullptr;
    return ptr_to_release;
  }

  /**
   * @brief convert `ClonedPtr` to `bool`, `false` if there is a managed object
   * is else converts to `true`.
   **/
  explicit operator bool() const noexcept { return ptr != nullptr; }

  /**
   * @brief Dereference the smart pointer.
   * @return `const` reference to manged object.
   **/
  const T& operator*() const { return *ptr; }

  /**
   * @brief Dereference the smart pointer.
   * @return Reference to managed object.
   **/
  T& operator*() { return *ptr; }

  /**
   * @brief Struct dereference to access members of the managed object.
   * @return The `const` raw pointer to the managed object.
   **/
  const T* operator->() const { return ptr; }

  /**
   * @brief Struct dereference to access members of the managed object.
   * @return The raw pointer to the managed object.
   **/
  T* operator->() { return ptr; }

  /**
   * @brief Get the pointer managed by `ClonedPtr` object.
   * @return The `const` raw pointer to the managed object.
   **/
  const T* get() const { return ptr; }

  /**
   * @brief Get the pointer managed by `ClonedPtr` object.
   * @return The raw pointer to the managed object.
   **/
  T* get() { return ptr; }

  /**
   * @brief Equal to another `ClonedPtr` object.
   * @param other `ClonedPtr` object to compare to.
   * @return `true` if this `ClonedPtr` manages the same object as `other`
   * else `false`.
   **/
  bool operator==(const ClonedPtr& other) const { return ptr == other.ptr; }

  /**
   * @brief Not equal to another `ClonedPtr` object.
   * @param other `ClonedPtr` object to compare to.
   * @return `true` if this `ClonedPtr` does not manage the same object as
   * `other` else `false`.
   **/
  bool operator!=(const ClonedPtr& other) const { return ptr != other.ptr; }

  /**
   * @brief Less than another `ClonedPtr` object.
   * @param other `ClonedPtr` object to compare to.
   * @return `true` if the address of the managed object is less than that of
   * the object managed by `other` else `false`.
   **/
  bool operator<(const ClonedPtr& other) const { return ptr < other.ptr; }

  /**
   * @brief Greater than another `ClonedPtr` object.
   * @param other `ClonedPtr` object to compare to.
   * @return `true` if the address of the managed object is greater than that of
   * the object managed by `other` else `false`.
   **/
  bool operator>(const ClonedPtr& other) const { return ptr > other.ptr; }

  /**
   * @brief Less than or equal to another `ClonedPtr` object.
   * @param other `ClonedPtr` object to compare to.
   * @return `true` if the address of the managed object is less than or equal
   * to that of the object managed by `other` else `false`.
   **/
  bool operator<=(const ClonedPtr& other) const { return ptr <= other.ptr; }

  /**
   * @brief Greater than or equal to another `ClonedPtr` object.
   * @param other `ClonedPtr` object to compare to.
   * @return `true` if the address of the managed object is greater than or
   *equal to that of the object managed by `other` else `false`.
   **/
  bool operator>=(const ClonedPtr& other) const { return ptr >= other.ptr; }

  /**
   * @brief Swap the managed objects of the this `ClonedPtr` object and
   * another.
   * @param other `ClonedPtr` object to swap managed object with.
   **/
  void swap(ClonedPtr& other) { std::swap(ptr, other.ptr); }

private:
  // Copy function applies cloning policy.
  T* copy(T* p) { return (p ? Cloner::apply(p) : nullptr); }

  // The raw pointer to the managed object.
  T* ptr = nullptr;
};

/**
 * @brief `CopiedPtr` is an alias to `ClonedPtr` but with the cloning policy
 * set to `shallow` copy.
 * @tparam T Type of the object pointed to by `CopiedPtr`.
 **/
template <typename T>
using CopiedPtr = ClonedPtr<T, ShallowCopy<T>>;

/**
 * @brief Construct an object of type `T` based on arguments passed and point
 * to it with a `ClonedPtr` object.
 * @param args The arguments to be forwarded to create an object that will be
 * pointed to by `ClonedPtr`.
 * @tparam T The type of object that will be constructed and pointed to by the
 * `ClonedPtr` object.
 * @tparam Args Packed type for arguments to construct a `ClonedPtr`.
 * @return A `ClonedPtr` object pointing to the newly constructed object.
 **/
template <typename T, typename... Args>
ClonedPtr<T> makeClonedPtr(Args&&... args)
{
  return ClonedPtr<T>(new T(std::forward<Args>(args)...));
}

} // namespace helib

#endif // HELIB_CLONEDPTR_H
