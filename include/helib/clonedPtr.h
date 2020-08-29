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
/**
 * @file clonedPtr.h
 * @brief Implementation of smart pointers with "deep cloning" semantics
 *
 * Based (loosely) on code from
 *
 *  http://github.com/yonat/smart_ptr/blob/master/cloned_ptr.h
 **/
/*
 * This seems like a very useful template class for smart
 * pointers with "deep cloning" semantics and class hierarchies.
 *
 * To work, a class needs a "clone" method that creates a new copy.
 *
 * NOTES:
 *
 * this did not compile...had to change the definition of the class clone<X>
 *
 * added boolean operators
 *
 * added assign operation set_ptr(X* p) removed get(), introduced public
 * method get_ptr()... these names chosen to make it easier to hunt down
 * "loopholes"
 *
 * added noclone<X> template
 *
 * Made "const"-ness apply to both the pointer and the object pointed to.
 * For this kind of "deep copy" semantics, this makes the most sense.
 *
 * Made a corresponding class copied_ptr, with "shallow" cloning...
 * To bad this can't really be done with templates... C++11 defines
 * the concept of "template aliases", which would do the trick...
 * unfortunately, gcc does not yet implement this...
 *
 * Changed the name "clone" to avoid conflicts...
 *****/

/* For ANSI-challenged compilers, you may want to
 * #define NO_MEMBER_TEMPLATES or explicit */

namespace helib {

/**
 * @class deep_clone
 * @brief Deep copy: initialize with clone
 * @tparam X The class to which this points
 **/
template <typename X>
class deep_clone
{
public:
  static X* apply(const X* x) { return x->clone(); }
};

/**
 * @class shallow_clone
 * @brief Shallow copy: initialize with copy constructor
 * @tparam X The class to which this points
 **/
template <typename X>
class shallow_clone
{
public:
  static X* apply(const X* x) { return new X(*x); }
};

#ifndef NO_MEMBER_TEMPLATES

#define CLONED_PTR_TEMPLATE_MEMBERS(CLONED_PTR_TYPE)                           \
                                                                               \
  template <typename Y>                                                        \
  CLONED_PTR_TYPE(const CLONED_PTR_TYPE<Y>& r)                                 \
  {                                                                            \
    copy(r.ptr);                                                               \
  }                                                                            \
  template <typename Y>                                                        \
  CLONED_PTR_TYPE& operator=(const CLONED_PTR_TYPE<Y>& r)                      \
  {                                                                            \
    if (this != &r) {                                                          \
      delete ptr;                                                              \
      copy(r.ptr);                                                             \
    }                                                                          \
    return *this;                                                              \
  }

#else

#define CLONED_PTR_TEMPLATE_MEMBERS(CLONED_PTR_TYPE)

#endif

#define CLONED_PTR_DECLARE(CLONED_PTR_TYPE, CLONED_PTR_INIT)                   \
                                                                               \
  template <typename X, typename Cloner = CLONED_PTR_INIT<X>>                  \
  class CLONED_PTR_TYPE                                                        \
  {                                                                            \
  public:                                                                      \
    typedef X element_type;                                                    \
                                                                               \
    explicit CLONED_PTR_TYPE(X* p = 0) : ptr(p) {}                             \
    ~CLONED_PTR_TYPE() { delete ptr; }                                         \
    CLONED_PTR_TYPE(const CLONED_PTR_TYPE& r) { copy(r.ptr); }                 \
                                                                               \
    CLONED_PTR_TYPE& operator=(const CLONED_PTR_TYPE& r)                       \
    {                                                                          \
      if (this != &r) {                                                        \
        delete ptr;                                                            \
        copy(r.ptr);                                                           \
      }                                                                        \
      return *this;                                                            \
    }                                                                          \
                                                                               \
    void set_ptr(X* p)                                                         \
    {                                                                          \
      if (ptr != p) {                                                          \
        delete ptr;                                                            \
        ptr = p;                                                               \
      }                                                                        \
    }                                                                          \
                                                                               \
    CLONED_PTR_TEMPLATE_MEMBERS(CLONED_PTR_TYPE)                               \
                                                                               \
    const X& operator*() const { return *ptr; }                                \
    X& operator*() { return *ptr; }                                            \
                                                                               \
    const X* operator->() const { return ptr; }                                \
    X* operator->() { return ptr; }                                            \
                                                                               \
    bool null() const { return ptr == 0; }                                     \
                                                                               \
    const X* get_ptr() const { return ptr; }                                   \
    X* get_ptr() { return ptr; }                                               \
                                                                               \
    void swap(CLONED_PTR_TYPE& r)                                              \
    {                                                                          \
      X* tmp;                                                                  \
      tmp = r.ptr;                                                             \
      r.ptr = ptr;                                                             \
      ptr = tmp;                                                               \
    }                                                                          \
                                                                               \
  private:                                                                     \
    X* ptr;                                                                    \
                                                                               \
    void copy(X* p) { ptr = (p ? Cloner::apply(p) : 0); }                      \
  };

// declare the template class cloned_ptr<X>
CLONED_PTR_DECLARE(cloned_ptr, deep_clone)

// declare the template class copied_ptr<X>
CLONED_PTR_DECLARE(copied_ptr, shallow_clone)

// template <class X> using copied_ptr = cloned_ptr<X, noclone<X> >;

template <typename X, typename Cloner>
void swap(cloned_ptr<X, Cloner>& x, cloned_ptr<X, Cloner>& y)
{
  x.swap(y);
}

template <typename X, typename Cloner>
void swap(copied_ptr<X, Cloner>& x, copied_ptr<X, Cloner>& y)
{
  x.swap(y);
}

} // namespace helib

#endif // ifndef HELIB_CLONEDPTR_H
