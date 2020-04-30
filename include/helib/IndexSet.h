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
#ifndef HELIB_INDEXSET_H
#define HELIB_INDEXSET_H
/**
 * @file IndexSet.h
 * @brief A dynamic set of integers
 **/

#include <helib/NumbTh.h>

namespace helib {

//! @brief A dynamic set of non-negative integers.
//!
//! You can iterate through a set as follows:
//! \code
//!    for (long i = s.first(); i <= s.last(); i = s.next(i)) ...
//!    for (long i = s.last(); i >= s.first(); i = s.prev(i)) ...
//! \endcode
class IndexSet
{

  std::vector<bool> rep;
  // NOTE: modern versions of C++ are supposed
  // to implement this efficiently as a "specialized template class".
  // Older versions of C++ define the equivalent class bit_std::vector.

  long _first, _last, _card;

  // Invariant: if _card == 0, then _first = 0, _last = -1;
  // otherwise, _first (resp. _last) is the lowest (resp. highest)
  // index in the set.
  // In any case, the std::vector rep always defines the characteristic
  // function of the set.

  // private helper function
  void intervalConstructor(long low, long high);

public:
  /*** constructors ***/

  // @brief No-argument constructor, creates empty set
  IndexSet() : _first(0), _last(-1), _card(0) {}

  // @brief Constructs an interval, low to high
  IndexSet(long low, long high) { intervalConstructor(low, high); }

  // @brief Constructs a singleton set
  explicit IndexSet(long j) { intervalConstructor(j, j); }

  // copy constructor: use the built-in copy constructor

  /*** assignment ***/

  // assignment: use the built-in assignment operator

  //! @brief Returns the first element, 0 if the set is empty
  long first() const { return _first; }

  //! @brief Returns the last element, -1 if the set is empty
  long last() const { return _last; }

  //! @brief Returns the next element after j, if any; otherwise j+1
  long next(long j) const;

  // @brief Returns the previous element before j, if any; otherwise j-1
  long prev(long j) const;

  //! @brief The cardinality of the set
  long card() const { return _card; }

  //! @brief Returns true iff the set contains j
  bool contains(long j) const;

  //! @brief Returns true iff the set contains s
  bool contains(const IndexSet& s) const;

  //! @brief Returns true iff the set is disjoint from s
  bool disjointFrom(const IndexSet& s) const;

  /*** comparison ***/

  bool operator==(const IndexSet& s) const;

  bool operator!=(const IndexSet& s) const { return !(*this == s); }

  /*** update methods ***/

  //! @brief Set to the empty set
  void clear();

  //! @brief Add j to the set
  void insert(long j);

  //! @brief Remove j from the set
  void remove(long j);

  //! @brief Add s to the set (union)
  void insert(const IndexSet& s);

  //! @brief Remove s from the set (set minus)
  void remove(const IndexSet& s);

  //! @brief Retains only those elements that are also in s (intersection)
  void retain(const IndexSet& s);

  //! @brief Read-only access to an empty set
  static const IndexSet& emptySet();

  //! @brief Is this set a contiguous interval?
  bool isInterval() const { return (_card == (1 + _last - _first)); }

  /*** raw IO ***/
  void read(std::istream& str);
  void write(std::ostream& str) const;

  /*** code to allow one to write "for (long i: set)" ***/

  class iterator
  {
    friend class IndexSet;

  public:
    long operator*() const { return i_; }
    iterator& operator++()
    {
      i_ = s_.next(i_);
      return *this;
    }

    bool operator==(const iterator& other) const
    {
      return &s_ == &other.s_ && i_ == other.i_;
    }

    bool operator!=(const iterator& other) const { return !(*this == other); }

  protected:
    iterator(const IndexSet& s, long i) : s_(s), i_(i) {}

  private:
    const IndexSet& s_;
    long i_;
  };

  iterator begin() const { return iterator(*this, this->first()); }
  iterator end() const { return iterator(*this, this->last() + 1); }
};

// some high-level convenience methods...not very efficient...
// not sure if we really need these

//! @brief union
IndexSet operator|(const IndexSet& s, const IndexSet& t);

//! @brief intersection
IndexSet operator&(const IndexSet& s, const IndexSet& t);

//! @brief exclusive-or
IndexSet operator^(const IndexSet& s, const IndexSet& t);

//! @brief set minus
IndexSet operator/(const IndexSet& s, const IndexSet& t);

// I/O operator
std::ostream& operator<<(std::ostream& str, const IndexSet& set);
std::istream& operator>>(std::istream& str, IndexSet& set);

//! @brief Functional cardinality
long card(const IndexSet& s);

inline bool empty(const IndexSet& s) { return s.card() == 0; }

//! @brief Is s1 subset or equal to s2
bool operator<=(const IndexSet& s1, const IndexSet& s2);

//! @brief Is s1 strict subset of s2
bool operator<(const IndexSet& s1, const IndexSet& s2);

//! @brief Is s2 subset or equal to s2
bool operator>=(const IndexSet& s1, const IndexSet& s2);

//! @brief Is s2 strict subset of s1
bool operator>(const IndexSet& s1, const IndexSet& s2);

//! @brief Functional disjoint
inline bool disjoint(const IndexSet& s1, const IndexSet& s2)
{
  return s1.disjointFrom(s2);
}

} // namespace helib

#endif // ifndef HELIB_INDEXSET_H
