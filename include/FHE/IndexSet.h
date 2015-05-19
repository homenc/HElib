/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#ifndef _IndexSet
#define _IndexSet
/**
 * @file IndexSet.h
 * @brief A dynamic set of integers
 **/

#include "NumbTh.h"

//! @brief A dynamic set of non-negative integers.
//!
//! You can iterate through a set as follows:
//! \code 
//!    for (long i = s.first(); i <= s.last(); i = s.next(i)) ...
//!    for (long i = s.last(); i >= s.first(); i = s.prev(i)) ...
//! \endcode
class IndexSet {

  vector<bool> rep;
  // NOTE: modern versions of C++ are supposed
  // to implement this efficiently as a "specialized template class".
  // Older versions of C++ define the equivalent class bit_vector.

  long _first, _last, _card;

  // Invariant: if _card == 0, then _first = 0, _last = -1;
  // otherwise, _first (resp. _last) is the lowest (resp. highest)
  // index in the set.
  // In any case, the vector rep always defines the characterstic
  // function of the set.

  // private helper function
  void intervalConstructor(long low, long high);

public:

  /*** constructors ***/

  // @brief No-argument constructor, creates empty set
  IndexSet() {
    _first = 0;  _last = -1; _card = 0;
  }

  // @brief Constructs an interval, low to high
  IndexSet(long low, long high) {
    intervalConstructor(low, high);
  }

  // @brief Constructs a singleton set
  explicit
  IndexSet(long j) {
    intervalConstructor(j, j);
  }

  // copy constructor: use the built-in copy constructor

  /*** asignment ***/

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

  bool operator!=(const IndexSet& s) const {
    return !(*this == s);
  }

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
  bool isInterval() const {return (_card==(1+_last-_first));}
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
ostream& operator << (ostream& str, const IndexSet& set);
istream& operator >> (istream& str, IndexSet& set);

//! @brief Functional cardinality
long card(const IndexSet& s);

inline bool empty(const IndexSet& s) { return s.card()==0; }

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
{ return s1.disjointFrom(s2); }

#endif
