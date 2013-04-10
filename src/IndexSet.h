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

#include <vector>
#include <iostream>
#include <cassert>

using namespace std;

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

  // no-arg constructor - creates empty set
  IndexSet() {
    _first = 0;  _last = -1; _card = 0;
  }

  // constructs an interval, low to high
  IndexSet(long low, long high) {
    intervalConstructor(low, high);
  }

  // constructs a singleton set
  explicit
  IndexSet(long j) {
    intervalConstructor(j, j);
  }

  // copy constructor: use the built-in copy constructor

  /*** asignment ***/

  // assignment: use the built-in assignment operator

  /*** access methods ***/

  long first() const { return _first; }
  long last() const { return _last; }
  long card() const { return _card; }

  // return the next element after j, if any;
  // otherwise j+1

  long next(long j) const;

  // return the previous element before j, if any;
  // otherwise j-1
  long prev(long j) const;

  // NOTE: you can iterate through a set as follows:
  //    for (long i = s.first(); i <= s.last(); i = s.next(i)) ...

  // returns true iff the set contains j
  bool contains(long j) const;

  // returns true iff the set contains s
  bool contains(const IndexSet& s) const;

  // returns true iff the set is disjoint from s
  bool disjointFrom(const IndexSet& s) const;

  /*** comparison ***/

  bool operator==(const IndexSet& s) const;

  bool operator!=(const IndexSet& s) const {
    return !(*this == s);
  }

  /*** update methods ***/

  // set to the empty set
  void clear();

  // add j to the set
  void insert(long j);

  // remove j from the set
  void remove(long j);

  // add s to the set (union)
  void insert(const IndexSet& s);

  // remove s from the set (set minus)
  void remove(const IndexSet& s);

  // retain only those elements that are also in s (intersection)
  void retain(const IndexSet& s);

  // read-only access to an empty set
  static const IndexSet& emptySet();
};

// some high-level convenience methods...not very efficient...
// not sure if we really need these

// union
IndexSet operator|(const IndexSet& s, const IndexSet& t);

// intersection
IndexSet operator&(const IndexSet& s, const IndexSet& t);

// exclusive-or
IndexSet operator^(const IndexSet& s, const IndexSet& t);

// set minus
IndexSet operator/(const IndexSet& s, const IndexSet& t);

// I/O operator
ostream& operator << (ostream& str, const IndexSet& set);
istream& operator >> (istream& str, IndexSet& set);

// functional card
long card(const IndexSet& s);

inline bool empty(const IndexSet& s) { return s.card()==0; }

// functional "contains"
bool operator<=(const IndexSet& s1, const IndexSet& s2);

bool operator<(const IndexSet& s1, const IndexSet& s2);

bool operator>=(const IndexSet& s1, const IndexSet& s2);

bool operator>(const IndexSet& s1, const IndexSet& s2);

// functional disjoint
inline bool disjoint(const IndexSet& s1, const IndexSet& s2)
{ return s1.disjointFrom(s2); }

#endif
