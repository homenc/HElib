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
#ifndef HELIB_INDEXMAP_H
#define HELIB_INDEXMAP_H
/**
 * @file IndexMap.h
 * @brief Implementation of a map indexed by a dynamic set of integers.
 **/

#include <unordered_map>
#include <helib/IndexSet.h>
#include <helib/clonedPtr.h>

namespace helib {

//! @brief Initializing elements in an IndexMap
template <typename T>
class IndexMapInit
{
public:
  //! @brief Initialization function, override with initialization code
  virtual void init(T&) = 0;

  //! @brief Cloning a pointer, override with code to create a fresh copy
  virtual IndexMapInit<T>* clone() const = 0;
  virtual ~IndexMapInit() {} // ensure that derived destructor is called
};

//! @brief IndexMap<T> implements a generic map indexed by a dynamic index set.
//!
//! Additionally, it allows new elements of the map to be initialized in a
//! flexible manner.
template <typename T>
class IndexMap
{

  std::unordered_map<long, T> map;

  IndexSet indexSet;
  cloned_ptr<IndexMapInit<T>> init;

public:
  //! @brief The empty map
  IndexMap();

  //! @brief A map with an initialization object.
  //! This associates a method for initializing new elements in the map.
  //! When a new index j is added to the index set, an object t of type T is
  //! created using the default constructor for T, after which the function
  //! _init->init(t) is called (t is passed by reference). To use this
  //! feature, you need to derive a subclass of IndexMapInit<T> that defines
  //! the init function. This "helper object" should be created using
  //! operator new, and the pointer is "exclusively owned" by the map object.
  explicit IndexMap(IndexMapInit<T>* _init) : init(_init) {}

  //! @brief Get the underlying index set
  const IndexSet& getIndexSet() const { return indexSet; }

  //! @brief Access functions: will raise an error
  //! if j does not belong to the current index set
  T& operator[](long j)
  {
    assertTrue(indexSet.contains(j), "Key not found");
    return map[j];
  }
  const T& operator[](long j) const
  {
    assertTrue(indexSet.contains(j), "Key not found");
    // unordered_map does not support a const [] operator,
    // so we have to artificially strip away the const-ness here
    std::unordered_map<long, T>& map1 =
        const_cast<std::unordered_map<long, T>&>(map);
    return map1[j];
  }

  //! @brief Insert indexes to the IndexSet.
  //! Insertion will cause new T objects to be created, using the default
  //! constructor, and possibly initialized via the IndexMapInit<T> pointer.
  void insert(long j)
  {
    if (!indexSet.contains(j)) {
      indexSet.insert(j);
      if (!init.null())
        init->init(map[j]);
    }
  }
  void insert(const IndexSet& s)
  {
    for (long i = s.first(); i <= s.last(); i = s.next(i))
      insert(i);
  }

  //! @brief Delete indexes from IndexSet, may cause objects to be destroyed.
  void remove(long j)
  {
    indexSet.remove(j);
    map.erase(j);
  }
  void remove(const IndexSet& s)
  {
    for (long i = s.first(); i <= s.last(); i = s.next(i))
      map.erase(i);
    indexSet.remove(s);
  }

  void clear()
  {
    map.clear();
    indexSet.clear();
  }
};

//! @brief Comparing maps, by comparing all the elements
template <typename T>
bool operator==(const IndexMap<T>& map1, const IndexMap<T>& map2)
{
  if (map1.getIndexSet() != map2.getIndexSet())
    return false;
  const IndexSet& s = map1.getIndexSet();
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    if (map1[i] == map2[i])
      continue;
    return false;
  }
  return true;
}

template <typename T>
bool operator!=(const IndexMap<T>& map1, const IndexMap<T>& map2)
{
  return !(map1 == map2);
}

} // namespace helib

#endif // ifndef HELIB_INDEXMAP_H
