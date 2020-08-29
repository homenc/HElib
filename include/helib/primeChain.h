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
#ifndef HELIB_PRIMECHAIN_H
#define HELIB_PRIMECHAIN_H
/**
 * @file primeChain.h
 * @brief handling the chain of moduli
 */
#include <vector>
#include <helib/IndexSet.h>

namespace helib {

class Cmodulus;
class Context;

//! A helper class to map required modulo-sizes to primeSets
class ModuliSizes
{
public:
  typedef std::pair<double, IndexSet> Entry;
  // each table entry is a pair<double,IndexSet>=(size, set-of-primes)

  //! initialize helper table for a given chain
  void init(const Context& context);

  //! Find a suitable IndexSet of primes whose total size is in the
  //! target interval [low,high], trying to minimize the number of
  //! primes dropped from fromSet.
  //! If no IndexSet exists that fits in the target interval, returns
  //! the IndexSet that gives the largest value smaller than low
  //! (or the smallest value greater than low if reverse flag is set).
  IndexSet getSet4Size(double low,
                       double high,
                       const IndexSet& fromSet,
                       bool reverse) const;

  //! Find a suitable IndexSet of primes whose total size is in the
  //! target interval [low,high], trying to minimize the total number
  //! of primes dropped from both from1, from2.
  //! If no IndexSet exists that fits in the target interval, returns
  //! the IndexSet that gives the largest value smaller than low.
  //! (or the smallest value greater than low if reverse flag is set).
  IndexSet getSet4Size(double low,
                       double high,
                       const IndexSet& from1,
                       const IndexSet& from2,
                       bool reverse) const;

  // ASCII I/O
  friend std::istream& operator>>(std::istream& s, ModuliSizes& szs);
  friend std::ostream& operator<<(std::ostream& s, const ModuliSizes& szs);
  // Raw I/O
  void read(std::istream& str);
  void write(std::ostream& str) const;

private:
  std::vector<Entry> sizes;
  long iFFT_cost = -1;
};

std::ostream& operator<<(std::ostream& s, const ModuliSizes::Entry& e);
std::istream& operator>>(std::istream& s, ModuliSizes::Entry& e);
void write(std::ostream& s, const ModuliSizes::Entry& e);
void read(std::istream& s, ModuliSizes::Entry& e);

} // namespace helib

#endif // ifndef HELIB_PRIMECHAIN_H
