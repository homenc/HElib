/* Copyright (C) 2019 IBM Corp.
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

#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <vector>

#ifndef HELIB_POLYMODRING_H
#define HELIB_POLYMODRING_H

/**
 * @file PolyModRing.h
 * @brief Definition of the plaintext slot algebraic ring.
 **/
namespace helib {

/**
 * @struct PolyModRing
 * @brief Lightweight type for describing the structure of a single slot of the
 * plaintext space.
 *
 * A single slot of the plaintext space is isomorphic to
 * \f$\mathbb{Z}[X]/(G(x),p^r)\f$ for some irreducible factor G of
 * \f$\Phi_m(X)\f$, so the main useful members of this `struct` are `p`, `r`,
 * `G`, and `p2r`.
 *
 * The fields of this `struct` are all `const`, so they should be determined
 * at the time of construction.
 *
 * @note This `struct` aggregates this often-useful information into a single
 * placeholder for convenience.
 **/
struct PolyModRing
{
  /**
   * @brief The characteristic of the plaintext space.  This should be prime.
   **/
  const long p;
  /**
   * @brief The power of p used in the plaintext space coefficient modulus.
   **/
  const long r;
  /**
   * @brief The irreducible factor of Phi_m(X) used for the algebra of the
   * individual slots.
   **/
  const NTL::ZZX G;
  /**
   * @brief The plaintext space coefficient modulus, equal to p^r.
   **/
  const long p2r;

  // Delete the default constructor.
  PolyModRing() = delete;

  // Delete the assignment operators.
  PolyModRing& operator=(const PolyModRing&) = delete;
  PolyModRing& operator=(PolyModRing&&) = delete;

  /**
   * @brief Copy constructor.
   **/
  PolyModRing(const PolyModRing& other) = default;

  /**
   * @brief Move constructor.
   **/
  PolyModRing(PolyModRing&& other) = default;

  /**
   * @brief Destructor.
   **/
  ~PolyModRing() = default;

  /**
   * @brief Constructor.
   * @param p The characteristic of the plaintext space.
   * @param r The power of p used in the plaintext space coefficient modulus.
   * @param G The irreducible factor of Phi_m(X) used for the algebra of the
   * individual slots.
   *
   * p^r will be calculated automatically.
   *
   * @note p should be a prime number.
   **/
  PolyModRing(long p, long r, const NTL::ZZX& G);

  /**
   * @brief Equality operator.
   **/
  bool operator==(const PolyModRing& rhs) const noexcept;

  /**
   * @brief Not-equals operator.
   **/
  bool operator!=(const PolyModRing& rhs) const noexcept;

  /**
   * @brief Output shift operator.
   * @param os Output `std::ostream`.
   * @param ring `PolyModRing` object to be written.
   * @return Input `std::ostream` post writing.
   **/
  friend std::ostream& operator<<(std::ostream& os, const PolyModRing& ring);
};
} // namespace helib

#endif // HELIB_POLYMODRING_H
