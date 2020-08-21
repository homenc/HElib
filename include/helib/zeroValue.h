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
#ifndef HELIB_ZEROVALUE_H
#define HELIB_ZEROVALUE_H

#include <helib/Ctxt.h>
#include <helib/Ptxt.h>

// FIXME: There must be a better way to get a zero object.

namespace helib {

/**
 * @brief Given an object `x` return a zero object of the same type.
 * @tparam T Type of object to return.
 * @param x The object to use for returning a zero object of type `T`.
 * @return A zero object of type `T`.
 **/
template <typename T>
inline T zeroValue(const T& x)
{
  // x is unused by design, so cast to avoid error
  static_cast<void>(x);
  // type T must be able to convert int
  return T(0);
}

/**
 * @brief Given a `Ctxt` return a zero object of the same type.
 * @param x The `Ctxt` to use for returning a zero object of type `Ctxt`.
 * @return A zero object of type `Ctxt`.
 **/
template <>
inline Ctxt zeroValue<Ctxt>(const Ctxt& x)
{
  return Ctxt(ZeroCtxtLike, x);
}

/**
 * @brief Given a `Ptxt<BGV>` return a zero object of the same type.
 * @param x The object to use for returning a zero object of type `Ptxt<BGV>`.
 * @return A zero object of type `Ptxt<BGV>`.
 **/
template <>
inline Ptxt<BGV> zeroValue<Ptxt<BGV>>(const Ptxt<BGV>& x)
{
  return Ptxt<BGV>(x.getContext());
}

/**
 * @brief Given a `Ptxt<CKKS>` return a zero object of the same type.
 * @param x The object to use for returning a zero object of type `Ptxt<CKKS>`.
 * @return A zero object of type `Ptxt<CKKS>`.
 **/
template <>
inline Ptxt<CKKS> zeroValue<Ptxt<CKKS>>(const Ptxt<CKKS>& x)
{
  return Ptxt<CKKS>(x.getContext());
}

} // namespace helib

#endif // HELIB_ZEROVALUE_H
