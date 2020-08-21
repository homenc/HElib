
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

#include <helib/PolyModRing.h>

namespace helib {

PolyModRing::PolyModRing(long p, long r, const NTL::ZZX& G) :
    p(p), r(r), G(G), p2r(pow(p, r))
{}

bool PolyModRing::operator==(const PolyModRing& rhs) const noexcept
{
  return p == rhs.p && r == rhs.r && G == rhs.G && p2r == rhs.p2r;
}

bool PolyModRing::operator!=(const PolyModRing& rhs) const noexcept
{
  return !operator==(rhs);
}

std::ostream& operator<<(std::ostream& os, const PolyModRing& ring)
{
  os << "(p: " << ring.p << ", r: " << ring.r << ", p^r: " << ring.p2r
     << ", G: " << ring.G << ")";
  return os;
}

} // namespace helib
