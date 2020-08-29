/* Copyright (C) 2019-2020 IBM Corp.
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

#include <helib/PolyMod.h>
#include <helib/exceptions.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <vector>
#include <helib/NumbTh.h>

namespace helib {

PolyMod::PolyMod() : ringDescriptor(nullptr) {}
PolyMod::PolyMod(const std::shared_ptr<PolyModRing>& ringDescriptor) :
    PolyMod(NTL::ZZX(0), ringDescriptor)
{}
PolyMod::PolyMod(long input,
                 const std::shared_ptr<PolyModRing>& ringDescriptor) :
    PolyMod(NTL::ZZX(input), ringDescriptor)
{}
PolyMod::PolyMod(const std::vector<long>& input,
                 const std::shared_ptr<PolyModRing>& ringDescriptor) :
    PolyMod(ringDescriptor)
{
  *this = input;
}
PolyMod::PolyMod(const NTL::ZZX& input,
                 const std::shared_ptr<PolyModRing>& ringDescriptor) :
    ringDescriptor(ringDescriptor), data(input)
{
  this->modularReduce();
}

PolyMod& PolyMod::operator=(long input)
{
  assertValidity(*this);
  this->data = NTL::ZZX(input);
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator=(const std::vector<long>& input)
{
  assertValidity(*this);
  NTL::clear(data); // Make sure higher-degree terms don't remain
  for (std::size_t i = 0; i < input.size(); ++i)
    NTL::SetCoeff(data, i, input[i]);
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator=(const std::initializer_list<long>& input)
{
  assertValidity(*this);
  *this = std::vector<long>(input);
  return *this;
}

PolyMod& PolyMod::operator=(const NTL::ZZX& input)
{
  assertValidity(*this);
  this->data = input;
  this->modularReduce();
  return *this;
}

PolyMod::operator long() const
{
  assertValidity(*this);
  long ret;
  NTL::conv(ret, NTL::ConstTerm(this->data));
  return ret;
}

PolyMod::operator std::vector<long>() const
{
  assertValidity(*this);
  std::vector<long> ret(NTL::deg(ringDescriptor->G));
  for (std::size_t i = 0; i < ret.size(); ++i)
    NTL::conv(ret[i], NTL::coeff(data, i));
  return ret;
}

PolyMod::operator NTL::ZZX() const
{
  assertValidity(*this);
  return getData();
}

bool PolyMod::isValid() const { return ringDescriptor != nullptr; }

long PolyMod::getp2r() const { return ringDescriptor->p2r; }

NTL::ZZX PolyMod::getG() const { return ringDescriptor->G; }

const NTL::ZZX& PolyMod::getData() const
{
  assertValidity(*this);
  return this->data;
}

bool PolyMod::operator==(const PolyMod& rhs) const
{
  if (!isValid() && !rhs.isValid())
    return true;
  else
    return isValid() && rhs.isValid() &&
           *ringDescriptor == *(rhs.ringDescriptor) && data == rhs.data;
}

bool PolyMod::operator==(long rhs) const { return *this == NTL::ZZX(rhs); }

bool PolyMod::operator==(const std::vector<long>& rhs) const
{
  if (!this->isValid()) {
    return false;
  } else {
    PolyMod other(rhs, ringDescriptor);
    return *this == other;
  }
}

bool PolyMod::operator==(const NTL::ZZX& rhs) const
{
  if (!this->isValid()) {
    return false;
  } else {
    PolyMod copy(*this);
    // Using subtraction to ensure modularReduce is called.
    // We are checking for divisibility of difference by G in Z_p rather than
    // direct equality
    copy -= rhs;
    return copy.data == NTL::ZZX(0);
  }
}

PolyMod& PolyMod::negate()
{
  assertValidity(*this);
  *this *= -1;
  return *this;
}

PolyMod PolyMod::operator-() const
{
  assertValidity(*this);
  PolyMod poly(*this);
  poly.negate();
  return poly;
}

PolyMod PolyMod::operator*(const PolyMod& rhs) const
{
  assertInterop(*this, rhs);
  PolyMod result = *this;
  result *= rhs;
  return result;
}

PolyMod PolyMod::operator*(long rhs) const { return operator*(NTL::ZZX{rhs}); }

PolyMod PolyMod::operator*(const NTL::ZZX& rhs) const
{
  PolyMod result(*this);
  PolyMod multiplier(*this);
  multiplier = rhs;
  return result * multiplier;
}

PolyMod PolyMod::operator+(const PolyMod& rhs) const
{
  assertInterop(*this, rhs);
  PolyMod result(*this);
  result += rhs;
  return result;
}

PolyMod PolyMod::operator+(long rhs) const { return operator+(NTL::ZZX{rhs}); }

PolyMod PolyMod::operator+(const NTL::ZZX& rhs) const
{
  PolyMod result(*this);
  PolyMod addend(result);
  addend = rhs;
  return result + addend;
}

PolyMod PolyMod::operator-(const PolyMod& rhs) const
{
  assertInterop(*this, rhs);
  PolyMod result = *this;
  result -= rhs;
  return result;
}

PolyMod PolyMod::operator-(long rhs) const { return operator-(NTL::ZZX{rhs}); }

PolyMod PolyMod::operator-(const NTL::ZZX& rhs) const
{
  PolyMod result(*this);
  PolyMod subtrahend(result);
  subtrahend = rhs;
  return result - subtrahend;
}

PolyMod& PolyMod::operator*=(const PolyMod& otherPoly)
{
  assertInterop(*this, otherPoly);
  this->data *= otherPoly.data;
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator*=(long scalar)
{
  assertValidity(*this);
  this->data *= NTL::ZZX(scalar);
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator*=(const NTL::ZZX& otherPoly)
{
  assertValidity(*this);
  this->data *= otherPoly;
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator+=(const PolyMod& otherPoly)
{
  assertInterop(*this, otherPoly);
  this->data += otherPoly.data;
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator+=(long scalar)
{
  assertValidity(*this);
  this->data += NTL::ZZX(scalar);
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator+=(const NTL::ZZX& otherPoly)
{
  assertValidity(*this);
  this->data += otherPoly;
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator-=(const PolyMod& otherPoly)
{
  assertInterop(*this, otherPoly);
  this->data -= otherPoly.data;
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator-=(long scalar)
{
  assertValidity(*this);
  this->data -= NTL::ZZX(scalar);
  this->modularReduce();
  return *this;
}

PolyMod& PolyMod::operator-=(const NTL::ZZX& otherPoly)
{
  assertValidity(*this);
  this->data -= otherPoly;
  this->modularReduce();
  return *this;
}

void deserialize(std::istream& is, PolyMod& poly)
{
  PolyMod::assertValidity(poly);
  std::vector<std::stringstream> parts =
      extractTokenizeRegion(is, '[', ']', ',');

  long degree = NTL::deg(poly.ringDescriptor->G);
  if (lsize(parts) > degree) {
    // Too many elements. Raising an error.
    std::stringstream err_msg;
    err_msg << "Cannot deserialize to PolyMod: Degree is too small.  "
            << "Trying to deserialize " << parts.size() << " coefficients.  "
            << "Degree is " << degree << ".";
    throw IOError(err_msg.str());
  }

  // Actual parsing and setup
  NTL::clear(poly.data); // Make sure higher-degree terms don't remain
  for (std::size_t i = 0; i < parts.size(); ++i) {
    long tmp;
    parts[i] >> tmp;
    NTL::SetCoeff(poly.data, i, tmp);
  }

  // Normalization (removal of leading zeros) is done by modularReduce.
  poly.modularReduce();
}

void serialize(std::ostream& os, const PolyMod& poly)
{
  PolyMod::assertValidity(poly);
  if (poly.data == NTL::ZZX::zero()) {
    // Avoid string "[]" for zero ZZX.
    os << "[0]";
    return;
  }

  // TODO: Add stream modifier option for separator
  std::string sep = ", ";
  os << "[";
  for (auto ite = poly.data.rep.begin(); ite != poly.data.rep.end(); ++ite) {
    os << *ite;
    if (ite + 1 != poly.data.rep.end()) {
      os << ", ";
    }
  }
  os << "]";
}

std::istream& operator>>(std::istream& is, PolyMod& poly)
{
  PolyMod::assertValidity(poly);

  deserialize(is, poly);
  return is;
}

std::ostream& operator<<(std::ostream& os, const PolyMod& poly)
{
  PolyMod::assertValidity(poly);

  serialize(os, poly);
  return os;
}

void PolyMod::modularReduce()
{
  NTL::ZZ_pContext pContext;
  pContext.save();
  NTL::ZZ_p::init(NTL::ZZ(ringDescriptor->p2r));
  NTL::ZZ_pX poly_mod_p2r;
  NTL::conv(poly_mod_p2r, this->data);
  NTL::ZZ_pX G_mod_p2r;
  NTL::conv(G_mod_p2r, ringDescriptor->G);
  poly_mod_p2r %= G_mod_p2r;
  NTL::conv(this->data, poly_mod_p2r);
  pContext.restore();
  this->data.normalize();
}

void PolyMod::assertValidity(const PolyMod& poly)
{
  if (!poly.isValid()) {
    throw LogicError("Cannot operate on invalid (default constructed) PolyMod");
  }
}

void PolyMod::assertInterop(const PolyMod& lhs, const PolyMod& rhs)
{
  assertValidity(lhs);
  assertValidity(rhs);
  if (*(lhs.ringDescriptor) != *(rhs.ringDescriptor))
    throw LogicError("Ring descriptors are not equal between PolyMod objects");
}

} // namespace helib
