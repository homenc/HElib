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

#include <random>
#include <helib/Ptxt.h>
#include <helib/apiAttributes.h>

#include "io.h"

namespace helib {

template <>
PolyMod Ptxt<BGV>::convertToSlot(const Context& context, long slot)
{
  PolyMod data(NTL::ZZX(slot), context.getSlotRing());
  return data;
}

template <>
std::complex<double> Ptxt<CKKS>::convertToSlot(const Context&, long slot)
{
  return {static_cast<double>(slot), 0};
}

template <typename Scheme>
Ptxt<Scheme>::Ptxt() : context(nullptr)
{}

template <typename Scheme>
Ptxt<Scheme>::Ptxt(const Context& context) :
    context(&context),
    slots(context.getEA().size(),
          SlotType{Ptxt<Scheme>::convertToSlot(*(this->context), 0L)})
{}

template <typename Scheme>
Ptxt<Scheme>::Ptxt(const Context& context, const SlotType& value) :
    context(std::addressof(context)),
    slots(context.getEA().size(),
          SlotType{Ptxt<Scheme>::convertToSlot(*(this->context), 0L)})
{
  setData(value);
}

template <>
template <>
void Ptxt<BGV>::setData(const NTL::ZZX& value)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call setData on default-constructed Ptxt");
  PolyMod poly(value, context->getSlotRing());
  std::vector<PolyMod> poly_vec(context->getEA().size(), poly);
  setData(poly_vec);
}

template <>
template <>
Ptxt<BGV>::Ptxt(const Context& context, const NTL::ZZX& value) :
    context(&context),
    slots(context.getEA().size(),
          SlotType{Ptxt<BGV>::convertToSlot(*(this->context), 0L)})
{
  setData(value);
}

template <typename Scheme>
Ptxt<Scheme>::Ptxt(const Context& context, const std::vector<SlotType>& data) :
    context(std::addressof(context)),
    slots(context.getEA().size(),
          SlotType{Ptxt<Scheme>::convertToSlot(*(this->context), 0L)})
{
  setData(data);
}

template <typename Scheme>
bool Ptxt<Scheme>::isValid() const
{
  return context != nullptr;
}

template <typename Scheme>
size_t Ptxt<Scheme>::size() const
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call size on default-constructed Ptxt");
  return this->slots.size();
}

template <typename Scheme>
long Ptxt<Scheme>::lsize() const
{
  assertTrue<LogicError>(isValid(),
                         "Cannot call lsize on default-constructed Ptxt");
  return this->size();
}

template <typename Scheme>
void Ptxt<Scheme>::setData(const std::vector<SlotType>& data)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call setData on default-constructed Ptxt");
  assertTrue<RuntimeError>(helib::lsize(data) <= context->getEA().size(),
                           "Cannot setData to Ptxt: not enough slots");

  // Need to verify that they all match
  assertSlotsCompatible(data);

  slots = data;
  if (helib::lsize(slots) < context->getEA().size()) {
    slots.resize(context->getEA().size(),
                 SlotType{Ptxt<Scheme>::convertToSlot(*(this->context), 0L)});
  }
}

template <typename Scheme>
void Ptxt<Scheme>::setData(const SlotType& value)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call setData on default-constructed Ptxt");
  setData(std::vector<SlotType>(context->getEA().size(), value));
}

template <>
template <>
void Ptxt<BGV>::decodeSetData(const NTL::ZZX& data)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call decodeSetData on default-constructed Ptxt");
  PolyMod poly(context->getSlotRing());
  std::vector<PolyMod> poly_vec(context->getEA().size(), poly);
  std::vector<NTL::ZZX> ptxt(context->getEA().size());
  context->getEA().decode(ptxt, data);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    poly_vec[i] = ptxt[i];
  }
  setData(poly_vec);
}

template <typename Scheme>
void Ptxt<Scheme>::clear()
{
  for (auto& slot : this->slots) {
    slot = 0;
  }
}

template <typename Scheme>
typename Scheme::SlotType randomSlot(const Context& context);

template <>
BGV::SlotType randomSlot<BGV>(const Context& context)
{
  std::vector<long> coeffs(context.getOrdP());
  NTL::VectorRandomBnd(coeffs.size(),
                       coeffs.data(),
                       context.getSlotRing()->p2r);
  return PolyMod(coeffs, context.getSlotRing());
}

template <>
CKKS::SlotType randomSlot<CKKS>(UNUSED const Context& context)
{
  // VJS-FIXME: see new function RandomComplex in NumbTh.h
  // for an alternative
  std::mt19937 gen{std::random_device{}()};
  std::uniform_real_distribution<> dist{-1e10, 1e10};

  return {dist(gen), dist(gen)};
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::random()
{
  for (auto& slot : slots)
    slot = randomSlot<Scheme>(*context);
  return *this;
}

template <typename Scheme>
const std::vector<typename Ptxt<Scheme>::SlotType>& Ptxt<Scheme>::getSlotRepr()
    const
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call getSlotRepr on default-constructed Ptxt");
  return slots;
}

/**
 * @brief BGV specialisation of the `getPolyRepr` function.
 * @return Single encoded polynomial.
 * @note Only enabled for the `BGV` scheme.
 **/

template <>
NTL::ZZX Ptxt<BGV>::getPolyRepr() const
{
  assertTrue<LogicError>(isValid(),
                         "Cannot call getPolyRepr on default-constructed Ptxt");
  NTL::ZZX repr;
  std::vector<NTL::ZZX> slots_data(context->getEA().size());
  for (std::size_t i = 0; i < slots_data.size(); ++i) {
    slots_data[i] = slots[i].getData();
  }
  context->getEA().encode(repr, slots_data);
  return repr;
}

template <>
void Ptxt<BGV>::encode(EncodedPtxt& eptxt, double mag, OptLong prec) const
{
  assertTrue<LogicError>(isValid(),
                         "Cannot call encide on default-constructed Ptxt");

  assertTrue<LogicError>(mag < 0 && !prec.isDefined(),
                         "mag,prec must be defaulted for BGV");

  std::vector<NTL::ZZX> slots_data(context->getEA().size());
  for (std::size_t i = 0; i < slots_data.size(); ++i) {
    slots_data[i] = slots[i].getData();
  }

  context->getEA().encode(eptxt, slots_data);
}

/**
 * @brief CKKS specialisation of the `getPolyRepr` function.
 * @return Single encoded polynomial.
 * @note Only enabled for the `CKKS` scheme.
 **/
template <>
void Ptxt<CKKS>::encode(EncodedPtxt& eptxt, double mag, OptLong prec) const
{
  assertTrue<LogicError>(isValid(),
                         "Cannot call encode on default-constructed Ptxt");
  context->getEA().encode(eptxt, slots, mag, prec);
}

template <typename Scheme>
typename Ptxt<Scheme>::SlotType& Ptxt<Scheme>::operator[](long i)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot access elements of default-constructed Ptxt");
  return this->slots[i];
}

template <typename Scheme>
typename Ptxt<Scheme>::SlotType Ptxt<Scheme>::operator[](long i) const
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot access elements of default-constructed Ptxt");
  return this->slots[i];
}

template <typename Scheme>
typename Ptxt<Scheme>::SlotType& Ptxt<Scheme>::at(long i)
{
  assertInRange(i, 0l, lsize(), "Index out of range");
  return operator[](i);
}

template <typename Scheme>
typename Ptxt<Scheme>::SlotType Ptxt<Scheme>::at(long i) const
{
  assertInRange(i, 0l, lsize(), "Index out of range");
  return operator[](i);
}

template <typename Scheme>
bool Ptxt<Scheme>::operator==(const Ptxt<Scheme>& other) const
{
  return (!isValid() && !other.isValid()) ||
         (this->slots == other.slots && *(this->context) == *(other.context));
}

template <typename Scheme>
bool Ptxt<Scheme>::operator!=(const Ptxt<Scheme>& other) const
{
  return !(*this == other);
}

template <typename Scheme>
Ptxt<Scheme> Ptxt<Scheme>::operator*(const Ptxt<Scheme>& rhs) const
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call operator* on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      rhs.isValid(),
      "Cannot call operator* with a default-constructed Ptxt as the right "
      "operand");
  assertEq<LogicError>(*context,
                       *(rhs.context),
                       "Ptxts must have matching contexts");
  Ptxt<Scheme> result = *this;
  result *= rhs;
  return result;
}

template <typename Scheme>
Ptxt<Scheme> Ptxt<Scheme>::operator+(const Ptxt<Scheme>& rhs) const
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call operator+ on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      rhs.isValid(),
      "Cannot call operator+ with a default-constructed Ptxt as the right "
      "operand");
  assertEq<LogicError>(*context,
                       *(rhs.context),
                       "Ptxts must have matching contexts");
  Ptxt<Scheme> result = *this;
  result += rhs;
  return result;
}

template <typename Scheme>
Ptxt<Scheme> Ptxt<Scheme>::operator-(const Ptxt<Scheme>& rhs) const
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call operator- on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      rhs.isValid(),
      "Cannot call operator- with a default-constructed Ptxt as the right "
      "operand");
  assertEq<LogicError>(*context,
                       *(rhs.context),
                       "Ptxts must have matching contexts");
  Ptxt<Scheme> result = *this;
  result -= rhs;
  return result;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::operator*=(const Ptxt<Scheme>& otherPtxt)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call operator*= on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      otherPtxt.isValid(),
      "Cannot call operator*= with a default-constructed Ptxt as the right "
      "operand");
  assertEq<LogicError>(*context,
                       *(otherPtxt.context),
                       "Ptxts must have matching contexts");
  for (unsigned i = 0; i < this->slots.size(); i++) {
    this->slots[i] *= otherPtxt.slots[i];
  }
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::operator*=(const SlotType& scalar)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call operator*= on default-constructed Ptxt");
  for (auto& x : this->slots)
    x *= scalar;
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::operator+=(const Ptxt<Scheme>& otherPtxt)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call operator+= on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      otherPtxt.isValid(),
      "Cannot call operator+= with a default-constructed Ptxt as the right "
      "operand");
  assertEq<LogicError>(*context,
                       *(otherPtxt.context),
                       "Ptxts must have matching contexts");
  for (unsigned i = 0; i < this->slots.size(); i++) {
    this->slots[i] += otherPtxt.slots[i];
  }
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::operator+=(const SlotType& scalar)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call operator+= on default-constructed Ptxt");
  for (auto& x : this->slots)
    x += scalar;
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::operator-=(const Ptxt<Scheme>& otherPtxt)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call operator-= on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      otherPtxt.isValid(),
      "Cannot call operator-= with a default-constructed Ptxt as the right "
      "operand");
  assertEq<LogicError>(*context,
                       *(otherPtxt.context),
                       "Ptxts must have matching contexts");
  for (unsigned i = 0; i < this->slots.size(); i++) {
    this->slots[i] -= otherPtxt.slots[i];
  }
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::operator-=(const SlotType& scalar)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call operator-= on default-constructed Ptxt");
  for (auto& x : this->slots)
    x -= scalar;
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::negate()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call negate on default-constructed Ptxt");
  for (auto& slot : slots) {
    slot = -slot;
  }
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::multiplyBy(const Ptxt<Scheme>& otherPtxt)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call multiplyBy on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      otherPtxt.isValid(),
      "Cannot call multiplyBy with default-constructed Ptxt as argument");
  assertEq<LogicError>(*context,
                       *(otherPtxt.context),
                       "Ptxts must have matching contexts");
  if (size() != otherPtxt.size()) {
    throw RuntimeError("Cannot multiply by plaintext of different size");
  }
  return *this *= otherPtxt;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::multiplyBy2(const Ptxt& otherPtxt1,
                                        const Ptxt& otherPtxt2)
{
  assertTrue<RuntimeError>(
      isValid(),
      "Cannot call multiplyBy2 on default-constructed Ptxt");
  assertTrue<RuntimeError>(
      otherPtxt1.isValid(),
      "Cannot call multiplyBy2 with default-constructed Ptxt as first "
      "argument");
  assertTrue<RuntimeError>(
      otherPtxt2.isValid(),
      "Cannot call multiplyBy2 with default-constructed Ptxt as second "
      "argument");
  assertEq<LogicError>(*context,
                       *(otherPtxt1.context),
                       "Ptxts must have matching contexts");
  assertEq<LogicError>(*context,
                       *(otherPtxt2.context),
                       "Ptxts must have matching contexts");
  assertEq<RuntimeError>(size(),
                         otherPtxt1.size(),
                         "Cannot multiply by plaintext of different size - "
                         "first argument has wrong size");
  assertEq<RuntimeError>(size(),
                         otherPtxt2.size(),
                         "Cannot multiply by plaintext of different size - "
                         "second argument has wrong size");
  for (std::size_t i = 0; i < size(); ++i)
    slots[i] *= otherPtxt1.slots[i] * otherPtxt2.slots[i];

  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::square()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call square on default-constructed Ptxt");
  return multiplyBy(*this);
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::cube()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call cube on default-constructed Ptxt");
  return multiplyBy2(*this, *this);
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::power(long e)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call power on default-constructed Ptxt");
  if (e < 1) {
    throw InvalidArgument("Cannot raise a Ptxt to a non positive "
                          "exponent");
  } else if (e > 1) {
    // exponentiation through squaring.
    std::vector<SlotType> multiplier(slots);
    std::vector<SlotType> result(
        multiplier.size(),
        Ptxt<Scheme>::convertToSlot(*(this->context), 1l));
    for (; e; e >>= 1u) {
      if (e & 1u)
        for (std::size_t i = 0; i < size(); ++i)
          result[i] *= multiplier[i];
      for (auto& slot : multiplier)
        slot *= slot;
    }
    slots = std::move(result);
  }
  // If e == 1 do nothing
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::rotate(long amount)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call rotate on default-constructed Ptxt");
  amount = mcMod(amount, size());
  if (amount == 0)
    return *this;
  std::vector<SlotType> rotated_slots(size());
  for (long i = 0; i < lsize(); ++i) {
    rotated_slots[i] = slots[mcMod(i - amount, size())];
  }
  slots = std::move(rotated_slots);
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::rotate1D(long dim, long amount)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call rotate1D on default-constructed Ptxt");
  if (slots.size() == 1)
    return *this; // Nothing to do (only one slot)
  const PAlgebra& zMStar = context->getZMStar();
  long num_gens = zMStar.numOfGens();
  assertInRange<LogicError>(dim,
                            0l,
                            num_gens,
                            "Dimension must be between 0 and "
                            "number of generators");
  // Copying in slots to avoid default PolyMod issues.
  std::vector<SlotType> new_slots(slots);
  long ord = context->getEA().sizeOfDimension(dim);
  amount = mcMod(amount, ord); // Make amount smallest positive integer < ord
  if (amount == 0)
    return *this; // Nothing to do

  // This for loop iterates over the slots of a flat array structure and
  // converts the index of each slot to its equivalent coordinate
  // representation with respect to the generators of the quotient group
  // Zm*/<p> where coordinate i represents the i'th generator.
  // After the conversion the relevant generator (specified by dim) is
  // incremented by amount.  This new set of coordinates is then converted back
  // to the new index of the slot.
  for (long index = 0; index < lsize(); ++index) {
    // Vector to hold the coordinate representation of the current index.
    std::vector<long> coord(indexToCoord(index));
    // Increments the coordinate of the specific dimension (dim) by amount
    // and reduces it modulo the order of the dimension.
    coord[dim] = mcMod(amount + coord[dim], ord);
    // Convert the new coordinates post rotation into the correct index.
    long new_index = coordToIndex(coord);
    // Set the new index of the current slot.
    new_slots[new_index] = slots[index];
  }
  slots = std::move(new_slots);
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::shift(long amount)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call shift on default-constructed Ptxt");
  if (amount == 0)
    return *this;
  if (std::abs(amount) >= lsize()) {
    clear();
    return *this;
  }

  rotate(amount);
  for (long i = 0; i < lsize(); ++i)
    if ((i - amount) < 0 || (i - amount) >= lsize())
      slots[i] = 0;

  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::shift1D(long dim, long amount)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call shift1D on default-constructed Ptxt");
  if (amount == 0)
    return *this;
  if (slots.size() == 1 ||
      std::abs(amount) >= context->getEA().sizeOfDimension(dim)) {
    clear();
    return *this;
  }
  // NOTE: There is some code duplication here and in rotate1D
  const PAlgebra& zMStar = context->getZMStar();
  long num_gens = zMStar.numOfGens();
  assertInRange<LogicError>(dim,
                            0l,
                            num_gens,
                            "Dimension must be between 0 and "
                            "number of generators");
  // Copying in slots to avoid default PolyMod issues.
  std::vector<SlotType> new_slots(slots);
  long ord = context->getEA().sizeOfDimension(dim);

  // This for loop performs similar logic in rotate1D to obtain the new index
  // post shift via the conversion to and from the corresponding coordinate
  // representation.
  // An extra check is then performed to see if the shift operation caused an
  // element to wrap around in which case it is replaced with 0.
  for (long new_index = 0; new_index < lsize(); ++new_index) {
    // Vector to hold the coordinate representation of the current index.
    std::vector<long> coord(indexToCoord(new_index));
    // Perform the shift backwards to obtain the locations of the 0 values
    // first.
    coord[dim] -= amount;
    // If the coordinate exceeds the bounds of the order of the generator then
    // set the new slot to 0.
    if (coord[dim] < 0 || coord[dim] >= ord) {
      new_slots[new_index] = 0;
    } else { // Otherwise set the old index to the value of the new index.
      long old_index = coordToIndex(coord);
      new_slots[new_index] = slots[old_index];
    }
  }
  slots = std::move(new_slots);
  return *this;
}

template <>
Ptxt<BGV>& Ptxt<BGV>::automorph(long k)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call automorph on default-constructed Ptxt");
  assertTrue<RuntimeError>(context->getZMStar().inZmStar(k),
                           "k must be an element in Zm*");
  NTL::ZZX poly;
  switch (context->getEA().getTag()) {
  case PA_GF2_tag: {
    decodeSetData(automorph_internal<PA_GF2>(k));
    break;
  }
  case PA_zz_p_tag: {
    decodeSetData(automorph_internal<PA_zz_p>(k));
    break;
  }
  default:
    throw LogicError("Could not find valid tag for BGV automorphism");
  }
  return *this;
}

template <>
Ptxt<CKKS>& Ptxt<CKKS>::automorph(long k)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call automorph on default-constructed Ptxt");
  assertTrue<RuntimeError>(context->getZMStar().inZmStar(k),
                           "k must be an element in Zm*");
  return rotate(context->getZMStar().indexOfRep(k));
}

template <>
template <>
Ptxt<BGV>& Ptxt<BGV>::frobeniusAutomorph(long j)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call frobeniusAutomorph on "
                           "default-constructed Ptxt");
  long d = context->getOrdP();
  if (d == 1)
    return *this; // Nothing to do.
  long exponent = NTL::PowerMod(context->getSlotRing()->p,
                                mcMod(j, d),
                                context->getZMStar().getM());
  return automorph(exponent);
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::replicate(long pos)
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call replicate on default-constructed Ptxt");
  for (auto& slot : slots)
    slot = slots[pos];
  return *this;
}

template <typename Scheme>
std::vector<Ptxt<Scheme>> Ptxt<Scheme>::replicateAll() const
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call replicateAll on "
                           "default-constructed Ptxt");
  std::vector<Ptxt<Scheme>> ret(size(), *this);
  for (std::size_t i = 0; i < size(); ++i) {
    ret[i].replicate(i);
  }
  return ret;
}

template <>
template <>
Ptxt<CKKS>& Ptxt<CKKS>::complexConj()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call complexConj on "
                           "default-constructed Ptxt");
  for (auto& x : this->slots)
    x = std::conj(x);
  return *this;
}

template <>
template <>
Ptxt<CKKS> Ptxt<CKKS>::real() const
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call real on default-constructed Ptxt");
  Ptxt<CKKS> real(*this);
  for (auto& x : real.slots)
    x = x.real();
  return real;
}

template <>
template <>
Ptxt<CKKS> Ptxt<CKKS>::imag() const
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call imag on default-constructed Ptxt");
  Ptxt<CKKS> imag(*this);
  for (auto& x : imag.slots)
    x = x.imag();
  return imag;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::runningSums()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call runningSums on "
                           "default-constructed Ptxt");
  for (std::size_t i = 1; i < size(); ++i)
    slots[i] += slots[i - 1];
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::totalSums()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call totalSums on "
                           "default-constructed Ptxt");
  SlotType sum = slots[0];
  for (std::size_t i = 1; i < size(); ++i)
    sum += slots[i];
  setData(sum);
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::incrementalProduct()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call incrementalProduct on "
                           "default-constructed Ptxt");
  for (std::size_t i = 1; i < size(); ++i)
    slots[i] *= slots[i - 1];
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::totalProduct()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call totalProduct on "
                           "default-constructed Ptxt");
  SlotType product = slots[0];
  for (std::size_t i = 1; i < size(); ++i)
    product *= slots[i];
  setData(product);
  return *this;
}

template <typename Scheme>
Ptxt<Scheme>& Ptxt<Scheme>::mapTo01()
{
  assertTrue<RuntimeError>(isValid(),
                           "Cannot call mapTo01 on default-constructed Ptxt");
  for (auto& slot : slots)
    if (slot != Ptxt<Scheme>::convertToSlot(*context, 0l))
      slot = 1;
  return *this;
}

// VJS-FIXME: all of this logic for coordToIndex and indexToCoord
// duplicates logic already implemented in EncryptedArrayBase.

template <typename Scheme>
long Ptxt<Scheme>::coordToIndex(const std::vector<long>& coords)
{
  const PAlgebra& zMStar = context->getZMStar();
  assertEq<LogicError>(coords.size(),
                       static_cast<std::size_t>(zMStar.numOfGens()),
                       "Coord must have same size as hypercube structure");
  long index = 0;
  // Convert the coordinates into its corresponding index by computing the
  // expression Sum_{i} (g_i Prod_{j} e_j), where g is the generator of the
  // Zm*<p> group and e is the respective order of each generator.
  for (long i = coords.size() - 1; i >= 0; --i) {
    long product = 1;
    for (std::size_t j = i + 1; j <= coords.size() - 1; ++j) {
      product *= zMStar.OrderOf(j);
    }
    index += coords.at(i) * product;
  }
  return index;
}

// TODO: Refactor this out into a free function, also see the Ctxt version for
// refactoring.
template <typename Scheme>
std::vector<long> Ptxt<Scheme>::indexToCoord(long index)
{
  const PAlgebra& zMStar = context->getZMStar();
  long num_gens = zMStar.numOfGens();
  assertInRange<LogicError>(index, 0l, lsize(), "Index out of range");
  std::vector<long> coords(num_gens);
  long coord = 0;
  long product = 1;
  // Convert an index into its corresponding coordinate representation by
  // computing the equation c_i = floor(index / Prod_{j>i+1} e_j) where c_i is
  // the i'th coordinate, e is the order of each generator and Index is the
  // input index to convert.
  for (long j = 1; j < num_gens; ++j)
    product *= zMStar.OrderOf(j);
  for (long i = 0; i < num_gens; ++i) {
    coord = index / product;
    index -= coord * product;
    if (i < num_gens - 1)
      product /= zMStar.OrderOf(i + 1);
    coords[i] = coord;
  }
  return coords;
}

template <>
template <>
PA_GF2::RX Ptxt<BGV>::slotsToRX<PA_GF2>() const
{
  assertEq<LogicError>(context->getAlMod().getPPowR(),
                       2l,
                       "Plaintext modulus p^r must be equal to 2^1");
  return NTL::conv<NTL::GF2X>(getPolyRepr());
}

template <>
template <>
PA_zz_p::RX Ptxt<BGV>::slotsToRX<PA_zz_p>() const
{
  assertNeq<LogicError>(context->getAlMod().getPPowR(),
                        2l,
                        "Plaintext modulus p^r must not be equal to 2^1");
  return NTL::conv<NTL::zz_pX>(getPolyRepr());
}

template <>
void Ptxt<BGV>::assertSlotsCompatible(const std::vector<SlotType>& slots) const
{
  for (const auto& slot : slots) {
    if (slot.getp2r() != context->getSlotRing()->p2r)
      throw RuntimeError("Mismatching p^r found");
    if (slot.getG() != context->getSlotRing()->G)
      throw RuntimeError("Mismatching G found");
  }
}

template <>
void Ptxt<CKKS>::assertSlotsCompatible(const std::vector<SlotType>&) const
{
  // Do nothing.  A vector of std::complex<double> is always valid.
}

template <>
template <typename type>
NTL::ZZX Ptxt<BGV>::automorph_internal(long k)
{
  NTL::zz_pContext pContext;
  pContext.save();
  NTL::zz_p::init(context->getSlotRing()->p2r);
  long m = context->getM();
  auto old_slots = slotsToRX<type>();
  decltype(old_slots) new_slots;
  plaintextAutomorph(new_slots,
                     old_slots,
                     k,
                     m,
                     context->getAlMod().getDerived(type()).getPhimXMod());
  NTL::ZZX ret = NTL::conv<NTL::ZZX>(new_slots);
  pContext.restore();
  return ret;
}

template <typename Scheme>
void Ptxt<Scheme>::writeToJSON(std::ostream& os) const
{
  // TODO: add JSON try-catch wrapper
  assertTrue<RuntimeError>(this->isValid(),
                           "Cannot operate on invalid "
                           "(default constructed) Ptxt");

  executeRedirectJsonError<void>([&]() { os << this->writeToJSON(); });
}

template <typename Scheme>
JsonWrapper Ptxt<Scheme>::writeToJSON() const
{
  assertTrue<RuntimeError>(this->isValid(),
                           "Cannot operate on invalid "
                           "(default constructed) Ptxt");
  auto body = [this]() {
    json jslots;

    if constexpr (std::is_same_v<Scheme, CKKS>) {
      jslots = slots;
    } else {
      jslots = writeVectorToJSON(slots);
    }

    json j{{"scheme", Scheme::schemeName}, {"slots", jslots}};

    return wrap(toTypedJson<Ptxt<Scheme>>(j));
  };

  return executeRedirectJsonError<JsonWrapper>(body);
}

template <typename Scheme>
Ptxt<Scheme> Ptxt<Scheme>::readFromJSON(std::istream& is,
                                        const Context& context)
{
  Ptxt<Scheme> ret{context};
  ret.readJSON(is);
  return ret;
}

template <typename Scheme>
Ptxt<Scheme> Ptxt<Scheme>::readFromJSON(const JsonWrapper& tjw,
                                        const Context& context)
{
  Ptxt<Scheme> ret{context};
  ret.readJSON(tjw);
  return ret;
}

template <typename Scheme>
void Ptxt<Scheme>::readJSON(std::istream& is)
{
  assertTrue<RuntimeError>(this->isValid(),
                           "Cannot operate on invalid "
                           "(default constructed) Ptxt");
  executeRedirectJsonError<void>([&]() {
    json j;
    is >> j;
    this->readJSON(wrap(j));
  });
}

template <typename Scheme>
void Ptxt<Scheme>::readJSON(const JsonWrapper& tjw)
{
  assertTrue<RuntimeError>(this->isValid(),
                           "Cannot operate on invalid "
                           "(default constructed) Ptxt");
  auto body = [&]() {
    json tj = unwrap(tjw);
    json jslots;
    // if the input is just an array short-circuit to slot deserialization
    // (assuming there is no type-header).
    if (tj.is_array()) {
      jslots = tj;

    } else {
      json j = fromTypedJson<Ptxt<Scheme>>(tj);

      std::string expected_scheme{j.at("scheme").get<std::string>()};
      assertTrue<IOError>(
          Scheme::schemeName == expected_scheme,
          "Scheme mismatch in deserialization.\nExpected: " + expected_scheme +
              ", actual: " + std::string(Scheme::schemeName) + ".");

      jslots = j.at("slots");

      if (!jslots.is_array()) {
        throw IOError("Slot content is not a JSON array");
      }
    }

    if (static_cast<long>(jslots.size()) > this->context->getEA().size()) {
      std::stringstream err_msg;
      err_msg << "Cannot deserialize to Ptxt: not enough slots.  "
              << "Trying to deserialize " << jslots.size() << " elements.  "
              << "Got " << this->context->getEA().size() << " slots.";
      throw IOError(err_msg.str());
    }

    if constexpr (std::is_same_v<Scheme, CKKS>) {
      // Scheme is CKKS
      this->setData(jslots.get<std::vector<std::complex<double>>>());
    } else {
      // Scheme is BGV
      this->setData(
          readVectorFromJSON<PolyMod>(jslots, context->getSlotRing()));
    }
  };

  executeRedirectJsonError<void>(body);
}

template <typename Scheme>
std::istream& operator>>(std::istream& is, Ptxt<Scheme>& ptxt)
{
  assertTrue<RuntimeError>(ptxt.isValid(),
                           "Cannot operate on invalid "
                           "(default constructed) Ptxt");
  ptxt.readJSON(is);
  return is;
}

// Explicit operator << instantiation
template std::istream& operator>><BGV>(std::istream& is, Ptxt<BGV>& ptxt);
template std::istream& operator>><CKKS>(std::istream& is, Ptxt<CKKS>& ptxt);

template <typename Scheme>
std::ostream& operator<<(std::ostream& os, const Ptxt<Scheme>& ptxt)
{
  assertTrue<RuntimeError>(ptxt.isValid(),
                           "Cannot operate on invalid "
                           "(default constructed) Ptxt");

  ptxt.writeToJSON(os);
  return os;
}

// Explicit operator << instantiation
template std::ostream& operator<<<BGV>(std::ostream& os, const Ptxt<BGV>& ptxt);
template std::ostream& operator<<<CKKS>(std::ostream& os,
                                        const Ptxt<CKKS>& ptxt);

// Explicit class instantiation.
template class Ptxt<BGV>;
template class Ptxt<CKKS>;

} // namespace helib
