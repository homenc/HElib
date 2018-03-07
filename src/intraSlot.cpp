/* Copyright (C) 2012-2017 IBM Corp.
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
/* intraSlot.cpp - Packing/unpacking of mod-p integers in GF(p^d) slots.
 *
 * The normal basis and its inverse are computed when on the 1st call
 * to either ea.getNormalBasisMatrixInverse() or ea.getNormalBasisMatrix().
 */
#include <memory>
#include "replicate.h"
#include "intraSlot.h"

// Implementation classes for unpacking:
// buildUnpackSlotEncoding_pa_impl prepares the constants for the linear
// transformation, and unpack_pa_impl uses them to do the actual uppacking.

//! \cond FALSE (make doxygen ignore this code)
template<class type>
class buildUnpackSlotEncoding_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, 
                    std::vector<zzX>& unpackSlotEncoding)
  {
    FHE_NTIMER_START(buildUnpackSlotEncoding);
    RBak bak; bak.save(); ea.restoreContext();  // the NTL context for mod p^r

    long nslots = ea.size(); // how many slots
    long d = ea.getDegree(); // size of each slot

    const Mat<R>& CBi=ea.getNormalBasisMatrixInverse();
    // CBi contains a description of the normal-basis inverse transformation

    std::vector<RX> LM(d);
    for (long i = 0; i < d; i++) // prepare the linear polynomial
      LM[i] = CBi[i][0];

    std::vector<RX> C; 
    ea.buildLinPolyCoeffs(C, LM); // "build" the linear polynomial

    unpackSlotEncoding.resize(d);  // encode the coefficients
    for (long j = 0; j < d; j++) {
      std::vector<RX> v(nslots, C[j]);
      ea.encode(unpackSlotEncoding[j], v);
    }
  }
};
// A wrapper function, calls the apply method of the class above
void buildUnpackSlotEncoding(std::vector<zzX>& unpackSlotEncoding,
                             const EncryptedArray& ea)
{
  ea.dispatch<buildUnpackSlotEncoding_pa_impl>(unpackSlotEncoding);
}

template<class type>
class unpack_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, const CtPtrs& unpacked,
                    const Ctxt&ctxt,const std::vector<zzX>& unpackSlotEncoding)
  {
    long d = ea.getDegree(); // size of each slot

    //ctxt.cleanUp();

    // Convert the unpack constants to doubleCRT
    std::vector< std::shared_ptr<DoubleCRT> > coeff_vector(d);
    for (long i = 0; i < d; i++) {
      coeff_vector[i] = std::shared_ptr<DoubleCRT>(new
        DoubleCRT(unpackSlotEncoding[i],ctxt.getContext(),ctxt.getPrimeSet()));
    }
    // Compute the d Frobenius automorphisms of ctxt (use multi-threading)
    std::vector<Ctxt> frob(d, Ctxt(ZeroCtxtLike, ctxt));
    NTL_EXEC_RANGE(d, first, last)
	for (long j = first; j < last; j++) { // process jth Frobenius 
	  frob[j] = ctxt;
	  frob[j].frobeniusAutomorph(j);
	  frob[j].cleanUp();
          // NOTE: Why do we apply cleanup after the Frobenus?
	}
    NTL_EXEC_RANGE_END

    // compute the unpacked ciphertexts: the j'th slot of unpacked[i]
    // contains the i'th coefficient from the j'th clot of ctxt
    Ctxt tmp1(ZeroCtxtLike, ctxt);
    for (long i = 0; i < unpacked.size(); i++) {
      *(unpacked[i]) = frob[0];
      unpacked[i]->multByConstant(*coeff_vector[i]);
      for (long j = 1; j < d; j++) {
	tmp1 = frob[j];
	tmp1.multByConstant(*coeff_vector[mcMod(i+j, d)]);
	*(unpacked[i]) += tmp1;
      }
    } // NOTE: why aren't we using multi-threading here?
  }
};
//! \endcond

// A wrapper function, calls the apply method of the class above
void unpack(const CtPtrs& unpacked, const Ctxt& packed, 
            const EncryptedArray& ea,
            const std::vector<zzX>& unpackSlotEncoding)
// void unpack(const EncryptedArray& ea, 
// 	    std::vector<Ctxt*>& unpacked,
// 	    const Ctxt& ctxt, 
// 	    const std::vector<zzX>& unpackSlotEncoding)
{
  ea.dispatch<unpack_pa_impl>(unpacked, packed, unpackSlotEncoding);
}

// unpack many ciphertexts, returns the number of unpacked ciphertexts
long unpack(const CtPtrs& unpacked, const CtPtrs& packed,
            const EncryptedArray& ea, 
            const std::vector<zzX>& unpackSlotEncoding)
{
  long d = ea.getDegree(); // size of each slot
  long num2unpack = unpacked.size();
  assert(packed.size()*d >= num2unpack); // we must have enough ciphertexts
  long offset = 0;
  long idx = 0;
  while (num2unpack > 0) {
    if (num2unpack < d) d = num2unpack;
    const CtPtrs_slice nextSlice(unpacked, offset, d);
    ea.dispatch<unpack_pa_impl>(nextSlice,
                                *(packed[idx++]), unpackSlotEncoding);
    num2unpack -= d;
    offset += d;
  }
  return idx;
}

// An implementation classes for (re)packing.

//! \cond FALSE (make doxygen ignore this code)
template<class type>
class repack_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt,
                    const CtPtrs& unpacked)
  {
    RBak bak; bak.save(); ea.restoreContext();  // the NTL context for mod p^r
    long nslots = ea.size(); // how many slots

    const Mat<R>& CB=ea.getNormalBasisMatrix();
    // CB contains a description of the normal-basis transformation

    RX pow;
    zzX powInSlots;
    std::vector<RX> powVec(nslots);
    ctxt.clear();
    for (long i=0; i<unpacked.size(); i++) {
      conv(pow, CB[i]); // convert CB[i] from Vec<R> to RX
      for (long j=0; j < nslots; j++) powVec[j] = pow;
      ea.encode(powInSlots, powVec); // a constant with X^{p^i} in all slots
      Ctxt tmp(*(unpacked[i]));
      tmp.multByConstant(powInSlots);// accumulate unpacked[i] * X^{p^i}
      ctxt += tmp;
    }
  }
};
//! \endcond

// A wrapper function, calls the apply method of the class above
void repack(Ctxt& packed, const CtPtrs& unpacked, const EncryptedArray& ea)
{
  ea.dispatch<repack_pa_impl>(packed, unpacked);
}

// pack many ciphertexts, returns the number of packed ciphertexts
long repack(const CtPtrs& packed, const CtPtrs& unpacked, const EncryptedArray& ea)
{
  long d = ea.getDegree(); // size of each slot
  long num2pack = unpacked.size();
  assert(packed.size()*d >= num2pack); // we must have enough ciphertexts
  long offset = 0;
  long idx = 0;
  while (num2pack > 0) {
    if (num2pack < d) d = num2pack;
    const CtPtrs_slice nextSlice(unpacked, offset, d);
    ea.dispatch<repack_pa_impl>(*(packed[idx++]), nextSlice);
    num2pack -= d;
    offset += d;
  }
  return idx;
}


//! \cond FALSE (make doxygen ignore this code)
template<class type>
class packConstant_pa_impl {

public:
  PA_INJECT(type)

  // encode bottom nbits bits of data as a polynomial. It is assumed that
  // the NTL constant was set correctly before calling this function.
  static void int2Poly(RX& result, const EncryptedArrayDerived<type>& ea,
                unsigned long data, long nbits)
  {
    long d = ea.getDegree(); // size of each slot

    assert(nbits >= 0 && nbits <= d);

    const Mat<R>& CB=ea.getNormalBasisMatrix();
    // CB contains a description of the normal-basis transformation

    Vec<R> acc;
    acc.SetLength(d);
    clear(acc); // clear all d bits in acc, length stays =d

    for (long i = 0; i < nbits; i++) {  // compute sum_i bit_i * X^{p^i}
      if ((data >> i) & 1) 
        add(acc, acc, CB[i]); // add each CB[i][j] to acc[j]
    }
    conv(result, acc); // convert from vec<R> to a polynomial RX
  }

  // encode data in all slots of result
  static void apply(const EncryptedArrayDerived<type>& ea, 
		    unsigned long data,
		    long nbits,
		    zzX& result)
  {
    RBak bak; bak.save(); ea.restoreContext();  // the NTL context for mod p^r
    RX acc_poly;
    int2Poly(acc_poly, ea, data, nbits);     // endoce data as a polynomial RX
    long nslots = ea.size(); // how many slots
    vector<RX> acc_poly_vec(nslots, acc_poly);
    ea.encode(result, acc_poly_vec);
  }

  // encode different integers in the different slots
  static void apply(const EncryptedArrayDerived<type>& ea,
		    const std::vector<unsigned long>& data, long nbits,
                    zzX& result)
  {
    long nslots = ea.size(); // how many slots
    assert(lsize(data)==nslots);
    RBak bak; bak.save(); ea.restoreContext(); // the NTL context for mod p^r
    vector<RX> vec(nslots, RX::zero());
    RX acc_poly;
    for (long i=0; i<nslots; i++)
      int2Poly(vec[i], ea, data[i], nbits); // endoce data as a polynomial RX

    ea.encode(result, vec);
  }
};
//! \endcond

// this packs the low-order nbits of data into each slot,
// where the bits get mapped to coefficients on the normal basis
void packConstant(zzX& result, unsigned long data, long nbits,
                  const EncryptedArray& ea)
{
  ea.dispatch<packConstant_pa_impl>(data, nbits, result);
}


// this packs the low-order nbits of the entries of data into slots,
// where the bits get mapped to coefficients on the normal basis
void packConstants(zzX& result, const std::vector<unsigned long>& data,
                   long nbits, const EncryptedArray& ea)
{
  ea.dispatch<packConstant_pa_impl>(data, nbits, result);
}

//! \cond FALSE (make doxygen ignore this code)
template<class type>
class unpackSlots_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, 
                    NewPlaintextArray& pa, std::vector<unsigned long>& value)
  {
    PA_BOILER
    const Mat<R>& CBi=ea.getNormalBasisMatrixInverse();

    value.resize(n);
    for (long i = 0; i < n; i++) {
       Vec<R> v, w;
       VectorCopy(v, data[i], d);
       mul(w, v, CBi);
       unsigned long res = 0;
       for (long j = 0; j < d; j++) {
          if (w[j] != 0) res += (1UL << j);
       }
       value[i] = res;
    }
  }
};
//! \endcond


// Sets value to be a vector of size nslots.
// The bits of value[i] are equal to the coeffs of pa[i],
// as represented on the normal basis
void unpackSlots(std::vector<unsigned long>& value,
                 NewPlaintextArray& pa, const EncryptedArray& ea)
{
   ea.dispatch<unpackSlots_pa_impl>(pa, value);
}
