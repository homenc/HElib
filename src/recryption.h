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
#ifndef _RECRYPTION_H_
#define _RECRYPTION_H_
/** @file recryption.h
 *  @brief Define some data structures to hold recryption data
 */

#include "NumbTh.h"

class  PAlgebraMod;
class  EncryptedArray;
class  EvalMap;
class  PowerfulDCRT;
class  FHEcontext;

//! @class RecryptData
//! @brief A structure to hold recryption-related data inside the FHEcontext
class RecryptData {
public:
  //! default Hamming weight of recryption key
  static const long defSkHwt=56;

  //! Some data members that are only used for I/O
  Vec<long> mvec;     //! partition of m into co-prime factors
  long hwt;           //! Hamming weight of recryption secret-key
  bool conservative;  //! flag for choosing more conservatice parameters

  //! skey encrypted wrt space p^{e-e'+r}
  long e, ePrime;

  //! Hamming weight of recryption secret key
  long skHwt;

  //! an optimization parameter
  double alpha;

  //! for plaintext space p^{e-e'+r}
  PAlgebraMod *alMod;

  //! for plaintext space p^{e-e'+r}
  EncryptedArray *ea;

  //! linear maps
  EvalMap *firstMap, *secondMap;

  //! conversion between ZZX and Powerful
  PowerfulDCRT *p2dConv;

  //! linPolys for uppacking the slots
  vector<ZZX> unpackSlotEncoding;

  RecryptData() {
    hwt=0; conservative=false; e=ePrime=0; alpha=0.0;
    alMod=NULL; ea=NULL; firstMap=secondMap=NULL;p2dConv=NULL;
  }
  ~RecryptData();

  //! Initialize the recryption data in the context
  void init(const FHEcontext& context, const Vec<long>& mvec_,
	    long t=0/*min Hwt for sk*/, bool consFlag=false);

  bool operator==(const RecryptData& other) const;
  bool operator!=(const RecryptData& other) const {
    return !(operator==(other));
  }
};

#endif /* _RECRYPTION_H_ */
