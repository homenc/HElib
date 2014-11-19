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


#include "NumbTh.h"

class  PAlgebraMod;
class  EncryptedArray;
class  EvalMap;
class  PowerfulDCRT;
class  FHEcontext;

class RecryptData {
public:
  //! Default Hamming weight of recryption key
  static const long defSkHwt=56;

  long e, ePrime; // skey encrypted wrt space p^{e-e'+r}
  long skHwt;     // Hamming weight of recryption secret key
  double alpha;   // an optimization parameter

  PAlgebraMod *alMod; // for plaintext space p^{e-e'+r}
  EncryptedArray *ea; // for plaintext space p^{e-e'+r}
  EvalMap *firstMap, *secondMap; // linear maps
  PowerfulDCRT *p2dConv; // conversion between ZZX and Powerful

  vector<ZZX> unpackSlotEncoding;// linPolys for uppacking the slots

  RecryptData() { alMod=NULL; ea=NULL; firstMap=secondMap=NULL;p2dConv=NULL; }
  ~RecryptData();
  void init(const FHEcontext& context, const Vec<long>& mvec,
	    long t=0/*min Hwt for sk*/, bool conservative=false);
};

#endif /* _RECRYPTION_H_ */
