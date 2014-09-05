#ifndef _ALT_EVAL_MAP_H_
#define _ALT_EVAL_MAP_H_
/* Copyright (C) 2012-2014 IBM Corp.
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
/* AltEvalMap.h - Implementing the reCrypt linear transformations
 */
#include "FHE.h"
#include "Ctxt.h"
#include "EncryptedArray.h"
#include "permutations.h"

class AltEvalMap {
private:
  const EncryptedArray& ea;
  bool invert; 

  long nfactors;

  shared_ptr<PlaintextBlockMatrixBaseInterface> mat1;

  Vec< shared_ptr<PlaintextMatrixBaseInterface> > matvec;

  
public:

  AltEvalMap(const EncryptedArray& _ea, const Vec<long>& mvec, bool _invert);

  void apply(Ctxt& ctxt) const;
};

#endif
