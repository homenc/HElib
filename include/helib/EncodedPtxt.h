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
#ifndef HELIB_ENCODED_PTXT_H
#define HELIB_ENCODED_PTXT_H

#include <helib/DoubleCRT.h>

namespace helib {





class EncodedPtxt_BGV {

private:
  long ptxtSpace = -1; // ptxtSpace < 0 signifies an "uninitialized"
  zzX poly;


public:
  bool initialized() const { return ptxtSpace >= 0; }

  long getPtxtSpace() const { return ptxtSpace; }
  const zzX& getPoly() const { return poly; }


  friend class EncryptedArrayBase;
};

class EncodedPtxt_CKKS {

private:
  double mag = -1.0, scale = -1.0, err = -1.0;
  // mag < 0 signifies an "uninitialized"

  zzX poly;


public:

  bool initialized() const { return mag >= 0; }
  

  double getMag() const { return mag; }
  double getScale() const { return scale; }
  double getErr() const { return err; }

  const zzX& getPoly() const { return poly; }


  friend class EncryptedArrayCx;
};


/*

The following is a lot of "boilerplate" code that 
gives a type-safe union of the above two EncodedPtxt types.

The idea is that an EncodedPtxt object eptxt contains
either an EncodedPtxt_BGV object, an EncodedPtxt_CKKS object,
or nothing at all.

* eptxt.isBGV() return true if it contains an EncodedPtxt_BGV object 

* eptxt.isCKKS() return true if it contains an EncodedPtxt_CKKS object 

* eptxt.getBGV() returns a read-only reference to the contained
     EncodedPtxt_BGV object (or throws a std::bad_cast exception)

* eptxt.getCKKS() returns a read-only reference to the contained
     EncodedPtxt_CKKS object (or throws a std::bad_cast exception)

* eptxt.resetBGV() replaces the contents of eptxt with a new
     EncodedPtxt_BGV object and returns a reference to that object

* eptxt.resetCKKS() replaces the contents of eptxt with a new
     EncodedPtxt_CKKS object and returns a reference to that object

The default contructor creates an empty container.
Copy contructors and assignmemnt operators are available,
   and provide "deep" copy semantics

*/


class EncodedPtxt_base {
public:

  virtual ~EncodedPtxt_base() {}

  virtual EncodedPtxt_base* clone() const = 0;
  // makes this usable with cloned_ptr

  virtual bool isBGV() const { return false; }
  virtual bool isCKKS() const { return false; }

  virtual const EncodedPtxt_BGV& getBGV() const { throw std::bad_cast(); }
  virtual const EncodedPtxt_CKKS& getCKKS() const { throw std::bad_cast(); }

};

class EncodedPtxt_derived_BGV : public EncodedPtxt_base, 
                                public EncodedPtxt_BGV {
public:

  virtual EncodedPtxt_base* clone() const override 
  { return new EncodedPtxt_derived_BGV(*this); }

  virtual bool isBGV() const override { return true; }

  virtual const EncodedPtxt_BGV& getBGV() const override { return *this; }
};


class EncodedPtxt_derived_CKKS : public EncodedPtxt_base, 
                                public EncodedPtxt_CKKS {
public:

  virtual EncodedPtxt_base* clone() const override 
  { return new EncodedPtxt_derived_CKKS(*this); }

  virtual bool isCKKS() const override { return true; }

  virtual const EncodedPtxt_CKKS& getCKKS() const override { return *this; }
};


class EncodedPtxt {
  cloned_ptr<EncodedPtxt_base> rep;

public:

  bool isBGV() const { return !rep.null() && rep->isBGV(); }
  bool isCKKS() const { return !rep.null() && rep->isCKKS(); }

  const EncodedPtxt_BGV& getBGV() const 
  { if (rep.null()) throw std::bad_cast(); return rep->getBGV(); }

  const EncodedPtxt_CKKS& getCKKS() const 
  { if (rep.null()) throw std::bad_cast(); return rep->getCKKS(); }

  EncodedPtxt_BGV& resetBGV() 
  { 
    EncodedPtxt_derived_BGV *p = new EncodedPtxt_derived_BGV();
    rep.set_ptr(p);
    return *p;
  }

  EncodedPtxt_CKKS& resetCKKS() 
  { 
    EncodedPtxt_derived_CKKS *p = new EncodedPtxt_derived_CKKS();
    rep.set_ptr(p);
    return *p;
  }

};


}

#endif
