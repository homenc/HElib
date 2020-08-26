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
  zzX poly;
  long ptxtSpace; 


public:

  const zzX& getPoly() const { return poly; }
  long getPtxtSpace() const { return ptxtSpace; }

  EncodedPtxt_BGV(const zzX& poly_, long ptxtSpace_)
    : poly(poly_), ptxtSpace(ptxtSpace_) { }

};

class EncodedPtxt_CKKS {

private:
  zzX poly;
  double mag, scale, err;

public:

  const zzX& getPoly() const { return poly; }
  double getMag() const { return mag; }
  double getScale() const { return scale; }
  double getErr() const { return err; }

  EncodedPtxt_CKKS(const zzX& poly_, double mag_, double scale_, double err_)
    : poly(poly_), mag(mag_), scale(scale_), err(err_) { }

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

* eptxt.resetBGV(poly, ptxtSpace) replaces the contents of eptxt with a new
     EncodedPtxt_BGV object

* eptxt.resetCKKS(poly, mag, scale, err) replaces the contents of eptxt 
     with a new EncodedPtxt_CKKS object

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

  using EncodedPtxt_BGV::EncodedPtxt_BGV;
};


class EncodedPtxt_derived_CKKS : public EncodedPtxt_base, 
                                public EncodedPtxt_CKKS {
public:

  virtual EncodedPtxt_base* clone() const override 
  { return new EncodedPtxt_derived_CKKS(*this); }

  virtual bool isCKKS() const override { return true; }

  virtual const EncodedPtxt_CKKS& getCKKS() const override { return *this; }

  using EncodedPtxt_CKKS::EncodedPtxt_CKKS;
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

  void resetBGV(const zzX& poly, long ptxtSpace) 
  { 
    rep.set_ptr(new EncodedPtxt_derived_BGV(poly, ptxtSpace));
  }

  void resetCKKS(const zzX& poly, double mag, double scale, double err) 
  { 
    rep.set_ptr(new EncodedPtxt_derived_CKKS(poly, mag, scale, err));
  }

};


}

#endif
