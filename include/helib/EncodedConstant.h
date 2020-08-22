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
#ifndef HELIB_ENCODED_CONSTANT_H
#define HELIB_ENCODED_CONSTANT_H

#include <helib/DoubleCRT.h>

namespace helib {





class EncodedConstant {

private:
  long ptxtSpace;
  zzX poly;


public:
  EncodedConstant() : ptxtSpace(-1) { }
  // ptxtSpace == -1 signifies an "uninitialzed constant"

  bool initialized() const { return ptxtSpace != -1; }

  long getPtxtSpace() const { return ptxtSpace; }
  const zzX& getPoly() const { return poly; }


  friend class EncryptedArray;
};



}

#endif
