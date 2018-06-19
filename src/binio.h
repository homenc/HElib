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

#ifndef  _BINIO_H_
#define  _BINIO_H_

#include <type_traits>
#include "FHE.h"

#define BINIO_32BIT 4
#define BINIO_48BIT 6
#define BINIO_64BIT 8

#define BINIO_EYE_SIZE 4

#define BINIO_EYE_CONTEXTBASE_BEGIN "|BS["
#define BINIO_EYE_CONTEXTBASE_END   "]BS|"
#define BINIO_EYE_CONTEXT_BEGIN     "|CN["
#define BINIO_EYE_CONTEXT_END       "]CN|"
#define BINIO_EYE_CTXT_BEGIN        "|CX["
#define BINIO_EYE_CTXT_END          "]CX|"
#define BINIO_EYE_PK_BEGIN          "|PK["
#define BINIO_EYE_PK_END            "]PK|"
#define BINIO_EYE_SK_BEGIN          "|SK["
#define BINIO_EYE_SK_END            "]SK|"
#define BINIO_EYE_SKM_BEGIN         "|KM["
#define BINIO_EYE_SKM_END           "]KM|"

/* This struct (or similar) is a nice to have not used at the moment. */
//struct BinaryHeader {
//  uint8_t structId[4];
//  uint8_t version[4] = {0, 0, 0, 1};
//  uint64_t id;
//  uint64_t payloadSize;
//};

/* Some utility functions for binary IO */

int readEyeCatcher(istream& str, const char * expect);
void writeEyeCatcher(ostream& str, const char* eye);

void write_ntl_vec_long(ostream& str, const vec_long& vl, long intSize=BINIO_64BIT);
void read_ntl_vec_long(istream& str, vec_long& vl);

long read_raw_int(istream& str, long intSize=BINIO_64BIT);
void write_raw_int(ostream& str, long num, long intSize=BINIO_64BIT);

void write_raw_xdouble(ostream& str, const xdouble xd);
xdouble read_raw_xdouble(istream& str);

void write_raw_ZZ(ostream& str, const ZZ& zz);
void read_raw_ZZ(istream& str, ZZ& zz);

template<typename T> void write_raw_vector(ostream& str, const vector<T>& v)
{
  long sz = v.size(); 
  write_raw_int(str, v.size()); 

  for(auto n: v){
    n.write(str);
  }
};
// vector<long> has adifferent implementation, since long.write does not work
template<> void write_raw_vector<long>(ostream& str, const vector<long>& v);

template<typename T> void read_raw_vector(istream& str, vector<T>& v, T& init)
{
  long sz = read_raw_int(str);
  v.resize(sz, init); // Make space in vector

  for(auto& n: v){
    n.read(str);
  }
};

template<typename T> void read_raw_vector(istream& str, vector<T>& v)
{
  read_raw_vector<T>(str, v, T());
}
// vector<long> has adifferent implementation, since long.read does not work
template<> void read_raw_vector<long>(istream& str, vector<long>& v);

// KeySwitch::read(...) (in FHE.cpp) requires the context.
template<typename T> void read_raw_vector(istream& str, vector<T>& v, const FHEcontext& context)
{ 
  long sz = read_raw_int(str);
  v.resize(sz); // Make space in vector

  for(auto& n: v){
    n.read(str, context);
  }
}

#endif // ifndef _BINIO_H_
