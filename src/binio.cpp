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

#include <cstring>
#include "binio.h"

/* Some utility functions for binary IO */

int readEyeCatcher(istream& str, const char * expect)
{
  char eye[BINIO_EYE_SIZE];
  str.read(eye, BINIO_EYE_SIZE); 
  return memcmp(eye, expect, BINIO_EYE_SIZE);
}

void writeEyeCatcher(ostream& str, const char* eyeStr)
{
  char eye[BINIO_EYE_SIZE];
  memcpy(eye, eyeStr, BINIO_EYE_SIZE);
  str.write(eye, BINIO_EYE_SIZE);
}

// compile only 64-bit (-m64) therefore long must be at least 64-bit
long read_raw_int(istream& str, long intSize)
{
  long result = 0;
  char byte;

  for(long i=0; i<intSize; i++){
    str.read(&byte, 1); // read a byte
    result |= (static_cast<long>(byte)&0xff) << i*8; // must be in little endian
  }

  return result;
}

// compile only 64-bit (-m64) therefore long must be at least 64-bit
void write_raw_int(ostream& str, long num, long intSize) 
{
  char byte;

  for(long i=0; i<intSize; i++){
    byte = num >> 8*i; // serializing in little endian
    str.write(&byte, 1);  // write byte out
  }
}

void write_ntl_vec_long(ostream& str, const vec_long& vl, long intSize)
{
  write_raw_int(str, vl.length(), BINIO_32BIT); 
  write_raw_int(str, intSize, BINIO_32BIT); 

  for(long i=0, tmp; i<vl.length(); i++){
    write_raw_int(str, vl[i], intSize); 
  }
}

void read_ntl_vec_long(istream& str, vec_long& vl)
{
  long sizeOfVL = read_raw_int(str, BINIO_32BIT);
  long intSize  = read_raw_int(str, BINIO_32BIT);

  // Remember to check and increase Vec before trying to fill it.
  if(vl.length() < sizeOfVL){
    vl.SetLength(sizeOfVL);
  }

  for(long i=0, tmp; i<sizeOfVL; i++){
    vl[i] = read_raw_int(str, intSize);
  }
}


void write_raw_xdouble(ostream& str, const xdouble xd)
{
  
  double m = xd.mantissa();
  long e = xd.exponent();

  long *pm = reinterpret_cast<long*>(&m); 
  write_raw_int(str, *pm);
  write_raw_int(str, e);
}

xdouble read_raw_xdouble(istream& str)
{

  long m = read_raw_int(str);
  long e = read_raw_int(str);
  double* pm = reinterpret_cast<double*>(&m);

  return xdouble(*pm,e);

}

void write_raw_ZZ(ostream& str, const ZZ& zz)
{
  long noBytes = NumBytes(zz);
  assert(noBytes > 0);
  unsigned char zzBytes[noBytes];
  BytesFromZZ(zzBytes, zz, noBytes); // From ZZ.h
  write_raw_int(str, noBytes); 
  // TODO - ZZ appears to be endian agnostic
  str.write(reinterpret_cast<char*>(zzBytes), noBytes); 
}

void read_raw_ZZ(istream& str, ZZ& zz)
{
  long noBytes = read_raw_int(str);
  assert(noBytes > 0);
  unsigned char zzBytes[noBytes];
  // TODO - ZZ appears to be endian agnostic
  str.read(reinterpret_cast<char*>(zzBytes), noBytes); 
  zz = ZZFromBytes(zzBytes, noBytes);
}

template<> void read_raw_vector<long>(istream& str, vector<long>& v)
{

  long sz = read_raw_int(str);
  v.resize(sz); // Make space in vector

  for(long i=0; i<sz; i++)
    v[i] = read_raw_int(str); 
    
};

template<> void write_raw_vector<long>(ostream& str, const vector<long>& v)
{
  long sz = v.size();  
  write_raw_int(str, sz); 

  for(long n: v)
    write_raw_int(str, n); 
};
