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

#include "binio.h"

/* Some utility functions for binary IO */

void write_ntl_vec_long(ostream& str, const vec_long& vl)
{
  
  long sizeOfVL = vl.length();
  str.write(reinterpret_cast<char*>(&sizeOfVL), sizeof(sizeOfVL)); 

  for(long i=0, tmp=0; i<sizeOfVL; i++){ 
    tmp=vl[i];
    str.write(reinterpret_cast<char*>(&tmp), sizeof(tmp)); 
//    cerr << "[write ntl vec long] value:" << tmp << endl;
  }
}

void read_ntl_vec_long(istream& str, vec_long& vl)
{
  long sizeOfVL;
  str.read(reinterpret_cast<char*>(&sizeOfVL), sizeof(sizeOfVL)); 
//  cerr << "[read ntl vec long] size of vec long " << sizeOfVL << endl;

  // Remeber to check and increase Vec before trying to fill it.
  if(vl.length() < sizeOfVL){
    vl.SetLength(sizeOfVL);
  }

  for(long i=0, tmp=0; i<sizeOfVL; i++){
    str.read(reinterpret_cast<char*>(&tmp), sizeof(tmp)); 
    vl[i] = tmp;
//    cerr << "[read ntl vec long] value:" << vl[i] << endl;
  }
}

void write_raw_long(ostream& str, long n)
{  
  str.write(reinterpret_cast<char*>(&n), sizeof(n)); 
//  cerr << "Look at what I am writting out:" << n << endl;
}

void read_raw_long(istream& str, long& n)
{  
  str.read(reinterpret_cast<char*>(&n), sizeof(n)); 
//  cerr << "Look at what I am reading in:" << n << endl;
}


void write_raw_xdouble(ostream& str, const xdouble xd)
{
  
  double m = xd.mantissa();
  long e = xd. exponent();

  str.write(reinterpret_cast<char*>(&m), sizeof(m)); 
  str.write(reinterpret_cast<char*>(&e), sizeof(e)); 
  
}

void read_raw_xdouble(istream& str, xdouble& xd)
{

  double m;
  long e;

  str.read(reinterpret_cast<char*>(&m), sizeof(m)); 
  str.read(reinterpret_cast<char*>(&e), sizeof(e)); 

  xd = xdouble(m,e);

}

void write_raw_ZZ(ostream& str, const ZZ& zz)
{
  long noBytes = NumBytes(zz);
//  cerr << "Number of bytes: " << noBytes << endl;
  assert(noBytes > 0);
  unsigned char zzBytes[noBytes];
  // From ZZ.h
  BytesFromZZ(zzBytes, zz, noBytes);
  str.write(reinterpret_cast<char*>(&noBytes), sizeof(noBytes)); 
  str.write(reinterpret_cast<char*>(zzBytes), noBytes); 
  
//  cerr << "[PORK] Write Complete\n";
}

void read_raw_ZZ(istream& str, ZZ& zz)
{
  long noBytes = 0;
//  cerr << "Number of bytes: " << noBytes << endl;
  str.read(reinterpret_cast<char*>(&noBytes), sizeof(noBytes)); 
  assert(noBytes > 0);
  unsigned char zzBytes[noBytes];
  str.read(reinterpret_cast<char*>(zzBytes), noBytes); 
  zz = ZZFromBytes(zzBytes, noBytes);

//  cerr << "[PORK] Read Complete\n";
}


template<> void read_raw_vector<long>(istream& str, vector<long>& v)
{

  long sz; 
  str.read(reinterpret_cast<char*>(&sz), sizeof(sz)); 
  v.resize(sz); // Make space in vector

  for(long i=0, n=0; i<sz; v[i]=n, i++)
    str.read(reinterpret_cast<char*>(&n), sizeof(n)); 
    
};

template<> void write_raw_vector<long>(ostream& str, const vector<long>& v)
{
  long sz = v.size();  
  str.write(reinterpret_cast<char*>(&sz), sizeof(sz)); 

  for(long n: v)
    str.write(reinterpret_cast<char*>(&n), sizeof(n)); 
};
