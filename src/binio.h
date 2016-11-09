#ifndef  _BINIO_H_
#define  _BINIO_H_


#include <type_traits>
#include "FHE.h"


/* Some utility functions for binary IO */

void write_ntl_vec_long(ostream& str, const vec_long& vl);

void read_ntl_vec_long(istream& str, vec_long& vl);

void write_raw_long(ostream& str, long n);

void read_raw_long(istream& str, long& n);

void write_raw_xdouble(ostream& str, const xdouble xd);

void read_raw_xdouble(istream& str, xdouble& xd);

void write_raw_ZZ(ostream& str, const ZZ& zz);

void read_raw_ZZ(istream& str, ZZ& zz);


template<typename T> void write_raw_vector(ostream& str, const vector<T>& v)
{
  long sz = v.size();  
  str.write(reinterpret_cast<char*>(&sz), sizeof(sz)); 
  cerr << "[write_raw_vector] vector size:" << sz << endl;

  for(auto n: v){
    n.write(str);
    cerr << "[write_raw_vector] value:" << n << endl;
  }
};

template<> void write_raw_vector<long>(ostream& str, const vector<long>& v);

template<typename T> void read_raw_vector(istream& str, vector<T>& v, T& init)
{

  long sz; 
  str.read(reinterpret_cast<char*>(&sz), sizeof(sz)); 
  cerr << "[read_raw_vector] resizing vector to:" << sz << endl;
  v.resize(sz, init); // Make space in vector

  for(auto& n: v){
    n.read(str);
    cerr << "[read_raw_vector] value read:" << n << endl;
  }   

};

template<typename T> void read_raw_vector(istream& str, vector<T>& v)
{
  T init = T(); 
  read_raw_vector<T>(str, v, init);  
}

template<> void read_raw_vector<long>(istream& str, vector<long>& v);

// At least KeySwitch requires the context.
template<typename T> void read_raw_vector(istream& str, vector<T>& v, const FHEcontext& context)
{
 
  long sz; 
  str.read(reinterpret_cast<char*>(&sz), sizeof(sz)); 
  cerr << "[read_raw_vector] resizing vector to:" << sz << endl;
  v.resize(sz); // Make space in vector

  for(auto& n: v){
    n.read(str, context);
    cerr << "[read_raw_vector] value read:" << n << endl;
  }   
 
}

#endif // ifndef _BINIO_H_
