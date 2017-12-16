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
#include <iostream>
#include <vector>
#include <iomanip>
#include <cassert>

using namespace std;


/**
 * @class Cube
 * @brief Indexing into a hypercube
 **/
class Cube {
private:
   vector<int> dims;  // dims[i] is the size along the i'th diemnsion
   vector<int> prods; // prods[i] = \prod_{j=i}^{n-1} dims[i]
   vector<int> data;
   long size;

   Cube();

public:
   Cube(const vector<int>& idims) {
      dims = idims;
      long n = dims.size();
      assert(n > 0);

      
      prods.resize(n+1);
      prods[n] = 1;
      for (long i = n-1; i >= 0; i--) {
         assert(dims[i] > 0);
         prods[i] = dims[i]*prods[i+1];
      }


      size = prods[0];
      data.resize(size);
      for (long i = 0; i < size; i++) data[i] = 0; 
   }

   bool operator==(const Cube& that) const {
      return dims == that.dims && data == that.data;
   }

   bool operator!=(const Cube& that) const {
      return !(*this == that);
   }

   const vector<int>& getDims() const { return dims; }

   long getSize() const { return size; }

   long getNumDims() const { return dims.size(); }

   long getDim(long d) const { return dims.at(d); }

   long getCoord(long i, long d) const {
      assert(i >= 0 && i < size);
   
      return (i % prods.at(d)) / prods.at(d+1); 
   }

   long addCoord(long i, long d, long offset) const {
      assert(i >= 0 && i < size);
      
      offset = offset % dims.at(d);
      if (offset < 0) offset += dims.at(d);

      long i_d = getCoord(i, d);
      long i_d1 = (i_d + offset) % dims.at(d);

      long i1 = i + (i_d1 - i_d) * prods.at(d+1);

      return i1;
   }

   int& at(long i) { return data.at(i); }

   const int& at(long i) const { return data.at(i); }

};

Cube simpleRotate(const Cube& c, long offset) {

   assert(offset >= 0);
   Cube c1(c.getDims());

   long size = c.getSize();

   for (long i = 0; i < size; i++) {
      c1.at((i+offset)%size) = c.at(i);
   }

   return c1;
}

Cube rotate1D(const Cube& c, long d, long offset) {

   assert(offset >= 0);

   Cube c1(c.getDims());

   long size = c.getSize();

   for (long i = 0; i < size; i++) {
      c1.at(c.addCoord(i, d, offset)) = c.at(i);
   }

   return c1;
}

Cube operator+(const Cube& c1, const Cube& c2) {
   const vector<int>& dims1 = c1.getDims();
   const vector<int>& dims2 = c2.getDims();
   assert(dims1 == dims2);

   Cube c(dims1);
   long size = c.getSize();

   for (long i = 0; i < size; i++) 
      c.at(i) = c1.at(i) + c2.at(i); 
   
   return c;
}

Cube operator*(const Cube& c1, const Cube& c2) {
   const vector<int>& dims1 = c1.getDims();
   const vector<int>& dims2 = c2.getDims();
   assert(dims1 == dims2);

   Cube c(dims1);
   long size = c.getSize();

   for (long i = 0; i < size; i++) 
      c.at(i) = c1.at(i) * c2.at(i); 
   
   return c;
}

Cube operator!(const Cube& c1) {
   const vector<int>& dims1 = c1.getDims();

   Cube c(dims1);
   long size = c.getSize();

   for (long i = 0; i < size; i++) 
      c.at(i) = c1.at(i) == 0 ? 1 : 0;
   
   return c;
}

Cube computeMask(const vector<int>& dims, long d, long k) {
   Cube c(dims);
   long size = c.getSize();

   for (long i = 0; i < size; i++) {
      if (c.getCoord(i, d) >= k)
         c.at(i) = 1;
   }

   return c;
}

void print3D(const Cube& c) {
   const vector<int>& dims = c.getDims();
   assert(dims.size() == 3);

   long size = c.getSize();

   for (long i = 0; i < size; i++) {
      cout << setw(3) << c.at(i);
      if ((i+1) % dims.at(2) == 0) cout << endl;
      if ((i+1) % (dims.at(1)*dims.at(2)) == 0) cout << endl;
   }
}

Cube fancyRotate(const Cube& c, long offset) {

   const vector<int>& dims = c.getDims();

   assert(offset >= 0);
   offset = offset % c.getSize();

   Cube c1 = c;
   Cube mask = !Cube(dims); // the all-1 cube

   for (long d = dims.size()-1; d >= 0; d--) {
      long k = c.getCoord(offset, d);

      c1 = (rotate1D(c1, d, k)*mask) + (rotate1D(c1, d, k+1)*(!mask));
      mask = computeMask(dims, d, k)*mask + computeMask(dims, d, k+1)*(!mask); 
   }

   return c1;
}

int main()
{
   vector<int> dims(5);
   dims[0] = 5;
   dims[1] = 4;
   dims[2] = 3;
   dims[3] = 3;
   dims[4] = 4;

   Cube c(dims);
   long size = c.getSize();
   for (long i = 0; i < size; i++) c.at(i) = i;

   bool fail=false;
   for (long offset = 1; offset < size; offset++) {
      Cube c1 = simpleRotate(c, offset);
      Cube c2 = fancyRotate(c, offset);
      if (c1 != c2) { 
	cout << offset << ": BOO HOO\n";
	fail = true;
      }
   }
   if (!fail) cout << "YEE HA\n";
}

