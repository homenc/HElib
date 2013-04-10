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
#include <NTL/vec_long.h>
#include "IndexMap.h"

NTL_CLIENT

// typedef IndexMapInit<vec_long> HelperBase;

class Helper : public IndexMapInit<vec_long> {
public:

  long val;
  Helper(long _val) { val = _val; }

  virtual void init(vec_long& v) { cerr << "YES\n"; v.SetLength(val); }
};

  



int main() 
{

  IndexMap<vec_long> map(new Helper(10));
  IndexMap<vec_long> map2;

  vec_long v1;
  vec_long v2;
  v1.SetLength(10);
  v2.SetLength(5);
  v1[0]=1;
  v2[0]=2;

  map.insert(7);
  map.remove(8);
  map.insert(2);
  map[7] = v1;
  map[2] = v2; // this should fail, I think


  map2.insert(3);
  map2.insert(5);
  map2[3] = v1;
  map2[5] = v2;

  cout << "map[2]=" << map[2] << "\n";
  cout << "map[7]=" << map[7] << "\n";

  cout << "map2[3]=" << map2[3] << "\n";
  cout << "map2[5]=" << map2[5] << "\n";

  return 0;

}



