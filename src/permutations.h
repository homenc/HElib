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
/**
 * @file permutations.h
 * @brief Permutations over Hypercubes and their slices.
 **/
#ifndef _PERMUTATIONS_H_
#define _PERMUTATIONS_H_

#include "matching.h"
#include "hypercube.h"
#include <iostream>

using namespace std;
using namespace NTL;

//! A simple permutation is just a vector with p[i]=\pi_i
typedef Vec<long> Permut;

//! @brief A random size-n permutation
void randomPerm(Permut& perm, long n);


class SlicePerm; // a "row permutation"
class ColPerm;   // a "column permutation"

/**
 * @class SlicePerm
 * @brief Permuting subcubes (rows/slices) of a hypercube.
 *
 * SlicePerm is derived from a HyperCube<long>, and it uses the cube object
 * to store the actual permutation data. The interpretation of this data,
 * however, depends on the data member dim.
 * 
 * A SlicePerm is breaking the cube into consecutive subcubes and permuting
 * each subcube separately. Specifically, the cube is partitioned into 
 * n=getSize()/getProd(dim) subcubes, so we have n consequtive permutations,
 * all over the domain [0,getProd(dim)-1]. A SlicePerm with dim=0 is just a
 * regular permutation over the entire cube.
 *
 * For example, for a 2x3x2 cube and dim=1 (i.e., the 2nd dimension), the data
 * vector could be something like [1 0 2 5 4 3 3 0 5 4 1 2]. This means the
 * permutations [1 0 2 5 4 3] over the 1st 2x3 subcube and [3 0 5 4 1 2] over
 * the 2nd subcube. The same permutation written explicitly would be
 * [1 0 2 5 4 3 9 6 11 10 7 8].
 **/
class SlicePerm : public HyperCube<long> {
private:
  long dim;
  SlicePerm(); // disabled

public:
  friend class ColPerm;

  explicit SlicePerm(const CubeSignature& _sig): HyperCube<long>(_sig)
  { dim = -1; }
  SlicePerm(const CubeSignature& _sig, const Permut& p): HyperCube<long>(_sig)
  { *this = p; } // use assignment operator to copy data

  long getPermDim() const {return dim;}

  void makeExplicit(Permut& out) const;   //! Write the permutation explicitly
  void extractSlice(Permut& out, long i) const;//! The perm over the ith subcube

  //! When dim=0 (the first dimension), SlicePerm and Permut are equivalent
  SlicePerm& operator=(const Permut &other) {
    dim = 0;           // Set the dimension to 0
    getData() = other; // Copy the raw data
    return *this;
  }

  //! When dim=m (the last dimension), ColPerm and SlicePerm are equivalent
  SlicePerm& operator=(const ColPerm &other);

  //! @brief Break a row permutation into column-row-column format.
  //! The input permutes each dimension-i subcube, and in the output rho1,rho3
  //! permute only along the i'th dimension and rho2 permutes each dimension-i+1
  //! subcube. This method cannot permute in-place, it is assumed that rho2 and
  //! *this point to disjoint objects.
  void breakPermTo3(ColPerm& rho1, SlicePerm& rho2, ColPerm& rho3) const;

 //! A test/debugging method
  void printout(ostream& s) { // a test/debugging method
    s << "Cube signature: " << getSig() << endl;
    s << "  dim="<<dim<<endl;
    s << "  data="<<getData()<<endl;
  }
};

/**
 * @class ColPerm
 * @brief Permuting a single dimension (column) of a hypercube
 * 
 * ColPerm is derived from a HyperCube<long>, and it uses the cube object
 * to store the actual permutation data. The interpretation of this data,
 * however, depends on the data member dim.
 * 
 * The cube is partitioned into n=getSize()/getDim(dim) subcubes, and each of
 * them is permuted separatly. Hence pi consists of n interleaved permutations,
 * each over the domain [0,getDim(dim)-1], and the "step" between elements of
 * any single permutation inside the data vector is getProd(dim+1). Permuting
 * along the last dimension is the same as SlicePerm for that dimension.
 * 
 * For example, permuting a 2x3x2 cube along dim=1 (the 2nd dimention), we
 * could have the data vector as [ 1  1  0  2  2  0  2  0  1  1  0  2 ]. (In
 * this example we have n=4 and step=2.) This means the four subcubes are
 * permuted using the permutation [1 0 2],[1 2 0],[2 1 0],[0 1 2]. Written
 * explicitly, this permutation would be [2 3 0 5 4 1 10 7 8 9 11].
 **/
class ColPerm : public HyperCube<long> {
private:
  long dim;
  ColPerm(); // disabled

public:
  friend class SlicePerm;

  explicit ColPerm(const CubeSignature& _sig): HyperCube<long>(_sig)
  { dim = -1; }

  long getPermDim() const {return dim;}

  void makeExplicit(Permut& out) const;      // Write the permutation explicitly
  void extractSlice(Permut& out, long i) const; // The perm over the ith subcube

  //! When dim=m (the last dimension), ColPerm and SlicePerm are equivalent
  ColPerm& operator=(const SlicePerm &other);

 //! A test/debugging method
  void printout(ostream& s) { // a test/debugging method
    s << "Cube signature: " << getSig() << endl;
    s << "  dim="<<dim<<endl;
    s << "  data="<<getData()<<endl;
  }
};

/**
 * @brief Takes a permutation pi over m-dimensional cube C=Z_{n1} x...x Z_{nm}
 * and expresses pi as a product pi = rho_{2m-1} o ... o rho_2 o rho_1 where
 * each rho_i is a column permutation along one dimension. Specifically for
 * i<m, the permutations rho_i and rho_{2(m-1)-i} permute the i'th dimension
 **/
void breakPermByDim(vector<ColPerm>& out, 
		    const Permut &pi, const CubeSignature& sig);

#endif /* ifndef _PERMUTATIONS_H_ */
