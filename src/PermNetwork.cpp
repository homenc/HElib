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
#include <NTL/ZZ.h>
#include "Ctxt.h"
#include "permutations.h"
#include "EncryptedArray.h"

void PermNetwork::BuildNetwork(const Permut& pi, const GeneratorTrees& trees)
{
  Vec<long> dims;
  trees.getCubeDims(dims);

  std::cerr << "\ndims="<<dims<<endl;
  std::cerr << "pi =      "<<pi<<endl;
  std::cerr << "map2cube ="<<trees.Map2Cube()<<endl;
  std::cerr << "map2array="<<trees.Map2Array()<<endl;

  // Compute the permutation on the cube, rho = map2cube o pi o map2array
  Permut rho;
  ApplyPermsToFunc(rho, trees.Map2Cube(), pi, trees.Map2Array());
  std::cerr << "rho =     "<<rho<<endl;

  // Break rho along the different dimensions
  CubeSignature sig(dims); // make a cube-signature object
  vector<ColPerm> perms;
  breakPermByDim(perms, rho, sig);

  for (long i=0; i<(long)perms.size(); i++) { // debugging printouts
    Permut tmp;
    perms[i].makeExplicit(tmp);
    std::cerr << " prems["<<i<<"]="<<tmp<<endl;
  }

  // Compute the total number of levels in the target permutation network
  long nLvls=0;
  for (long g=0; g<trees.length(); g++) {
    const OneGeneratorTree &T = trees[g];
    for (long i=0, leaf=T.FirstLeaf(); leaf>=0; i++, leaf=T.NextLeaf(leaf)) {
      long lvls=1; // how many levels to add for this leaf

      const SubDimension& leafData = T[leaf].getData();
      if (leafData.benes != NULL) // a benes leaf
	lvls = 3; /* FIXME: leafData.benes->getNLvls(); */

      if (g<trees.length()-1 || T.NextLeaf(leaf)>=0) 
	lvls *= 2; // not last leaf in last tree, so not middle level
      nLvls += lvls;
    }
  }
  layers.SetLength(nLvls); // allocate space

  // Go over the different permutations and build the corresponding levels
  Vec<long> shifts;
  long dimIdx=0, lyrIdx=0;
  for (long g=0; g<trees.length(); g++) {
    const OneGeneratorTree &T = trees[g];
    for (long leaf=T.FirstLeaf(); leaf>=0; leaf=T.NextLeaf(leaf)) {
      const SubDimension& leafData = T[leaf].getData();

      // This leaf corresponds to permutations idx, n-idx
      const ColPerm& p1 = perms[dimIdx];
      const ColPerm& p2 = perms[perms.size()-dimIdx-1];

      if (leafData.benes != NULL) { // A Benes leaf
	// FIXME: do something here..
	lyrIdx += 3; /* FIXME: leafData.benes->getNLvls(); */
	continue;
      }

      // A naive permutation leaf, corresponding to 1 or 2 levels
      PermNetLayer& l1 = layers[lyrIdx];
      l1.isID = !p1.getShiftAmounts(shifts);
      if (!l1.isID) {
	std::cerr << "layer "<<lyrIdx<<": "<<shifts<<endl;
	if (leafData.good) // For good leaves, shift by -x is the same as size-x
	  for (long k=0; k<shifts.length(); k++)
	    if (shifts[k]<0) shifts[k] += leafData.size;
	// li.shifts[i] := shifts[map2cube[i]]
	ApplyPermToFunc(l1.shifts, shifts, trees.Map2Cube());
	std::cerr << "       : "<<l1.shifts<<endl;
      }
      else std::cerr << "layer "<<lyrIdx<<"= identity\n";

      if (lyrIdx == layers.length()-lyrIdx-1) break;// middle layer (last leaf)
      lyrIdx = layers.length()-lyrIdx-1; // the other layers for this leaf

      PermNetLayer& l2 = layers[lyrIdx];
      l2.isID = !p2.getShiftAmounts(shifts);
      if (!l2.isID) {
	std::cerr << "layer "<<lyrIdx<<": "<<shifts<<endl;
	if (leafData.good) // For good leaves, shift by -x is the same as size-x
	  for (long k=0; k<shifts.length(); k++)
	    if (shifts[k]<0) shifts[k] += leafData.size;
	// li.shifts[i] := shifts[map2cube[i]]
	ApplyPermToFunc(l2.shifts, shifts, trees.Map2Cube());
	std::cerr << "       : "<<l2.shifts<<endl;
      }
      else std::cerr << "layer "<<lyrIdx<<"= identity\n";

      lyrIdx = layers.length()-lyrIdx;
    }
  }
}

void PermNetwork::ApplyToArray(HyperCube<long>& v)
{
}

void PermNetwork::ApplyToPtxt(ZZX& p, const EncryptedArray& ea)
{
}


void PermNetwork::ApplyToCtxt(Ctxt& c)
{
#if 0
  const PAlgebra& al = c.getContext().zMStar;
  EncryptedArray ea(c.getContext()); // use G(X)=X for this ez object

  // Apply the layers, one at a time
  for (long i=0; i<this->length(); i++) {
    PermNetLayer& lyr = layers[i];
    long g = PowerMod(al.ZmStarGen(lyr.genIdx), lyr.e, al.getM());

    Ctxt sum(c.getPubKey(), c.getPtxtSpace()); // an empty ciphertext
    for (long j=0; j<(long)lyr.shiftAmountUsed.size(); j++) {
      long shft = lyr.shiftAmountUsed[j];

      // Prepare the mask for this shift amount
      vector<long> mask = lyr.shifts;
      ZZX maskPoly;
      for (long k=0; k<(long) mask.size(); k++)
	mask[k] = (mask[k]==shft)? 1 : 0;
      ea.encode(maskPoly, mask);

      // Then mask and rotate
      Ctxt tmp = c;
      tmp.multByConstant(maskPoly);
      if (shft!=0) tmp.smartAutomorph(PowerMod(g, shft, al.getM()));
      if (j==0) sum = tmp;
      else sum += tmp;
    }
    c = sum; // Done with this layer, update the cipehrtext c
  }
#endif
}
