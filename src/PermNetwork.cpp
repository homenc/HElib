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
#include <NTL/ZZ.h>
#include <helib/Context.h>
#include <helib/Ctxt.h>
#include <helib/permutations.h>
#include <helib/EncryptedArray.h>

namespace helib {

std::ostream& operator<<(std::ostream& s, const PermNetwork& net)
{
  s << "[";
  for (long i = 0; i < net.layers.length(); i++) {
    const PermNetLayer& lyr = net.layers[i];
    s << "[" << lyr.genIdx << " " << lyr.e << " " << lyr.isID << " "
      << lyr.shifts << "]\n";
  }
  return s << "]";
}

// Compute one or more layers corresponding to one network of one leaf
void PermNetwork::setLayers4Leaf(long lyrIdx,
                                 const ColPerm& p,
                                 const NTL::Vec<long>& benesLvls,
                                 long gIdx,
                                 const SubDimension& leafData,
                                 const Permut& map2cube)
{
#ifdef HELIB_DEBUG
  std::cerr << "Layer " << lyrIdx << ", column-permutation=" << p << std::endl;
#endif
  // Compute the shift amounts for all the layers in this network
  NTL::Vec<bool> isID;
  NTL::Vec<Permut> shifts;
  if (benesLvls.length() == 1) { // Special case for a "trivial" 1-layer network
    shifts.SetLength(1);
    isID.SetLength(1);
    isID[0] = !p.getShiftAmounts(shifts[0]);
  } else // The general case of a multi-layer Benes network
    p.getBenesShiftAmounts(shifts, isID, benesLvls);

  // Copy the shift amounts to the right place in the bigger network,
  // renaming the slots from a linear array to the hyper cube
  for (long i = 0; i < benesLvls.length(); i++) {
    PermNetLayer& lyr = layers[lyrIdx + i];
    lyr.genIdx = gIdx;
    lyr.isID = isID[i];
    lyr.e = leafData.e;
    if (!lyr.isID) {
#ifdef HELIB_DEBUG
      std::cerr << "layer " << lyrIdx + i << ": " << shifts[i] << std::endl;
#endif
      if (leafData.good) // For good leaves, shift by -x is the same as size-x
        for (long k = 0; k < shifts[i].length(); k++)
          if (shifts[i][k] < 0)
            shifts[i][k] += leafData.size;
      applyPermToVec(lyr.shifts, shifts[i], map2cube); // do the renaming
#ifdef HELIB_DEBUG
      std::cerr << "       : " << lyr.shifts << std::endl;
#endif
    }
    //    else std::cerr << "layer "<<lyrIdx+i<<"= identity\n";
  }
}

// Build a full permutation network
void PermNetwork::buildNetwork(const Permut& pi, const GeneratorTrees& trees)
{
  if (trees.numTrees() == 0) { // the identity permutation, nothing to do
    layers.SetLength(0);
    return;
  }

  NTL::Vec<long> dims;
  trees.getCubeSubDims(dims);

  //  std::cerr << "pi =      "<<pi<<std::endl;
  //  std::cerr << "map2cube ="<<trees.mapToCube()<<std::endl;
  //  std::cerr << "map2array="<<trees.mapToArray()<<std::endl;

  // Compute the permutation on the cube, rho = map2cube o pi o map2array

  Permut rho;
  applyPermsToVec(rho, trees.mapToCube(), pi, trees.mapToArray());
  //  std::cerr << "rho =     "<<rho<<std::endl;

  // Break rho along the different dimensions
  CubeSignature sig(dims); // make a cube-signature object
  std::vector<ColPerm> perms;
  breakPermByDim(perms, rho, sig);

  //  for (long i=0; i<(long)perms.size(); i++) { // debugging printouts
  //    Permut tmp;
  //    perms[i].makeExplicit(tmp);
  //    std::cerr << " prems["<<i<<"]="<<tmp<<std::endl;
  //  }

  layers.SetLength(trees.numLayers()); // allocate space

  // Go over the different permutations and build the corresponding layers
  long dimIdx = 0;
  long frntLyr = 0, backLyr = layers.length();
  // go over all the generators/trees
  for (long g = 0; g < trees.numTrees(); g++) {
    const OneGeneratorTree& T = trees[g];

    // In each tree, go over all the leaves
    for (long leaf = T.firstLeaf(); leaf >= 0; leaf = T.nextLeaf(leaf)) {
      const SubDimension& leafData = T[leaf].getData();

      // This leaf determines layers frntLyr...frntLyr+frst.length()-1, and
      // if it isn't the middle then also backLyr-scnd.length()...backLyr-1

      // handle the first Benes network
      setLayers4Leaf(/*1st-layer-index=*/frntLyr,
                     /*permutation    =*/perms[dimIdx],
                     /*Benes levels   =*/leafData.frstBenes,
                     /*generator index=*/T.getAuxKey(),
                     /*(size,good,e)  =*/leafData,
                     /*hypercube renaming permutation=*/trees.mapToCube());
      frntLyr += leafData.frstBenes.length(); // how many layers were used
      dimIdx++;

      if (leafData.scndBenes.length() > 0) {  // Also a second Benes network
        long dimIdx2 = perms.size() - dimIdx; // dimIdx was incremented above
        backLyr -= leafData.scndBenes.length();
        setLayers4Leaf(/*1st-layer-index=*/backLyr,
                       /*permutation    =*/perms[dimIdx2],
                       /*Benes levels   =*/leafData.scndBenes,
                       /*generator index=*/T.getAuxKey(),
                       /*(size,good,e)  =*/leafData,
                       /*hypercube renaming permutation=*/trees.mapToCube());
      }
    }
  }
}

// Apply a permutation network to a hypercube, used mostly for debugging
void PermNetwork::applyToCube(HyperCube<long>& cube) const
{
  if (layers.length() == 0)
    return;
  long n = cube.getSize();
  NTL::Vec<long> tmp(NTL::INIT_SIZE, n); // temporary vector

  // Apply the layers, one at a time
  for (long i = 0; i < layers.length(); i++) {
    const PermNetLayer& lyr = layers[i];
    if (lyr.isID)
      continue; // this layer is the identity permutation

    assertEq(lyr.shifts.length(), n, "layer has incorrect size");

    // This layer shift elements along the dimension lyr.genIdx
    long dim = lyr.genIdx;

    // Move elements as dictated by this layer
    for (long j = 0; j < n; j++) {
      long shamt = lyr.e * lyr.shifts[j]; // how much to shift this slot
      if (shamt < 0)
        shamt += cube.getDim(dim);            // addCoord expects shamt>=0
      long j2 = cube.addCoord(j, dim, shamt); // new index for this slot
      tmp[j2] = cube[j];
    }
    // Copy back to cube
    for (long j = 0; j < n; j++)
      cube[j] = tmp[j];
#ifdef HELIB_DEBUG
    std::cerr << " after layer " << i << ", cube=" << cube.getData()
              << std::endl;
#endif
  }
}

// TODO: Implement this method.
// void PermNetwork::applyToPtxt(NTL::ZZX& p, const EncryptedArray& ea) const
// {
//   throw LogicError("PermNetwork::applyToPtxt is not implemented");
// }

// Upon return, mask[i]=1 if haystack[i]=needle, 0 otherwise.
// Also set to 0 all the entries in haystack where mask[i]=1.
// Return the index of the first nonzero entry in haystack at the end
// of the pass (-1 if they are all zero). Also return a flag saying if
// any entries of the mask are nonzero.
static std::pair<long, bool> makeMask(std::vector<long>& mask,
                                      NTL::Vec<long>& haystack,
                                      long needle)
{
  bool found = false;
  long fstNonZeroIdx = -1;
  for (long i = 0; i < (long)mask.size(); i++) {
    if (haystack[i] == needle) { // found a needle
      found = true;
      mask[i] = 1;
      haystack[i] = 0; // remove this needle from haystack
    } else {           // no needle here
      mask[i] = 0;
      if (haystack[i] != 0 && fstNonZeroIdx < 0)
        fstNonZeroIdx = i; // first nonzero entry in haystack
    }
  }
  return std::make_pair(fstNonZeroIdx, found);
}

// Apply a permutation network to a ciphertext
void PermNetwork::applyToCtxt(Ctxt& c, const EncryptedArray& ea) const
{
  const PAlgebra& al = ea.getPAlgebra();

  // Apply the layers, one at a time
  for (long i = 0; i < layers.length(); i++) {
    const PermNetLayer& lyr = layers[i];
    if (lyr.isID)
      continue; // this layer is the identity permutation

    // This layer is shifted via powers of g^e mod m
    long g2e = NTL::PowerMod(al.ZmStarGen(lyr.genIdx), lyr.e, al.getM());

    NTL::Vec<long> unused = lyr.shifts;          // copy to a new vector
    std::vector<long> mask(lyr.shifts.length()); // buffer to hold masks
    Ctxt sum(c.getPubKey(), c.getPtxtSpace());   // an empty ciphertext

    long shamt = 0;
    bool frst = true;
    while (true) {
      std::pair<long, bool> ret = makeMask(mask, unused, shamt); // compute mask
      if (ret.second) { // non-empty mask
        Ctxt tmp = c;
        NTL::ZZX maskPoly;
        ea.encode(maskPoly, mask);    // encode mask as polynomial
        tmp.multByConstant(maskPoly); // multiply by mask
        if (shamt != 0)               // rotate if the shift amount is nonzero
          tmp.smartAutomorph(NTL::PowerMod(g2e, shamt, al.getM()));
        if (frst) {
          sum = tmp;
          frst = false;
        } else
          sum += tmp;
      }
      if (ret.first >= 0)
        shamt = unused[ret.first]; // next shift amount to use

      else
        break; // unused is all-zero, done with this layer
    }
    c = sum; // update the ciphertext c before the next layer
  }
}

} // namespace helib
