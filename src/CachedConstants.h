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
#ifndef FHE_CACHED_CONST_H_
#define FHE_CACHED_CONST_H_
/**
 * @file CachedConstants.h
 * @brief Cached constants for linear transformations
 */

// SH: This is written using "old style" coding, just because I can
// never remember how to do it using inheritance.

class CachedConstants {
public:
  enum CacheTag { tagEmpty, tagZero, tagZZX, tagDCRT };

  struct Zero_type {}; // used to select a constructor
  static Zero_type Zero;

private:
  struct PolyPtr {
    CacheTag tag;
    std::unique_ptr<NTL::ZZX> zzx;
    std::unique_ptr<DoubleCRT> dcrt;

    PolyPtr(): tag(tagEmpty) {}
    PolyPtr(Zero_type zero): tag(tagZero) {}
    PolyPtr(NTL::ZZX* zzxPtr): tag(tagZZX), zzx(zzxPtr) {}
    PolyPtr(DoubleCRT* dcrtPtr): tag(tagDCRT), dcrt(dcrtPtr) {}
  };


  std::vector< PolyPtr > polys;

public:
  long size() const { return polys.size(); }
  void resize(long sz) { polys.resize(sz); }
  void clear() { polys.clear(); }

  bool isZero(int i) const {

<<<<<<< HEAD
<<<<<<< 32300783d9179c5fe7b2982d743180ab99546561
<<<<<<< e657e294e21b2b33e4a1c1a8c8a91ddcb953dc92
=======
<<<<<<< febe625a45df829493261ac83d2abdb14e3fff05
>>>>>>> .
=======
<<<<<<< febe625a45df829493261ac83d2abdb14e3fff05
>>>>>>> 16d4870b680f4ba4b493467999cc848e6b1a4696
=======
#if 0
    // This is not thread safe

    if (polys.at(i).tag == tagZZX) {
      if ( NTL::IsZero(*(polys[i].zzx)) ) {
	polys[i].tag = tagZero;
	polys[i].zzx.reset();
      }
    }
#endif

>>>>>>> .
    // We may want to do the same for DoubleCRT, but checking for zero
    // would be more expensive
    return (polys[i].tag == tagZero);
  }

  void setZero(long i) {
    polys.at(i).tag = tagZero;
    polys[i].zzx.reset();
    polys[i].dcrt.reset();
  }

  bool isEmpty(int i) const { return (polys.at(i).tag == tagEmpty); }
  bool isZZX(int i) const { return (polys.at(i).tag == tagZZX); }
  bool isDCRT(int i) const { return (polys.at(i).tag == tagDCRT); }

  NTL::ZZX* const setAt(long i, NTL::ZZX* zzxPtr) {
    if (polys.at(i).tag != tagEmpty && polys[i].tag != tagZZX)
      throw std::logic_error("Cannot assign ZZX pointer to another cache type");
    polys[i].zzx.reset(zzxPtr);
    polys[i].tag = tagZZX;
    polys[i].dcrt.reset(); // just for good measure, this should not matter
    return zzxPtr;
  }

  NTL::ZZX* const getZZX(long i) const {
    if (polys.at(i).tag != tagEmpty && polys[i].tag != tagZZX)
      throw std::logic_error("Cannot return ZZX pointer from another cache type");
    return polys[i].zzx.get();
  }

  DoubleCRT* const setAt(long i, DoubleCRT* const dcrtPtr) {
    if (polys.at(i).tag == tagZero)
      throw std::logic_error("Cannot assign DCRT pointer to cache type zero");
    polys[i].dcrt.reset(dcrtPtr);
    polys[i].tag = tagDCRT;
    polys[i].zzx.reset(); // if "upgraded" from ZZX to DoubleCRT
    return dcrtPtr;
  }

  DoubleCRT* const getDCRT(long i) const {
    if (polys.at(i).tag != tagEmpty && polys[i].tag != tagDCRT)
      throw std::logic_error("Cannot return DCRT pointer from another cache type");
    return polys.at(i).dcrt.get();
  }
};




inline bool
allConstsAvailable(const CachedConstants& cache, long idx, long blockSize)
{
  long offset = idx*blockSize;
  for (long i=0; i<blockSize; i++)
    if (cache.isEmpty(i+offset)) return false;
  return true;
}

inline bool
allConstsZero(const CachedConstants& cache, long idx, long blockSize) 
{
  long offset = idx*blockSize;
  for (long i=0; i<blockSize; i++)
    if (!cache.isZero(i+offset)) return false;
  return true;
}


#endif // FHE_CACHED_CONST_H_
