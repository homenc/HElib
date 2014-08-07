#ifndef _EVAL_MAP_H_
#define _EVAL_MAP_H_

#include "FHE.h"
#include "Ctxt.h"
#include "EncryptedArray.h"
#include "permutations.h"

class Step2aShuffleBase;
class TowerBase;



class EvalMap {
private:
  const EncryptedArray& ea;
  bool invert; 

  bool easy;  // easy => d1 == d, 
              // !ease => d1 != d (but we d1 * d2 == d)

  long nfactors;

  shared_ptr<PlaintextBlockMatrixBaseInterface> mat1;
    // use for both easy and !easy

  shared_ptr<Step2aShuffleBase> shuffle;
  shared_ptr<TowerBase> tower;
    // use only in the !easy case

  Vec< shared_ptr<PlaintextMatrixBaseInterface> > matvec;

  shared_ptr<PermNetwork> net;
  Vec<long> slot_rotate;
    // used for the initial/final inter- and intra-slot rotations

  
public:

  // mvec: the factorization of m
  // width: a bound on the width used for the permutation network
  // invert: false => forward eval [ coeffs packed in slots to coeffs ]
  //         true  => reverse eval [ coeffs to coeffs packed in slots ]

  // so for bootstrapping, we want the reverse eval, followed by digit
  // extraction, followed by forward eval

  EvalMap(const EncryptedArray& _ea, const Vec<long>& mvec, long width, 
          bool _invert);

  void apply(Ctxt& ctxt) const;
};




#endif
