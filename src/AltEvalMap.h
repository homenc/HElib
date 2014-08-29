#ifndef _ALT_EVAL_MAP_H_
#define _ALT_EVAL_MAP_H_

#include "FHE.h"
#include "Ctxt.h"
#include "EncryptedArray.h"
#include "permutations.h"



class AltEvalMap {
private:
  const EncryptedArray& ea;
  bool invert; 

  long nfactors;

  shared_ptr<PlaintextBlockMatrixBaseInterface> mat1;

  Vec< shared_ptr<PlaintextMatrixBaseInterface> > matvec;

  
public:

  AltEvalMap(const EncryptedArray& _ea, const Vec<long>& mvec, bool _invert);

  void apply(Ctxt& ctxt) const;
};




#endif
