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
#ifndef _MULTIAUTOMORPH_H
#define _MULTIAUTOMORPH_H
#include <memory>
#include <vector>
#include <unordered_map>

/* The automorphism tree is represented by a list of the form:
 *
 *   fromNode1: list of toNodes
 *   fromNode2: list of toNodes
 *   [...]
 *
 * where each node is represented by its representative in Zm*, and
 * having an edge a->b in the tree implies that we have a key-switching
 * matrix for k = b * a^{-1} mod m. The root of the tree is always 1
 * (which in particular means that autGraph[1] must exist).
 *
 * The "high level" multiAutomorph implementation below gets as input a
 * ciphertext and a tree as above, and generates all the automorphisms
 * in a DFS order on this tree. Note that the order is well defined
 * by the order of the std:vector's, even though we use unordered_map
 * for the different vectors. (A cycle in the underlying graph may
 * cause the same rotation amount to be returned more than once.)
 */
class AutGraph: public std::unordered_map< long, std::vector<long> >
{
  void dfsOrder(std::vector<long>& vec, long from) const
  {
    for (long v : this->at(from)) {
      vec.push_back(v);
      if (this->count(v) > 0) // internal node
        dfsOrder(vec,v);
    }
  }
public:
  void getDFSorder(std::vector<long>& vec)
  {
    vec.clear();
    vec.push_back(1);
    dfsOrder(vec,1);
  }
};


// Interface #1: using call-backs:
//--------------------------------
// The AutomorphHandler class below is a virtual class for call-backs
// to get the output of the multiAutomorph implementation. This is called
// as ctxt.multiAutomorph(tree, handler), where tree is a tree as above
// and handler is derived from AutomorphHandler.
class Ctxt; // forward decleration
class AutomorphHandler {
public:
  virtual bool handle(std::unique_ptr<Ctxt>& ctxt, long amt) {return true;}
          // returns false to stop processing, true otherwise
  virtual ~AutomorphHandler() {}
};
/* Applications can derive from this class a handler that actually
 * does something with the "rotated" cipehrtexts. The AutomorphHandler
 * class above can also be used by itself as a do-nothing processor
 * for debugging, or to calculate the required automorphisms (see
 * automorphVals in numbTh.h)
 *
 * The application can take ownership of Ctxt object, e.g., by calling
 * ctxt.swap(...) or std::move(ctxt), and then the application must
 * make sure that it frees the object when it no longer needs it.
 */

void multiAutomorph(Ctxt& ctxt, const AutGraph& tree,
                    AutomorphHandler& handler);


// Interface #2: using an iterator:
//---------------------------------
// The AutoIterator::next(ctxt) method will be called repeatedly by the
// application, each time returning in ctxt the next automorphed ctxt
// and returning the automorphism amount in its return value as an
// elemnt in Zm*. Returns 0 When no more automorphisms are available.
class AutoIterator {
public:
  static AutoIterator* build(Ctxt& c, const AutGraph& tree); // factory
  virtual ~AutoIterator() {} // virtual destructor
  virtual long next(Ctxt& c) =0;
};
/* Applications will use code similar to this:
 *
 *     const AutGraph& tree = publicKey.getTree4dim(i);
 *     std::unique_ptr<AutoIterator> it(AutoIterator::build(ctxt,tree));
 *     Ctxt tmp(ZeroCtxtLike, ctxt);
 *     while (long val = it->next(tmp)) {
 *       ... do something with val and tmp ...
 *     }
 */
#endif // _MULTIAUTOMORPH_H
