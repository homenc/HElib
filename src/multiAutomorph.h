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

class Ctxt; // forward decleration

// A virtual class to handle call-backs to get the output of the "low
// level" multiAutomorph implementation. This "low level" is called as
// ctxt.multiAutomorph(vec, handler), where vec is a vector of values
// to use for automorphism and handler is derived from AutomorphHandler.
class AutomorphHandler {
public:
  virtual bool handle(std::unique_ptr<Ctxt>& ctxt, long amt) {return true;}
          // returns false to stop processing, true otherwise
  virtual ~AutomorphHandler() {}
};
/* Applications will derive from this class a handler that actually
 * does something with the "rotated" cipehrtexts. The AutomorphHandler
 * class above can also be used by itself as a do-nothing processor
 * for debugging, or to calculate the required automorphisms (see
 * automorphVals in numbTh.h)
 *
 * Upon calling handle, the multiAutomorph routine has ownership of
 * Ctxt object. The application can take ownership, e.g., by calling
 * ctxt.swap(...) or std::move(ctxt), and then the application must
 * make sure that it frees the object when it no longer needs it.
 */


typedef std::unordered_map< long, std::vector<long> > AutGraph;
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
void
multiAutomorph(Ctxt& ctxt, const AutGraph& tree, AutomorphHandler& handler);

#endif // _MULTIAUTOMORPH_H
