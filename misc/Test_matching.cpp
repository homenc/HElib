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
#include "matching.h"
#include <queue>

/*********** Simple networks to test the code from matching.cpp ***********/
/******* a directed flow graph *******
         /-------3-\   /-------3-\
        /           \ /           \
A -3-> B -2-> C -3-> D -4-> E -3-> F
 \           / \           /
  \-------3-/   \-------2-/
**************************************/
int main()
{
  // initializing in C++ is such a pain
  FlowGraph fg;
  { // Node 0 with edges to 1 (capacity 3) and 2 (capacity 3)
    FNeighborList l;
    l.insert(pair<long,FlowEdge>(1,FlowEdge(3)));
    l.insert(pair<long,FlowEdge>(2,FlowEdge(3)));
    fg.push_back(l);
  }
  { // Node 1 with edges to 2 (capacity 2) and 3 (capacity 3)
    FNeighborList l;
    l.insert(pair<long,FlowEdge>(2,FlowEdge(2)));
    l.insert(pair<long,FlowEdge>(3,FlowEdge(3)));
    fg.push_back(l);
  }
  { // Node 2 with edges to 3 (capacity 3) and 4 (capacity 2)
    FNeighborList l;
    l.insert(pair<long,FlowEdge>(3,FlowEdge(3)));
    l.insert(pair<long,FlowEdge>(4,FlowEdge(2)));
    fg.push_back(l);
  }
  { // Node 3 with edges to 4 (capacity 4) and 5 (capacity 3)
    FNeighborList l;
    l.insert(pair<long,FlowEdge>(4,FlowEdge(4)));
    l.insert(pair<long,FlowEdge>(5,FlowEdge(3)));
    fg.push_back(l);
  }
  { // Node 4 with an edge to 5 (capacity 3)
    FNeighborList l;
    l.insert(pair<long,FlowEdge>(5,FlowEdge(3)));
    fg.push_back(l);
  }

  long f = maximum_flow(fg, 0, 5); // compute max-flow
  cout << "Max-flow value = " << f << endl;
  printFlow(fg);

/*********** an undirected bipartite graph ***********
 0 -> 1,3,4    2 -> 0,2,3    4 -> 1,3,4
 1 -> 0,1,4    3 -> 0,2,2
*****************************************************/
  BipartitleGraph bg;
  bg.addEdge(0,1,1);  bg.addEdge(0,3,2);  bg.addEdge(0,4,3);
  bg.addEdge(1,0,4);  bg.addEdge(1,1,5);  bg.addEdge(1,4,6);
  bg.addEdge(2,0,7);  bg.addEdge(2,2,8);  bg.addEdge(2,3,9);
  bg.addEdge(3,0,10);  bg.addEdge(3,2,11);  bg.addEdge(3,2,12);
  bg.addEdge(4,1,13);  bg.addEdge(4,3,14);  bg.addEdge(4,4,15);
  bg.partitionToMatchings();
  bg.printout();
}
/*********************************************************************/
