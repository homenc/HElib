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
/**
 * @file matching.h
 * @brief Classes and functions for max-flow in a generic graph
 **/
#ifndef HELIB_MATCHING_H
#define HELIB_MATCHING_H

#include <unordered_map>
#include <helib/NumbTh.h>

namespace helib {

//! An edge in a flow graph
class FlowEdge
{
public:
  long capacity, flow;
  explicit FlowEdge(long c = 0, long f = 0)
  {
    capacity = c;
    flow = f;
  }
};

typedef std::unordered_map<long, FlowEdge> FNeighborList;
// FNeighborList[i]=edge-to-node-i

typedef std::vector<FNeighborList> FlowGraph;
// FlowGraph[i][j] is the edge i->j

long maximum_flow(FlowGraph& fg, long src, long sink);
// Use the Edmonds-Karp max-flow algorithm

//! @class LabeledEdge
//! @brief A generic directed edge in a graph with some labels
class LabeledEdge
{
public:
  long from, to;
  long label;
  long color;
  LabeledEdge(long f, long t, long l = 0, long c = 0)
  {
    from = f;
    to = t;
    label = l;
    color = c;
  }
};

typedef std::unordered_multimap<long, LabeledEdge> LNeighborList;

//! @class LabeledVertex
//! @brief A generic node in a graph with some labels
class LabeledVertex
{
public:
  long name;
  long label; // We don't really use the label, but why not..
  LNeighborList neighbors;
  explicit LabeledVertex(long n, long l = 0) { name = n, label = l; }

  void addEdge(long nn, long l = 0, long c = 0)
  { // allow parallel edges
    neighbors.insert(
        std::pair<long, LabeledEdge>(nn, LabeledEdge(name, nn, l, c)));
  }
  void addNeighbor(long nn, long l = 0, long c = 0)
  { // dont insert a parallel edge
    if (neighbors.count(nn) == 0)
      neighbors.insert(
          std::pair<long, LabeledEdge>(nn, LabeledEdge(name, nn, l, c)));
  }
};

//! A bipartite flow graph
class BipartitleGraph
{
  // Construct a flow graph corresponding to this bipartite graph
  void buildFlowGraph(FlowGraph& fg);

public:
  std::vector<LabeledVertex> left; //  the right side is implicit

  void addEdge(long from, long to, long label, long color = 0)
  {
    for (long sz = left.size(); sz <= from; sz++) // insert nodes if needed
      left.push_back(LabeledVertex(sz));
    left.at(from).addEdge(to, label, color); // insert the edge itself
  }

  // Partition the graph into maximum matchings, color by i=1,2,3,...
  // the edges of the i'th matching.
  void partitionToMatchings();
  void printout();
};

} // namespace helib

#endif // ifndef HELIB_MATCHING_H
