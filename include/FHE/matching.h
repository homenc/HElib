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
 * @file matching.h
 * @brief Classes and functions for max-flow in a generic graph
 **/
#ifndef _MATCHING_H_
#define _MATCHING_H_

#include "NumbTh.h"


//! An edge in a flow graph
class FlowEdge {
public:
  long capacity, flow;
  explicit FlowEdge(long c=0, long f=0){capacity=c; flow=f;}
};

typedef unordered_map<long,FlowEdge> FNeighborList;
     // FNeighborList[i]=edge-to-node-i

typedef vector<FNeighborList> FlowGraph;
     // FlowGraph[i][j] is the edge i->j

long maximum_flow(FlowGraph& fg, long src, long sink);
     // Use the Edmonds-Karp max-flow algorithm

//! @class LabeledEdge
//! @brief A generic directed edge in a graph with some labels
class LabeledEdge {
public:
  long from, to;
  long label;
  long color;
  LabeledEdge(long f, long t, long l=0, long c=0)
  {from=f; to=t; label=l; color=c;}
};

typedef unordered_multimap<long,LabeledEdge> LNeighborList;

//! @class LabeledVertex
//! @brief A generic node in a graph with some labels
class LabeledVertex {
public:
  long name;
  long label; // We don't really use the label, but why not..
  LNeighborList neighbors;
  explicit LabeledVertex(long n, long l=0) {name=n, label=l;}

  void addEdge(long nn, long l=0, long c=0) { // allow parallel edges
    neighbors.insert(pair<long,LabeledEdge>(nn,LabeledEdge(name,nn,l,c)));
  }
  void addNeighbor(long nn, long l=0, long c=0){ // dont insert a parallel edge
    if (neighbors.count(nn)==0)
      neighbors.insert(pair<long,LabeledEdge>(nn,LabeledEdge(name,nn,l,c)));
  }
};

//! A bipartite flow graph
class BipartitleGraph {
  // Construct a flow graph corresponding to this bipartite graph
  void buildFlowGraph(FlowGraph& fg);

public:
  vector<LabeledVertex> left; //  the right side is implicit

  void addEdge(long from, long to, long label, long color=0) {
    for (long sz = left.size(); sz <= from; sz++) // insert nodes if needed
      left.push_back(LabeledVertex(sz));
    left.at(from).addEdge(to, label, color); // insert the edge itself
  }

  // Partition the graph into maximum matchings, color by i=1,2,3,...
  // the edges of the i'th matching.
  void partitionToMatchings();
  void printout();
};

#endif /* ifdef _MATCHING_H_ */
