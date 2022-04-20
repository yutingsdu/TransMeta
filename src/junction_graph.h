#ifndef JUNCTIONGRAPH_H
#define JUNCTIONGRAPH_H

#include <iostream>
#include "simplify_graph.h"

using namespace std;
typedef vector<node_idx_t> triplet_t;
class Combing:public SimplifyGraph
{
private:
  int h;
  map<triplet_t,double> triplet_map;
  map<node_idx_t,vector<triplet_t> > node_triplet_map;
  map<pair<edge_t,edge_t>,int> packing_map;
public:
  //From  pathsearch/get_junction_graph.h
  vector<int> Edges_left;// left node of an edge in junction graph
  vector<int> Edges_right;// right node of an edge in junction graph
  vector<double> Weights;//weights of an edge in junction graph
public:
  Combing(){}
  void get_triplet_map(vector<path_t>& PairPath, vector<double>& PairPath_cov);
  void packing();
  void packing(node_idx_t n, vector< pair<node_idx_t,double> > children, vector< pair<node_idx_t,double> > parents);
  //check current_triplets is ok or not
  bool check_triplet(vector<triplet_t> current_triplets, vector<triplet_t> pair_triplets,
		     vector< pair<node_idx_t,double> > children,vector< pair<node_idx_t,double> > parents);
  double solve_qp(vector<triplet_t>& current_triplets,
            vector< pair<node_idx_t,double> > children,vector< pair<node_idx_t,double> >parents);
  void trivial_packing(node_idx_t n, vector< pair<node_idx_t,double> > children,  vector< pair<node_idx_t,double> > parents);

  void show_graph();
  void get_packing_result(vector<vector<int> >Vec_edges);

  void show_triplets(vector<triplet_t> vT);
  void show_triplet(triplet_t T);
};

#endif
