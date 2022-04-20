#ifndef SIMPLFYGRAPH_H
#define SIMPLFYGRAPH_H

#include <iostream> 
#include <map>
#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <sstream>
using namespace std;

typedef int node_idx_t;
typedef vector<node_idx_t> path_t;
typedef pair<node_idx_t,node_idx_t> edge_t;
typedef pair<edge_t,double> edge_with_cov_t;


class SimplifyGraph{
public:
  
   class Node{
     public:
	double coverage_;
	int length_;
	vector<int> sequence;
	vector< pair<node_idx_t,double> > children;
	vector< pair<node_idx_t,double> > parents;


	path_t path;//some node merge to one
     public:
  	void add_child(node_idx_t child,double weight){

	  pair<node_idx_t,double> p = make_pair(child,weight);
	  this->children.push_back(p);
	  return;
	}
	void add_parent(node_idx_t parent,double weight) 
	{
	  pair<node_idx_t,double> p = make_pair(parent,weight);
          this->parents.push_back(p);
          return;
	}
	void delete_parent(node_idx_t n)
	{
	    for(size_t i=0;i<parents.size();)
            {
                if(n == parents[i].first) parents.erase(parents.begin() + i);
                else i++;
            }
            return;
	}
	void delete_child(node_idx_t n)
	{
	    for(size_t i=0;i<children.size();)
	    {
		if(n == children[i].first) children.erase(children.begin() + i);
		else i++;
	    }
	    return;
	}
	void add_path_node(node_idx_t node) {
	  path.push_back(node);
	  return;
	}
	double get_child_coverage(node_idx_t child)
     	{
	    for(size_t i=0;i<children.size();i++)
	    {
		if(child == children[i].first) return children[i].second;
	    }
	    return -1;
    	}
	string path_str()
	{
	    if(path.empty()) return "*; ";
	    stringstream ss;
	    for(size_t i=0;i<path.size();i++) ss<<path[i]<<"-";
	    return ss.str();
	}
	string str()
	{
	    if(sequence.empty()) return "*; ";
	    stringstream ss;
	  
	    for(size_t i=0;i<sequence.size()-1;){
		ss<<sequence[i]<<" "<<sequence[i+1]<<"; ";
		i = i+2;
	    }
	    return ss.str();
	}
	void addseq(vector<int> vec)
	{
	    for(size_t i=0;i<vec.size();i++) sequence.push_back(vec[i]);
	    return;
	}
	double coverage()
	{
	    return coverage_;
	}
	int length()
	{
	    return length_;
	}
	void clear()
	{
	    children.clear();
	    parents.clear();
	    path.clear();
	    sequence.clear();
	    coverage_ = 0;
	    length_ = -1;
	}
	
   };
public:
   vector<Node> node_set;
   string Strand, Chr;
   int Graph_size;
   int RawSize;
   int size_;
   map<node_idx_t, node_idx_t> rawNode_eNode_map;

   map<edge_t,edge_t > pair_Left_to_Right_edge_map;

   map<edge_t,edge_t > pair_Right_to_Left_edge_map;

   map<edge_t,vector<edge_with_cov_t> > Left_to_MultiRight_edge_map;
   map<edge_t,vector<edge_with_cov_t> > Right_to_MultiLeft_edge_map; 
   typedef map<edge_t,edge_t >::iterator iter;
   vector< pair<edge_t,double> > edge_descendcov;

   map<edge_t,double> EdgeCov_map;
   vector< edge_t >  reserved_junc;

public:
  void get_reserved_junc(vector<path_t> PairPath);
  void get_CovInfo();
  void outread_empirical();
  void build_node_set_and_simplify(vector<int>& exon_l, vector<int>& exon_r, vector<double>& exon_cov,
		      vector<node_idx_t>& ,vector<node_idx_t>&, vector<double>&, 
		      vector<int>& Single_node,vector<path_t>& PairPath,vector<int>& invalid_nodes);

  void remove_all_partial(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  bool keep_edge(edge_t e);


  void remove_edges_by_average_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_unbalance_edges( vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_partial_junction(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);



  void remove_small_exons(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  bool partial(node_idx_t n1, node_idx_t n2);
  void remove_intron_contamination(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_partial_end_by_edge_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);

  void remove_edges_onemulti(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_lowcov_edges_of_bifurcation_nodes(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);

  void delete_children_edges(vector<edge_t> );
  void delete_parents_edges(vector<edge_t> );
  
  void show_graph(); 
  void contract_graph();
  void get_graph(vector<vector<int> >& vecEdges,vector<double>& vecCov,vector<path_t>& PairPath, vector<double>& PairPath_cov);

  double InCov(node_idx_t);
  double OutCov(node_idx_t);

  double average_coverage();
};

#endif
