#ifndef PAIRPATH_GRAPH_H
#define PAIRPATH_GRAPH_H


#include<iostream>
#include<algorithm>
#include<map>
#include<vector>

#include <boost/unordered_map.hpp>
//#include"Find_junc_last_map.h"

using namespace std;

// pairpath graph
typedef pair<int,int> edge_t;
 bool sorter(const pair<edge_t,int>& p1, const pair<edge_t,int>& p2)
  {
        return( (p1.second > p2.second) ||
                 (p1.second == p2.second && p1.first < p2.first) );
  }
  bool sorter_(const pair<int,double>& p1, const pair<int,double>& p2)
  {
	return( (p1.second > p2.second) ||
		(p1.second == p2.second && p1.first < p2.first) );
  }

class PairPath_Graph
{

public:
  typedef int node_idx_t;
  map< path_t, vector<int> > prefix_subpath_map;
  map< path_t, vector<int> > suffix_subpath_map;

  typedef map<path_t, vector<int> >::iterator iter;

  typedef pair<int,int> edge_t;

  class Node   
  {
    public:
	path_t node_ppath;
	double node_coverage;
	double raw_node_coverage;
        vector< pair<node_idx_t,int> > children;
	vector< pair<node_idx_t,int> > parents;
    public:

	Node(path_t path,double cov):node_ppath(path),node_coverage(cov) {}

	void add_child(node_idx_t child, int cov)
 	{
	    pair<node_idx_t,int> p = make_pair(child,cov);
	    children.push_back(p);
	}
	void add_parent(node_idx_t parent, int cov)
	{
	    pair<node_idx_t,int> p = make_pair(parent,cov);
	    parents.push_back(p);
	}
  };
public: vector<path_t> Ppath;
	vector<Node> node_set;
	vector<path_t> final_pair_path;

	vector< pair<edge_t,int> > edge_with_descend_cov;
	vector< pair<int,double> > node_with_descend_cov;
public:

  PairPath_Graph(vector<path_t>& pair_path)
  {
	Ppath = pair_path;
  }
  void subpath_hash()
  {
    for(int i=0;i<Ppath.size();i++)
    {
        path_t path = Ppath[i];
        path.pop_back();
        while(!path.empty())
        {
            if(prefix_subpath_map.find( path) != prefix_subpath_map.end())
                prefix_subpath_map[path].push_back(i);
            else
            {
                vector<int> temp(1,i);
                prefix_subpath_map[path] = temp;
            }
            path.pop_back();
        }
        path = Ppath[i];
        path.erase(path.begin());
        while(!path.empty())
        {
            if(suffix_subpath_map.find(path) != suffix_subpath_map.end())
                suffix_subpath_map[path].push_back(i);
            else
            {
                vector<int> temp(1,i);
                suffix_subpath_map[path] = temp;
            }
            path.erase(path.begin());
        }
    }
  }
  void show_graph()
  {
    cerr<<"Pair-Path Graph ..."<<endl;
    cerr<<"** Edge:"<<endl;
    for(int i = 0;i<node_set.size();i++)
    {
	for(int j=0;j<node_set[i].children.size();j++) 
	{
	    cerr<<i<<"->"<<node_set[i].children[j].first<<": "<<node_set[i].children[j].second<<endl;
	}
    }
    cerr<<"** Node:"<<endl;
    for(int i = 0;i<node_set.size();i++)
    {
	cerr<<i<<" ";
	for(int j=0;j<node_set[i].node_ppath.size();j++) cerr<<node_set[i].node_ppath[j]<<"->";
	cerr<<"  "<<node_set[i].node_coverage;
	cerr<<endl;
    }
    cerr<<"-------------------------------------------------"<<endl;
  }
  double get_minimum_cov(path_t path, vector<int> junc_l,vector<int> junc_r,vector<double> junc_cov,
			vector<int> exon_l,vector<int> exon_r)
  {
    double cov = 1000000;
    for(int i = 0;i<path.size()-1;i++)
    {
	int jl = exon_r[ path[i] ] + 1;
	int jr = exon_l[ path[i+1] ];

	for(int j =0;j<junc_l.size();j++)
	{
	    if(jl == junc_l[j] && jr == junc_r[j])
	    {
		if(cov > junc_cov[j])
		    cov = junc_cov[j];
		continue;
	    }
	}
    }
    return cov;
  }
  void build_pairpath_graph(vector<int> junc_l,vector<int> junc_r,vector<double> junc_cov,
			    vector<int> exon_l,vector<int> exon_r)
  {
    subpath_hash();

    for(int i=0;i<Ppath.size();i++)
    {

	path_t temp = Ppath[i];
	double cov = get_minimum_cov(temp,junc_l,junc_r,junc_cov,exon_l,exon_r);

	Node node(temp,cov);
	node.raw_node_coverage = cov;
	pair<int,double> p = make_pair(i,cov);
	node_with_descend_cov.push_back(p);
	//add children to current path
	temp.erase(temp.begin());
	//while( !temp.empty() )
	while(temp.size() > 1)
	{
	    int cov = temp.size();

	    iter it = prefix_subpath_map.find(temp);

	    if(it != prefix_subpath_map.end())
	    {
		vector<int> nodes = it -> second;
		for(int j = 0;j<nodes.size();j++)
		{
		   node.add_child(nodes[j],cov);

		  edge_t edge = make_pair(i,nodes[j]);
		
		  pair<edge_t,int> p = make_pair(edge,cov); 
		  edge_with_descend_cov.push_back(p);
		}
	    }

	    temp.erase(temp.begin());
	}

	//add parents to current path
	temp = Ppath[i];
	temp.pop_back();
	//while( !temp.empty() )
	while(temp.size() > 1)
	{
	    int cov = temp.size();
	    iter it = suffix_subpath_map.find(temp);
	
	    if(it != suffix_subpath_map.end())
	    {
		vector<int> nodes = it -> second;
	 	for(int j = 0;j<nodes.size();j++) node.add_parent(nodes[j],cov);
	    }

	    temp.pop_back();
	}

	node_set.push_back(node);
    }
    //show_graph();
  }

  bool get_reverse_extend_node(int check_node, int& extend_node, vector<bool>& node_used_flag)
  {
    bool unused_flag = false;//record if there exists unused edge;
    bool used_flag = false;

    int used_maxcov_node,unused_maxcov_node;
    double used_maxcov = 0,unused_maxcov = 0;
    vector< pair<node_idx_t,int> > parents = node_set[check_node].parents;

    for(int i = 0;i < parents.size();i++)
    {
	//edge_t edge = make_pair(parents[i].first,check_node);
	//if(edge_used_flag.find(edge) == edge_used_flag.end()) // NOT used
	/*
	if( !node_used_flag[parents[i].first] ) // NOT used
	{
	    unused_flag = true;
	    if(unused_maxcov < node_set[parents[i].first].node_coverage)
	    {
		unused_maxcov = node_set[parents[i].first].node_coverage;
		unused_maxcov_node = parents[i].first;
	    }
	}
	else//used
	{
	    used_flag = true;
	    if(used_maxcov < node_set[parents[i].first].node_coverage)
	    {
		used_maxcov = node_set[parents[i].first].node_coverage;
		used_maxcov_node = parents[i].first;
	    }
	}
	*/
	if(node_set[parents[i].first].node_coverage == 0) continue;
	if(unused_maxcov < node_set[parents[i].first].node_coverage)
	{
	    unused_flag = true;
	    unused_maxcov = node_set[parents[i].first].node_coverage;
	    unused_maxcov_node = parents[i].first;
	}
    }
    if(unused_flag) extend_node = unused_maxcov_node;
    else return false;
/*
    if(unused_flag)
	extend_node = unused_maxcov_node;
    else if(used_flag)
	extend_node = used_maxcov_node;
    else return false;
*/
    return true;

  }
  bool get_forward_extend_node(int check_node, int& extend_node, vector<bool>& node_used_flag)
  {
    bool unused_flag = false;//record if there exists unused edge;
    bool used_flag = false;

    int used_maxcov_node,unused_maxcov_node;
    double used_maxcov = 0,unused_maxcov = 0;
    vector< pair<node_idx_t,int> > children = node_set[check_node].children;

    for(int i = 0;i < children.size();i++)
    {
	//edge_t edge = make_pair(check_node,children[i].first);
	//if(edge_used_flag.find(edge) == edge_used_flag.end()) // NOT used
	/*
	if( !node_used_flag[children[i].first] ) // NOT used
	{
	    unused_flag = true;
	    if(unused_maxcov < node_set[children[i].first].node_coverage)
	    {
		unused_maxcov = node_set[children[i].first].node_coverage;
		unused_maxcov_node = children[i].first;
	    }
	}
	else//used
	{
	    used_flag = true;
	    if(used_maxcov < node_set[children[i].first].node_coverage)
	    {
		used_maxcov = node_set[children[i].first].node_coverage;
		used_maxcov_node = children[i].first;
	    }
	}
	*/
	if(node_set[children[i].first].node_coverage == 0) continue;
	if(unused_maxcov < node_set[children[i].first].node_coverage)
	{
	    unused_flag = true;
	    unused_maxcov = node_set[children[i].first].node_coverage;
	    unused_maxcov_node = children[i].first;
	}
    }
    if(unused_flag) extend_node = unused_maxcov_node;
    else return false;
/*
    if(unused_flag) 
	extend_node = unused_maxcov_node;
    else if(used_flag)
	extend_node = used_maxcov_node;
    else return false;
*/
    return true;
  }
  void reverse_extend(path_t& path_in_ppgraph, path_t& current_extended_pair_path, 
		      vector<bool>& node_used_flag)
  {
     int check_node = path_in_ppgraph.front();
     int extend_node = -1;

     while(get_reverse_extend_node(check_node,extend_node,node_used_flag))
     {
	//edge_t edge = make_pair(extend_node,check_node);
	//edge_used_flag[edge] = true;

	node_used_flag[extend_node] = true;

	path_in_ppgraph.insert(path_in_ppgraph.begin(),extend_node);

	path_t path_ = node_set[extend_node].node_ppath;
	path_t::iterator it = find(path_.begin(),path_.end(),current_extended_pair_path.front());
	
	it--;

	for(;it != path_.begin();it--) 
	   current_extended_pair_path.insert(current_extended_pair_path.begin(),*it);

	current_extended_pair_path.insert(current_extended_pair_path.begin(),path_.front());

	check_node = path_in_ppgraph.front();
	extend_node = -1;
     }
     return;
  }

  void forward_extend(path_t& path_in_ppgraph, path_t& current_extended_pair_path,
		      vector<bool>& node_used_flag)
  {
    int check_node = path_in_ppgraph.back();
    int extend_node = -1;
 
    while(get_forward_extend_node(check_node,extend_node,node_used_flag))
    {
	//edge_t edge = make_pair(check_node,extend_node);
	//edge_used_flag[edge] = true;

	node_used_flag[extend_node] = true;

	path_in_ppgraph.push_back(extend_node);

	path_t path_ = node_set[extend_node].node_ppath;
	path_t::iterator it = find(path_.begin(),path_.end(),current_extended_pair_path.back());
	it++;
	for(;it != path_.end();it++) 

	    current_extended_pair_path.push_back(*it);

	check_node = path_in_ppgraph.back();
	extend_node = -1;
    }
    return;
  }

  void modify_coverage(path_t path) //path in pairpath graph
  {
    double min_cov = 1000000;
    for(int i=0;i<path.size();i++)
    {
	int node = path[i];
	if(node_set[node].node_coverage <min_cov) min_cov = node_set[node].node_coverage;
    }

    for(int i=0;i<path.size();i++)
    {
	int node = path[i];
	//if( double(node_set[node].node_coverage - min_cov) / node_set[node].raw_node_coverage <= 0.1) node_set[node].node_coverage = 0;
	//else node_set[node].node_coverage = node_set[node].node_coverage - min_cov;
	//
	node_set[node].node_coverage = node_set[node].node_coverage - min_cov;
    }
    return;
  }
  int get_maxcov_node()
  {
    double max_cov = 0;
    int maxcov_node = -1;
    for(int i=0;i<node_set.size();i++)
    {
	if(max_cov <= node_set[i].node_coverage)
	{
	    max_cov = node_set[i].node_coverage;
	    maxcov_node = i;
	}
    }
    return maxcov_node;
  }
  void search_path()
  {
/*
    for(int i=0;i<node_set.size();i++)
    {
	if(node_set[i].parents.empty() && node_set[i].children.empty())
	{
	    final_pair_path.push_back(node_set[i].node_ppath);
	}
    }
*/
    //sort(edge_with_descend_cov.begin(),edge_with_descend_cov.end(),sorter);
    //sort(node_with_descend_cov.begin(),node_with_descend_cov.end(),sorter_);

    vector<bool> node_used_flag(node_set.size(),false);

    map<edge_t,bool> edge_used_flag;

    //for(int i = 0;i<node_with_descend_cov.size();i++)
    if(node_set.empty()) return;

    while(1)
    {
	//sort(node_with_descend_cov.begin(),node_with_descend_cov.end(),sorter_);
	//cerr<<"node: "<<node_with_descend_cov.front().first<<" -- cov:"<<node_with_descend_cov.front().second<<endl;
	//if(node_with_descend_cov.front().second == 0) break;

	//int node = node_with_descend_cov.front().first;

	int node = get_maxcov_node();

	//cerr<<"node: "<<node<<" current_cov: "<<node_set[node].node_coverage<<" raw_cov: "<<node_set[node].raw_node_coverage<<endl;

	if(node_set[node].node_coverage == 0) break;
	
        path_t path_in_ppgraph; //path in pair path graph;
	path_t current_extended_pair_path;//path in RAW SPLICING GRAPH;

	path_in_ppgraph.push_back(node);

	current_extended_pair_path = Ppath[ node ];

	reverse_extend(path_in_ppgraph,current_extended_pair_path,node_used_flag);
	forward_extend(path_in_ppgraph,current_extended_pair_path,node_used_flag);

	modify_coverage(path_in_ppgraph);

	final_pair_path.push_back(current_extended_pair_path);
/*
	cerr<<"  current_extended_pair_path: "<<endl;
 	cerr<<"  ";
 	for(int j=0;j<current_extended_pair_path.size();j++) cerr<<current_extended_pair_path[j]<<"==>>";
	cerr<<endl;
	cerr<<"  ";
 	for(int j = 0;j<path_in_ppgraph.size();j++) cerr<<path_in_ppgraph[j]<<"->";
	cerr<<endl;
*/	
	node_used_flag[node] = true;

    }

    return;
  }
};

#endif
