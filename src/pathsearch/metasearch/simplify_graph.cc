
#include <iostream> 
#include <map>
#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>
#include "simplify_graph.h"

using namespace std;
extern double SampleSize;
extern bool unstranded;
extern double SEED_Filter;
extern string mode;//R|I|U  reference/ipac/combine
extern int pack_graph_num;
extern int unpack_graph_num;
extern bool PackingFlag;
extern bool MyFlag;
extern double AVERAGE_REMOVE_RATE;
extern double UNBALANCE_RATE ;
extern bool SFlag;
extern int rg_index;
extern int trans_id;
extern string out_name;
extern ofstream out_gtf;
extern ofstream out_info;

typedef vector<int> triple_t;
bool sorter(pair<edge_t,double> p1, pair<edge_t,double>p2)
{
    return ((p1.second>p2.second) || (p1.second == p2.second && p1.first < p2.first));
}
  void SimplifyGraph::get_reserved_junc(vector<path_t> PairPath)
  {
	for(size_t i=0;i<PairPath.size();i++)
	{
	    path_t p = PairPath[i];
	    for(size_t j=0;j<p.size()-1;j++)
	    {
		edge_t e = make_pair(p[j],p[j+1]);
		reserved_junc.push_back(e);
	    }
	}
	sort(reserved_junc.begin(),reserved_junc.end());
	reserved_junc.erase(unique(reserved_junc.begin(),reserved_junc.end()),reserved_junc.end());

	return;
  }
  bool SimplifyGraph::keep_edge(edge_t e)
  {
	if(find(reserved_junc.begin(),reserved_junc.end(),e) == reserved_junc.end()) 
		return false;

	return true;
  }
  void SimplifyGraph::get_CovInfo()
  {
    //cout<<"Info: "<<Graph_size<<" ";
    int good=0,bad=0;
    for(int i=0;i<Graph_size;i++)
    {
	double Ccov=0,Pcov=0;
	if(node_set[i].children.size()>0 && node_set[i].parents.size()>0)
	{
	    if(node_set[i].children.size()>1 || node_set[i].parents.size()>1)
	    {
		for(size_t j=0;j<node_set[i].children.size();j++) Ccov += node_set[i].children[j].second;
		for(size_t j=0;j<node_set[i].parents.size();j++) Pcov += node_set[i].parents[j].second;
		double min = Ccov<Pcov?Ccov:Pcov;
		double max = Ccov<Pcov?Pcov:Ccov;
		double rate = min/max;
		//cout<<rate<<" ";
		if(rate>0.7) good++;
		else bad++;
	    }
	}
    }
    //cout<<"; "<<good<<"/"<<bad;
    //cout<<endl;
    if(good > 5*bad) PackingFlag = true;
    if(PackingFlag) pack_graph_num++;
    else unpack_graph_num++;
  }
  void SimplifyGraph::build_node_set_and_simplify(vector<int>& exon_l, vector<int>& exon_r, vector<double>& exon_cov,
			       vector<node_idx_t>&edge_out ,vector<node_idx_t>&edge_in, vector<double>& edge_weight,
				vector<int>& Single_node,vector<path_t>& PairPath , vector<int>& invalid_nodes)
  {
    Graph_size = exon_l.size();
    size_ = Graph_size;
    RawSize = Graph_size;
  
    for(int i=0;i < Graph_size;i++)
    {
	Node node;
	node.coverage_ = exon_cov[i];
        node.length_ = exon_r[i] - exon_l[i] + 1;

	vector<int> vec;
        vec.push_back(exon_l[i]); vec.push_back(exon_r[i]);
        node.addseq(vec);

	node_set.push_back(node);
    }
    for(size_t i=0;i<edge_out.size();i++)
    {
	node_set[edge_out[i]].add_child(edge_in[i],edge_weight[i]);
 	node_set[edge_in[i]].add_parent(edge_out[i],edge_weight[i]);

    }
    for(int i=0;i<Graph_size;i++) {
	 if(node_set[i].children.empty() && node_set[i].parents.empty()) Single_node.push_back(i);
    }

    vector<edge_t> delete_edges_child,delete_edges_parent;

    remove_partial_junction(delete_edges_child,delete_edges_parent);

    remove_partial_end_by_edge_coverage(delete_edges_child,delete_edges_parent);

    AVERAGE_REMOVE_RATE = 0.1/SampleSize;
    remove_edges_by_average_coverage(delete_edges_child,delete_edges_parent);

    remove_intron_contamination(delete_edges_child,delete_edges_parent);

    //remove_all_partial(delete_edges_child,delete_edges_parent);
    remove_small_exons(delete_edges_child,delete_edges_parent);

    delete_children_edges(delete_edges_child);
    delete_parents_edges(delete_edges_parent);
    return;
  }
void SimplifyGraph::Statistic(string strand, string chr)
{
    //STAT
    for(int j = 0;j<Graph_size;j++)
    {
	
	if(node_set[j].children.size() == 2 && node_set[j].parents.size() == 2)
	{
	    int p1 = node_set[j].parents.front().first;
	    double covp1 = node_set[j].parents.front().second;

	    int p2 = node_set[j].parents.back().first;
	    double covp2 = node_set[j].parents.back().second;

	    int c1 = node_set[j].children.front().first;
	    double covc1 = node_set[j].children.front().second;

	    int c2 = node_set[j].children.back().first;
	    double covc2 = node_set[j].children.back().second;

	    if(partial(p1,j) || partial(p1,j) || partial(j,c1) || partial(j,c2)) break;

	    pair<int,int> j1 = make_pair(node_set[p1].sequence.back(),node_set[j].sequence.front());
	    vector<double> vw1 = get_cov_distribution(j1);
	    pair<int,int> j2 = make_pair(node_set[p2].sequence.back(),node_set[j].sequence.front());
	    vector<double> vw2 = get_cov_distribution(j2);

	    pair<int,int> j1_=make_pair(node_set[j].sequence.back(),node_set[c1].sequence.front());
	    vector<double> vw1_ = get_cov_distribution(j1_);
	    pair<int,int> j2_=make_pair(node_set[j].sequence.back(),node_set[c2].sequence.front());
	    vector<double> vw2_ = get_cov_distribution(j2_);

//	    cout<<"ok: "<<strand<<" "<<chr<<" "<<rg_index<<" "<<j<<endl;
	    cout<<j1.first<<"->"<<j1.second<<"  "<<j1_.first<<"->"<<j1_.second<<endl;
	    cout<<j2.first<<"->"<<j2.second<<"  "<<j2_.first<<"->"<<j2_.second<<endl;
	    /*
	    if(covp1>covp2 && covc1>covc2){
	        cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
		cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
	    }
	    else if(covp1>covp2 && covc1<covc2){
	        cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
		cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
	    }
	    else if(covp1<covp2 && covc1>covc2){
	        cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
		cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
	    }
	    else if(covp1<covp2 && covc1<covc2){
	        cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
		cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
	    }
	    */
	   /* 
	    cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
	    cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
	    cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
	    cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
	   */
	    double cos11 = cosine(vw1,vw1_);
	    double cos12 = cosine(vw1,vw2_);

	    double cos21 = cosine(vw2,vw1_);
	    double cos22 = cosine(vw2,vw2_);
	    cout<<cos11<<" "<<cos12<<" "<<cos21<<" "<<cos22<<endl;
	    
	    if(cos11*cos22 > cos12*cos21){// && cos11>0.5 && cos22 >0.5){
	    //if(cos12<0.5 && cos21<0.5 && cos11>0.5 && cos22 >0.5){
	    //if(cos11> cos12 && cos11>cos21 && cos22>cos12 &&cos22>cos21){// && cos11*cos22 > cos12*cos21){
	    //if(0&& cos11>0.5 && cos22 >0.5){
	    
	        cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
		cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j2_.first<<" "<<j2_.second<<endl;

	    } 
	    else if(cos11*cos22 < cos12*cos21){// && cos12>0.5 && cos21>0.5){
	    //else if(cos11<0.5 && cos22 < 0.5 && cos12>0.5 && cos21>0.5){
	    //else if(cos11 < cos12 && cos11 < cos21 && cos22<cos12 &&cos22<cos21){// && cos11*cos22 < cos12*cos21){
	    //if(0&& cos12>0.5 && cos21>0.5){
	    	cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
		cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j1_.first<<" "<<j1_.second<<endl;

	    }
	    
	    
	      if(0 && cos11>0.5 && cos12<=0.5 && cos21 <= 0.5 && cos22>0.5){
		      cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
		      cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
		      //cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
		      //cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j1_.first<<" "<<j1_.second<<endl;

	      }
	      vector<edge_t> delete_edges;
	      /*
	      if(cos11 < 0.1 && cos12 < 0.1) delete_edges.push_back(make_pair(p1,j));
	      if(cos21 < 0.1 && cos22 < 0.1) delete_edges.push_back(make_pair(p2,j));
	      delete_children_edges(delete_edges);
	      */
	      /*
	      if(cos11>0.5)cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
	      if(cos12>0.5) cout<<strand<<" "<<chr<<" "<<j1.first<<" "<<j1.second<<" "<<j2_.first<<" "<<j2_.second<<endl;

	      if(cos21>0.5) cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j1_.first<<" "<<j1_.second<<endl;
	      if(cos22>0.5) cout<<strand<<" "<<chr<<" "<<j2.first<<" "<<j2.second<<" "<<j2_.first<<" "<<j2_.second<<endl;
	      */
	    
/*
   	    for(size_t i=0;i<node_set[current_edge.second].children.size() ;i++)

		  edge_t e = make_pair(current_edge.second, node_set[current_edge.second].children[i].first);
		  double c = node_set[current_edge.second].children[i].second;

		  pair<int,int> junc = make_pair(node_set[e.first].sequence.back(), node_set[e.second].sequence.front());
		  vector<double> junc_cov_distribution = get_cov_distribution(junc);

		  double d = distance(current_junc_cov_distribution,junc_cov_distribution);
		  //edges_cov.push_back(make_pair(e,double(1.0/d)));
		  
		  double p = pearson(current_junc_cov_distribution,junc_cov_distribution);

		  double cos = cosine(current_junc_cov_distribution,junc_cov_distribution);
		  double adcos = adcosine(current_junc_cov_distribution,junc_cov_distribution);
	 	  double dt = dot(current_junc_cov_distribution,junc_cov_distribution);
		  edges_cov.push_back(make_pair(e,cos*c));
*/		  //edges_cov.push_back(make_pair(e,node_set[current_edge.second].children[i].second));

	    }
	}
    return;
  }
  void SimplifyGraph::remove_all_partial(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    for(int i=0;i<Graph_size;i++)
    {
	for(size_t j = 0;j<node_set[i].children.size();j++)
	{
	    node_idx_t c = node_set[i].children[j].first;
	    if(MultiJunction(i,c)) continue;
	    if(partial(i,c) && node_set[i].children[j].second < 0.05)
	    {
		delete_edges_child.push_back(make_pair(i,c));
	    }
	}
    }
    return;
  }

  void SimplifyGraph::remove_edges_by_average_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	double all_coverage = 0;
	int edges_number = 0;
	map<node_idx_t, bool> used_nodes;
	vector<vector<int> > Comps;

	for(int i=0;i<Graph_size;i++)
	{
	    if(used_nodes[i]) continue;
	    used_nodes[i] = true;
	    vector<node_idx_t> comp_nodes;
	    comp_nodes.push_back(i);
	    for(size_t j=0;j<comp_nodes.size();j++)
	    {
		for(size_t k=0;k<node_set[comp_nodes[j]].children.size();k++)
		{
		    node_idx_t c = node_set[comp_nodes[j]].children[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), c) == comp_nodes.end())
		    {
			//cerr<<"  "<<c<<endl;
			comp_nodes.push_back(c);
			used_nodes[c] = true;
		    }
		}
		for(size_t k=0;k<node_set[comp_nodes[j]].parents.size();k++)
		{
		    node_idx_t p = node_set[comp_nodes[j]].parents[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), p) == comp_nodes.end())
		    {
		 	comp_nodes.push_back(p);
			used_nodes[p] = true;
		    }
		}
	    }
	    Comps.push_back(comp_nodes);
 	}
	map<node_idx_t,double> nodes_average_coverage_map;
	//cout<<"hahaha: "<<Comps.size()<<endl;
	for(size_t i=0;i<Comps.size();i++)
	{
	    vector<int> comp = Comps[i];
	    double all_cov = 0;
	    double edges_number = 0;
	    for(size_t j=0;j<comp.size();j++)
	    {
		node_idx_t n = comp[j];
		for(size_t k=0;k<node_set[n].children.size();k++)
		{
		    if(partial(n,node_set[n].children[k].first)) continue;
		    all_cov += node_set[n].children[k].second;
		    edges_number++;
		}
	    }
	    double ave = 1;
	    if(edges_number != 0) ave = all_cov/edges_number;
	    if(ave < 4.0/SampleSize) 
	    {
		for(size_t j=0;j<comp.size();j++) node_set[comp[j]].clear();
	    }
		
	    for(size_t j=0;j<comp.size();j++)
		nodes_average_coverage_map[comp[j]] = ave;
	}

	
	for(int i=0;i<Graph_size;i++)
	{
	    for(size_t j=0;j<node_set[i].children.size();j++)
	    {
		double rate = node_set[i].children[j].second / nodes_average_coverage_map[i];

		node_idx_t c = node_set[i].children[j].first;
		if(MultiJunction(i,c)) continue;
		if(rate < AVERAGE_REMOVE_RATE)
	 	{
		    delete_edges_child.push_back(make_pair(i,node_set[i].children[j].first));
		}
	    }
 	}
	return;
  }
  bool SimplifyGraph::MultiJunction(node_idx_t i, node_idx_t c)
  {
      if(node_set[i].sequence.empty() || node_set[c].sequence.empty()) return false;
      pair<int,int> junc = make_pair(node_set[i].sequence.back(), node_set[c].sequence.front());
      vector<pair<int,int> >::iterator it = find(junctions.begin(),junctions.end(),junc);
      if(it == junctions.end()) return false;
      int k = it - junctions.begin();
      if(junctions_MappingInfo[k].size() <= 1) return false;
      if(junctions_MappingInfo[k].size() <= 0.05*SampleSize) return false;//1->0.05*SampleSize ipsc

      return true;
  }
  void SimplifyGraph::remove_partial_junction(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	double R=0;
	if(!PackingFlag) R = 0.12;
	else R=0.08;
	if(SFlag) R=0.5;//3.1 spk

	R=0.5;
    	for(int i=1;i<Graph_size - 1;i++)
        {
	    if(partial(i,i+1))
	    {
	  	double partial_cov = node_set[i].get_child_coverage(i+1);
		if(partial_cov != -1)
		{
		    for(size_t j=0;j<node_set[i].children.size();j++)
	   	    {
			node_idx_t c = node_set[i].children[j].first;
			if(MultiJunction(i,c)) continue;
		 	if(node_set[i].children[j].second/partial_cov < R) 
			    delete_edges_child.push_back(make_pair(i,c));
		    }
		}
	    }
	    if(partial(i-1,i))
	    {
		double partial_cov = node_set[i-1].get_child_coverage(i);
		if(partial_cov != -1)
		{
		    for(size_t j=0;j<node_set[i].parents.size();j++)
		    {
			node_idx_t p = node_set[i].parents[j].first;
			if(MultiJunction(p,i)) continue;
			if(node_set[i].parents[j].second / partial_cov < R)
			    delete_edges_parent.push_back(make_pair(p,i));
		    }
		}
	    }
	}
	return;
  }
  void SimplifyGraph::remove_unbalance_edges( vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	for(int i=0;i<Graph_size;i++)
	{
	    if( !node_set[i].parents.empty() && node_set[i].children.size() == 1)
	    {
		double Incov = InCov(i);
		node_idx_t c = node_set[i].children.front().first;
		double rate = node_set[i].children.front().second/Incov;
		if(node_set[i].children.front().second > 5 && keep_edge(make_pair(i,c)) ) continue;
		if(rate<UNBALANCE_RATE) delete_edges_child.push_back(make_pair(i,c));
	    }
	    if(node_set[i].parents.size() == 1 && !node_set[i].children.empty())
	    {
		double outcov = OutCov(i);
		node_idx_t p = node_set[i].parents.front().first;
		double rate = node_set[i].parents.front().second/outcov;
		if(node_set[i].parents.front().second > 5 && keep_edge(make_pair(p,i)) ) continue;
		if(rate<UNBALANCE_RATE) delete_edges_parent.push_back(make_pair(p,i));
	    }
	}
	return;
  }
  void SimplifyGraph::remove_small_exons(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	for(int i=0;i<Graph_size;i++)
	{
	    //one out
	    if(node_set[i].parents.empty() && node_set[i].children.size() == 1)
	    {
		node_idx_t c =  node_set[i].children.front().first;
		double ecov = node_set[i].get_child_coverage(c);
		double r1 = double((1.0*ecov)/node_set[c].coverage());
		double r2 = double((1.0*node_set[i].coverage())/node_set[c].coverage());
		if(MultiJunction(i,c)) continue;
		//if(node_set[i].coverage() <= 2  && ecov < 3 && r1<0.2 &&  node_set[c].parents.size() > 1)
		if(node_set[i].coverage() < 5.0/SampleSize  && ecov < 5.0/SampleSize && r1<0.2 && r2<0.2 && node_set[c].parents.size() > 1 )
		//if(node_set[i].coverage() < 1 && ecov < 1 && r1<0.8 && r2<0.8 && node_set[c].parents.size() > 1 )
		{
		    delete_edges_child.push_back(make_pair(i,c));
		}
	    }
	    //one in
	    if(node_set[i].children.empty() && node_set[i].parents.size() == 1)
	    {
		node_idx_t p = node_set[i].parents.front().first;
		double ecov = node_set[p].get_child_coverage(i);

		double r1 = double(1.0*ecov/node_set[p].coverage());
		double r2 = double((1.0*node_set[i].coverage())/node_set[p].coverage());
		if(MultiJunction(p,i)) continue;
		//if(node_set[i].coverage() <=2 && ecov < 3 && r1<0.15 && node_set[p].children.size() > 1)
		if(node_set[i].coverage() < 5.0/SampleSize && ecov < 5.0/SampleSize && r1<0.2 && r2<0.2 && node_set[p].children.size() > 1)
		//if(node_set[i].coverage() < 1&& ecov < 1 && r1<0.8 && r2<0.8 && node_set[p].children.size() > 1)
		{
		    delete_edges_parent.push_back(make_pair(p,i));
		}
	    }
	}
        return;
  }
  bool SimplifyGraph::partial(node_idx_t n1, node_idx_t n2)
  {
	if(node_set[n1].sequence.empty() || node_set[n2].sequence.empty() ) return false;

	if(node_set[n1].sequence.back() + 1 == node_set[n2].sequence.front()) return true;
	return false;
  }
  void SimplifyGraph::remove_intron_contamination(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	for(int i=0;i<Graph_size;i++)
	{
	    if(node_set[i].children.size() == 1 && node_set[i].parents.size() == 1)
	    {
		bool partial_flag = false;
		node_idx_t p = node_set[i].parents.front().first;
		node_idx_t c = node_set[i].children.front().first;
		double pc_cov = node_set[p].get_child_coverage(c);
		if(partial(p,i) && partial(i,c) && pc_cov != -1)
		{
		    if( (node_set[i].coverage() < 10.0/SampleSize && node_set[i].coverage()<pc_cov) )
		      
		    {
			delete_edges_child.push_back(make_pair(i,c));
			delete_edges_parent.push_back(make_pair(p,i));
		    } 
		}

	    }
	}
	return;
  }

  void SimplifyGraph::remove_partial_end_by_edge_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	double R=1;//1
	for(int i=0;i<Graph_size;i++)
	{
	    if(node_set[i].children.size() > 1 || node_set[i].parents.size() > 1)
	    {
		bool partial_flag = false;
		node_idx_t partial_node;
		double max_cov = 0,partial_cov = 0;

		for(size_t j = 0;j<node_set[i].children.size();j++)
		{
		    if(partial(i,node_set[i].children[j].first)){
			partial_node = node_set[i].children[j].first;
			partial_cov = node_set[i].children[j].second;
			partial_flag = true;
		    }
		    if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
		}
	    	if(partial_flag)
		{
		    double r = double(1.0*partial_cov)/max_cov; //0.2
		    if(node_set[partial_node].children.empty() && r <= 0.2 && node_set[partial_node].length() < 10000 && node_set[partial_node].parents.size() == 1)
			delete_edges_child.push_back(make_pair(i,partial_node));
			size_t len = node_set[partial_node].length();
                        if( (node_set[partial_node].children.empty() && node_set[partial_node].parents.size() == 1)
                           && r<=R && len < 200) 
                           //10.11-noref
                            delete_edges_child.push_back(make_pair(i,partial_node));
		}

		partial_flag = false; max_cov = 0;
		for(size_t j=0;j<node_set[i].parents.size();j++)
		{
		    if(partial(node_set[i].parents[j].first,i)){
			partial_node = node_set[i].parents[j].first;
			partial_cov = node_set[i].parents[j].second;
			partial_flag = true;
		    }
		    if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
		}
	  	if(partial_flag)
		{
		    double r = double(1.0*partial_cov)/max_cov; //0.2
		    if(node_set[partial_node].parents.empty() && r <= 0.2&& node_set[partial_node].length() < 10000 && node_set[partial_node].children.size() == 1)
			delete_edges_parent.push_back(make_pair(partial_node,i));

			size_t len = node_set[partial_node].length();
                        if( (node_set[partial_node].parents.empty() && node_set[partial_node].children.size() == 1)
                           && r<=R && len < 200) //1,20
                         //10.11-noref
                        delete_edges_parent.push_back(make_pair(partial_node,i));
		}
	    } 
	}
	return;
  }
  void SimplifyGraph::remove_edges_onemulti(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    for(int i=0;i<Graph_size;i++)
    {
	if(node_set[i].children.size() > 1 && node_set[i].parents.size() <=1)//one to multi
	{
	    double max_cov = 0;
	    for(size_t j = 0;j<node_set[i].children.size();j++)
		if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
	    for(size_t j = 0;j<node_set[i].children.size();j++)
	    {
		node_idx_t c = node_set[i].children[j].first;
		double cov = node_set[i].children[j].second;
		if(MultiJunction(i,c)) continue;
//		if(cov <= 1.0/SampleSize  ||(cov <= 10.0/SampleSize && double(cov/max_cov) <= 0.5))
		if(double(cov/max_cov) <= 0.01)
	 	{
		  if( node_set[c].children.empty())// || node_set[c].parents.size() > 1)
		  {
		    node_set[c].length_ = -1;
		    delete_edges_child.push_back(make_pair(i,c));
		  }
		}
	    }
	}
	if(node_set[i].children.size() <= 1 && node_set[i].parents.size() > 1)//multi to one
	{
	    double max_cov = 0;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    	if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
	    	node_idx_t p = node_set[i].parents[j].first;
	    	double cov = node_set[i].parents[j].second;
		if(MultiJunction(p,i)) continue;
	    	//if(cov <= 1.0/SampleSize  ||(cov <= 10.0/SampleSize && double(cov/max_cov) <= 0.5))
		if(double(cov/max_cov) <= 0.01)
	    	{
		  if(node_set[p].parents.empty())// || node_set[p].children.size() > 1)
		  {
		    node_set[p].length_ = -1;
		    delete_edges_parent.push_back(make_pair(p,i));
		  }
	    	}
	    }
	}
    }
    return;
  }
  void SimplifyGraph::remove_lowcov_edges_of_bifurcation_nodes2(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    for(int i=0;i<Graph_size;i++)
    {
	if(node_set[i].children.size() > 1 )
	{
	    double max_cov = 0;
	    for(size_t j = 0;j<node_set[i].children.size();j++) 
	    {
		if( !partial(i,node_set[i].children[j].first)) continue;

		if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
	    }

	    for(size_t j = 0;j<node_set[i].children.size();j++)
	    {
		node_idx_t c = node_set[i].children[j].first;
	   	double cov = node_set[i].children[j].second;
		double r = double(cov/max_cov);
		if(MultiJunction(i,c)) continue;
		if( !partial(i,c) && cov < max_cov)
		{
			node_set[c].length_ = -1;
			delete_edges_child.push_back(make_pair(i,c));
		}
	    }
	    
	    max_cov = 0;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		if( !partial(node_set[i].parents[j].first,i)) continue;
		if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
	    }

	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		node_idx_t p = node_set[i].parents[j].first;
		double cov = node_set[i].parents[j].second;
		double r = double(cov/max_cov);
		if(MultiJunction(p,i)) continue;
		if(!partial(p,i) &&  cov < max_cov)//5 0.18
		{
			node_set[p].length_ = -1;
			delete_edges_parent.push_back(make_pair(p,i));
		}
	    }
	    //cerr<<endl;
	}	
    }
    return;
  }
  void SimplifyGraph::remove_lowcov_edges_of_bifurcation_nodes(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    double COV1 = 3.0/SampleSize;
    double R = 0.2;
    double COV2 = 5.0/SampleSize;
    //if(PackingFlag) R = 0.1;
    for(int i=0;i<Graph_size;i++)
    {
	//cerr<<"node: "<<i<<endl;
	if(node_set[i].children.size() > 1 && node_set[i].parents.size() > 1)
	{
	    //cerr<<"child"<<endl;
	    double max_cov = 0;
	    for(size_t j = 0;j<node_set[i].children.size();j++) 
	    {
		if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
	    }

	    for(size_t j = 0;j<node_set[i].children.size();j++)
	    {
		node_idx_t c = node_set[i].children[j].first;
	   	double cov = node_set[i].children[j].second;
		if(MultiJunction(i,c)) continue;
		//cerr<<c<<": "<<cov<<"; ";
		double r = double(cov/max_cov);
		if(cov < COV1 || (cov <= COV2 && r <= R)) //5,0.18
		{
		    if(1 || node_set[c].children.empty() || node_set[c].parents.size() > 1)
		    {
			node_set[c].length_ = -1;
			delete_edges_child.push_back(make_pair(i,c));
		    }
		}
	    }
	    //cerr<<endl;
	    //cerr<<"parent"<<endl;
	    
	    max_cov = 0;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
	    }

	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		node_idx_t p = node_set[i].parents[j].first;
		double cov = node_set[i].parents[j].second;
		if(MultiJunction(p,i)) continue;
		//cerr<<p<<": "<<cov<<"; ";
		double r = double(cov/max_cov);
		if(cov < COV1 || (cov <= COV2 && r <= R))//5 0.18
		{
		    if(1 || node_set[p].parents.empty() || node_set[p].children.size() > 1)
		    {
			node_set[p].length_ = -1;
			delete_edges_parent.push_back(make_pair(p,i));
		    }
		}
	    }
	    //cerr<<endl;
	}	
    }
    return;
  }
  void SimplifyGraph::delete_children_edges( vector<edge_t> delete_edges_child)
  {
    sort(delete_edges_child.begin(),delete_edges_child.end());
    delete_edges_child.erase(unique(delete_edges_child.begin(),delete_edges_child.end()),delete_edges_child.end());
    for(size_t i = 0;i<delete_edges_child.size();i++)
    {
	node_idx_t n1 = delete_edges_child[i].first, n2 = delete_edges_child[i].second;
    	if( !node_set[n1].sequence.empty() && !node_set[n2].sequence.empty() )
        {
            double cov = node_set[n1].get_child_coverage(n2);

            pair<int,int> junction = make_pair(node_set[n1].sequence.back(),node_set[n2].sequence.front());
	}
	node_set[n1].delete_child(n2);
	node_set[n2].delete_parent(n1);
    }
    /*
    for(size_t i = 0;i<delete_edges_child.size();i++)
    {
	node_idx_t c = delete_edges_child[i].second;
	vector< pair<node_idx_t ,double> > c_children = node_set[c].children;
	while(1)
	{
	    if(c_children.empty() ) break;
	    if(c_children.size() >= 2) break;

	    node_idx_t c_child = c_children.front().first;
	    if(c_children.front().second > 1.5) break;
	    node_set[c].delete_child(c_child);
	    node_set[c_child].delete_parent(c);

	    if(!node_set[c_child].parents.empty()) break;//other in_edges

	    c_children = node_set[c_child].children;
	}
    }
    */
    
    return;
  }    
  void SimplifyGraph::delete_parents_edges( vector<edge_t> delete_edges_parent)
  {
    sort(delete_edges_parent.begin(),delete_edges_parent.end());
    delete_edges_parent.erase(unique(delete_edges_parent.begin(),delete_edges_parent.end()), delete_edges_parent.end());
    for(size_t i = 0;i<delete_edges_parent.size();i++)
    {
	node_idx_t n1 = delete_edges_parent[i].first, n2 = delete_edges_parent[i].second;
	if( !node_set[n1].sequence.empty() && !node_set[n2].sequence.empty() )
        {
             double cov = node_set[n1].get_child_coverage(n2);
 
             pair<int,int> junction = make_pair(node_set[n1].sequence.back(),node_set[n2].sequence.front());
        }
	node_set[n1].delete_child(n2);
	node_set[n2].delete_parent(n1);
    }
    /*
    for(size_t i = 0;i<delete_edges_parent.size();i++)
    {
	node_idx_t p = delete_edges_parent[i].first;
	vector< pair<node_idx_t ,double> > p_parents = node_set[p].parents;
	while(1)
	{
	    if(p_parents.empty()) break;
	    if(p_parents.size() >= 2) break;
	    node_idx_t p_parent = p_parents.front().first;
	    if(p_parents.front().second > 1.5) break;
	    node_set[p_parent].delete_child(p);
	    node_set[p].delete_parent(p_parent);

	    if(!node_set[p_parent].children.empty()) break; //other out_edges;
	    p_parents = node_set[p_parent].parents;
	}
    }
    */
    return;
  }
  void SimplifyGraph::show_graph()
  {
    for(int i=0;i<Graph_size;i++)
    {
        Node n = node_set[i];
        for(size_t j=0;j<n.children.size();j++)
        {
            cout<<i<<"->"<<n.children[j].first<<": "<<n.children[j].second<<endl;
        }
    }
    for(int i=0;i<Graph_size;i++)
    {
	cout<<i<<": "<<node_set[i].path_str()<<endl;
    }
    return;
  }
  void SimplifyGraph::contract_graph()
  {
    map<node_idx_t,bool> used_nodes;
    int S = Graph_size;
    for(int i=0;i<S;i++)
    {
	if(used_nodes[i]) continue;
	used_nodes[i] = true;
	path_t P(1,i);
	while(1)
	{
	    node_idx_t n = P.back();
	    bool eflag=false;
	    if(node_set[n].children.size() == 1) 
	    {
		node_idx_t c = node_set[n].children.front().first;
		if(node_set[c].parents.size() == 1){
		   P.push_back(c);
		   used_nodes[c] = true;
		   eflag = true;
		}
	    }
	    if(!eflag) break;
	}
	while(1)
	{
	    node_idx_t n = P.front();
	    bool eflag=false;
	    if(node_set[n].parents.size() == 1)
	    {
		node_idx_t p = node_set[n].parents.front().first;
		if(node_set[p].children.size() == 1){
		    P.insert(P.begin(),p);
		    used_nodes[p] = true;
		    eflag = true;
		}
	    }
	    if(!eflag) break;
	}
	if(P.size() == 1){
	    rawNode_eNode_map[i] = i;
	    continue;
	}
	Node N;
	N.coverage_ = 1;
	N.length_ = 1;
	N.path = P;
	node_set.push_back(N);
	Graph_size++;
	node_idx_t Nid = Graph_size - 1;
	node_set[Nid].children = node_set[P.back()].children;
	node_set[Nid].parents = node_set[P.front()].parents;
	for(size_t i=0;i<node_set[P.back()].children.size();i++){
	    node_idx_t c = node_set[P.back()].children[i].first;
	    double cov = node_set[P.back()].children[i].second;
	    node_set[ c ].delete_parent(P.back());
	    node_set[ c ].add_parent(Nid,cov);
	}
	for(size_t i=0;i<node_set[P.front()].parents.size();i++){
	    node_idx_t p = node_set[P.front()].parents[i].first;
	    double cov = node_set[P.front()].parents[i].second;
	    node_set[ p ].delete_child(P.front());
	    node_set[ p ].add_child(Nid,cov);
	}
	for(size_t i=0;i<P.size();i++){ 
	    rawNode_eNode_map[P[i]] = Nid;
	    node_set[P[i]].clear();
	}
    }
    return;
  }
  void SimplifyGraph::get_graph(vector<vector<int> >& vecEdges,vector<double>& vecCov, 
			  vector<path_t>& PairPath, vector<double>& PairPath_cov)
  {
	//cout<<"GetGraph: "<<endl;
        for(size_t i=0;i<PairPath.size();i++)
     	{
	    if(rawNode_eNode_map.empty()) break;
	    path_t pp = PairPath[i];
	    //for(size_t j=0;j<pp.size();j++) cout<<pp[j]<<"* ";
	    //cout<<endl;
	    for(size_t j=0;j<pp.size();j++) pp[j] = rawNode_eNode_map[pp[j]];
	    //for(size_t j=0;j<pp.size();j++) cout<<pp[j]<<"& ";
 	    //cout<<endl;
	    path_t::iterator it = pp.begin(),it1;
	    it++;
 	    for(;it != pp.end();)
  	    {
		it1=find(pp.begin(),it,*it);
	 	if(it1 != it) pp.erase(it);
		else it++;
	    }
	    //for(size_t j=0;j<pp.size();j++) cout<<pp[j]<<"+ ";
	    //cout<<endl;
	    PairPath[i] = pp;
	}
	for(int i=0;i<Graph_size;i++)
	{
	    for(size_t j=0;j<node_set[i].children.size();j++)
	    {
		vector<int> edge;
		double cov;
		edge.push_back(i);edge.push_back(node_set[i].children[j].first);
	   	cov = node_set[i].children[j].second;

		vecEdges.push_back(edge); 
		if(cov == 0) cov = 0.1;
		vecCov.push_back(cov);
	    }
	}
	double R = 4;//3.1
	if(SFlag) R = 5000000000;
	for(size_t i=0;i<PairPath.size();)
	{
	    path_t p = PairPath[i];
	    //cerr<<"check2: "<<endl;
	    //for(size_t j=0;j<p.size();j++) cerr<<p[j]<<" ";
	    //cerr<<": "<<PairPath_cov[i]<<endl;
	    bool flag = false;
	    if(p.size() < 3) flag = true;//NEW
	    for(size_t j=0;j<p.size() - 1;j++)
	    {
	        if(node_set[p[j]].get_child_coverage(p[j+1]) == -1)
		{
	 	    flag = true;
		    break;
		}
	    }
	    
	    if(PairPath_cov[i]<R) {//3.1 spk
		flag = true;
	    }
	    
	    if(flag)
	    {
	    	PairPath.erase(PairPath.begin() + i);
	    	PairPath_cov.erase(PairPath_cov.begin() + i);
	    }
	    else i++;
	}
	pairpath = PairPath;
	pairpath_cov = PairPath_cov;
	//partial pair
	//bool flag = false;
	/*
	for(int i=0;i<Graph_size;i++)
        {
	    if(node_set[i].parents.size()>1 && node_set[i].children.size()>1)
	    {
		cerr<<"bifurcation node: "<<i<<endl;
		cerr<<"parents "<<endl;
	 	for(size_t j=0;j<node_set[i].parents.size();j++)
		{
		    cerr<<node_set[i].parents[j].first<<" "<<node_set[i].parents[j].second<<endl;
		}
		cerr<<"children: "<<endl;
		for(size_t j=0;j<node_set[i].children.size();j++)
		{
		    cerr<<node_set[i].children[j].first<<" "<<node_set[i].children[j].second<<endl;
		}
	    }
	}
	*/
	return;
  }
  double SimplifyGraph::InCov(node_idx_t node)
  {
    double cov = 0;
    for(size_t i=0;i<node_set[node].parents.size();i++)
    {
	cov += node_set[node].parents[i].second;
    }
    return cov;
  }

  double SimplifyGraph::OutCov(node_idx_t node)
  {
    double cov = 0;
    for(size_t i=0;i<node_set[node].children.size();i++)
    {
	cov += node_set[node].children[i].second;
    }
    return cov;
  }
  double SimplifyGraph::average_coverage()
  {
	double cov = 0;
	int number = 0;
	for(int i=0;i<Graph_size;i++)
	{
	    for(size_t j=0;j<node_set[i].children.size();j++)
	    {
		cov += node_set[i].children[j].second;
		number++;
	    }
	}
	if(number == 0) return 1;
	return (cov/number);
  }
  void SimplifyGraph::get_packing_result_new(vector<vector<int> >Vec_edges, vector<int> Edges_legt, vector<int>Edges_right, vector<double> Weights)
  {
//      packing_map
	packing_map.clear();      
	for(size_t i=0;i<Weights.size();i++)
	{
	    if(Weights[i] == 1) continue;
	    int i1 = Edges_legt[i], i2 = Edges_right[i];
	    vector<int> E1 = Vec_edges[i1], E2 = Vec_edges[i2];

	    pair<int,int> e1 = make_pair(E1[0],E1[1]), e2 = make_pair(E2[0],E2[1]);
	    packing_map[make_pair(e1,e2)] = 0;
	}
	return;
  }
  vector<double> SimplifyGraph::get_cov_distribution(pair<int,int> junc)
  {
	
	vector<double> vdistrib;
	for(size_t i=0;i<SampleSize;i++) vdistrib.push_back(0.0);
	vector<pair<int,int> >::iterator it = find(junctions.begin(),junctions.end(),junc);
        if(it == junctions.end()) return vdistrib;
        int k = it - junctions.begin();
	vector<pair<int,double> > junc_info = junctions_MappingInfo[k];
	for(size_t i=0;i<junc_info.size();i++)
	{
	    int sample_id = junc_info[i].first;
	    double cov = junc_info[i].second;
	    vdistrib[sample_id] = cov;
	}

	return vdistrib;

  }

  double SimplifyGraph::distance(vector<double> junc1, vector<double> junc2)
  {
	double d;
	for(size_t i=0;i<SampleSize;i++)
	{
	//    if(junc1[i] == 0 || junc2[i]==0) continue;

	    d += (junc1[i] - junc2[i])*(junc1[i] - junc2[i]);
	}
	
	if(d>0) return sqrt(d);

	return 0.001;
  }
  double SimplifyGraph::pearson(vector<double> junc1, vector<double> junc2)
  {
	double cov1 = 0.0, cov2 = 0.0;  
	for(size_t i=0;i<SampleSize;i++){
	    cov1 += junc1[i];
	    cov2 += junc2[i];
	}
	cov1 = cov1/(1.0*SampleSize);
	cov2 = cov2/(1.0*SampleSize);

	double d1 = 0.0, d2 = 0.0, d3 = 0.0;
	for(size_t i=0;i<SampleSize;i++){
	    d1 += (junc1[i] - cov1)*(junc2[i] - cov2);
	    d2 += (junc1[i] - cov1)*(junc1[i] - cov1);
	    d3 += (junc2[i] - cov2)*(junc2[i] - cov2);
	}
	double d4 = double((sqrt(d2))*(sqrt(d3)));
	return d1/d4;
  }

  double SimplifyGraph::adcosine(vector<double> junc1, vector<double> junc2)
  {
      double c1 = 0, c2 = 0;
      for(size_t i=0;i<SampleSize;i++)
      {
          c1 += junc1[i];
	  c2 += junc2[i];
      }
      c1 = c1/(1.0*SampleSize);
      c2 = c2/(1.0*SampleSize);
      double d = 0.0, norm1 = 0.1, norm2 = 0.0;
      for(size_t i=0;i<SampleSize;i++)
      {
 /*          d += (junc1[i]-c1)*(junc2[i]-c2);
	   norm1 += (junc1[i]-c1)*(junc1[i]-c1);
	   norm2 += (junc2[i]-c2)*(junc2[i]-c2);
*/
	      //sqrt
	   d += sqrt(junc1[i]*junc2[i]);
	   norm1 += junc1[i];
	   norm2 += junc2[i];
      }

      return (d/(sqrt(norm1)*sqrt(norm2)));
  }
  double SimplifyGraph::cosine(vector<double> junc1, vector<double> junc2){
      double d = 0.0, norm1 = 0.1, norm2 = 0.0;
      for(size_t i=0;i<SampleSize;i++)
      {
           d += junc1[i]*junc2[i];
	   norm1 += junc1[i]*junc1[i];
	   norm2 += junc2[i]*junc2[i];

      }
      double cos = (d/(sqrt(norm1)*sqrt(norm2)));
      if(cos == 0) return 0.001;
      return cos;
      //return (d/(sqrt(norm1)*sqrt(norm2)));
  }
  double SimplifyGraph::dot(vector<double> junc1, vector<double> junc2){
      double d = 0.0, norm1 = 0.1, norm2 = 0.0;
      for(size_t i=0;i<SampleSize;i++)
	      d += junc1[i]*junc2[i];

      return d;
  }
 
  void SimplifyGraph::forward_extend(edge_t e, vector<path_t>& paths)
  {
      path_t path_temp;
      path_temp.push_back(e.first);
      path_temp.push_back(e.second);
      vector<path_t> paths_temp;
      paths_temp.push_back(path_temp);
      while(1)
      {
	  if(paths_temp.empty()) break;

	  path_t current_path = paths_temp.front();
	  paths_temp.erase(paths_temp.begin());

          edge_t current_edge = make_pair(current_path[current_path.size() - 2], current_path[current_path.size() - 1]);
	  if(node_set[current_edge.second].children.size() == 0)
	  {
	      paths.push_back(current_path);    
	  }
	  else
	  {
	      vector<pair<edge_t,double> > edges_cov;
	      pair<int,int> current_junc =  make_pair(node_set[current_edge.first].sequence.back(), node_set[current_edge.second].sequence.front());
	      vector<double> current_junc_cov_distribution = get_cov_distribution(current_junc);

	      for(size_t i=0;i<node_set[current_edge.second].children.size() ;i++)
	      {

		  edge_t e = make_pair(current_edge.second, node_set[current_edge.second].children[i].first);
		  double c = node_set[current_edge.second].children[i].second;

		  pair<int,int> junc = make_pair(node_set[e.first].sequence.back(), node_set[e.second].sequence.front());
		  vector<double> junc_cov_distribution = get_cov_distribution(junc);

		  //double d = distance(current_junc_cov_distribution,junc_cov_distribution);
		  
		  //double p = pearson(current_junc_cov_distribution,junc_cov_distribution);

		  double cos = cosine(current_junc_cov_distribution,junc_cov_distribution);
		  //double adcos = adcosine(current_junc_cov_distribution,junc_cov_distribution);
	 	  //double dt = dot(current_junc_cov_distribution,junc_cov_distribution);
		  edges_cov.push_back(make_pair(e,c*cos));
		  //edges_cov.push_back(make_pair(e,node_set[current_edge.second].children[i].second));

		  
	      }
	      sort(edges_cov.begin(),edges_cov.end(),sorter); //revision

	      if(edges_cov.size() == 1)
	      {
	          edge_t e = edges_cov[0].first;
		  current_path.push_back(e.second);
		  paths_temp.push_back(current_path);
	      }
	      else {
	        for(size_t i=0;i<edges_cov.size();i++)
	        {
	          edge_t e = edges_cov[i].first;
		  //if(partial(e.first,e.second)) continue;
	          pair<edge_t ,edge_t> ee = make_pair(current_edge, e);
		  if(1||packing_map.find(ee) != packing_map.end())
		  {
		      current_path.push_back(e.second);
		      paths_temp.push_back(current_path);
		      break;
		  }
		}
	      }
	  }
      }
      return;

  }
  void SimplifyGraph::reverse_extend(edge_t e, vector<path_t>& paths)
  {
  
      path_t path_temp;
      path_temp.push_back(e.first);
      path_temp.push_back(e.second);
      vector<path_t> paths_temp;
      paths_temp.push_back(path_temp);
      while(1)
      {
	  if(paths_temp.empty()) break;

	  path_t current_path = paths_temp.front();
	  paths_temp.erase(paths_temp.begin());

          edge_t current_edge = make_pair(current_path[0], current_path[1]);
	  if(node_set[current_edge.first].parents.size() == 0)
	  {
	      paths.push_back(current_path);    
	  }
	  else
	  {
	      vector<pair<edge_t,double> > edges_cov;
	      pair<int,int> current_junc =  make_pair(node_set[current_edge.first].sequence.back(), node_set[current_edge.second].sequence.front());
	      vector<double> current_junc_cov_distribution = get_cov_distribution(current_junc);
	      /*
	      cout<<endl<<"--------------------------------------------------------------"<<endl;
	      for(size_t i = 0;i<current_junc_cov_distribution.size();i++)
	          cout<<current_junc_cov_distribution[i]<<" ";
	      cout<<endl;
	      */
	      for(size_t i=0;i<node_set[current_edge.first].parents.size() ;i++)
	      {
		  edge_t e = make_pair(node_set[current_edge.first].parents[i].first,current_edge.first);
		  double c = node_set[current_edge.first].parents[i].second;

		  pair<int,int> junc = make_pair(node_set[e.first].sequence.back(), node_set[e.second].sequence.front());
		  vector<double> junc_cov_distribution = get_cov_distribution(junc);

		  //double d = distance(current_junc_cov_distribution,junc_cov_distribution);
		  
		  //double p = pearson(current_junc_cov_distribution,junc_cov_distribution);

		  double cos = cosine(current_junc_cov_distribution,junc_cov_distribution);
		  //double adcos = adcosine(current_junc_cov_distribution,junc_cov_distribution);
		  //double dt = dot(current_junc_cov_distribution,junc_cov_distribution);


		  /*
		  for(size_t j = 0;j<junc_cov_distribution.size();j++)
			  cout<<junc_cov_distribution[j]<<" ";
		  cout<<"; cov:"<<c<<"; cos:"<<cos<<"; dot:"<<dt<<endl;
		  */
		  edges_cov.push_back(make_pair(e,c*cos));

		  //edges_cov.push_back(make_pair(e,node_set[current_edge.first].parents[i].second));
	      }
	      sort(edges_cov.begin(),edges_cov.end(),sorter); //revision
	      if(edges_cov.size() == 1)
	      {
	          edge_t e = edges_cov[0].first;
		  current_path.insert(current_path.begin(), e.first);
		  paths_temp.push_back(current_path);
	      }
	      else 
	      {
	        for(size_t i=0;i<edges_cov.size();i++)
	        {
	          edge_t e = edges_cov[i].first;
		  //if(partial(e.first,e.second)) continue;
	          pair<edge_t ,edge_t> ee = make_pair(e,current_edge);
		  if(1 || packing_map.find(ee) != packing_map.end())
		  {
		      //current_path.push_back(e.second);
		      current_path.insert(current_path.begin(), e.first);
		      paths_temp.push_back(current_path);
		      break;
		  }
		}
	      }
	  }
      }
      return;
  }

  bool SimplifyGraph::check_overlap_of_2paths(path_t p1_, path_t p2_)
  {
    path_t p1,p2;
    if(p1_.size() <= p2_.size()) { p1 = p1_; p2= p2_;}
    else {p1 = p2_; p2 = p1_;}

    stringstream ss1,ss2;
    string s_p1, s_p2;
    for(size_t i=0;i<p1.size();i++) ss1<<p1[i]<<"_";
    s_p1 = ss1.str();
    for(size_t i=0;i<p2.size();i++) ss2<<p2[i]<<"_";
    s_p2 = ss2.str();

    if(s_p2.find(s_p1) != string::npos) return true;
    return false;   
  }
  void SimplifyGraph::get_components()
  {
        map<node_idx_t, bool> used_nodes;
	for(int i=0;i<Graph_size;i++)
	{
	    if(used_nodes[i]) continue;
	    used_nodes[i] = true;
	    vector<node_idx_t> comp_nodes;
	    comp_nodes.push_back(i);
	    for(size_t j=0;j<comp_nodes.size();j++)
	    {
		for(size_t k=0;k<node_set[comp_nodes[j]].children.size();k++)
		{
		    node_idx_t c = node_set[comp_nodes[j]].children[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), c) == comp_nodes.end())
		    {
			//cerr<<"  "<<c<<endl;
			comp_nodes.push_back(c);
			used_nodes[c] = true;
		    }
		}
		for(size_t k=0;k<node_set[comp_nodes[j]].parents.size();k++)
		{
		    node_idx_t p = node_set[comp_nodes[j]].parents[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), p) == comp_nodes.end())
		    {
		 	comp_nodes.push_back(p);
			used_nodes[p] = true;
		    }
		}
	    }
	    //Comps.push_back(comp_nodes);
	    nodes_of_components.push_back(comp_nodes);
 	}
  }
  void SimplifyGraph::update_graph(path_t path)
  {
      double mincov = 10000000;
      for(size_t i=0;i<path.size() -  1;i++)
      {
          int n1=path[i],n2=path[i+1];
	  double cov = node_set[n1].get_child_coverage(n2);
	  if(cov < mincov) mincov=cov;
      }
      for(size_t i=0;i<path.size() -  1;i++)
      {
          int n1=path[i],n2=path[i+1];
	  node_set[n1].reduced_child_coverage(n2,mincov/2);
	  node_set[n2].reduced_parent_coverage(n1,mincov/2);
      }
      return;
  }
  void SimplifyGraph::path_search(string strand, string chr)
  {
  get_components();
  //cout<<"components size: "<<nodes_of_components.size()<<endl;

  for(size_t I=0;I<nodes_of_components.size();I++)
  {
    vector<int> nodes = nodes_of_components[I];
    //vector<int> nodes;
    //for(int i=0;i<Graph_size;i++) nodes.push_back(i);
    vector<pair<int,int> > unused_junctions;
    vector<double> unused_junctions_cov;
    vector<int> unused_junctions_sample_number;
    for(int k = 0;k<nodes.size();k++)
    {
	int i = nodes[k];
        for(size_t j=0;j<node_set[i].children.size();j++)
	{
	    if(partial(i,node_set[i].children[j].first)) continue;

	    int c = node_set[i].children[j].first;
	    double cov = node_set[i].children[j].second;

	    pair<int,int> junc = make_pair(node_set[i].sequence.back(), node_set[c].sequence.front());
	    vector<pair<int,int> >::iterator it = find(junctions.begin(),junctions.end(),junc);
	    if(it == junctions.end()) continue;
	    int k = it - junctions.begin();
	    if(SampleSize > 10)
	    {
	      if(junctions_MappingInfo[k].size() == 1)
		    if(cov < 100) continue;
	    }
	    if(SampleSize > 20)
	    {
	      if(junctions_MappingInfo[k].size() <=2 ) if(cov < 50) continue;
	    }
	    if(junctions_MappingInfo[k].size() <=3 ){
	//    	    if(cov < 0.5) continue;
	    }
	    /*
	    if(junctions_MappingInfo[k].size() == 3)
	    {
	    	    if(cov < 0.6) continue;
	    }
	    */
	    pair<int,int> junc_ = make_pair(i,node_set[i].children[j].first);
	    unused_junctions.push_back(junc_);
	    unused_junctions_cov.push_back(node_set[i].children[j].second);
	    unused_junctions_sample_number.push_back(junctions_MappingInfo[k].size());

	}

    }
    if(unused_junctions.empty())
    {
        for(int k=0;k<nodes.size();k++)
	{
	    int i = nodes[k];
	    path_t p(1,i);
	    continue;
	    final_paths.push_back(p);
	}
    }
    //get pair of current components

    vector<path_t> pairpath_current;
    vector<double> pairpath_current_cov;
    for(size_t k=0;k<pairpath.size();k++)
    {
	path_t p = pairpath[k];
       	bool flag = false;
	for(size_t m=0;m<nodes.size();m++)
	{
	    if(find(p.begin(), p.end(),nodes[m]) != p.end())
	    {
		flag = true;
		break;
	    }
  	}
	if(1 || SampleSize > 1)
	{
	  int normal_num = 0, partial_num = 0;
 	  if(flag)
	  {
	    for(size_t m=0;m<p.size() - 1;m++)
	    {
		if(partial(p[m],p[m+1]) ) 
		    partial_num++;
		else normal_num++;
	    }
	    if(normal_num > 2 && partial_num < 2) flag = true;
	    else flag = false;
 	  }

	  if(flag){
	    for(size_t m=0;m<p.size() - 1;m++){
		if(!MultiJunction(p[m],p[m+1]) ){
		    flag = false;
		    break;
		}
	    }
	  }
	}
	//flag = false;
	if( flag){
	    pairpath_current.push_back(p);
	    pairpath_current_cov.push_back(pairpath_cov[k]);
	}
    }
    map<edge_t,bool> used_edges;
    while(0 && !pairpath_current.empty()) //pair as seed
    {
	vector<double>::iterator biggest = max_element(pairpath_current_cov.begin(), pairpath_current_cov.end());
	int k = biggest - pairpath_current_cov.begin();
	path_t pp = pairpath_current[k];

	if(pairpath_current_cov[k] == 0) break;
	else pairpath_current_cov[k] = 0;
        vector<path_t> forward_paths, reverse_paths;

	edge_t e = make_pair(pp[pp.size() - 2], pp[pp.size() - 1]);
	forward_extend(e,forward_paths);
	e = make_pair(pp[0],pp[1]);
	reverse_extend(e,reverse_paths);
	path_t final_path = reverse_paths.front();
	for(size_t i=2;i<pp.size();i++) final_path.push_back(pp[i]);
	for(size_t i=2;i<forward_paths.front().size();i++) final_path.push_back(forward_paths.front()[i]);
	for(size_t i=0;i<final_path.size()-1;i++){
	    edge_t edge = make_pair(final_path[i],final_path[i+1]);
	    used_edges[edge] = true;
	}

	final_paths.push_back(final_path);
//	path_info pi={seed_cov,normal_num,partial_num,max,min,ave,mj,Nmj,sample_number};
	path_info pi={1,1,1,1,1,1,1,1,1};
	final_paths_info.push_back(pi);
	for(size_t j=0;j<pairpath_current.size();){
	    if(check_overlap_of_2paths(pairpath_current[j],final_path))
	    {
	        pairpath_current.erase(pairpath_current.begin() + j);
		pairpath_current_cov.erase(pairpath_current_cov.begin() + j);
	    }
	    else j++;
	}

	
    }
    double seed_filter = SEED_Filter/SampleSize;
    //ofstream out("final.path.txt",ios::app);
    int path_number = 0;
    while(1 && !unused_junctions.empty())//junction as seed
    {
        vector<double>::iterator biggest = max_element(unused_junctions_cov.begin(),unused_junctions_cov.end());
	int k = biggest - unused_junctions_cov.begin();
	double seed_cov = unused_junctions_cov[k];
	edge_t e = unused_junctions[k];
	int sample_number = unused_junctions_sample_number[k];
	if(unused_junctions_cov[k] < seed_filter) break;
	else unused_junctions_cov[k] = 0;	
	
	if(path_number >= 4) seed_filter = 5*SEED_Filter/SampleSize;
	if(path_number >= 5 ) seed_filter = 10*SEED_Filter/SampleSize;
	if(path_number >= 6 ) seed_filter = 15*SEED_Filter/SampleSize;
	if(path_number >= 7 ) seed_filter = 20*SEED_Filter/SampleSize;
	if(path_number >= 8 ) seed_filter = 50*SEED_Filter/SampleSize;
	
	
	if(used_edges.find(e) != used_edges.end()) continue;

	vector<path_t> forward_paths, reverse_paths;
	forward_extend(e,forward_paths);
	reverse_extend(e,reverse_paths);
	path_t final_path = reverse_paths.front();
	for(size_t i=2;i<forward_paths.front().size();i++) final_path.push_back(forward_paths.front()[i]);
	for(size_t i=0;i<final_path.size()-1;i++){
	    edge_t edge = make_pair(final_path[i],final_path[i+1]);
	    used_edges[edge] = true;
	}
	int partial_num = 0, normal_num = 0;
	double max = 0.0,min = 10000000.0,AllCov=0;
	int Nmj = 0,mj=0;
	for(size_t i=0;i<final_path.size()-1;i++)
	{
	    int n1 = final_path[i], n2 = final_path[i+1];
	    if(partial(n1,n2)){
		    partial_num++;
		    continue;
	    }

	    if(!MultiJunction(n1,n2)) Nmj++;
	    else mj++;

	    normal_num++;
	    double cov = node_set[n1].get_child_coverage(n2);
	    if(max < cov ) max = cov;
	    if(min > cov) min = cov;
	    AllCov += cov;
	}
	//if(M > 1) continue;
	//if(path_number >= 5 && max>20*min) continue;
	//if(path_number >= 5 && partial_num>5) continue;
	/*
	out<<chr<<"TRANSMETA."<<rg_index<<" "<<seed_cov<<" "<<final_paths.size()+1<<": ";
	for(size_t i=0;i<final_path.size()-1;i++)
	{
	    int n1 = final_path[i], n2 = final_path[i+1];
	    if(partial(n1,n2)) out<<"(P) ";
	    else out<<"(nP) ";
	    out<<n1<<"->"<<n2<<" "<<node_set[n1].get_child_coverage(n2)<<"; ";
	    //out<<node_set[n1].get_child_coverage(n2)<<"; ";
	}
	out<<endl;
	*/
	double ave = double(AllCov/(1.0*normal_num));
	path_info pi={seed_cov,normal_num,partial_num,max,min,ave,mj,Nmj,sample_number};
	final_paths.push_back(final_path);
	final_paths_info.push_back(pi);
	//update_graph(final_path);
	path_number++;
    }
    //out<<endl;
    output(strand,chr);
    final_paths.clear();
    final_paths_info.clear();
 }
    return;
  }
  void SimplifyGraph::output(string strand, string chr)
  {
    //L = 800, COV1 = 20, COV2=30;//ipac2
    int L = 800;
    double COV1 = 20 ,COV2 = 30;
    //double COV1 = 20 ,COV2 = 25;
    rg_index++;
    trans_id = 1;
    for(size_t i=0;i<final_paths.size();i++)
    {
	//if(final_paths.size() == 1) continue;
	//if(i != 0) continue;
	//if(final_paths.size() >9 && i>=2) break;
	if(i >= 6) break;  
      	path_t p = final_paths[i];
	path_info pi = final_paths_info[i];
	if(final_paths.size() == 1)
	{
		//if(pi.seed_cov == pi.min_cov && pi.max_cov > pi.min_cov*100.0)
		if(pi.seed_cov < 0.5)
			continue;;
	}
	else if(final_paths.size() == 2){
		//if(pi.seed_cov < 0.1) continue;
	}
	int length = 0;
	for(size_t j=0;j<p.size();j++) length += node_set[p[j]].length();

	double cov = 0;
	for(size_t j=0;j<p.size();j++) cov += node_set[p[j]].coverage()*node_set[p[j]].length();
	cov = double(cov/length);

	//if(cov < 2 && final_paths.size() == 1) continue;
	if(unstranded && p.size() == 1 && strand != ".") continue;//single
	if(p.size() == 1 && ( !(length>L && cov>COV1) && !(cov>COV2) )) continue; //single

	double single_flag = true;
	for(int j=0;j<int(p.size()) - 1;j++){
	    int intron_len = node_set[p[j+1]].sequence.front() - node_set[p[j]].sequence.back();
	    if(intron_len > 1){
	        single_flag = false;
		break;
	    }
	}
	if(unstranded && single_flag && cov <= 10) continue; //single
	if(single_flag && cov <= 10) continue; //single Merge

	if(length >500 || (final_paths.size() == 1 && cov > 10))
	{
      	  out_gtf<<chr<<"	"<<"TRANSMETA"<<"	"
	      	<<"transcript"<<"	"<<node_set[p.front()].sequence.front()<<"	"<<node_set[p.back()].sequence.back()<<"	"
	      	<<1000<<"	"<<strand<<"	"
	      	<<".	gene_id "<<"\""<<chr<<"."<<rg_index<<"\""<<"; "
	      	<<"transcript_id "<<"\""<<chr<<"."<<rg_index<<"."<<trans_id<<"\";"<<endl;
      	  for(size_t j=0;j<p.size();j++)
	  {
             out_gtf<<chr<<"	"<<"TRANSMETA"<<"	"
		    <<"exon"<<"	"<<node_set[p[j]].sequence.front()<<"	"<<node_set[p[j]].sequence.back()<<"	"
		    <<1000<<"	"<<strand<<"	"
		    <<".	gene_id "<<"\""<<chr<<"."<<rg_index<<"\""<<"; "
		    <<"transcript_id "<<"\""<<chr<<"."<<rg_index<<"."<<trans_id<<"\"; "
		    <<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
      	  }

	  out_info<<"Info: \""<<chr<<"."<<rg_index<<"."<<trans_id<<"\"; "
		  <<pi.seed_cov<<" "<<pi.normal_edge_number<<" "<<pi.partial_edge_number<<" "
		  <<pi.max_cov<<" "<<pi.min_cov<<" "<<pi.ave_cov<<" "
		  <<pi.mj_number<<" "<<pi.Nmj_number<<" "<<pi.seed_sample_number<<" "
		  <<" :  "<<final_paths.size()<<endl;
	  trans_id++;
	}
    }

    return;
  }
