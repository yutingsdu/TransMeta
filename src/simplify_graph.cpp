
#include <iostream> 
#include <map>
#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <fstream>
#include "simplify_graph.h"

using namespace std;
extern ofstream out_junction_read;
extern map<pair<int,int>, vector<string> > plus_junction_readid_map;
extern map<pair<int,int>, vector<string> > minus_junction_readid_map;
extern map<pair<int,int>, vector<string> > unstranded_junction_readid_map;
extern string strand_specific;
extern int pack_graph_num;
extern int unpack_graph_num;
extern bool PackingFlag;
extern bool MyFlag;
extern double AVERAGE_REMOVE_RATE;
extern double UNBALANCE_RATE ;
extern bool SFlag;
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
  void SimplifyGraph::outread_empirical()
  {
    return;
    for(int i=0;i<Graph_size;i++)
    {
	node_idx_t n1 = i;
	for(size_t j=0;j<node_set[i].children.size();j++)
	{
	  node_idx_t n2 = node_set[i].children[j].first;
	  //if(node_set[i].children[j].second > 3) continue;
	  if(!node_set[n1].sequence.empty() && !node_set[n2].sequence.empty() )
	  {
	    pair<int,int> junc = make_pair(node_set[n1].sequence.back(), node_set[n2].sequence.front());
	    if(Strand == "+"){
		//out_junction_read<<">"<<Strand<<" "<<Chr<<" "<<junc.first<<" "<<junc.second<<'\n';
		map<pair<int,int>, vector<string> >::iterator it = plus_junction_readid_map.find(junc);
		if(it != plus_junction_readid_map.end()){
		   for(size_t j=0;j<it->second.size();j++){
			//out_junction_read<<it->second[j]<<" ";
			

		   }
		}
		//out_junction_read<<'\n';
	    }
	    else if(Strand == "-"){
		//out_junction_read<<">"<<Strand<<" "<<Chr<<" "<<junc.first<<" "<<junc.second<<'\n';
		map<pair<int,int>, vector<string> >::iterator it = minus_junction_readid_map.find(junc);
		if(it != minus_junction_readid_map.end()){
		    for(size_t j=0;j<it->second.size();j++){
			    //out_junction_read<<it->second[j]<<" ";
		    }
		}
		//out_junction_read<<'\n';
	    }
	  }
	}
    }
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
//return;
    int partial_number = 0;
    for(int i=0;i<Graph_size-1;i++)
    {
	if(partial(i,i+1)) partial_number++;
    }
    if(partial_number < 10) get_reserved_junc(PairPath);
    /*
    cout<<"size: "<<minus_junction_readid_map.size()<<endl;
    for(map<pair<int,int>, vector<string> >::iterator it = minus_junction_readid_map.begin(); it != minus_junction_readid_map.end();it++){
        cout<<it->first.first<<" "<<it->first.second<<" : ";
	for(size_t i=0;i<it->second.size();i++) cout<<it->second[i]<<",";
	cout<<endl;
    }
    */
    if(SFlag) reserved_junc.clear(); //3.1 spkied
    //AVERAGE_REMOVE_RATE= 0.04;
    //if(PackingFlag) AVERAGE_REMOVE_RATE = 0.01;
    //if(partial_number == 0) AVERAGE_REMOVE_RATE = 0.04;
    //if(partial_number > 1 && Graph_size > 30) AVERAGE_REMOVE_RATE = 2;
    vector<edge_t> delete_edges_child,delete_edges_parent;
    get_CovInfo();
    
    if(SFlag) PackingFlag = true;//3.1 spkied
    //PackingFlag = false;
	remove_unbalance_edges(delete_edges_child,delete_edges_parent);

	remove_partial_junction(delete_edges_child,delete_edges_parent);

    delete_children_edges(delete_edges_child);
    delete_parents_edges(delete_edges_parent);
    {
      if(PackingFlag) AVERAGE_REMOVE_RATE = 0.02;
      else AVERAGE_REMOVE_RATE = 0.04;
    }
    if(strand_specific == "unstranded") AVERAGE_REMOVE_RATE = 0.02; //for simulation star transref 
    if(SFlag) AVERAGE_REMOVE_RATE = 0.015;//3.1 spkied
    //if(!PackingFlag) 

   if(!MyFlag) remove_edges_by_average_coverage(delete_edges_child,delete_edges_parent);

    delete_children_edges(delete_edges_child);
    delete_parents_edges(delete_edges_parent);
	remove_partial_end_by_edge_coverage(delete_edges_child,delete_edges_parent);

   if(!PackingFlag) 
    remove_lowcov_edges_of_bifurcation_nodes(delete_edges_child,delete_edges_parent);//NewNote

      remove_intron_contamination(delete_edges_child,delete_edges_parent);

//   remove_edges_onemulti(delete_edges_child,delete_edges_parent);
   
   //remove_small_exons(delete_edges_child,delete_edges_parent);
   
    delete_children_edges(delete_edges_child);
    delete_parents_edges(delete_edges_parent);
    delete_edges_child.clear(); delete_edges_parent.clear();
    if(partial_number > 1)// && Graph_size > 30){
    {
     	remove_small_exons(delete_edges_child,delete_edges_parent);
	//remove_edges_onemulti(delete_edges_child,delete_edges_parent);	
        delete_children_edges(delete_edges_child);
        delete_parents_edges(delete_edges_parent);
    }
    outread_empirical();
 //   if(Strand == "+") plus_junction_readid_map.clear();
 //   else if(Strand == "-") minus_junction_readid_map.clear();
    //get_CovInfo();
    //show_graph();
    //cout<<"Contract: "<<endl;
    if(0 && PackingFlag) //ipac2
	    contract_graph();
    //show_graph();
    for(int i=0;i<Graph_size;i++) {
        if(i>=RawSize && node_set[i].children.empty() && node_set[i].parents.empty()) Single_node.push_back(i);
    }
    //ipac2
    for(int i=0;i<Graph_size;i++) {
	if(node_set[i].children.empty() && node_set[i].parents.empty())
	{
	    if(find(Single_node.begin(),Single_node.end(),i) == Single_node.end())
	    {
		invalid_nodes.push_back(i);
	    }
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
	    if(partial(i,c) && node_set[i].children[j].second < 10)
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
	    if(ave < 2) 
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
		if(node_set[i].children[j].second > 5 && keep_edge(make_pair(i,c)) )continue; 
		if(rate < AVERAGE_REMOVE_RATE)
	 	{
		    delete_edges_child.push_back(make_pair(i,node_set[i].children[j].first));
		}
	    }
 	}
	return;
  }
  void SimplifyGraph::remove_partial_junction(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	double R=0;
	if(!PackingFlag) R = 0.12;
	else R=0.08;
	if(SFlag) R=0.5;//3.1 spk
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
			if(node_set[i].children[j].second > 5 && keep_edge(make_pair(i,c)) ) continue;
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
			if(node_set[i].parents[j].second > 5 && keep_edge(make_pair(p,i)) ) continue;
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
		//if(node_set[i].coverage() <= 2  && ecov < 3 && r1<0.2 &&  node_set[c].parents.size() > 1)
		if(node_set[i].coverage() < 5  && ecov < 5 && r1<0.2 && r2<0.2 && node_set[c].parents.size() > 1 )
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
		//if(node_set[i].coverage() <=2 && ecov < 3 && r1<0.15 && node_set[p].children.size() > 1)
		if(node_set[i].coverage() < 5 && ecov < 5 && r1<0.2 && r2<0.2 && node_set[p].children.size() > 1)
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
		    if( (node_set[i].coverage() < 10 && node_set[i].coverage()<pc_cov) )
		      
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
                           && r<=R && len < 20) 
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
                           && r<=R && len < 20) //1,20
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
		if(cov <= 1  ||(cov <= 10 && double(cov/max_cov) <= 0.5))
	 	{
		  //if( node_set[c].children.empty() || node_set[c].parents.size() > 1)
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
	    	if(cov <= 1  ||(cov <= 10 && double(cov/max_cov) <= 0.5))
	    	{
		  //if(node_set[p].parents.empty() || node_set[p].children.size() > 1)
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
  
  void SimplifyGraph::remove_lowcov_edges_of_bifurcation_nodes(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    double R = 0.18;
    double COV = 5;
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
		//cerr<<c<<": "<<cov<<"; ";
		double r = double(cov/max_cov);
		if(cov < 2 || (cov <= COV && r <= R)) //5,0.18
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
		//cerr<<p<<": "<<cov<<"; ";
		double r = double(cov/max_cov);
		if(cov < 2 || (cov <= COV && r <= R))//5 0.18
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
	node_set[n1].delete_child(n2);
	node_set[n2].delete_parent(n1);
	/*
	if(!node_set[n1].sequence.empty() && !node_set[n2].sequence.empty() )
	{
	    pair<int,int> junc = make_pair(node_set[n1].sequence.back(), node_set[n2].sequence.front());
	    if(Strand == "+"){
		map<pair<int,int>, vector<string> >::iterator it = plus_junction_readid_map.find(junc);
		if(it != plus_junction_readid_map.end()){
		    for(size_t j=0;j<it->second.size();j++) out_junction_read<<">"<<(it->second)[j]<<endl;
		}
	    }
	    else if(Strand == "-"){
		map<pair<int,int>, vector<string> >::iterator it = minus_junction_readid_map.find(junc);
		if(it != minus_junction_readid_map.end()){
		    for(size_t j=0;j<it->second.size();j++) out_junction_read<<">"<<it->second[j]<<endl;
		}
	    }
	}
	*/
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
	node_set[n1].delete_child(n2);
	node_set[n2].delete_parent(n1);
	/*
	if(!node_set[n1].sequence.empty() && !node_set[n2].sequence.empty() )
	{
	    pair<int,int> junc = make_pair(node_set[n1].sequence.back(), node_set[n2].sequence.front());
	    if(Strand == "+"){
		map<pair<int,int>, vector<string> >::iterator it = plus_junction_readid_map.find(junc);
		if(it != plus_junction_readid_map.end()){
		   for(size_t j=0;j<it->second.size();j++) out_junction_read<<">"<<it->second[j]<<endl;
		}
	    }
	    else if(Strand == "-"){
		map<pair<int,int>, vector<string> >::iterator it = minus_junction_readid_map.find(junc);
		if(it != minus_junction_readid_map.end()){
		    for(size_t j=0;j<it->second.size();j++) out_junction_read<<">"<<it->second[j]<<endl;
		}
	    }
	}
	*/
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
