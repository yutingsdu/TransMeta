#ifndef PAIRPATH_H
#define PAIRPATH_H


#include<iostream>
#include<algorithm>
#include<map>
#include<vector>

#include <boost/unordered_map.hpp>
//#include"Find_junc_last_map.h"

using namespace std;
extern int rg_index;
int check_graph = 379;
typedef vector<int> path_t;
typedef pair<path_t,path_t> pair_path_t;

//extern multimap<string, vector< pair<int,int> > > readid_cigar_map;

typedef  vector<pair<int,int> > cigar_t;
typedef  pair<int,int> junc_type;
//extern map<string, pair<cigar_t,cigar_t> > readid_pairInfo_map;

class PairPath
{
private:
  class Node
  {
     public: 
	pair<int,int> boundary;
	int length;
	vector<int> children,parents;

    public:
	Node(int i,int j)
 	{
	    boundary = make_pair(i,j);
	    length = j - i + 1;
	}
	void add_child(int child) 
	{
	    children.push_back(child);
	}
	void add_parent(int parent)
	{
	    parents.push_back(parent);
	}
     
  };
public:
  vector<Node> node_set;
  size_t size;
  vector<int>exon_l,exon_r,junc_l,junc_r;
  vector<pair<int,int> > junc;

  map<pair_path_t,int> pairpath_map;

  vector<path_t> final_pair_path;


  vector<pair<int,int> > reserved_junc; //the junction that should NOT delete in following step

  boost::unordered_map<vector<int>, int> triplet_map;
  boost::unordered_map<vector<int>, int>::iterator triplet_map_iter;

public:
  PairPath(vector<int>& exonl,vector<int>& exonr, vector<int>& juncl,vector<int>& juncr){
	exon_l = exonl; exon_r = exonr;
	junc_l = juncl; junc_r = juncr; 
	for(int i=0;i<juncl.size();i++)
	{
	    pair<int,int> p = make_pair(juncl[i],juncr[i]);
	    junc.push_back(p);
	}
  };
  void show_graph()
  {
    cerr<<"raw_graph:"<<endl;
    cerr<<"edge:"<<endl;
     for(int i=0;i<size;i++)
    {
	for(int j=0;j<node_set[i].children.size();j++)
	{
	    cerr<<i<<"->"<<node_set[i].children[j]<<endl;
 	}
    }
    cerr<<"node:"<<endl;
    for(int i=0;i<size;i++) cerr<<i<<" "<<node_set[i].boundary.first<<" "<<node_set[i].boundary.second<<endl;
    cerr<<"====================================="<<endl;
  }
  void build()
  {
    size = exon_l.size();
    vector<int> temp;
    vector< vector<int> > out_jcidx_of_exon(size,temp);

    for(int i=0;i<exon_l.size();i++)
    {
        for(int j=0;j<junc_l.size();j++)
        {
            if(exon_r[i] == junc_l[j] - 1) {
                out_jcidx_of_exon[i].push_back(j);
            }
        }
    }

    for(int i=0;i<exon_l.size();i++)
    {
 	Node node(exon_l[i],exon_r[i]);

	vector<int> out_jcidx = out_jcidx_of_exon[i];
	for(int j = 0;j<out_jcidx.size();j++)
	{
	    for(int k = i+1;k<exon_l.size();k++){
		if(exon_l[k] == junc_r[out_jcidx[j]]) node.add_child(k);
	    }
	}
	node_set.push_back(node);
    }

    return;
  }
  bool exist_junc(pair<int,int> j)
  {
    if(find(junc.begin(),junc.end(),j) == junc.end()) return false;
    else return true;
  }

  bool get_supported_node_of_oneseg(pair<int,int> seg,vector<int>& subpath)
  {
    //cerr<<"seg: "<<seg.first<<" "<<seg.second<<endl;
    for(int i=0;i<size;i++)
    {
	if(exon_l[i] > seg.first) return false;

	if(exon_l[i]<= seg.first && exon_r[i] > seg.first)
	{
	
	    subpath.push_back(i);
	    //cerr<<"exon "<<exon_l[i]<<" "<<exon_r[i]<<endl;
	    if(exon_r[i] >= seg.second) return true;  //segment all in one exon
	    else
	    {
		for(int k = i;k<size-1;k++)
		{
		    if(exon_l[k+1] - exon_r[k] >1) return false;
		    pair<int,int> junc = make_pair(exon_r[k] + 1,exon_l[k+1]);
	//	    if(rg_index == check_graph ) cerr<<"hhhh "<<k<<" "<<k+1<<" "<<exon_r[k] + 1<<" "<<exon_l[k+1]<<endl;
		    if(exist_junc(junc))
		    {
			subpath.push_back(k+1);
			if(exon_r[k+1] >= seg.second) return true;
		    }
		    else return false;
		    
		}
	    }

	}
    }//for each exon finish
    return false;
  }

  bool find_supported_path_of_oneread(vector<pair<int,int> > read_cigar,vector<int>& path)
  {
    for(int i=0;i<read_cigar.size();i++)
    {
	vector<int> subpath;
	if(get_supported_node_of_oneseg(read_cigar[i],subpath))
	{
	    if( !path.empty())
	    {	
		pair<int,int> junc = make_pair(exon_r[path.back()] + 1,exon_l[subpath.front()]);
//		if(rg_index == check_graph ) cerr<<"HERE "<<path.back()<<" "<<subpath.front()<<" -- "<<junc.first<<" "<<junc.second<<endl;	
		if( !exist_junc(junc) ) return false; 
	    }

	    for(int k=0;k<subpath.size();k++) path.push_back(subpath[k]);
	}
	else return false;
	
    }
    return true;
  }
  void show_path(vector<int> v,vector<int> v2)
 { 
    for(int i=0;i<v.size();i++) cerr<<v[i]<<"->";
    cerr<<"...";
    for(int i=0;i<v2.size();i++) cerr<<v2[i]<<"->";
    cerr<<endl;
  }

  void find_pair_path(boost::unordered_map<string, pair<cigar_t,cigar_t> >& readid_pairInfo_map)
  {
    build();

    //typedef multimap<string, vector< pair<int,int> > >::iterator iter;
    //map<string, pair<cigar_t,cigar_t> > readid_pairInfo_map;
    typedef boost::unordered_map<string, pair<cigar_t,cigar_t> >::iterator iter;

    for(iter i = readid_pairInfo_map.begin();i != readid_pairInfo_map.end();i++)
    {
	//cerr<<i->first<<endl;

	vector<int> left_path,right_path;
	if(i -> second.first.empty() || i->second.second.empty()) continue;

	vector< pair<int,int> >  left = i->second.first;
/*
	if(readid_pairInfo_map.size() == 424) 
	{
	  for(int j=0;j<left.size();j++) cerr<<left[j].first<<" L ->"<<left[j].second<<endl;
	  cerr<<endl;
	}
*/
	vector< pair<int,int> >  right = i->second.second;
/*
	if(readid_pairInfo_map.size() == 424) 
	{
	  for(int j=0;j<right.size();j++) cerr<<right[j].first<<" R -> "<<right[j].second<<endl;
	  cerr<<endl;
	}
*/	
	if(left.front().first > right.front().first) continue;

	bool left_flag = find_supported_path_of_oneread(left,left_path);
	bool right_flag = find_supported_path_of_oneread(right,right_path);
	
	//if(rg_index == check_graph )   show_path(left_path,right_path);
	//if(rg_index == check_graph )   cerr<<"&&&"<<endl<<endl;

	if(left_flag && right_flag)
	{
	    pair_path_t p = make_pair(left_path,right_path);
	    if(pairpath_map.find(p) == pairpath_map.end())
	       pairpath_map[p] = 1;
	    else pairpath_map[p] += 1;

	    i -> second.first.clear(); i->second.second.clear();
	}

    }

  }
  int path_length(path_t p)
  {
    int length = 0;
    for(int i=1;i<p.size();i++)
	length += node_set[p[i]].length;

    return length;
  }

  bool get_path_between_2node(int n1,int n2, path_t& path)
  {
    path_t temp(1,n1);
    vector<path_t> vec_path(1,temp);

    vector<path_t> vec_path_;

    while( !vec_path.empty())
    {
 	path_t curr_path = vec_path.front();
	//for(int i = 0;i<curr_path.size();i++)   cerr<<curr_path[i]<<" --- ";
	//cerr<<endl;
	int q = curr_path.back();

	vec_path.erase(vec_path.begin());

 	vector<int> q_children = node_set[q].children;

	for(int i=0;i<q_children.size();i++)
	{
	    if(q_children[i] == n2)
	    {
		//curr_path.push_back(n2);
		//vec_path_.push_back(curr_path);
		path_t path_temp = curr_path;//NEW1
		path_temp.push_back(n2);//NEW1
		vec_path_.push_back(path_temp);//NEW1

		if(vec_path_.size() >= 2) return false;
	
		continue;
	    }
	    if(q_children[i] > n2) continue;//NEW1 nodes sorted, if q_children[i]>n2 then never find a path to n2

	    path_t path_temp = curr_path;
	    path_temp.push_back(q_children[i]);

	    if( path_length(path_temp) > 250) continue;//note path_length(); //NEW8.14

	    vec_path.push_back(path_temp);
	    
	}
	if(vec_path.size() > 1000) return false;
    }
    if(vec_path_.size() == 1) 
    {
	path = vec_path_.front();
	//if(path.size() == 2) // NEW8.14
	    return true;
    }
    return false;
  }
  void connect_pair_path(path_t& left,path_t& right, 
			 map<pair<int,int>,path_t>& node_which_has_find_path)
  {
   if(right.front() <= left.back())
   {

	path_t::iterator it = find(left.begin(),left.end(),right.front());

	if( it!= left.end())
	{

	    int i = it - left.begin();
	    int k = 0;
	    for(;k<left.size() - i;k++)
	    {
		if( left[i + k] == right[k] ) continue;
		else 
		{
		   if(left.size() > 2) final_pair_path.push_back(left);
		   if(right.size() > 2) final_pair_path.push_back(right);
		   return;
		}
	    }
	    for(;k<right.size();k++) left.push_back(right[k]);
	    if(left.size() > 2) final_pair_path.push_back(left);
	}
    } 
    else {
	//??? what
	//return;//NEW8.14
        path_t path_between_2node;

	pair<int,int> p = make_pair(left.back(),right.front());
	if(node_which_has_find_path.find(p) != node_which_has_find_path.end())
	{
	    path_between_2node = node_which_has_find_path[p];
	 

	    for(int i = 1;i<path_between_2node.size()-1;i++) left.push_back(path_between_2node[i]);
	    for(int i=0;i<right.size();i++) left.push_back(right[i]);

	    if(left.size() > 2)  final_pair_path.push_back(left);
	}
	else if(get_path_between_2node(left.back(),right.front(),path_between_2node))
	{
	  node_which_has_find_path[p] = path_between_2node;

	  for(int i = 1;i<path_between_2node.size()-1;i++) left.push_back(path_between_2node[i]);
	  for(int i=0;i<right.size();i++) left.push_back(right[i]);
	  if(left.size() > 2) final_pair_path.push_back(left);
	}
	
    }
    return;
  }
  void show_final_path()
  {
    cerr<<"final path -----"<<endl;
    for(int i=0;i<final_pair_path.size();i++) 
    {
	path_t p = final_pair_path[i];
	for(int j=0;j<p.size();j++) cerr<<p[j]<<" ";
	cerr<<endl;
    }
  }
  void get_triplet_hash()
  {
	for(int i=0;i<final_pair_path.size();i++)
	{
	    path_t path = final_pair_path[i];
	    for(int j=0;j<=path.size()-3;j++)
	    {
		vector<int> triplet;
		triplet.push_back(path[j]); triplet.push_back(path[j+1]); triplet.push_back(path[j+2]);

		if(triplet_map.find(triplet) != triplet_map.end())
		{
		    triplet_map[triplet] ++;
		}
	  	else triplet_map[triplet] = 1;
	    }

	}
  }
  void get_credible_pair_path_by_tripletnum()
  {
	//get_strict_triplet_hash();

	vector<path_t> final_ppath_temp = final_pair_path;
	final_pair_path.clear();

 	while(!final_ppath_temp.empty())
	{
	    path_t temp = final_ppath_temp.front();
	    vector<int> flag(temp.size(),1);
	    //remove incredible triplet;
	    //
	    for(int i = 0;i <= temp.size() - 3;i++)
	    {
		vector<int> triplet;
		triplet.push_back(temp[i]); triplet.push_back(temp[i+1]); triplet.push_back(temp[i+2]);

		//if(triplet_map[triplet] < 2) flag[i+1] = 0;
		if(triplet_map.find(triplet) == triplet_map.end()) flag[i+1] = 0;

	    }
	    //1 2 3 4 5 6 7 8 9
	    //1 0 1 1 1 0 0 1 1
	    path_t c_path;
	    for(int i = 0;i < flag.size();i++)
	    {
		if(flag[i] == 1) c_path.push_back(temp[i]);
		else{
		    c_path.push_back(temp[i]);
		    if(c_path.size() > 2) final_pair_path.push_back(c_path);

		    c_path.clear();
		    c_path.push_back(temp[i]);
		}
	    }
	    if(c_path.size() > 2) final_pair_path.push_back(c_path);
	    final_ppath_temp.erase(final_ppath_temp.begin());
	}
  }
  void get_credible_pair_path_by_pathlength()
  {
    for(int i=0;i<final_pair_path.size();)
    {
	vector<int> current_path = final_pair_path[i];
	current_path.erase(current_path.begin());
	current_path.pop_back();
	if(path_length(current_path) > 200)
	{
	    final_pair_path.erase(final_pair_path.begin() + i);
	}	
	else i++;
    }
  }
  void get_strict_triplet_hash(path_t left,path_t right) 
  {
    //boost::unordered_map<vector<int>, int> triplet_map;
    //boost::unordered_map<vector<int>, int>::iterator triplet_map_iter;
    if(left.empty() || right.empty()) return;

    if(right.front() <= left.back()){//left: 1-2-3 right: 2-3-4 
        path_t::iterator it = find(left.begin(),left.end(),right.front());

        if( it!= left.end())
        {
            int i = it - left.begin();
            int k = 0;
            for(;k<left.size() - i;k++)
            {
                if( left[i + k] == right[k] ) continue;
                else return;
            }
            for(;k<right.size();k++) left.push_back(right[k]);

            if(left.size() > 2) {
		for(int j=0;j<=left.size()-3;j++)
		{
		    vector<int> triplet;
		    triplet.push_back(left[j]); triplet.push_back(left[j+1]); triplet.push_back(left[j+2]);

		    if(triplet_map.find(triplet) != triplet_map.end()) triplet_map[triplet] ++;
		    else triplet_map[triplet] = 1;
		}
	    }
        }
 	return;
    }

    pair<int,int> check_junc = make_pair(exon_r[left.back()]+1,exon_l[right.front()]);

    vector< pair<int,int> >::iterator it = find(junc.begin(),junc.end(),check_junc);
    if(it == junc.end())//no direct junction from left to right left:1-2-3 .....right:4-5-6
    {
	if(left.size() > 2) {
		for(int j=0;j<=left.size()-3;j++)
		{
		    vector<int> triplet;
		    triplet.push_back(left[j]); triplet.push_back(left[j+1]); triplet.push_back(left[j+2]);

		    if(triplet_map.find(triplet) != triplet_map.end()) triplet_map[triplet] ++;
		    else triplet_map[triplet] = 1;
		}
	}
	if(right.size() > 2){
		for(int j=0;j<=right.size()-3;j++)
		{
		    vector<int> triplet;
		    triplet.push_back(right[j]); triplet.push_back(right[j+1]); triplet.push_back(right[j+2]);

		    if(triplet_map.find(triplet) != triplet_map.end()) triplet_map[triplet] ++;
		    else triplet_map[triplet] = 1;
		}
	}
    } 
    else//exist direct junction from left to right left:1-2-3 ..(junc)->..right:4-5-6
    {
	for(int i=0;i<right.size();i++) left.push_back(right[i]);
	if(left.size() > 2) {
		for(int j=0;j<=left.size()-3;j++)
		{
		    vector<int> triplet;
		    triplet.push_back(left[j]); triplet.push_back(left[j+1]); triplet.push_back(left[j+2]);

		    if(triplet_map.find(triplet) != triplet_map.end()) triplet_map[triplet] ++;
		    else triplet_map[triplet] = 1;
		}
	}
	
    }
     
  }
  void get_credible_pair_path_by_juncsupport()
  {
	
  }
  void connect_pair_path()
  {
    //cerr<<"H1"<<endl;
    map<pair<int,int>,path_t> node_which_has_find_path;
    typedef map<pair_path_t,int>::iterator iter;

    for(iter i = pairpath_map.begin();i!=pairpath_map.end();i++)
    {
	path_t left = i->first.first,right = i->first.second;

	if(left == right && left.size() <= 2) continue;

	get_strict_triplet_hash(left,right);

	connect_pair_path(left,right,node_which_has_find_path);
    }

    //get credible pair path according to triplet
    //get_credible_pair_path_by_pathlength();
    //cerr<<"get credible pair path according to triplet: "<<endl;
    
    get_credible_pair_path_by_tripletnum();

    sort(final_pair_path.begin(),final_pair_path.end());
    final_pair_path.erase( unique(final_pair_path.begin(),final_pair_path.end()), final_pair_path.end());

    //?????? what
    //get credible pair path by node number
    /*
    for(int i=0;i<final_pair_path.size();)
    {
	if(final_pair_path[i].size() > 4) final_pair_path.erase(final_pair_path.begin() + i);
	else i++;
    }
    */
    
    
    //get_credible_pair_path_by_pathlength();

    for(int i=0;i<final_pair_path.size();)
    {
	bool delete_flag = false;
	for(int j = 0;j<final_pair_path.size();j++)
	{
	    if(j == i) continue;
	    if( path_compatible(final_pair_path[i],final_pair_path[j]) )
	    {
		final_pair_path.erase(final_pair_path.begin() + i);
		delete_flag = true;
		break;
	    }
	}
	if( !delete_flag ) i++;
    }
    //show_final_path();
    return;
  }

  bool path_compatible(path_t P1,path_t P2)//p1.size() < p2.size() !!!
  {
        if(P1 == P2) return true;
	if(P1.size() >= P2.size()) return false;
        if(P1.empty() || P2.empty()) return false;


        path_t::iterator i1 = find(P2.begin(),P2.end(),P1.front());
        if(i1 == P2.end()) return false;
        path_t::iterator i2 = find(P2.begin(),P2.end(),P1.back());
        if(i2 == P2.end()) return false;

        if(i1>i2) return false;

        path_t::iterator i3=P1.begin();


        while( i1 <= i2 && i3!= P1.end())
        {
            if((*i1) == (*i3))
            {
                i1++;
                i3++;
            }
            else return false;
        }
	return true;
    }

  vector< vector<int> > get_final_path_exon()
  {
    vector< vector<int> > final_pair_path_exon;
    for(int i=0;i<final_pair_path.size();i++)
    {
	    vector< int > exon;
	    path_t path = final_pair_path[i];
	    for(int j=0;j<path.size();j++)
	    {
		int n = path[j];
		exon.push_back(exon_l[n]);
		exon.push_back(exon_r[n]);
	    }
	    
	    final_pair_path_exon.push_back(exon);
    }
    return final_pair_path_exon;
  }
  void get_reserved_junc()
  {
    for(int i = 0;i < final_pair_path.size(); i++)
    {
	path_t curr_ppath = final_pair_path[i];
 	for(int j = 0;j<curr_ppath.size() - 1;j++)
	{
	    int n1 = curr_ppath[j], n2 = curr_ppath[j+1];

	    pair<int,int> p = make_pair(exon_r[n1] + 1, exon_l[n2]);
	    reserved_junc.push_back(p);
	}
    }
    sort(reserved_junc.begin(),reserved_junc.end());
    reserved_junc.erase(unique(reserved_junc.begin(),reserved_junc.end()),reserved_junc.end());

    return;
  }

};
#endif
