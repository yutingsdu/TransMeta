#ifndef JUNCTION_simplifyTHS
#define JUNCTION_simplifyTHS

//junction_paths.h

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include<algorithm>
#include<cmath>
#include <map>
#include "../simplify_graph.h"
using namespace std;
extern bool PackingFlag;
typedef vector<int> path_t;
typedef pair<int,int> edge_t;
bool big_enough(vector<pair<int,double> > v,bool cov)
{
    for(size_t i=0;i<v.size();i++)
    {
	if(v[i].second < 20 ) return false;
    }
    return true;

}
class Get_Junction_Paths
{
	private:
	vector<int> Edges_left;// left node of an edge in junction graph
        vector<int> Edges_right;// right node of an edge in junction graph
        vector<double> Weights;
	vector<int> DAG_left;
	vector<int> DAG_right;
	vector<double> DAG_weights;
	vector<int> Unused_junctions;
	vector<vector<int> > Unused_pair_paths;
	vector<double> Pair_path_weights;
        vector<int> Cons_left;
        vector<int> Cons_right;
        vector<int> S,T;
        vector<int> Single_nodes;
        int max_node_num;
	double SEED_Filter;
	public:
	Get_Junction_Paths(vector<int> edges_left,vector<int> edges_right,vector<double> weights,vector<int> dag_left,vector<int> dag_right,vector<double> dag_weights,vector<int> cons_left,vector<int> cons_right,vector<int> s,vector<int> t,vector<int> single_nodes,int max_node,double seed_filter,vector<int> unused_junctions,vector<vector<int> > unused_pair_paths,vector<double> pair_path_weights);
	~Get_Junction_Paths();
	void Search_Junction_Paths(SimplifyGraph& simplify);
	vector<double> Weights_Pair;
	vector<vector<int> > Junction_Path_cover;


 	bool check_overlap_of_2paths(vector<int> , vector<int> );
	bool CheckOverlap(vector<int> );//YU
	void show()
	{
	    for(size_t i=0;i<Edges_left.size();i++)
	    {
		cout<<DAG_left[ Edges_left[i] ]<<"->"<<DAG_right[Edges_left[i]]<<" ==>> "
		    <<DAG_left[Edges_right[i] ]<<"->"<<DAG_right[Edges_right[i]]<<": "<<DAG_weights[Edges_right[i]]<<"  "<<Weights[i]<<endl;
	    }
	    cout<<"*********"<<endl;
	    for(size_t i=0;i<DAG_left.size();i++)
	    {
		cout<<DAG_left[i]<<"++>>"<<DAG_right[i]<<endl;
	    }
	}
};
Get_Junction_Paths::Get_Junction_Paths(vector<int> edges_left,vector<int> edges_right,vector<double> weights,vector<int> dag_left,vector<int> dag_right,vector<double> dag_weights,vector<int> cons_left,vector<int> cons_right,vector<int> s,vector<int> t,vector<int> single_nodes,int max_node,double seed_filter,vector<int> unused_junctions,vector<vector<int> > unused_pair_paths,vector<double> pair_path_weights)
{
	
	Edges_left=edges_left;
	Edges_right=edges_right;
	Weights=weights;
	DAG_left=dag_left;
	DAG_right=dag_right;
	DAG_weights=dag_weights;
	Unused_junctions=unused_junctions;
	Unused_pair_paths=unused_pair_paths;
	Pair_path_weights=pair_path_weights;
	
	Cons_left=cons_left;
	Cons_right=cons_right;
	S=s;
        T=t;
	Single_nodes=single_nodes;
	max_node_num=max_node;
	
	SEED_Filter=seed_filter;
	
	return;
}
Get_Junction_Paths::~Get_Junction_Paths()
{
	return;
}
void Get_Junction_Paths::Search_Junction_Paths(SimplifyGraph& simplify)
{
	//show();
	size_t i,j;
	int seed,seed_temp,SEED;
	double max_weights;
	vector<double> Weights_Pair;
	vector<int>::iterator pos;
	vector<int> Extension_left;
        vector<int> Extension_right;
	vector<int> Extension_num_0, Extension_num_1;
	vector<int> Extension_num_unused;
        vector<int> Extension;
	//vector<int> Edges_new_L, Edges_new_R;	
	int Raw_size=DAG_weights.size();
	vector<vector<int> > Raw_Unused_pair_paths=Unused_pair_paths;
	int edge_id=DAG_weights.size();
    	vector<int> Unused_pair_ids;
	for (int p=0;p<Unused_pair_paths.size();p++)
	{
		int seed_left=-1;
		int seed_right=-1;
		vector<int> Curr_path=Unused_pair_paths[p];
		for (i=0;i<DAG_left.size();i++)
		{
			if (DAG_left[i]==Curr_path[0] && DAG_right[i]==Curr_path[1])
			{
				seed_left=i;
			}
			if (DAG_left[i]==Curr_path[Curr_path.size()-2] && DAG_right[i]==Curr_path[Curr_path.size()-1])
			{
				seed_right=i;
			}
			if (seed_left>=0 && seed_right>=0)
			{
				break;
			}
		}
		Edges_left.push_back(seed_left);
		Edges_right.push_back(edge_id);
		Weights.push_back(0);//(-1); //YU 0
		Edges_left.push_back(edge_id);
                Edges_right.push_back(seed_right);
		Weights.push_back(0);//(-1); //YU 0
		DAG_left.push_back(-1);
		DAG_right.push_back(-1);
		DAG_weights.push_back(Pair_path_weights[p] - 1); //-1 YU
		Unused_junctions.push_back(edge_id); 
		Unused_pair_ids.push_back(edge_id);
		edge_id++;
	}
	Weights_Pair=Weights; 
	//YU
	for (i=0;i<Weights_Pair.size();i++) {
		break;
		Weights_Pair[i] =0;
	}
	map<int,vector<pair<edge_t,edge_t> > >node_packing_map;
	for (i=0;i<Weights_Pair.size();i++) {
	    //break;
	    edge_t e1 = make_pair(DAG_left[ Edges_left[i] ],DAG_right[Edges_left[i]]);
	    edge_t e2 = make_pair(DAG_left[Edges_right[i] ],DAG_right[Edges_right[i]]);
	    if(e1.first == -1) break;

	    int packing_node = e1.second;
	    if(simplify.node_set[packing_node].parents.size() == 2 && simplify.node_set[packing_node].children.size() == 2 && Weights_Pair[i] == 0)
	    {
		pair<edge_t,edge_t> p = make_pair(e1,e2);
	   	if(node_packing_map.find(packing_node) == node_packing_map.end())
		{
		    vector<pair<edge_t,edge_t> > v(1,p);
		    node_packing_map[packing_node] = v;
		}
		else node_packing_map[packing_node].push_back(p);
	    }
	}
	for (i=0;i<Weights_Pair.size();i++) 
	{
	    //Weights_Pair[i] = 0;
	    //break;
	    if(!PackingFlag) Weights_Pair[i] = 0;
	    else
	    {
		PackingFlag= false;
		break;//New-borrow
	    }
	    edge_t e1 = make_pair(DAG_left[ Edges_left[i] ],DAG_right[Edges_left[i]]);
	    edge_t e2 = make_pair(DAG_left[Edges_right[i] ],DAG_right[Edges_right[i]]);
	    if(Weights_Pair[i] == 0) continue;
	    if(e1.first == -1) break;
	    int packing_node = e1.second;
	    bool flag = false;
	    if(node_packing_map.find(packing_node)!=node_packing_map.end())
	    {
		vector<pair<edge_t,edge_t> > packing_result = node_packing_map[packing_node];
		if(packing_result.size() == 2)
		{
		    //cerr<<"check "<<packing_node<<endl;
		
		    pair<edge_t,edge_t> pack1 = packing_result[0];
		    pair<edge_t,edge_t> pack2 = packing_result[1];

		    edge_t e1 = pack1.first, e2 = pack1.second;
		    edge_t e3 = pack2.first, e4 = pack2.second;


		    double cov1 = simplify.node_set[e1.first].get_child_coverage(e1.second);
		    double cov2 = simplify.node_set[e2.first].get_child_coverage(e2.second);
		    double cov3 = simplify.node_set[e3.first].get_child_coverage(e3.second);
		    double cov4 = simplify.node_set[e4.first].get_child_coverage(e4.second);

		    //cerr<<e1.first<<"->"<<e1.second<<" "<<e2.first<<"->"<<e2.second<<" "<<cov1<<" "<<cov2<<endl;
                    //cerr<<e3.first<<"->"<<e3.second<<" "<<e4.first<<"->"<<e4.second<<" "<<cov3<<" "<<cov4<<endl;
		  

		    double min1 = cov1<cov2?cov1:cov2;
		    double max1 = cov1>cov2?cov1:cov2;
		    double min2 = cov3<cov4?cov3:cov4;
		    double max2 = cov3>cov4?cov3:cov4;


		    double r1 = min1/max1, r2 = min2/max2;

		    //same side
		    double min3 = cov1<cov3?cov1:cov3;
		    double max3 = cov1>cov3?cov1:cov3;
		    double min4 = cov2<cov4?cov2:cov4;
		    double max4 = cov2>cov4?cov2:cov4;

		    double r3 = min3/max3;
		    double r4 = min4/max4;

		    //cerr<<r1<<" "<<r2<<endl;
		    //if(r1 > 0.6 && r2 > 0.6 && r1 + r2 >1.4) flag = true;
		    //if(r1 > 0.5 && r2 > 0.5 && r3<0.5 && r4<0.5) flag = true;
		    //
		    if(r1 > 0.49 && r2 > 0.49) flag = true;//0.49
	
		}
	    }

	    if(flag && simplify.node_set[packing_node].parents.size() > 1 && simplify.node_set[packing_node].children.size() > 1 )
	    {
		//cerr<<packing_node<<" size: "
		    //<<simplify.node_set[packing_node].parents.size()<<" "<<simplify.node_set[packing_node].children.size()<<endl;
		//keep packing
		continue;
	    }
	    else   Weights_Pair[i] = 0;//YU NEW
	    
	    
	    Weights_Pair[i] = 0;
	}
	/*
	for (i=0;i<Cons_left.size();i++)
	{
	    for (j=0;j<Weights.size();j++)
	    {
		if (Cons_left[i]==Edges_left[j] && Cons_right[i]==Edges_right[j])
		{
			Weights_Pair[j]=-1;
			break;
		}
	    }
	}
	*/
	
   // Search paths from Unused_junctions


    while (Unused_junctions.size()>0)
    {
	if ( Unused_pair_ids.size()>0) 
	{
		max_weights=0.0;
		seed=Unused_pair_ids[0];
		for (i=0;i<Unused_pair_ids.size();i++)
		{
			if (DAG_weights[Unused_pair_ids[i]]>max_weights)
			{
				max_weights=DAG_weights[Unused_pair_ids[i]];
				seed=Unused_pair_ids[i];
			}
		}
	}
	else
	{
		max_weights=0.0;
		seed=Unused_junctions[0];
		Extension_left.clear();
		Extension_right.clear();
		Extension_num_0.clear();
	        Extension_num_1.clear();
		for (i=0;i<Unused_junctions.size();i++)
		{
			if (DAG_weights[Unused_junctions[i]]>max_weights)
			{
				max_weights=DAG_weights[Unused_junctions[i]];
				seed=Unused_junctions[i];
			}
		}
		if (DAG_weights[seed]<SEED_Filter)
		{
			break;
		}
	}
	Unused_junctions.erase(remove(Unused_junctions.begin(),Unused_junctions.end(),seed),Unused_junctions.end());
	SEED=seed;
//cerr<<"RAW_SIZE: "<<Raw_size<<endl;
//cerr<<"seed: "<<SEED<<endl;
//if(SEED>Raw_size){
//cerr<<"here"<<endl;
//vector<int> Curr_path = Raw_Unused_pair_paths[SEED-Raw_size];
//for(int j=0;j<Curr_path.size()-1;j++) cerr<<Curr_path[j]<<" ** ";
//cerr<<endl;
//}
	for (i=0;i<Weights_Pair.size();i++)
	{
		if (seed==Edges_right[i] && Weights_Pair[i]==0)
		{
			Extension_num_0.push_back(Edges_left[i]);
		}
		else if (seed==Edges_right[i] && Weights_Pair[i]==-1)
		{
			Extension_num_1.push_back(Edges_left[i]);
		}
	}
	while (Extension_num_0.size()>0 || Extension_num_1.size()>0) //YU
	{
		max_weights=0.0;
		if (Extension_num_1.size()>0)
	   	//if(0)//YU
		{
			seed_temp=Extension_num_1[0];
			for (i=0;i<Extension_num_1.size();i++)
			{
				if (DAG_weights[Extension_num_1[i]]>max_weights)
				{
					max_weights=DAG_weights[Extension_num_1[i]];
					seed_temp=Extension_num_1[i];
				}
			}
		}
		else
		{
			seed_temp=Extension_num_0[0];
			for (i=0;i<Extension_num_0.size();i++)
                         {
                                 if (DAG_weights[Extension_num_0[i]]>max_weights)
                                 {
                                         max_weights=DAG_weights[Extension_num_0[i]];
                                         seed_temp=Extension_num_0[i];
                                 }
                         }
		}
		Extension_left.push_back(seed_temp);
		Unused_junctions.erase(remove(Unused_junctions.begin(),Unused_junctions.end(),seed_temp),Unused_junctions.end());
		seed=seed_temp;
		Extension_num_0.clear();
		Extension_num_1.clear();
		for (i=0;i<Weights_Pair.size();i++)
		{
			if (seed==Edges_right[i] && Weights_Pair[i]==0)
			{
				Extension_num_0.push_back(Edges_left[i]);
			}
			else if (seed==Edges_right[i] && Weights_Pair[i]==-1)
			{
				Extension_num_1.push_back(Edges_left[i]);
			}
		}
	}//while (Extension_num_0.size()>0 || Extension_num_1.size()>0)
	seed=SEED;
	Extension_num_0.clear();
	Extension_num_1.clear();
	for (i=0;i<Weights_Pair.size();i++)
	{
		if (seed==Edges_left[i] && Weights_Pair[i]==0)
		{
			Extension_num_0.push_back(Edges_right[i]);
		}
		else if (seed==Edges_left[i] && Weights_Pair[i]==-1)
		{
			Extension_num_1.push_back(Edges_right[i]);
		}
	}
	while (Extension_num_0.size()>0 || Extension_num_1.size()>0) //YU
	{
		max_weights=0.0;
		if (Extension_num_1.size()>0)
		//if(0) //YU
		{
			seed_temp=Extension_num_1[0];
			for (i=0;i<Extension_num_1.size();i++)
			{
				if (DAG_weights[Extension_num_1[i]]>max_weights)
				{
					max_weights=DAG_weights[Extension_num_1[i]];
					seed_temp=Extension_num_1[i];
				}
			}
		}
		else
		{
			seed_temp=Extension_num_0[0];
			for (i=0;i<Extension_num_0.size();i++)
                         {
                                 if (DAG_weights[Extension_num_0[i]]>max_weights)
                                 {
                                         max_weights=DAG_weights[Extension_num_0[i]];
                                         seed_temp=Extension_num_0[i];
                                 }
                         }
		}
		Extension_right.push_back(seed_temp);
		Unused_junctions.erase(remove(Unused_junctions.begin(),Unused_junctions.end(),seed_temp),Unused_junctions.end());
		seed=seed_temp;
		Extension_num_0.clear();
		Extension_num_1.clear();
		for (i=0;i<Weights_Pair.size();i++)
		{
			if (seed==Edges_left[i] && Weights_Pair[i]==0)
			{
				Extension_num_0.push_back(Edges_right[i]);
			}
			else if (seed==Edges_left[i] && Weights_Pair[i]==-1)
			{
				Extension_num_1.push_back(Edges_right[i]);
			}
		}
	}//while (Extension_num_0.size()>0 || Extension_num_1.size()>0)
	for (i=0;i<Extension_left.size();i++)
	{
		Extension.push_back(Extension_left[Extension_left.size()-1-i]);
	}
	Extension.push_back(SEED);
	for (i=0;i<Extension_right.size();i++)
	{
		Extension.push_back(Extension_right[i]);
	}
	//for(i=0;i<Extension.size();i++) cerr<<Extension[i]<<"-> ";
	//cerr<<endl;
	vector<int> Extension_to_DAG;
	vector<int> Extension_new;
	vector<int> vec_rm;
	vector<vector<int> > Unused_pair_paths_update;
	vector<int> Unused_pair_ids_update;
	vector<double> Pair_path_weights_update;
	for (i=0;i<Extension.size();i++)
	{
		if (Extension[i]<Raw_size)
		{
			if (Extension_new.size() == 0 || Extension_new[Extension_new.size()-1] != Extension[i])
			{
				Extension_new.push_back(Extension[i]);
			}
		}
		else
		{
			vector<int> Curr_path = Raw_Unused_pair_paths[Extension[i]-Raw_size];
			for(int j=0;j<Curr_path.size()-1;j++)
			{
				for (int k=0;k<DAG_left.size();k++)
				{
					if (DAG_left[k]==Curr_path[j] && DAG_right[k]==Curr_path[j+1])
					{
						if (Extension_new.size() == 0 || Extension_new[Extension_new.size()-1] != k)
						{
							Extension_new.push_back(k);
						}
						break;
					}
				}
			}
		}
	}
	Extension=Extension_new;
//	cerr<<"new: ";
//	for(i=0;i<Extension.size();i++) cerr<<Extension[i]<<"-> ";
//	cerr<<endl;
	//bool flag = CheckOverlap(Extension);//YU
	//if(!flag)
	Junction_Path_cover.push_back(Extension);
	//if(Junction_Path_cover.size() > 5) SEED_Filter = 5;
/*
	for (i=0;i<Unused_junctions.size();i++)
	{
		cerr<<"A: "<<Unused_junctions[i]<<endl;
	}
*/
	for (i=0;i<Extension.size();i++)
	{
		//cerr<<"B: "<<Extension[i]<<endl;
		Unused_junctions.erase(remove(Unused_junctions.begin(),Unused_junctions.end(),Extension[i]),Unused_junctions.end());
	}
	Extension_to_DAG.push_back(DAG_left[Extension[0]]);
	for (i=0;i<Extension.size();i++)
	{
		Extension_to_DAG.push_back(DAG_right[Extension[i]]);
	}
//cerr<<"**************"<<endl;
//for(i = 0;i<Extension_to_DAG.size();i++)
//cerr<<Extension_to_DAG[i]<<" =>> ";
//cerr<<endl;
	for (i=0;i<Unused_pair_paths.size();i++)
	{
		vector<int> vec_j,vec_k;
		int j_size=Unused_pair_paths[i].size();
		int k_size=Extension_to_DAG.size();
		if (j_size<=k_size)
		{
			vec_j=Unused_pair_paths[i];
			//for(int k = 0;k<vec_j.size();k++) cerr<<vec_j[k]<<" * ";
			//cerr<<endl;
			for (int t=0;t<k_size-j_size+1;t++)
			{
				for (int m=t;m<t+j_size;m++)
				{
					vec_k.push_back(Extension_to_DAG[m]);
				}
				//for(int k = 0;k<vec_k.size();k++) cerr<<vec_k[k]<<" - ";
				//cerr<<endl;
				if (vec_k == vec_j)
				{
					vec_rm.push_back(i);
					break;
				}
				vec_k.clear();
			}
		}
	}
	for (i=0;i<Unused_pair_paths.size();i++)
	{
		pos=find(vec_rm.begin(),vec_rm.end(),i);
		if (pos == vec_rm.end())
		{
			Unused_pair_paths_update.push_back(Unused_pair_paths[i]);
			Unused_pair_ids_update.push_back(Unused_pair_ids[i]);
			Pair_path_weights_update.push_back(Pair_path_weights[i]);
		}
		else
		{
			Unused_junctions.erase(remove(Unused_junctions.begin(),Unused_junctions.end(),Unused_pair_ids[i]),Unused_junctions.end());
		}
	}
	Unused_pair_paths=Unused_pair_paths_update;
	Unused_pair_ids=Unused_pair_ids_update;
	Pair_path_weights=Pair_path_weights_update;
//cerr<<Unused_pair_paths.size()<<" "<<Unused_pair_ids.size()<<" "<<Pair_path_weights.size()<<" "<<endl;
	Extension.clear();
	Extension_left.clear();
	Extension_right.clear();
    }//end of while (Unused_junctions.size()>0)

}
bool Get_Junction_Paths::check_overlap_of_2paths(path_t p1_, path_t p2_)
{
    path_t p1,p2;
    if(p1_.size() <= p2_.size()) { p1 = p1_; p2= p2_;}
    else {p1 = p2_; p2 = p1_;}

    if(p1.size() <= 3) return false;
/*
    cout<<endl;
    for(size_t i=0;i<p2.size();i++) cout<<p2[i]<<"->";
    cout<<endl;
    for(size_t i=0;i<p1.size();i++) cout<<p1[i]<<"->";
    cout<<endl;
*/

    stringstream ss,ss2;
    string s_p2;
    for(size_t i=0;i<p2.size();i++) ss<<p2[i]<<"_";
    s_p2 = ss.str();
    //cout<<s_p2<<endl;
    vector<path_t> subpath;
    int length = int (0.9 * p1.size()) + 1;
    for(int i=0;i<p1.size() - length;i++)
    {
	path_t sp;
	stringstream sstr;
	string sp_str;
	for(int j=i;j<i+length;j++) {
	    sp.push_back(p1[j]);
	    sstr<<p1[j]<<"_";
	    sp_str = sstr.str();
	}
	//cout<<" "<<sp_str<<endl;
	if(s_p2.find(sp_str) != string::npos) return true;
	//note: sp.size() = length;
	
	for(int j = i+length; j<p1.size();j++){
	    sp.push_back(p1[j]);
	    sstr<<p1[j]<<"_";
	    sp_str = sstr.str();
	    //cout<<" "<<sp_str<<endl;
	    if(s_p2.find(sp_str) != string::npos) return true;
	}
	
    }
    return false;   
}
//if(overlap) return true;
bool Get_Junction_Paths::CheckOverlap(path_t path)
{
    if(Junction_Path_cover.empty()) return false;
    int overlap_number = 0;
    for(size_t i=0;i<Junction_Path_cover.size();i++)
    {
	if(check_overlap_of_2paths(path,Junction_Path_cover[i]))
	    return true;
    }

    return false;
    
}

#endif







