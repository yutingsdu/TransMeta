#include "junction_graph.h"
#include "QuadProg++.h" 
using namespace std;
extern int rg_index;
void enumerate_subset(int* Set, int Set_length, int seq, vector<int>& subset)
{
    for(int i=0;i<Set_length;i++)
    {
        if( seq & (1<<i) )
            subset.push_back(Set[i]);
    }
}
void Combing::get_triplet_map(vector<path_t>& PairPath, vector<double>& PairPath_cov)
{
    for(size_t i=0;i<PairPath.size();i++)
    {
	path_t pp = PairPath[i];
	double cov = PairPath_cov[i];
	if(pp.size() < 3) continue;
	for(size_t j=0;j<pp.size()-2;j++)
	{
	    triplet_t trip;
	    trip.push_back(pp[j]);trip.push_back(pp[j+1]);trip.push_back(pp[j+2]);
	    if(triplet_map.find(trip) == triplet_map.end()) triplet_map[trip] = cov;
	    else triplet_map[trip] += cov;
	
	    if(node_triplet_map.find(trip[1]) == node_triplet_map.end())
	    {
		vector<triplet_t> v(1,trip);
		node_triplet_map[trip[1]] = v;
	    }
	    else node_triplet_map[trip[1]].push_back(trip);
	}
    }
    map<node_idx_t,vector<triplet_t> > ::iterator it = node_triplet_map.begin();
    for(;it!=node_triplet_map.end();it++)
    {
	sort(it->second.begin(),it->second.end());
	it->second.erase(unique(it->second.begin(),it->second.end()),it->second.end());
	/*
	cout<<it->first<<": "<<endl;
	vector<triplet_t> vT = it->second;
	for(int j=0;j<vT.size();j++) cout<<" "<<vT[j][0]<<" "<<vT[j][1]<<" "<<vT[j][2]<<endl;
	*/
	
    }

    return;
}
void Combing::packing()
{
    /*
    int Length = 5;
    int Set[Length];
    for(int i=0;i<Length;i++) Set[i] = i;
    for(int i=0;i<(1<<Length);i++){
        vector<int> subset;
        enumerate_subset(Set,Length,i,subset);
        //if(subset.size() == parents.size())  All_subset.push_back(subset);
        //for(size_t j=0;j<subset.size();j++) cout<<subset[j]<<" ";
	//cout<<endl;
    }
    */
    //cout<<"Graph: "<<rg_index<<endl;
    //show_graph();
    for(int i=0;i<size_;i++)
    {
	Node n = node_set[i];
	if(n.children.empty() || n.parents.empty()) continue;

	if(n.children.size() == 1 || n.parents.size() == 1)
	{
	    trivial_packing(i,n.children,n.parents);
	    continue;
	}
	packing(i,n.children,n.parents);
    }
 //   if(rg_index>30000) exit(0);
    return ;
}
bool Combing::check_triplet(vector<triplet_t> current_triplets, vector<triplet_t> pair_triplets,
				 vector< pair<node_idx_t,double> > children,vector< pair<node_idx_t,double> > parents)
{
    /*//obviousely
    for(size_t i=0;i<pair_triplets.size();i++)
    {
	if(find(current_triplets.begin(),current_triplets.end(),pair_triplets[i]) == current_triplets.end()) return false;
    }
    */
    vector<node_idx_t> vP,vC;
    for(size_t i=0;i<current_triplets.size();i++)
    {
	vP.push_back(current_triplets[i][0]);
	vC.push_back(current_triplets[i][2]);
    }
    for(size_t i=0;i<children.size();i++) 
	if(find(vC.begin(),vC.end(),children[i].first) == vC.end()) return false;
    for(size_t i=0;i<parents.size();i++)
	if(find(vP.begin(),vP.end(),parents[i].first) == vP.end()) return false;
/*
    cout<<"current: "<<endl;
    show_triplets(current_triplets);
    cout<<"Pair: "<<endl;
    show_triplets(pair_triplets);
    cout<<"*"<<endl;
*/
    return true;
   
}
/*
 *x-0-0 --> w0
 *x-1-0 --> w1
 *x-1-1 --> w2
 *x-2-0 --> w3
 *x-2-1 --> w4
 *x-3-2 --> w5
 * s0-w0(1,0,0,0,0,0) / s1-w1,w2(0,1,1,0,0,0) / s2-w3,w4(0,0,0,1,1,0) / s3-w5(0,0,0,0,0,1)
 * c0-w0,w1,w3(1,1,0,1,0,0) / c1-w2,w4(0,0,1,0,1,0) / c2-w5(0,0,0,0,0,1)
 * */
double Combing::solve_qp(vector<triplet_t>& current_triplets,
			 vector< pair<node_idx_t,double> > children,vector< pair<node_idx_t,double> >parents)
{
    //cout<<"****************** solve_qp ********************"<<endl;
    //cout<<"triplets: "<<endl;
    //show_triplets(current_triplets);
    int unknown_num = current_triplets.size();
    double w[unknown_num];
    vector<int> v_temp(unknown_num,0);

    //parents obj
    vector< vector<int> > SI(parents.size(),v_temp);
    for(size_t i=0;i<parents.size();i++)
    {
	for(size_t j=0;j<current_triplets.size();j++)
	{
	    if(parents[i].first == current_triplets[j].front())
		SI[i][j] = 1;
	}
    }
    
    //children obj
    vector< vector<int> > CI(children.size(),v_temp);
    for(size_t i=0;i<children.size();i++)
    {
	for(size_t j=0;j<current_triplets.size();j++)
	{
	    if(children[i].first == current_triplets[j].back())
		CI[i][j] = 1;
	}
    }
    double g[MATRIX_DIM][MATRIX_DIM]={0}, g0[MATRIX_DIM]={0}, ci[MATRIX_DIM][MATRIX_DIM]={0}, ci0[MATRIX_DIM]={0}, x[MATRIX_DIM];
    double ce[MATRIX_DIM][MATRIX_DIM]={0}, ce0[MATRIX_DIM]={0};
    // cout<<"MATRIX_DIM: "<<MATRIX_DIM<<endl;
    //print_matrix("G",g,unknown_num);
    for(size_t i_=0;i_<SI.size();i_++)//for each in_edge;
    {
	vector<int> v = SI[i_];	
	double cov = parents[i_].second;
	for(size_t i=0;i<v.size();i++)
	{
	    for(size_t j=0;j<v.size();j++)
	    {
		g[i][j] += 2*v[i]*v[j];
		//if(i == j)
		    //ci[i][j] += 0-2*cov*v[i]*v[j];
	    }
	   g0[i] += 0-2*cov*v[i];
	}
    }
    for(size_t i_=0;i_<CI.size();i_++)
    {
	vector<int> v = CI[i_];
	double cov = children[i_].second;
	for(size_t i=0;i<v.size();i++)
	{
	    for(size_t j=0;j<v.size();j++)
	    {
		g[i][j] += 2*v[i]*v[j];
		//if(i == j)
		    //ci[i][j] += 0-2*cov*v[i]*v[j];
	    }
	    g0[i] += 0-2*cov*v[i];
	}
    }
    for(int i=0;i<unknown_num;i++) g[i][i] += 0.0000001;
    for(int i=0;i<unknown_num;i++) ci[i][i] = 1;
    //print_matrix("G",g,unknown_num);
    //print_vector("g0",g0,unknown_num);
    double g_[MATRIX_DIM][MATRIX_DIM] ;
    for(int i=0;i<unknown_num;i++)
    {
	for(int j=0;j<unknown_num;j++) g_[i][j] = g[i][j];
    }
    double obj = solve_quadprog(g, g0, unknown_num, 
		   		ce, ce0, 0, 
		   		ci, ci0, unknown_num, x);
    //print_matrix("G",g,unknown_num);
   
    /*
    for(size_t i=0;i<SI.size();i++)
    {
	cout<<parents[i].first<<": ";
	vector<int> v = SI[i];
	for(size_t j=0;j<v.size();j++) cout<<v[j]<<" ";
	cout<<endl;
    }
    for(size_t i=0;i<CI.size();i++)
    {
	cout<<children[i].first<<": ";
	vector<int> v= CI[i];
	for(size_t j=0;j<v.size();j++) cout<<v[j]<<" ";
	cout<<endl;
    }
    */

    /*compute obj
    double wg[unknown_num];
    for(int i=0;i<unknown_num;i++)
    {
	wg[i] = 0;
	for(int j=0;j<unknown_num;j++)
	{
	   wg[i] += x[j]*g_[j][i];
	}
    }
    double wgw = 0;
    for(int i=0;i<unknown_num;i++) wgw += wg[i]*x[i];
    double g0w = 0;
    for(int i=0;i<unknown_num;i++) g0w += g0[i]*x[i];
    
    cout<<"obj_yu: "<<wgw<<" "<<g0w<<" "<<wgw/2 + g0w<<endl;
    */
    
    //for(int i=0;i<unknown_num;i++) cout<<x[i]<<" ";
    //cout<<endl<<"obj: "<<obj<<endl;
    //cout<<"*** Finish solve ***"<<endl<<endl;
    return obj;
}

void Combing::packing(node_idx_t n, vector< pair<node_idx_t,double> > children, vector< pair<node_idx_t,double> > parents)
{
    vector<vector<int> > All_subsets;
    vector<triplet_t> pair_triplets = node_triplet_map[n];
    vector<triplet_t> no_pair_triplets;
    //cout<<"packing node: "<<n<<"("<<parents.size()<<","<<children.size()<<")"<<endl;
    //cout<<"pair-triplet: "<<endl;
    //show_triplets(pair_triplets);
    for(size_t i=0;i<parents.size();i++)
    {
	node_idx_t p = parents[i].first;
	for(size_t j=0;j<children.size();j++)
	{
	    node_idx_t c = children[j].first;
	    triplet_t trip;
	    trip.push_back(p);trip.push_back(n);trip.push_back(c);
	    if(find(pair_triplets.begin(),pair_triplets.end(),trip) != pair_triplets.end()) continue;
	    no_pair_triplets.push_back(trip);
	}
    }
    if( no_pair_triplets.size()>2) return trivial_packing(n,children,parents);
    int max_size = children.size()>parents.size()?children.size():parents.size();

    int Length = no_pair_triplets.size();;
    int Set[Length];
    for(int i=0;i<Length;i++) Set[i] = i;
    for(int i=0;i<(1<<Length);i++){
	vector<int> subset;
	enumerate_subset(Set,Length,i,subset);
	if(subset.size() + pair_triplets.size() >= max_size)  All_subsets.push_back(subset);
    }
    int M = max_size>pair_triplets.size()?max_size:pair_triplets.size();
    int M0 = no_pair_triplets.size() + pair_triplets.size();
    //cout<<"M: "<<M<<endl;
    //M+=1;
    double Obj_for_current_M = 10000000;
    vector<triplet_t> final_triplets;
    while(1)
    {
	double Obj_for_current_subset = 10000000;
	vector<triplet_t> current_final_triplets;
	bool flag = false;
	for(size_t i=0;i<All_subsets.size();i++)
	{
	    vector<int> subset = All_subsets[i];
	    if(subset.size() + pair_triplets.size() != M ) continue;

	    vector<triplet_t> current_triplets;
	    for(size_t j=0;j<subset.size();j++)
	    {
		int index = subset[j];
		current_triplets.push_back(no_pair_triplets[index]);
	    }
	    for(size_t j=0;j<pair_triplets.size();j++)
		current_triplets.push_back(pair_triplets[j]);

	    flag = check_triplet(current_triplets,pair_triplets,children,parents);
	    if(flag)
	    {
		double obj_ = solve_qp(current_triplets,children,parents);
		if(obj_<Obj_for_current_subset ){
		    Obj_for_current_subset = obj_;	
		    current_final_triplets = current_triplets;
		}
	    }
	}
	if(Obj_for_current_subset != 10000000){
	    if(Obj_for_current_M == 10000000)
	    {
		final_triplets = current_final_triplets;
		Obj_for_current_M = Obj_for_current_subset;
		break;//NEW4.3
	    }
	    else if( Obj_for_current_subset<Obj_for_current_M )
	    {
		final_triplets = current_final_triplets;
		Obj_for_current_M = Obj_for_current_subset;
	    }
	    else break;
	}
	//if(!flag) M++;
	M++;
	if(M+1 > M0) break;
    }
    if(final_triplets.empty()) return trivial_packing(n,children,parents);
    for(int i=0;i<final_triplets.size();i++)
    {
	edge_t edge1 = make_pair(final_triplets[i][0],final_triplets[i][1]);
	edge_t edge2 = make_pair(final_triplets[i][1],final_triplets[i][2]);
	packing_map[make_pair(edge1,edge2)] = 0;
    }

}
void Combing::trivial_packing(node_idx_t n, vector< pair<node_idx_t,double> > children, vector< pair<node_idx_t,double> > parents)
{
    for(size_t i=0;i<children.size();i++)
    {
	edge_t edge1 = make_pair(n,children[i].first);
	for(size_t j=0;j<parents.size();j++)
	{
	    edge_t edge2 = make_pair(parents[j].first,n);
	    //cout<<edge2.first<<"->"<<edge2.second<<" --> "<<edge1.first<<"->"<<edge1.second<<endl;
	    packing_map[make_pair(edge2,edge1)] = 0;
	}
    }
}
void Combing::get_packing_result(vector<vector<int> >Vec_edges)
{
    vector<edge_t> Edges;
    for(size_t i=0;i<Vec_edges.size();i++)
    {
	edge_t e = make_pair(Vec_edges[i][0],Vec_edges[i][1]);
	Edges.push_back(e);
    }
    //map<pair<edge_t,edge_t>,int> packing_map;
    map<pair<edge_t,edge_t>,int>::iterator it = packing_map.begin();
    for(;it != packing_map.end();it++)
    {
	edge_t e1 = it->first.first,e2 = it->first.second;
	int w = it->second;
	vector<edge_t>::iterator i1 = find(Edges.begin(),Edges.end(),e1);
	vector<edge_t>::iterator i2 = find(Edges.begin(),Edges.end(),e2);
	int i1_= i1 - Edges.begin(),i2_ = i2 - Edges.begin();

	Edges_left.push_back(i1_); 
	Edges_right.push_back(i2_);
	Weights.push_back(w);
		
    }
    return;
}
void Combing::show_triplets(vector<triplet_t> vT)
{
    for(size_t i=0;i<vT.size();i++)
    {
	cout<<"w"<<i<<" - ";
	show_triplet(vT[i]);
    }
    return;
}
void Combing::show_triplet(triplet_t T)
{
    cout<<T[0]<<" "<<T[1]<<" "<<T[2]<<endl;
    return;
}
void Combing::show_graph()
{
    for(int i=0;i<size_;i++)
    {
        Node n = node_set[i];
	for(size_t j=0;j<n.children.size();j++)
	{
	    cout<<i<<"->"<<n.children[j].first<<": "<<n.children[j].second<<endl;
	}
    }   
    return;
}
