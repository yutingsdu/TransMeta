//#include"expression_level.h"
#include"Trim.h"
#include"Function.h"

//#include"get_junction_graph_new.h"
//#include"junction_paths_new.h"
//#include"recover_paths.h"
//#include"PairPacker.h"
#include "pair_path.h"
#include "PairPath_graph.h"
#include "path_search.h" 
#include <algorithm>
#include <boost/unordered_map.hpp>
using namespace std;
extern string suffix;
extern ofstream out_graph;
extern int Path_length;

extern int trim;

extern float SEED_Filter;

extern string out_name;
extern string out_dir;

extern int rg_index;
extern int trans_id;

typedef vector<int> path_t;
string base_name()
{
   stringstream idx ;
   idx << "comp" << rg_index ;
   rg_index++;
   trans_id=1;
   return idx.str();
}
/*
bool Psorter(path_t p1,path_t p2)
{
    return ((p1.size()>p2.size()) || 
	    (p1.size() == p2.size() && p1 < p2));
}
*/
void trim_node(int gene_l,vector<int>& exon_l,vector<int>& exon_r,vector<double>& v_exon_cov,vector<double>& v_Gene,
	      vector<int>& junc_l,vector<int>& junc_r)
{
    
    if(1 || trim == 1){ //YU-Borrow
	vector<int> exon_l_p,exon_r_p;
        exon_l_p=exon_l;exon_r_p=exon_r;
	if(exon_l_p.size()>1)
	{
	    for(int i=0;i<exon_l_p.size();i++)
	    {
	        if(!is_in(exon_l_p[i],junc_r))
		{
	            int l=exon_r_p[i]-exon_l_p[i]+1;
	            int site=exon_l_p[i]-gene_l;
	            vector<int> seg_node=find_trims(gene_l,exon_l_p[i],exon_r_p[i],v_Gene);
	            if(seg_node.size()>0 && seg_node[0]!=0){
	                 exon_l_p[i]=seg_node[0];
	            }
	            double cov=0.00000;
	            int l_new=exon_r_p[i]-exon_l_p[i]+1;
	            for(int j=exon_l_p[i]-gene_l;j<=l+site-1;j++) cov=cov+v_Gene[j];
	             v_exon_cov[i]=cov/l_new;
	        }
	    }
	    for(int i=0;i<exon_l_p.size();i++)
	    {
	        if(!is_in(exon_r_p[i]+1,junc_l))
		{
	            int l=exon_r_p[i]-exon_l_p[i]+1;
	            int site=exon_l_p[i]-gene_l;
	            vector<int> seg_node=find_trims(gene_l,exon_l_p[i],exon_r_p[i],v_Gene);
	            if(seg_node.size()>0 && seg_node[1]!=0){
	                exon_r_p[i]=seg_node[1];
	            }
	            double cov=0.00000;
	            int l_new=exon_r_p[i]-exon_l_p[i]+1;
	            for(int j=site;j<=exon_r_p[i]-gene_l;j++) cov=cov+v_Gene[j];
	            v_exon_cov[i]=cov/l_new;
	        }
	    }
	    exon_l=exon_l_p;exon_r=exon_r_p;

	} //if(exon_l_p.size>1);
     }
	return;
	
}
void get_exon_sequence( int gene_l, vector<int>& exon_l,vector<int>& exon_r,
			vector<char>& Gene_sequence, vector<string>& exon_sequence)
{
    for(int i = 0;i<exon_l.size();i++)
    {
	string s = "";
	for(int j = exon_l[i] - gene_l; j<=exon_r[i] - gene_l;j++) s.append(1,Gene_sequence[j]);
	exon_sequence.push_back(s);
    }

    return;
}
vector< vector<int> >get_pair_path(vector< vector<int> > final_pair_path_exon,
				   vector<int>& exon_l,vector<int>& exon_r,
				   vector<int>& junc_l,vector<int>& junc_r)
{
    vector<pair<int,int> > exon,junc;
    for(int i=0;i<exon_l.size();i++)
    {
	pair<int,int> p = make_pair(exon_l[i],exon_r[i]);
	exon.push_back(p);
    }
    for(int i = 0;i<junc_l.size();i++)
    {
	pair<int,int> p = make_pair(junc_l[i],junc_r[i]);
	junc.push_back(p);
    }
    vector<vector<int> > pair_path;
    for(int i = 0;i<final_pair_path_exon.size();i++)
    {
	vector<int> current_path_exon = final_pair_path_exon[i];
	bool nojunc_flag = false;
	for(int j = 1;j<current_path_exon.size()-1;)
	{
	    pair<int,int> junc_ = make_pair(current_path_exon[j] + 1,current_path_exon[j+1]);
	    if(find(junc.begin(),junc.end(),junc_) == junc.end())
	    {
		nojunc_flag = true;
		break;
	    }
	    j = j+2;
	}
	if(nojunc_flag) continue;

	vector<int> current_path;
	for(int j = 0;j<current_path_exon.size()-1;)
	{
	    int el = current_path_exon[j], er = current_path_exon[j+1];
	    pair<int,int> p = make_pair(el,er);

	    vector<pair<int,int> >::iterator it = find(exon.begin(),exon.end(),p);

	    if(it == exon.end()) continue;

	    int k = it - exon.begin();
	    current_path.push_back(k);
	    j = j+2;
	}
	pair_path.push_back(current_path);	
    }
    return pair_path;
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
void describe_graph(int gene_l,vector<int>& junc_l,vector<int>& junc_r,vector<double>& junc_cov,
		      vector<int>& exon_l,vector<int>& exon_r,vector<double>& v_exon_cov,
		      string ch,vector<double>& v_Gene,
		      int strand,int line,int XS_plus,int XS_minus,
		      boost::unordered_map<string, pair<cigar_t,cigar_t> >& tu_readid_pairInfo_map,
		      vector< vector<int> > final_pair_path_exon
		      )
{
    int count_partial=0;
    for(int i=1;i<exon_l.size();i++){
        if(exon_l[i]-exon_r[i-1]==1) count_partial++;
    }
    vector<vector<int> > Vec_edges;
    vector<double> Vec_weights;

    vector<vector<int> > Unused_pair_paths;

    vector<int> Unused_junctions;
 
    //PairPathGraph
    Unused_pair_paths = get_pair_path(final_pair_path_exon,exon_l,exon_r,junc_l,junc_r);
     
    PairPath_Graph ppath_graph(Unused_pair_paths);
    
    ppath_graph.build_pairpath_graph(junc_l,junc_r,junc_cov,exon_l,exon_r);

    ppath_graph.search_path();

    for(int i = 0;i<ppath_graph.final_pair_path.size();i++)
	Unused_pair_paths.push_back(ppath_graph.final_pair_path[i]);
    
    vector<double> Unused_pair_paths_cov;
    for(int i = 0; i<Unused_pair_paths.size(); )
    {
	double cov = get_minimum_cov(Unused_pair_paths[i], junc_l, junc_r, junc_cov, exon_l, exon_r);
	if(cov < 3) Unused_pair_paths.erase(Unused_pair_paths.begin() + i);
        else
        {
	    Unused_pair_paths_cov.push_back(cov);
	    i++;
        }
    }

    // ** trim	**
    trim_node(gene_l,exon_l, exon_r, v_exon_cov, v_Gene,junc_l,junc_r);

    
    vector<string> vec_edges,vec_pairs,Node_seq_1;
    vector<int>Node_seq_2,Node_seq_3;
    vector<double>Node_seq_4;

    // ** get edge **
    for(int i=0;i<junc_l.size();i++)
    {
	string edge;
	vector<int> edge_;
	stringstream ss;
	for(int j=0;j<exon_l.size();j++)
	{
	    if(junc_l[i]-1==exon_r[j]) {
		//ss<<j<<"->"; 
		edge_.push_back(j);
	    }
	    if(junc_r[i]==exon_l[j]) {
		//ss<<j<<":"<<junc_cov[i]<<";"; 

		edge_.push_back(j); 
		Vec_edges.push_back(edge_);
		Vec_weights.push_back(junc_cov[i]);		
	    }
	}
	//ss>>edge;
	vec_edges.push_back(edge);
    }
    vector<int> edge_out,edge_in;
    vector<double> edge_weight;
    for(int i=0;i<Vec_edges.size();i++)
    {
	edge_out.push_back(Vec_edges[i][0]);
	edge_in.push_back(Vec_edges[i][1]);
	edge_weight.push_back(Vec_weights[i]);
    }

    // ** get node infomation **
    for(int i=0;i<exon_l.size();i++)
    {
	Node_seq_1.push_back(ch);
	Node_seq_2.push_back(exon_l[i]);
	Node_seq_3.push_back(exon_r[i]);
	Node_seq_4.push_back(v_exon_cov[i]);
    }

    string strd = "";
    if(line == 2 && strand == 2) strd = "+";
    if(line == 2 && strand == 1) strd = "-";
    if(line == 1){
        if(XS_plus > XS_minus) strd = "+";
	if(XS_plus < XS_minus) strd = "-";
	if(XS_plus==0 && XS_minus==0) strd = ".";
/*
	if(XS_plus == XS_minus) //Borrow
	{
	    if(edge_out.size() > 0) strd = "+";
	    else strd = ".";
	}
*/
    }
    
    //string graph_name ="./" + out_dir + "/Graph";
    string graph_name;
  if(0) 
  {
    //graph_name = out_graphdir + "/MyGraph.frame";
    //char* graph_name_ =  const_cast<char*>(graph_name.c_str());
    //ofstream out_graph(graph_name_,ios::app);
    out_graph<<"Edges"<<endl;
    for(int i=0;i<Vec_edges.size();i++) out_graph<<Vec_edges[i][0]<<" -> "<<Vec_edges[i][1]<<" : "<<Vec_weights[i]<<endl;
    out_graph<<"Nodes"<<endl;
    for(int i=0;i<Node_seq_1.size();i++)out_graph<<strd<<" "<<Node_seq_1[i]<<" "<<Node_seq_2[i]<<" "<<Node_seq_3[i]<<" "<<Node_seq_4[i]<<endl;
    out_graph<<"Pair"<<endl;

    for(size_t i=0;i<Unused_pair_paths.size();i++) {
        vector<int> pp = Unused_pair_paths[i];
        for(size_t j=0;j<pp.size();j++) out_graph<<pp[j]<<" ";
        out_graph<<": "<<Unused_pair_paths_cov[i]<<endl;
    }
    out_graph<<"left-boundary: "<<gene_l<<endl;
    out_graph<<"vec-gene: "<<endl;
    for(size_t j=0;j<v_Gene.size();j++)
    {
       if(v_Gene[j] == 0 ) continue;
       out_graph<<gene_l+j<<" "<<v_Gene[j]<<" ; ";
    }
    out_graph<<'\n';
    /*
    size_t j=0;
    for(;j<v_Gene.size() - 100;) 
    {
	for(size_t k=j;k<j+100;k++)
	    out_graph<<v_Gene[k]<<" ";
	out_graph<<'\n';
	j += 100;
    }
    for(size_t k=j;k<v_Gene.size();k++) out_graph<<v_Gene[k]<<" ";
    out_graph<<'\n';
    */
    out_graph<<"#Graph "<<rg_index-1<<endl;
  }
  else if(0)
  {
     graph_name = out_dir + "/MyGraph.raw."+suffix +".g" ;
     char* graph_name_ =  const_cast<char*>(graph_name.c_str());
     ofstream Out_(graph_name_,ios::app);
     Out_<<"Edges"<<endl;
     for(int i=0;i<Vec_edges.size();i++) Out_<<Vec_edges[i][0]<<" -> "<<Vec_edges[i][1]<<" : "<<Vec_weights[i]<<endl;
     Out_<<"Nodes"<<endl;
     for(int i=0;i<Node_seq_1.size();i++)Out_<<strd<<" "<<Node_seq_1[i]<<" "<<Node_seq_2[i]<<" "<<Node_seq_3[i]<<" "<<Node_seq_4[i]<<endl;
    Out_<<"Pair"<<endl;

     for(size_t i=0;i<Unused_pair_paths.size();i++) {
         vector<int> pp = Unused_pair_paths[i];
         for(size_t j=0;j<pp.size();j++) Out_<<pp[j]<<" ";
         Out_<<": "<<Unused_pair_paths_cov[i]<<endl;
     }
     Out_<<"#Graph "<<rg_index-1<<endl;
   }
    
/*
    out_<<"Gene_"<<rg_index<<endl;
    out_<<"** Chr **"<<endl;
    out_<<ch+strd<<endl;
    out_<<"** Nodes **"<<endl;
    for(int i=0;i<Node_seq_1.size();i++)out_<<Node_seq_2[i]<<" "<<Node_seq_3[i]<<" "<<1<<" "<<2<<" "<<Node_seq_4[i]<<endl;
    out_<<"** Edges **"<<endl;
    for(int i=0;i<Vec_edges.size();i++) out_<<Vec_edges[i][0]<<"  "<<Vec_edges[i][1]<<"  "<<Vec_weights[i]<<endl;

    sort(Unused_pair_paths.begin(),Unused_pair_paths.end(),Psorter);
    out_<<"** Combine gene and pair paths **"<<endl;
    for(size_t i=0;i<Unused_pair_paths.size();i++) {
        vector<int> pp = Unused_pair_paths[i];
        for(size_t j=0;j<pp.size();j++) out_<<pp[j]<<" ";
	out_<<endl;
    }
    out_<<"** End **"<<endl;
    //return;
    //
    //
    */

/*
    out_<<"Edges"<<endl;
    for(int i=0;i<Vec_edges.size();i++) out_<<Vec_edges[i][0]<<" -> "<<Vec_edges[i][1]<<" : "<<Vec_weights[i]<<endl;
    out_<<"Nodes"<<endl;
    for(int i=0;i<Node_seq_1.size();i++)out_<<strd<<" "<<Node_seq_1[i]<<" "<<Node_seq_2[i]<<" "<<Node_seq_3[i]<<" "<<Node_seq_4[i]<<endl;
    out_<<"Pair"<<endl;

    for(size_t i=0;i<Unused_pair_paths.size();i++) {
        vector<int> pp = Unused_pair_paths[i];
        for(size_t j=0;j<pp.size();j++) out_<<pp[j]<<" ";
	out_<<": "<<Unused_pair_paths_cov[i]<<endl;
    }
    out_<<"#Graph "<<rg_index-1<<endl;
*/  
    /*
    for(size_t i=0;i<edge_weight.size();i++)
    {
	stringstream ss;
	ss<<edge_weight[i];
	string s = ss.str();
	edge_weight[i] = atof(s.c_str());
    }
    for(size_t i=0;i<Node_seq_4.size();i++)
    {
	stringstream ss;
	ss<<Node_seq_4[i];
	string s = ss.str();
	Node_seq_4[i] = atof(s.c_str());
    }
    for(size_t i=0;i<Unused_pair_paths_cov.size();i++) {
	stringstream ss;
	ss<<Unused_pair_paths_cov[i];
	string s = ss.str();
	Unused_pair_paths_cov[i] = atof(s.c_str());
    }
    */
    PathSearch ps(edge_out, edge_in, edge_weight,
		  Unused_pair_paths, Unused_pair_paths_cov,
		  Node_seq_1, Node_seq_2, Node_seq_3,Node_seq_4,
		  strd,ch);
 		  //line,strand,XS_plus,XS_minus);
     ps.path_search();

    // ** search paths ** 
    return;
}
