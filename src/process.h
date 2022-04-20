#include<iostream>
#include<string>
#include<algorithm>
#include<vector>
#include<fstream>
#include<sstream>
#include<map>
#include<cstring>
#include<stdlib.h>
#include "describe_graph_new.h"
//#include "Find_junc_last_map.h"
//#include "edge_node_pair.h"
#include  <boost/unordered_map.hpp> 

extern double UNBALANCE_RATE;
extern bool SingleEndFlag;

double SEP_EXON_RATE = 0.25;
int single_exon_cov = 3;
int single_exon_len = 250;
double partial_exon_min_cov = 5;//
int partial_legth = 20;

int flase_junc_cov = 3;// 

extern int min_delete_graph_cov;

extern multimap<string, vector< pair<int,int> > > readid_cigar_map;
//junction of two adjacent exon
vector<int> False_junc(vector<int> exon_l,vector<int> exon_r)
{
        vector<int> false_junc;
        for(int i=1;i<exon_l.size();i++)
	{
                if(exon_l[i]-exon_r[i-1]==1) false_junc.push_back(exon_r[i-1]);
        }
        return false_junc;
}

//coverage of junction of two adjacent exon
vector<double>  False_junc_cov(vector<int> false_junc,vector<int> seg_l,vector<int>seg_r,vector<int>seg_NH)
{
     vector<double> false_junc_cov(false_junc.size(),0);

     //for(int i=0;i<false_junc.size();i++){false_junc_cov.push_back(0);}

   for(int k=0;k<seg_l.size();k++)
   {
	if(seg_l[k] > false_junc.back() || seg_r[k]-1 < false_junc.front()) continue;

        for(int i=0;i<false_junc.size();i++){
	    if(seg_l[k] <= false_junc[i] && seg_r[k]-1 >= false_junc[i]) 
	    {
               double d = 1.0000 / seg_NH[k];
               false_junc_cov[i] = false_junc_cov[i] + d;
            }
	    //else if(seg_l[k] > false_junc[i]) break; //NEW
        }
     }
     return false_junc_cov;
}

double ave_cov(vector<int>& junc_l, vector<int>& junc_r, vector<double>& junc_cov)
{
    int n = 0;
    double sum_cov=0;
    for(int i = 0;i<junc_l.size();i++)
    {
	if(junc_l[i] == junc_r[i]) continue;

	sum_cov += junc_cov[i];
	n++;
    }

    if(n == 0) return 0;
    double ave = (sum_cov) / n;
    return ave;
}


void delete_junction_by_CovRate(vector<int>& junc_l, vector<int>& junc_r, vector<double>& junc_cov,vector<pair<int,int> > reserved_junc)
{
    vector<int> temp_l = junc_l,temp_r = junc_r;
    //sort(temp_l.begin(),temp_l.end()); sort(temp_r.begin(),temp_r.end());
    if(junc_cov.size()>0){
         //double max_cov=*max_element(junc_cov.begin(),junc_cov.end());
	 double ave = ave_cov(junc_l,junc_r,junc_cov);
	 if(ave == 0) return;
             for(int i=0;i<junc_cov.size();){
                if( junc_cov[i] / ave <= 0.04)//COV_RATE) //raw 0.05
    		{
		    pair<int,int> p = make_pair(junc_l[i],junc_r[i]);
		    if( find(reserved_junc.begin(),reserved_junc.end(),p) != reserved_junc.end() && junc_cov[i] > 5){
			i++;
			continue;//exist in reserved junction
		    }

                    junc_l.erase(junc_l.begin()+i);
                    junc_r.erase(junc_r.begin()+i);
                    junc_cov.erase(junc_cov.begin()+i);
                }
                else i++;
             }
    }
    return;
}

bool graph_delete(vector<int>& junc_l, vector<int>& junc_r, vector<double>& junc_cov)
{   
    double sum_cov = 0;
    int n = 0;
    for(int i = 0;i<junc_l.size();i++)
    {
	if(junc_l[i] == junc_r[i]) continue;
	sum_cov += junc_cov[i];
	n++;
    }
    if(n == 0) return false;
    double ave = (sum_cov) / n;

    if(ave < min_delete_graph_cov ) return true;
    else return false;

}
void seprate_exon( vector<int>& exon_l, vector<int>& exon_r,
		  vector<double>& v_Gene, int gene_l)
{
    int count_partial=0;
    for(int i=1;i<exon_l.size();i++){
        if(exon_l[i]-exon_r[i-1]==1) count_partial++;
    }
    if(exon_l.size()>1&& count_partial>=10)
    {
        for(int i=0;i<exon_l.size()-1;i++)
        {
            if(exon_r[i]-exon_l[i] + 1 >= 200)
            {
                    vector<double> exon_site;
                    for(int k=exon_l[i]-gene_l;k<=exon_r[i]-gene_l;k++)
		    {
                        exon_site.push_back(v_Gene[k]);
                    }
                    double cov_l=0;double cov_r=0;
                    double cov_min=0;int seq_site;
                    for(int j = 0;j < 50;j++) cov_l = cov_l + exon_site[j];
                    for(int j = exon_site.size()-1;j > exon_site.size() - 51;j--) cov_r=cov_r+exon_site[j];

                    if(cov_l < cov_r) {cov_min = cov_l;seq_site = 0;}
                    else {cov_min = cov_r;seq_site = exon_site.size() - 50;}

                    for(int j = 50;j <= exon_site.size()-102;j++)
		    {
                        double win_cov=0;
                        for(int l=j;l<j + 50;l++){
                            win_cov = win_cov + exon_site[l];
                        }
                        if(cov_min > win_cov) {cov_min = win_cov;seq_site = j;}
                    }
                    if(cov_min/cov_l < SEP_EXON_RATE && cov_min/cov_r < SEP_EXON_RATE )
		    {
                        for(int l = exon_l[i] - gene_l + seq_site;l < exon_l[i]-gene_l+seq_site+50;l++){ v_Gene[l]=0; }

                        exon_r.insert(exon_r.begin() + i,exon_l[i]+seq_site-1);
                        exon_l.insert(exon_l.begin()+i+1,exon_l[i]+seq_site+50);
                        i++;
                    }
            }

        }//for each exon
   }
}//seperate exon

void process_gene(vector<int> exon_l,vector<int> exon_r, vector<int> junc_l,vector<int>junc_r,
		  vector<int> seg_l,vector<int> seg_r,vector<int> seg_NH,
		  vector<double> exon_cov,vector<double> junc_cov,
		  vector<double> v_Gene,
		  int strand,int line,int XS_plus,int XS_minus,string Chr,
		  int gene_l,
		  boost::unordered_map<string, pair<cigar_t,cigar_t> >& tu_readid_pairInfo_map,
		  bool abnormal_flag = false)
{
    //seperate exon
    seprate_exon(exon_l,exon_r,v_Gene,gene_l);

    //get exon coverage
    for(int i1=0;i1<exon_l.size();i1++)
    {
          double  k=0.00000;
          for(int i2=exon_l[i1]-gene_l;i2<=exon_r[i1]-gene_l;i2++){
              k=k+v_Gene[i2];
           }
           double cov=k/(exon_r[i1]-exon_l[i1]+1);
           exon_cov.push_back(cov);
    }
	
    if( exon_l.size()==1 && (exon_cov[0] <= single_exon_cov || exon_r[0]-exon_l[0] < single_exon_len))  return;

    const string & basename = base_name();

    //delete partial exon
    if(exon_l.size()>1)
    {
        for(int i = 1;i < exon_l.size()-1;)
	{
           if( (exon_l[i]-exon_r[i-1] == 1|| exon_l[i+1] - exon_r[i] == 1) 
		&& !is_in(exon_l[i],junc_r) //exist NO in junction
		&& !is_in(exon_r[i]+1,junc_l) // exist NO out junction
		&& (exon_cov[i] <= partial_exon_min_cov) )
	   {
	        exon_l.erase(exon_l.begin() + i);
	        exon_r.erase(exon_r.begin() + i);
	        exon_cov.erase(exon_cov.begin() + i);
	   }
           else i++;
	}
    }
	
    //get partial junction 
    vector<int> false_junc;
    vector<double> false_junc_cov;
    false_junc = False_junc(exon_l,exon_r);
    if(false_junc.size() > 0){ false_junc_cov = False_junc_cov(false_junc,seg_l,seg_r,seg_NH); }

    for(int i=0;i<false_junc.size();){
	if(false_junc_cov[i] <= flase_junc_cov){
	    false_junc.erase(false_junc.begin()+i);
	    false_junc_cov.erase(false_junc_cov.begin()+i);
	}
	else i++;
    }
    //add partial junction to junction
    if(junc_l.size() > 0)
    {
	for(int i=0;i<false_junc.size();i++)
	{
	    for(int j=junc_l.size()-1;j>0;j--){
	        if(false_junc[i]+1<=junc_l[j] && false_junc[i]+1>junc_l[j-1]){
		    //cout<<false_junc[i]<<" %%"<<endl;
		    junc_l.insert(junc_l.begin()+j,false_junc[i]+1);
		    junc_r.insert(junc_r.begin()+j,false_junc[i]+1);
 		    junc_cov.insert(junc_cov.begin()+j,false_junc_cov[i]);
		    break;
		}
	    }
			   
	    if(false_junc[i]+1 > junc_l.back()){
		    junc_l.push_back(false_junc[i]+1);
		    junc_r.push_back(false_junc[i]+1);
		    junc_cov.push_back(false_junc_cov[i]);
	    }
	    if(false_junc[i]+1 <= junc_l[0]){
	        junc_l.insert(junc_l.begin(),false_junc[i]+1);
                junc_r.insert(junc_r.begin(),false_junc[i]+1);
                junc_cov.insert(junc_cov.begin(),false_junc_cov[i]);
	    }
	}
     }

    vector<int> partial_junc=false_junc;
    false_junc.clear();false_junc_cov.clear();
//*****************************************************************

   int count_partial=0;
   for(int i=1;i<exon_l.size();i++){
        if(exon_l[i]-exon_r[i-1]==1) count_partial++;
   }
   vector< pair<int,int> > reserved_junc;
   vector< vector<int> > final_pair_path_exon;

   if(junc_l.size()>3 && count_partial < 50)
   {

	PairPath pairpath(exon_l,exon_r,junc_l,junc_r);
	pairpath.find_pair_path(tu_readid_pairInfo_map);
	pairpath.connect_pair_path();

	final_pair_path_exon = pairpath.get_final_path_exon();

	if(count_partial < 10)//10
	{
	    pairpath.get_reserved_junc();
	    reserved_junc = pairpath.reserved_junc;
	}
   }
//**********************************

    //delete junction
    //delete_junction(junc_l,junc_r,junc_cov,exon_l,exon_r,exon_cov,false_junc,false_junc_cov,partial_junc,reserved_junc); //8.14


	
    vector<int> node_l,node_r;
    vector<double> edge_cov;

    for (int i=0;i<junc_l.size();i++)
    {
	    for(int j=0;j<exon_l.size();j++)
	    {
	       if(junc_l[i]-1==exon_r[j]) node_l.push_back(j);
	       if(junc_r[i]==exon_l[j]) { node_r.push_back(j); edge_cov.push_back(junc_cov[i]); }
	    }
    }

   if(exon_l.size()==1){
			//const string & basename = base_name();
	                describe_graph(gene_l,junc_l,junc_r,junc_cov,exon_l,exon_r,exon_cov,
				   Chr,v_Gene,
				   strand,line,XS_plus,XS_minus,tu_readid_pairInfo_map,final_pair_path_exon);
    }
    map<int,int> node_in_all_sub;
    typedef map<int,int>::iterator ite;
    while(node_l.size()>0)
    {
	    map<int,int> con_subgraph;
	    vector<int> exon_l_con,exon_r_con,junc_l_con,junc_r_con,false_junc_con;
	    vector<double>junc_cov_con,false_junc_cov_con,exon_cov_con;

	    con_subgraph.insert(pair<int,int> (node_l[0],1));
	    con_subgraph.insert(pair<int,int> (node_r[0],1));

	    node_l.erase(node_l.begin());node_r.erase(node_r.begin());

	    for(ite it=con_subgraph.begin();it!=con_subgraph.end();it++)
	    {
	   	   for(int i=0;i<node_l.size();)
		   {
			if(node_l[i] == it->first || node_r[i] == it->first){

			    con_subgraph.insert(pair<int,int> (node_l[i],1));
			    con_subgraph.insert(pair<int,int> (node_r[i],1));

			    node_l.erase(node_l.begin()+i);node_r.erase(node_r.begin()+i);
			    edge_cov.erase(edge_cov.begin()+i);

			    it=con_subgraph.begin();
			}
			else i++;
		   }
	    }
	    if(con_subgraph.size() == exon_l.size()){//
		    for(int it=0;it<exon_l.size();it++) node_in_all_sub.insert( pair<int,int> (it,1));
		    exon_l_con=exon_l; exon_r_con=exon_r;  exon_cov_con=exon_cov;
		    junc_l_con=junc_l; junc_r_con=junc_r;  junc_cov_con=junc_cov;
		    false_junc_con=false_junc;  false_junc_cov_con=false_junc_cov;
	    }
	    else{
		    for(ite it=con_subgraph.begin();it!=con_subgraph.end();it++)
		    {
		    	int k=it->first;
			node_in_all_sub.insert(pair<int,int> (k,1));
		    	exon_l_con.push_back(exon_l[k]);
		    	exon_r_con.push_back(exon_r[k]);
		    	exon_cov_con.push_back(exon_cov[k]);
		    	for(int i=0;i<junc_l.size();i++)
			{
			     if(junc_l[i]-1==exon_r[k]){junc_l_con.push_back(junc_l[i]);junc_r_con.push_back(junc_r[i]);junc_cov_con.push_back(junc_cov[i]);}
		 	     if(junc_l[i]-1>exon_r[k]) break;    
		        }	
		    }
	    }

	    //if(graph_delete(junc_l_con,junc_r_con,junc_cov_con)) continue;//NEW //8.14

	    //delete_junction_by_CovRate(junc_l_con,junc_r_con,junc_cov_con,reserved_junc);//NEW; //8.14

	    //const string & basename = base_name();
	    describe_graph(gene_l,junc_l_con,junc_r_con,junc_cov_con,exon_l_con,exon_r_con,exon_cov_con,
		       Chr,v_Gene,
		       strand,line,XS_plus,XS_minus,tu_readid_pairInfo_map,final_pair_path_exon);
    }

    if(exon_l.size() > 1){
	for(int i=0;i<exon_l.size();i++){
	  ite index;
	  index = node_in_all_sub.find(i);
	  if(index == node_in_all_sub.end() && exon_r[i] - exon_l[i] >= single_exon_len && exon_cov[i] > single_exon_cov)
	  {
                vector<int>exon_l_con,exon_r_con,junc_l_con,junc_r_con,false_junc_con;
                vector<double>exon_cov_con,junc_cov_con,false_junc_cov_con;

                exon_l_con.push_back(exon_l[i]); exon_r_con.push_back(exon_r[i]);exon_cov_con.push_back(exon_cov[i]);

                describe_graph(gene_l,junc_l_con,junc_r_con,junc_cov_con,exon_l_con,exon_r_con,exon_cov_con,
		       Chr,v_Gene,
		       strand,line,XS_plus,XS_minus,tu_readid_pairInfo_map,final_pair_path_exon);
	  }
	}
    }//exon_l.size()>1

    return;
}

	


void process_unstrand_gene(vector<int> exon_l,vector<int> exon_r,vector<int> junc_l,vector<int>junc_r,
			   vector<int> seg_l,vector<int> seg_r,vector<int> seg_NH,vector<double> exon_cov,vector<double> junc_cov,
			   vector<double> v_Gene,
			   int strand,int line,int XS_plus,int XS_minus,string Chr,int gene_l,
			   double rate,vector<int> same,
			   boost::unordered_map<string, pair<cigar_t,cigar_t> >& tu_readid_pairInfo_map)
{
    seprate_exon(exon_l,exon_r,v_Gene,gene_l);

 	for(int i1=0;i1<exon_l.size();i1++){
	                double  k=0.00000;
	                    for(int i2=exon_l[i1]-gene_l;i2<=exon_r[i1]-gene_l;i2++){
	                        k=k+v_Gene[i2];
	                    }
	                    double cov=k/(exon_r[i1]-exon_l[i1]+1);
	                    exon_cov.push_back(cov);
	}
    if( exon_l.size()==1 && (exon_cov[0] <= single_exon_cov || exon_r[0]-exon_l[0] < single_exon_len))  return;
    const string & basename = base_name();

    //delete partial
    if(exon_l.size()>1){
      for(int i=1;i<exon_l.size()-1;){
         if( (exon_l[i]-exon_r[i-1]==1|| exon_l[i+1]-exon_r[i]==1)
	     && !is_in(exon_l[i],junc_r) //exist NO in junction
	     && !is_in(exon_r[i]+1,junc_l) //exist NO out junction
	     && (exon_cov[i] <= partial_exon_min_cov) )
	 {
       	     exon_l.erase(exon_l.begin()+i);
       	     exon_r.erase(exon_r.begin()+i);
	     exon_cov.erase(exon_cov.begin()+i);
         }
         else i++;
       }
    }

    vector<int> false_junc;
    vector<double> false_junc_cov;
    false_junc=False_junc(exon_l,exon_r);
    if(false_junc.size()>0){false_junc_cov=False_junc_cov(false_junc,seg_l,seg_r,seg_NH);}

    for(int i=0;i<false_junc.size();){
	if(false_junc_cov[i] <= flase_junc_cov){
 	    false_junc.erase(false_junc.begin()+i);
	    false_junc_cov.erase(false_junc_cov.begin()+i);
	}
	else i++;
    }
    for(int i=0;i<false_junc.size();i++){
	if(is_in(false_junc[i],same)) false_junc_cov[i]=false_junc_cov[i]*rate;
    }

    if(junc_l.size()>0)
    {
		for(int i=0;i<false_junc.size();i++)
		{
			for(int j=junc_l.size()-1;j>0;j--){
			    if(false_junc[i]+1<=junc_l[j] && false_junc[i]+1>junc_l[j-1]){
				//cout<<false_junc[i]<<" %%"<<endl;
				junc_l.insert(junc_l.begin()+j,false_junc[i]+1);
				junc_r.insert(junc_r.begin()+j,false_junc[i]+1);
				junc_cov.insert(junc_cov.begin()+j,false_junc_cov[i]);
				break;
			    }
			}
		    
			if(false_junc[i]+1>junc_l.back()){
			    junc_l.push_back(false_junc[i]+1);
			    junc_r.push_back(false_junc[i]+1);
			    junc_cov.push_back(false_junc_cov[i]);
			}
			if(false_junc[i]+1<=junc_l[0]){
			//cout<<false_junc[i]<<" %%"<<endl;
				     junc_l.insert(junc_l.begin(),false_junc[i]+1);
	                             junc_r.insert(junc_r.begin(),false_junc[i]+1);
	                             junc_cov.insert(junc_cov.begin(),false_junc_cov[i]);
			}
		}
		    
    }
//****************************************************
   int count_partial=0;
   for(int i=1;i<exon_l.size();i++){
        if(exon_l[i]-exon_r[i-1]==1) count_partial++;
   }
   vector<int> partial_junc=false_junc;
   false_junc.clear();false_junc_cov.clear();

   vector<pair<int,int> > reserved_junc;
   vector<vector<int> > final_pair_path_exon;

   if(junc_l.size() > 3 && !SingleEndFlag)
   {

        PairPath pairpath(exon_l,exon_r,junc_l,junc_r);
        pairpath.find_pair_path(tu_readid_pairInfo_map);
        pairpath.connect_pair_path();

        final_pair_path_exon = pairpath.get_final_path_exon();

        if(count_partial < 10)//10
        {
            pairpath.get_reserved_junc();
            reserved_junc = pairpath.reserved_junc;
        }

   }
//***********************************
    //delete junction

   // delete_junction(junc_l,junc_r,junc_cov,exon_l,exon_r,exon_cov,false_junc,false_junc_cov,partial_junc,reserved_junc); //8.14

    //add_virtue_node(junc_l,junc_r,junc_cov,exon_l,exon_r,exon_cov);

    vector<int> node_l,node_r;
    vector<double> edge_cov;
    for (int i=0;i<junc_l.size();i++){
		 for(int j=0;j<exon_l.size();j++){
	            if(junc_l[i]-1==exon_r[j]) node_l.push_back(j);
	            if(junc_r[i]==exon_l[j]) {node_r.push_back(j);edge_cov.push_back(junc_cov[i]);}
	         }
    }
    for(int i=0;i<false_junc.size();i++){
	   	for(int j=0;j<exon_l.size();j++){
	           if(false_junc[i]==exon_r[j]&&false_junc_cov[i]>=2){node_l.push_back(j);node_r.push_back(j+1);edge_cov.push_back(false_junc_cov[i]);}
	 	}
    }
    if(exon_l.size()==1){
		//const string & basename = base_name();
	        describe_graph(gene_l,junc_l,junc_r,junc_cov,exon_l,exon_r,exon_cov,
			   Chr,v_Gene,
			   strand,line,XS_plus,XS_minus,tu_readid_pairInfo_map,final_pair_path_exon);
    }
    map<int,int> node_in_all_sub;
    typedef map<int,int>::iterator ite;
    while(node_l.size()>0)//
    {
		map<int,int> con_subgraph;
		vector<int> exon_l_con,exon_r_con,junc_l_con,junc_r_con,false_junc_con;
		vector<double>junc_cov_con,false_junc_cov_con,exon_cov_con;
		con_subgraph.insert(pair<int,int> (node_l[0],1));
		con_subgraph.insert(pair<int,int> (node_r[0],1));
		node_l.erase(node_l.begin());node_r.erase(node_r.begin());
		for(ite it=con_subgraph.begin();it!=con_subgraph.end();it++){
	   	   for(int i=0;i<node_l.size();){
			if(node_l[i]==it->first || node_r[i]==it->first){
			    con_subgraph.insert(pair<int,int> (node_l[i],1));
			    con_subgraph.insert(pair<int,int> (node_r[i],1));
			    node_l.erase(node_l.begin()+i);node_r.erase(node_r.begin()+i);
			    edge_cov.erase(edge_cov.begin()+i);
			    it=con_subgraph.begin();
			}
			else i++;
		   }
		}
		if(con_subgraph.size()==exon_l.size())//
		{
		    for(int it=0;it<exon_l.size();it++)node_in_all_sub.insert(pair<int,int> (it,1));
		    exon_l_con=exon_l;
		    exon_r_con=exon_r;
		    exon_cov_con=exon_cov;
		    junc_l_con=junc_l;
		    junc_r_con=junc_r;
		    junc_cov_con=junc_cov;
		    false_junc_con=false_junc;
		    false_junc_cov_con=false_junc_cov;
		}
		else
		{
	          for(ite it=con_subgraph.begin();it!=con_subgraph.end();it++){
		    	 int k=it->first;
			 node_in_all_sub.insert(pair<int,int> (k,1));
		    	 exon_l_con.push_back(exon_l[k]);
		    	 exon_r_con.push_back(exon_r[k]);
		    	 exon_cov_con.push_back(exon_cov[k]);
		    	 for(int i=0;i<junc_l.size();i++){
			     if(junc_l[i]-1==exon_r[k]){junc_l_con.push_back(junc_l[i]);junc_r_con.push_back(junc_r[i]);junc_cov_con.push_back(junc_cov[i]);}
		 	     if(junc_l[i]-1>exon_r[k]) break;    
		         }
		    }
		}
		
		//if(graph_delete(junc_l_con,junc_r_con,junc_cov_con)) continue;//NEW //8.14
		//delete junction 
		//delete_junction_by_CovRate(junc_l_con,junc_r_con,junc_cov_con,reserved_junc);//NEW;//8.14

		//const string & basename = base_name();
		describe_graph(gene_l,junc_l_con,junc_r_con,junc_cov_con,exon_l_con,exon_r_con,exon_cov_con,
			   Chr,v_Gene,
			   strand,line,XS_plus,XS_minus,tu_readid_pairInfo_map,final_pair_path_exon);
    }//for each connected subgraph
    if(exon_l.size()>1)
    {
		for(int i=0;i<exon_l.size();i++)
		{
		  ite index;
		  index=node_in_all_sub.find(i);
		  if(index==node_in_all_sub.end()&&exon_r[i]-exon_l[i]>=single_exon_len && exon_cov[i]>single_exon_cov){
//		    const string & basename = base_name();
	            vector<int>exon_l_con,exon_r_con,junc_l_con,junc_r_con,false_junc_con;
	            vector<double>exon_cov_con,junc_cov_con,false_junc_cov_con;
	            exon_l_con.push_back(exon_l[i]);
	            exon_r_con.push_back(exon_r[i]);
	            exon_cov_con.push_back(exon_cov[i]);
	            describe_graph(gene_l,junc_l_con,junc_r_con,junc_cov_con,exon_l_con,exon_r_con,exon_cov_con,
			       Chr,v_Gene,
			       strand,line,XS_plus,XS_minus,tu_readid_pairInfo_map,final_pair_path_exon);
		  }
		}
    }

    return;

}
