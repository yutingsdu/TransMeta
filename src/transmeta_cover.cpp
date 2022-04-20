#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include "transmeta_cover.h"
using namespace std;
typedef vector<int> path_t;
map<string, map<pair<int,int>, bool> > Chr_Junc_map;
map<string, map<vector<int>, bool> > intron_trans_map;
string strand="";
void load_ref(char* file)
{
    //cout<<file<<endl;
    ifstream in(file);
    istringstream istr;
    string s,temp;
    bool edge_flag = false, node_flag = false, pair_flag = false;
    string strd = "", chr = "";
    vector<double> gene_vec;
    int start_pos=0;
    vector<int> exon_l,exon_r,junc_l,junc_r,edge_out,edge_in;
    vector<double> exon_cov, junc_cov, edge_weight;
    vector<pair<int,int> > junction, exon;
    int I = 0;
    vector<path_t> ppaths;
    vector<double> ppaths_cov;
    while(getline(in,s))
    {
	I++;
        //if(I % 100000 == 0) cerr<<"Loading "<<I<<" lines"<<endl;
        if( s == "Edges"){ edge_flag = true; node_flag = false;pair_flag = false;continue;}
        if( s == "Nodes") { edge_flag = false; node_flag = true;pair_flag = false;continue;}
        if( s == "Pair") {edge_flag = false; node_flag = false;pair_flag = true;continue;}
        if( s[0] == '#') {
            istr.str(s);
	    bool invalid_flag = false;
	    
	    for(int i=0;i<edge_out.size();i++)
	    {
		//junction.push_back(make_pair(exon_r[edge_out[i]], exon_l[edge_in[i]]));
		pair<int,int> junc = make_pair(exon_r[edge_out[i]], exon_l[edge_in[i]]);
		if(junc.second - junc.first == 1) continue;
		//cout<<strd<<" "<<chr<<" "<<junc.first<<" "<<junc.second<<endl;
		if(edge_weight[i] <=1) continue; //Merge
		string cs = chr + strd;
		pair<int,int> J = make_pair(junc.first,junc.second);
		if(Chr_Junc_map.find(cs) == Chr_Junc_map.end())
		{
		    map<pair<int,int>, bool> m;
		    m[J] = true;
		    Chr_Junc_map[cs] = m;
		}
		else Chr_Junc_map[cs][J] = true;
	    }
	
	    istr.clear();
	    junc_l.clear();junc_r.clear();junc_cov.clear(); gene_vec.clear();
	    junction.clear();
	    exon_l.clear();exon_r.clear();exon_cov.clear();
	    exon.clear();
	    edge_out.clear();edge_in.clear();edge_weight.clear();
	    ppaths.clear(); ppaths_cov.clear();
            edge_flag = false; node_flag = false;pair_flag = false;continue;
        }

        if(edge_flag)
        {
            istr.str(s);
            int out,in;
            vector<int> v;
            double cov;
            istr>>out>>temp>>in>>temp>>cov;
            istr.clear();
            edge_out.push_back(out); edge_in.push_back(in);
            edge_weight.push_back(cov);
        }
        if(node_flag)
        {
 	    istr.str(s);
            int left, right;
            double cov;
            istr>>strd>>chr>>left>>right>>cov;
	    if( cov == 0) cov = 0.01;
            istr.clear();
            exon_l.push_back(left); exon_r.push_back(right);exon_cov.push_back(cov);
	    exon.push_back(make_pair(left,right));
	    if(gene_vec.empty()) start_pos = left;
	    int S = gene_vec.size();
	    for(int i=gene_vec.size(); i < left-start_pos;i++) gene_vec.push_back(0); //intron 
	    for(int i=gene_vec.size(); i <= right-start_pos;i++) gene_vec.push_back(cov); //exon 

        }  
	if(pair_flag)
        {
            istr.str(s);
            double cov;
            vector<int> v;
            while(istr>>temp)
            {
                if(temp == ":"){
                    istr>>cov;
                    break;
                }
                else
                    v.push_back(atoi(temp.c_str()) );
            }
            istr.clear();
        }
    }
    return;
}
//./exe ReferenceGraph >junction(- chr 12234 22345)
void load_gtf(char* file)
{
    ifstream in(file);
    istringstream istr;
    string s;

    string chr, strand, lable;
    string tranid;
    string temp;
    int exon_l,exon_r;
    vector<int> vecExon;
    vector<string> one_trans;

    getline(in,s);
    istr.str(s);

    istr>>chr>>temp>>lable>>exon_l>>exon_r>>temp>>strand;
    while(istr>>temp) if( temp == "transcript_id") istr>>tranid;
    if(lable == "exon"){ one_trans.push_back(s);vecExon.push_back(exon_l); vecExon.push_back(exon_r);}
    istr.clear();


    while(getline(in,s))
    {
	istr.str(s);
 	string current_id;
	string curr_chr, curr_strand;
	int curr_el, curr_er;

	istr>>curr_chr>>temp>>lable>>curr_el>>curr_er>>temp>>curr_strand;
	if(lable != "exon") continue;
//	cout<<s<<endl;
	while(istr>>temp) if( temp == "transcript_id") istr>>current_id;
	istr.clear();
//	cout<<tranid<<endl;
	if(current_id == tranid)//commom transcript
	{
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	    one_trans.push_back(s);
	}
	else 
	{
	    //if(vecExon.size() == 2) continue;
	    if(vecExon.size() != 2)
	    {
	      sort(vecExon.begin(),vecExon.end());
	      for(size_t j=1;j<vecExon.size()-1;)
	      {
	          if(vecExon[j+1] - vecExon[j] == 1)
		  {
		      vecExon.erase(vecExon.begin() +j );
		      vecExon.erase(vecExon.begin() +j );
		  }
		  else j+=2;
	      }
	      if(vecExon.size() != 2)
	      {
	          vecExon.erase(vecExon.begin()); vecExon.pop_back();
	          sort(vecExon.begin(),vecExon.end());
	          string cs = chr+strand;
	          map<pair<int,int>,bool> junc_map;
	          if(Chr_Junc_map.find(cs) != Chr_Junc_map.end()) junc_map = Chr_Junc_map[cs];
	          bool flag = true;
	          int notCoveredEdgeNum = 0;
	          for(size_t i=0;i<vecExon.size()-1;)
	          {
	            pair<int,int>J = make_pair(vecExon[i],vecExon[i+1]);
		    if( !junc_map[J] ){
			  //flag = false;
			  //break;
			  notCoveredEdgeNum++;
	 	    }
		    i += 2;
	          }
	          int edgeNumber = vecExon.size()/2;
	          if(notCoveredEdgeNum == edgeNumber) flag = false;
	          if(flag)
		  { 
		    if(intron_trans_map.find(cs) != intron_trans_map.end())
		    {
		      if(intron_trans_map[cs].find(vecExon) == intron_trans_map[cs].end())//not in transref
		      {
		    	for(size_t j=0;j<one_trans.size();j++) cout<<one_trans[j]<<'\n';
		      }
		      else cerr<<chr<<" "<<strand<<" "<<tranid<<endl;
		    }
		    else
		    {
		      for(size_t j=0;j<one_trans.size();j++)
		    	cout<<one_trans[j]<<'\n';
		    }      
	      
	          }
		 }
	    }
	    vecExon.clear();
	    one_trans.clear();
	    one_trans.push_back(s);
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	    chr = curr_chr; strand = curr_strand;
	    tranid = current_id;
	    
	}
    }

    if(vecExon.size() != 2) {
     sort(vecExon.begin(),vecExon.end());
     vecExon.erase(vecExon.begin()); vecExon.pop_back();
    }
    return;
}
 

//./exe completed_frame_graph reference.gtf > covered_gtf
int main(int argc,char* argv[])
{
    //load_transref(char* file,  map<string, map<vector<int>, bool> >& intron_trans_map) //map: chrom-> intron,trans
    strand = argv[3];
    //load_transref(argv[1],intron_trans_map,Chr_Junc_map);
    if(strand != "unstranded") Chr_Junc_map.clear();
    load_ref(argv[1]);
    load_gtf(argv[2]);
    return 0;
}
