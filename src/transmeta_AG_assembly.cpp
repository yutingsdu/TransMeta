#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include "transmeta_AG_assembly.h"
using namespace std;
typedef vector<int> path_t;
map<string, map<pair<int,int>, bool> > Chr_Junc_map;
map<string, map<vector<int>, bool> > intron_trans_map;
map<string, map<pair<int,int>, bool> > Chr_Range_map;
map<string, map<vector<int>, bool> > Chr_vecExon_map;
int gene_index = 0;
int tr_index = 1;
string strand="";
string chr="";
ofstream outgtf;
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
        if(I % 100000 == 0) cerr<<"Loading "<<I<<" lines"<<endl;
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
		if(edge_weight[i] <=0.5) continue;
		//if(edge_weight[i] <=1) continue;
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

	while(istr>>temp) if( temp == "transcript_id") istr>>current_id;
	istr.clear();
	if(current_id == tranid)//commom transcript
	{
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	    one_trans.push_back(s);
	}
	else 
	{
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
	      vecExon.erase(vecExon.begin()); vecExon.pop_back();
	      sort(vecExon.begin(),vecExon.end());
	      string cs = chr+strand;
	      bool flag = true;
	      if(flag){ 
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
	    vecExon.clear();
	    one_trans.clear();
	    one_trans.push_back(s);
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	    chr = curr_chr; strand = curr_strand;
	    tranid = current_id;
	    
	}
    }

    sort(vecExon.begin(),vecExon.end());
    vecExon.erase(vecExon.begin()); vecExon.pop_back();
    string cs = chr+strand;
     if(intron_trans_map.find(cs) != intron_trans_map.end())
	     for(size_t j=0;j<one_trans.size();j++) cout<<one_trans[j]<<'\n';
    return;
}
 

//./exe completed_frame_graph reference.gtf > covered_gtf

  bool compatible(path_t P1,path_t P2)//p1.size() < p2.size() !!!
  {
        //if(P1 == P2) return true;
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
bool Cover(vector<int> vecExon,string chr, string strand)
{
    string cs = chr+strand;
    //if(Chr_Junc_map.find(cs) == Chr_Junc_map.end()) return false;
    return true;
    for(size_t i=1;i<vecExon.size()-1;)
    {
	pair<int,int>junc=make_pair(vecExon[i],vecExon[i+1]);

	if(Chr_Junc_map[cs].find(junc) != Chr_Junc_map[cs].end()) 
		return true;
	i=i+2;
    }   
    return false;

}
void output(vector<int> vecExon,string chr, string strand,string tranid)
{
	     if(!Cover(vecExon,chr,strand)) return;

	      outgtf<<chr<<"	"<<"GINGKO"<<"	"
		  <<"transcript"<<"	"<<vecExon.front()<<"	"<<vecExon.back()<<"	"
		  <<1000<<"	"<<strand<<"	.	"
		  <<"gene_id "<<"\""<<tranid<<"."<<gene_index<<"\""<<"; "
		  <<"transcript_id "<<"\""<<tranid<<"."<<gene_index<<"."<<tr_index<<"\""<<"; "<<endl;
	      for(size_t j=0;j<vecExon.size();)
	      {
	      
	          outgtf<<chr<<"	"<<"GINGKO"<<"	"
		      <<"exon"<<"	"<<vecExon[j]<<"	"<<vecExon[j+1]<<"	"
		      <<1000<<"	"<<strand<<"	.	"
		      <<"gene_id "<<"\""<<tranid<<"."<<gene_index<<"\""<<"; "
		      <<"transcript_id "<<"\""<<tranid<<"."<<gene_index<<"."<<tr_index<<"\""<<"; "
		      <<"exon_number "<<"\""<<(j/2) + 1<<"\""<<";"<<endl;
		  j+=2;
	      }
}
void process_one_gene(vector<vector<int> > gene)
{
    vector<int> empty;
    map<pair<int,int>,bool> single_exon;
    for(size_t i=0;i<gene.size();)
    {
        if(gene[i].size() == 2)
	{
	    single_exon[make_pair(gene[i][0],gene[i][1])] = true;
	    gene.erase(gene.begin() + i);
	}
	else i++;
    }
    //cout<<"process_one_gene: "<<gene.size()<<endl;
    for(size_t i=0;i<gene.size();i++)
    {
	//cout<<i<<endl;
	//for(int k = 0;k<gene[i].size();k++) cout<<gene[i][k]<<" ";
	//cout<<endl;
	bool flag = true;
	if(gene[i].empty()) continue;
	for(size_t j=0;j<gene.size();j++)
	{
	    if(i == j) continue;
	    if(gene[j].size() <= 2) continue;
	    vector<int> intron1 = gene[i],intron2 = gene[j];
	    if(intron2.size() == 2) continue;
	    intron1.erase(intron1.begin());intron1.pop_back();
	    intron2.erase(intron2.begin());intron2.pop_back();
	    if(intron1 == intron2) {
		    gene[i]=empty;
		    flag = false;
		    break;
	    }
	    /*
	    if(compatible(intron1,intron2))
	    {
	        flag = false;
		break;
	    }
	    */
	}
	if(flag)
	{
	    output(gene[i],chr,strand,"Ting");
	    tr_index++;
	}
    }
    if(single_exon.size() == 0) return;
    int l = single_exon.begin()->first.first;
    int r = single_exon.begin()->first.second;
    for(map<pair<int,int>,bool>::iterator i = single_exon.begin();i!=single_exon.end();i++)
    {
	if(i->first.first > r)
	{
	    vector<int> v;
	    v.push_back(l);v.push_back(r);
	    output(v,chr,strand,"Ting");
	    tr_index++;
	}
	else if(i->first.second > r)
	    r = i->first.second;
    }
}
typedef map<string, map<pair<int,int>, bool> >::iterator iter1;
typedef map<string, map<vector<int>, bool> >::iterator iter2;
void simple_merge()
{
    for(iter1 i = Chr_Range_map.begin();i != Chr_Range_map.end();i++)
    {
	int l=-1,r=-1;
	map<pair<int,int>, bool> & Range_map = i->second;
	map<vector<int>, bool>  & vecExon_map = Chr_vecExon_map[i->first];
	map<vector<int>, bool> ::iterator k = vecExon_map.begin();
	l = Range_map.begin()->first.first;
	r = Range_map.begin()->first.second;
	chr = i->first.substr(0,i->first.length() - 1);
	strand = i->first.substr(i->first.length() - 1);
	cout<<chr<<" "<<strand<<endl;
        for(map<pair<int,int>, bool>  ::iterator j = Range_map.begin();j!=Range_map.end();j++)
	{
	    if(j->first.first > r)
	    {
		gene_index++;
		tr_index = 1;
		vector<vector<int> > gene;
		//cout<<"range: "<<l<<" "<<r<<endl;
	        for(;k != vecExon_map.end();k++)
		{
		    //cout<<k->first.front()<<" "<<k->first.back()<<endl;
		    if(k->first.front() >= l && k->first.back() <= r)
		    {
		        gene.push_back(k->first);
		    }
		    else{
			    //process a gene;
			    //process_one_gene(gene);
			    break;
		    }
		}
		process_one_gene(gene);
		l = j->first.first;
		r = j->first.second;
	    }
	    else if( r < j->first.second) 
		    r = j->first.second;
	}
	gene_index++;
	tr_index = 1;
	vector<vector<int> > gene;
	for(;k != vecExon_map.end();k++) gene.push_back(k->first);
	process_one_gene(gene);
    }

}
//./exe a.gtf b.gtf output.gtf //only remove repeat ones
int main(int argc,char* argv[])
{
    
    outgtf.open(argv[3]);
    load_GTF(true,argv[1],Chr_Junc_map,Chr_Range_map,Chr_vecExon_map); //GINgko assembly
   cout<<Chr_Junc_map.size()<<endl;

   map<string, map<pair<int,int>, bool> >::iterator i = Chr_Junc_map.begin();
   for(;i!=Chr_Junc_map.end();i++) 
	cout<<i->first<<" "<<i->second.size()<<endl;

    load_GTF(false,argv[2],Chr_Junc_map,Chr_Range_map,Chr_vecExon_map); //ANNOtation
 cout<<Chr_Junc_map.size()<<endl;
    simple_merge();
    outgtf.close();
    return 0;
}
