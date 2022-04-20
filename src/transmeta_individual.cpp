#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include "transmeta_individual.h"
using namespace std;
typedef vector<int> path_t;
ofstream outgtf;
int SampleSize;
double Coverage = 0;
double SEED = 0;
bool MyFlag = false;
int ZN=0, nZN=0;
struct info
{
    int number;
    double seed;
    int normal_edge_number;
    int partial_edge_number;
    int mj_number;
    int Nmj_number;
    int seed_sample_number;
    int path_number;
    int path_index;
};
map<string, vector<double> > id_cov_map;
map<string, info > id_info_map;

void load_info(char*file)
{
    ifstream in(file);
    string s;
    istringstream istr;
    while(getline(in,s))
    {
        istr.str(s);
	string temp,id;
	double seed;
	int normal_edge,partial_edge,mj,Nmj,seed_sample_number,path_number,path_index;
	istr>>temp>>id>>seed>>normal_edge>>partial_edge>>temp>>temp>>temp>>mj>>Nmj>>seed_sample_number>>temp>>path_number>>path_index;;
	istr.clear();
	int N = atoi(id.substr(id.length() - 3,1).c_str());
	info info_={N,seed,normal_edge,partial_edge,mj,Nmj,seed_sample_number,path_number,path_index};
	id_info_map[id] = info_;

	double cov = 0;
        if(id_cov_map.find(id) != id_cov_map.end()){
          for(size_t i=0;i<id_cov_map[id].size();i++)
                cov += id_cov_map[id][i];
        }
        if(cov == 0) ZN++;
        else nZN++;
    }
    //cout<<"id_info_map.size(): "<<id_info_map.size()<<endl;
}
void process(string tranid,vector<string> oneTrans)
{
    double r1,r2,r3,r4,r5;
    r1=0.0;r2=0.5;r3=1.5;r4=5.0;r5=10.0;
    if(MyFlag) {
	r1=0;r2=0.5;r3=1;r4=1.5;r5=4.0;
    }
    /*
    if(SampleSize>30){
        r1=0.1;r2=1.5;r3=5.0;r4=10.0;r5=20.0;
    }
    else {
        r1=2.5;r2=1.0;r3=5.0;r4=10.0;r5=20.0;
    }*/
    if(id_info_map.find(tranid) != id_info_map.end() )//&& id_cov_map.find(tran_id) != id_cov_map.end())
    {
	    double cov = 0;
	    for(size_t i=0;i<id_cov_map[tranid].size();i++)
		    cov += id_cov_map[tranid][i];

	    bool flag = false;

	    info ti = id_info_map[tranid];
	    if(MyFlag && SEED>=0 && ti.seed_sample_number <= 2) return;
            if(MyFlag && ti.Nmj_number >=2) return;
	    if(SEED>10) {r1=0.1;}
	//    if(SEED>5 && ti.seed_sample_number <= 0.15*SampleSize) {return;}
	//    if(SEED>10){r1=0.01;}
//	    if(SEED>=1 && ti.seed_sample_number < SampleSize) return;
//	    if(ti.Nmj_number >=1) return;
//	    if( ti.normal_edge_number<=1)return;
	    
	    if(id_info_map[tranid].number == 1)
	        if(cov >= r1*SEED) flag = true;

	    if(id_info_map[tranid].number == 2)
	        if(cov >= r2*SEED) flag = true;

	    if(id_info_map[tranid].number == 3)
	        if(cov > r3*SEED) flag = true;

	    if(id_info_map[tranid].number ==4 )
	        if(cov > r4*SEED) flag = true;

	    if(id_info_map[tranid].number >4 )
		    if(cov > r5*SEED) flag = true;
	    if(!MyFlag && cov == 0 && ZN > 1.2*nZN && (ti.seed <= SEED ||ti.seed<=5) ) flag = false;
            if(!MyFlag && cov == 0 &&  ti.seed<=1) flag = false;    
	    if(MyFlag && ti.path_number > 8 && oneTrans.size() - 1 > 25) flag = false;
	    if(flag)  
	    {
		    for(size_t i=0;i<oneTrans.size();i++)
		    outgtf<<oneTrans[i]<<endl;
	    }
     }
    return;
}
void get_final_results(char*file)
{

    ifstream in(file);
    istringstream istr;
    string s;

    string chr, strand, lable;
    string tranid;
    string temp, Cov_s;
    double Cov;
    int exon_l,exon_r;
    vector<int> vecExon;
    vector<string> oneTrans;
    getline(in,s);
    istr.str(s);
    while(istr>>temp)
	    if( temp == "transcript_id")
		    istr>>tranid;
    istr.clear();
    oneTrans.push_back(s);

    while(getline(in,s))
    {
	istr.str(s);
 	string current_id;
	string curr_chr, curr_strand;
	int curr_el, curr_er;

	istr>>curr_chr>>temp>>lable>>curr_el>>curr_er>>temp>>curr_strand;


	while(istr>>temp)
	{
	    if( temp == "transcript_id")
		istr>>current_id;
	}
	if(current_id ==  tranid)
	{
	  oneTrans.push_back(s);
	}
	else
	{
	  process(tranid,oneTrans);
	  tranid = current_id;
	  oneTrans.clear();
	  oneTrans.push_back(s);
	}

	istr.clear();
    }
    process(tranid,oneTrans);
    return;

}
//./exe a.gtf b.gtf .. N.gtf input.info input.gtf output.gtf DEED //only remove repeat ones
//
int main(int argc,char* argv[])
{
    //load_transref(argv[1],intron_trans_map,Chr_Junc_map);
    //cout<<argc<<endl;
    SampleSize = 1;

    string sFlag = argv[6];
    if(sFlag=="1") MyFlag=true;
 
    
    string S = argv[5];

    SEED = atof(S.c_str());
    if(MyFlag) SEED=SEED+0.5;
    load_transref(argv[1],id_cov_map);
    load_info(argv[2]);
    cout<<ZN<<" "<<nZN<<endl;
    outgtf.open(argv[4]);
  

    get_final_results(argv[3]);
    outgtf.close();
    
    return 0;
}
