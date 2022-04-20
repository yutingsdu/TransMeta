#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include "transmeta_merge.h"
using namespace std;
typedef vector<int> path_t;
ofstream outgtf;
int SampleSize;
double Coverage = 0;
double SEED = 0;
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
	istr>>temp>>id>>seed>>normal_edge>>partial_edge>>temp>>temp>>temp>>mj>>Nmj>>seed_sample_number>>temp>>path_number>>path_index;
	istr.clear();
	int N = atoi(id.substr(id.length() - 3,1).c_str());
	info info_={N,seed,normal_edge,partial_edge,mj,Nmj,seed_sample_number,path_number,path_index};
	id_info_map[id] = info_;
    }
    //cout<<"id_info_map.size(): "<<id_info_map.size()<<endl;
}
void process(string tranid,vector<string> oneTrans)
{
    double r1,r2,r3,r4,r5;
    //r1=1.0;r2=1.5;r3=5.0;r4=10.0;r5=20.0; //RAW
    bool flag1=false;
    if(SampleSize <10)
    {
	//r1=0.1;r2=0.15;r3=0.5;r4=1;r5=2;
	r1=0.3;r2=0.45;r3=1.5;r4=3;r5=6;
    }
    if(SampleSize >= 10 && SampleSize < 20)
    {
        r1=0.75;r2=1.5;r3=4;r4=7.5;r5=15;
	flag1=true;
    }
    if(SampleSize>=20 &&  SampleSize < 50)
    {   
        r1=1.0;r2=2.5;r3=7.5;r4=15.0;r5=30.0;
    }
    if(SampleSize >=50){
	r1=1.5;r2=4;r3=10;r4=25.0;r5=40.0;
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
	    cov = cov/(1.0*SampleSize);
	    bool flag = false;

	    info ti = id_info_map[tranid];
	    if(SampleSize>20 && SEED>=1 && SEED<2 && ti.seed_sample_number < 0.1*SampleSize) return;
	    if(SampleSize>20 && SEED>=2 && SEED<50 && ti.seed_sample_number < 0.15*SampleSize) return;
	    if(SampleSize>20 && SEED>=50 && ti.seed_sample_number < 0.2*SampleSize) return;    
	    if(flag1  && SEED>=2 && ti.seed_sample_number < 0.2*SampleSize) return;
	//    if(SampleSize>20 && SEED>=8 && ti.Nmj_number >=1) return;
	//    if(SampleSize>20 && SEED>=8 && ti.normal_edge_number<=1)return;
	    
	    if(id_info_map[tranid].number == 1)
	        if(cov > r1*SEED/(SampleSize)) flag = true;

	    if(id_info_map[tranid].number == 2)
	        if(cov > r2*SEED/(1.0*SampleSize)) flag = true;

	    if(id_info_map[tranid].number == 3)
	        if(cov > r3*SEED/(1.0*SampleSize)) flag = true;

	    if(id_info_map[tranid].number ==4)
	        if(cov > r4*SEED/(1.0*SampleSize)) flag = true;

	    if(id_info_map[tranid].number >4)
		    if(cov > r5*SEED/(1.0*SampleSize)) flag = true;
		    
	    if(ti.path_number > 10 && oneTrans.size() - 1 > 25) flag = false;
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
//./exe a.gtf b.gtf .. N.gtf input.info input.gtf output.gtf SEED sampleSize //only remove repeat ones
int main(int argc,char* argv[])
{
    //load_transref(argv[1],intron_trans_map,Chr_Junc_map);
    //cout<<argc<<endl;
    string SampleS=argv[argc-1];
    SampleSize = atoi(SampleS.c_str());
    string S = argv[argc-2];
    SEED = atof(S.c_str());
    //cout<<"filter: "<<SEED<<endl;
    //cout<<"SampleSize: "<<SampleSize<<endl;
    load_info(argv[argc-5]);

    outgtf.open(argv[argc-3]);

    for(int i=1;i<=argc-6;i++){
	//cout<<"load: "<<argv[i]<<" "<<i<<" sample..."<<endl;
        load_transref(argv[i],id_cov_map);
    
    }

    get_final_results(argv[argc-4]);
    outgtf.close();
    return 0;
}
