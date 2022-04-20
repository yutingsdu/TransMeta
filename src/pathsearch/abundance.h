#ifndef ABUNDANCE_H
#define ABUNDANCE_H

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<map>
#include<set>
#include<iterator>
#include<stdlib.h>

using namespace std;
extern bool MyFlag;
class Abundance
{
	private:
	const char * path_tsv;
	const char * path_gtf;
	const char * path_out;
	const char * min_cov;
	public:
	ofstream file;
	Abundance(const char * p1,const char * p2,const char * p3,const char * p4);
	void Process();

};
Abundance::Abundance(const char * p1,const char * p2,const char * p3,const char * p4)
{
	path_tsv=p1;
	path_gtf=p2;
	path_out=p3;
	min_cov=p4;
	file.open(path_out);
}
void Abundance::Process()
{
//1: abundance.tsv; 2: temp.gtf; 3: min_cov;
map<string,float > Trans_abund;
ifstream ifs1(path_tsv);
istringstream istr;
string temp;
string trans_id;
float len,eff_len,cov,tmp;
getline(ifs1,temp);
while (getline(ifs1,temp))
{
	istr.str(temp);
	istr>>trans_id>>len>>eff_len>>cov>>tmp;
	istr.clear();
	double min_cov_ = atof(min_cov);
	if(MyFlag) min_cov_ += 0.5;
	if (tmp>0 && tmp>=min_cov_)
	{
		Trans_abund[trans_id]=tmp;
	}
}
ifs1.close();
ifstream ifs2(path_gtf);
string chr,name,exon,start_pos,end_pos,score,plus,dot,gene,gene_id,trans;
while (getline(ifs2,temp))
{
	istr.str(temp);
	istr>>chr>>name>>exon>>start_pos>>end_pos>>score>>plus>>dot>>gene>>gene_id>>trans>>trans_id;
	istr.clear();
	if (Trans_abund.find(trans_id.substr(1,trans_id.size()-3)) != Trans_abund.end())
	{
		file<<chr<<"	"<<name<<"	"<<exon<<"	"<<start_pos<<"	"<<end_pos<<"	"<<score<<"	"<<plus<<"	"<<dot<<"	"<<gene<<" "<<gene_id<<" "<<trans<<" "<<trans_id<<" "<<"tpm"<<" "<<Trans_abund[trans_id.substr(1,trans_id.size()-3)]<<";"<<endl;
	}
}
ifs2.close();
return;
}

#endif

















