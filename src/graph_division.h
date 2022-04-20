#include"Function.h"

using namespace std;

class Graph_division
{
	private: vector<int>junc_plus_l,junc_plus_r,junc_minus_l,junc_minus_r;
		 vector<int> exon_plus_l,exon_plus_r,exon_minus_l,exon_minus_r;
		 vector<int> same_false_junc;
		 vector<double> junc_cov_p,junc_cov_m;

	public:
		void division(vector<int>exon_l,vector<int> exon_r,vector<int>junc_l,vector<int>junc_r,vector<double> junc_cov,vector<double> junc_cov_plus,vector<double>junc_cov_minus);
		vector<int> get_junc_plus_l();
		vector<int> get_junc_plus_r();
		vector<int> get_junc_minus_l();
		vector<int> get_junc_minus_r();
		vector<int> get_exon_plus_l();
                vector<int> get_exon_plus_r();
                vector<int> get_exon_minus_l();
                vector<int> get_exon_minus_r();
		vector<int> get_same_false_junc();
		vector<double> get_junc_cov_plus();
		vector<double> get_junc_cov_minus();
};
void Graph_division::division(vector<int>exon_l,vector<int> exon_r,vector<int>junc_l,vector<int>junc_r,vector<double> junc_cov,vector<double> junc_cov_plus,vector<double>junc_cov_minus){
	junc_plus_l=junc_l;
	junc_plus_r=junc_r;
	junc_minus_l=junc_l;
	junc_minus_r=junc_r;
	exon_plus_l=exon_l;
        exon_plus_r=exon_r;
        exon_minus_l=exon_l;
        exon_minus_r=exon_r;
	junc_cov_p=junc_cov_plus;
	junc_cov_m=junc_cov_minus;
//cout<<"here0"<<endl;
//cout<<junc_cov_p.size()<<" "<<junc_plus_l.size()<<" "<<junc_plus_r.size()<<" "<<junc_cov.size()<<endl;
	for(int i=0;i<junc_cov_p.size();){
		if(junc_cov_p[i]==0) {
			junc_plus_l.erase(junc_plus_l.begin()+i);
			junc_plus_r.erase(junc_plus_r.begin()+i);
			junc_cov_p.erase(junc_cov_p.begin()+i);
		}
		else i++;
	}
	for(int i=0;i<junc_cov_m.size();){
                if(junc_cov_m[i]==0) {
                        junc_minus_l.erase(junc_minus_l.begin()+i);
                        junc_minus_r.erase(junc_minus_r.begin()+i);
                        junc_cov_m.erase(junc_cov_m.begin()+i);
                }
                else i++;
        }
//cout<<"here1"<<endl;
/* //10.5noref
	for(int i=0;i<exon_plus_l.size();){
		if(!is_in(exon_plus_l[i],junc_plus_r) && !is_in(exon_plus_r[i]+1,junc_plus_l)){
			exon_plus_l.erase(exon_plus_l.begin()+i);
			exon_plus_r.erase(exon_plus_r.begin()+i);
		}
		else i++;
	}
	for(int i=0;i<exon_minus_l.size();){
                if(!is_in(exon_minus_l[i],junc_minus_r) && !is_in(exon_minus_r[i]+1,junc_minus_l)){
                        exon_minus_l.erase(exon_minus_l.begin()+i);
                        exon_minus_r.erase(exon_minus_r.begin()+i);
                }
                else i++;
        }
*/
//cout<<"here2"<<endl;
	vector<int> false_junc,false_junc_p,false_junc_m;
//cout<<exon_plus_l.size()<<endl;
	for(int i=0;i<(int)(exon_plus_l.size()-1);i++){
                if(exon_plus_l[i+1]-exon_plus_r[i]==1) false_junc_p.push_back(exon_plus_r[i]);
        }
//cout<<"h"<<endl;
	for(int i=0;i<(int)(exon_minus_l.size()-1);i++){
                if(exon_minus_l[i+1]-exon_minus_r[i]==1) false_junc_m.push_back(exon_minus_r[i]);
        }
	for(int i=0;i<false_junc_p.size();i++){
		for(int j=0;j<false_junc_m.size();j++){
			if(false_junc_m[j]==false_junc_p[i]) {same_false_junc.push_back(false_junc_p[i]);break;}
		}
	}	
}
vector<int> Graph_division::get_junc_plus_l(){return junc_plus_l;}
vector<int> Graph_division::get_junc_plus_r(){return junc_plus_r;}
vector<int> Graph_division::get_junc_minus_l(){return junc_minus_l;}
vector<int> Graph_division::get_junc_minus_r(){return junc_minus_r;}
vector<int> Graph_division::get_exon_plus_l(){return exon_plus_l;}
vector<int> Graph_division::get_exon_plus_r(){return exon_plus_r;}
vector<int> Graph_division::get_exon_minus_l(){return exon_minus_l;}
vector<int> Graph_division::get_exon_minus_r(){return exon_minus_r;}
vector<int> Graph_division::get_same_false_junc(){return same_false_junc;}
vector<double> Graph_division::get_junc_cov_plus(){return junc_cov_p;}
vector<double> Graph_division::get_junc_cov_minus(){return junc_cov_m;}
