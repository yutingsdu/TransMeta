#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>

using namespace std;
extern double Coverage;
void load_transref(char* file,  map<string, vector<double> >& id_cov_map) //transid->coverage
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

    while(getline(in,s))
    {
	istr.str(s);
 	string current_id;
	string curr_chr, curr_strand;
	int curr_el, curr_er;

	istr>>curr_chr>>temp>>lable>>curr_el>>curr_er>>temp>>curr_strand;

	if(lable == "transcript") 
	{
	    while(istr>>temp){
		    if( temp == "transcript_id")
			    istr>>current_id;
		    if( temp == "cov")
			    istr>>Cov_s;
	    }
	    Cov = atof(Cov_s.substr(1,Cov_s.length() - 3).c_str());

	    //cout<<current_id<<" "<<Cov<<endl;
	    //return;
	    if(id_cov_map.find(current_id) == id_cov_map.end())
	    {
	        vector<double> v(1,Cov);
		id_cov_map[current_id] = v;
	    }
	    else  id_cov_map[current_id].push_back(Cov);
	}
	istr.clear();
    }

    return;
}
 
/*
int main(int argc,char* argv[])
{
    map<string, map< vector<int>,string> >Chr_vecExon_map;
    map<string, vector<pair<int,int> > > Chr_single_map;
    typedef map<string, map< vector<int>,string> >::iterator iter;

    load_transref(argv[1], Chr_vecExon_map,Chr_single_map);
    cerr<<Chr_vecExon_map.size()<<endl;
    for(iter i = Chr_vecExon_map.begin();i!= Chr_vecExon_map.end();i++)
    {
	cerr<<i->first<<": "<<endl;
	map<vector<int>,string> m = i->second;

	for(map<vector<int>,string>::iterator j = m.begin();j != m.end();j++)
	{
	    vector<int> vec = j->first;
	    cerr<<j->second<<"  ";
	    for(int k = 0;k<vec.size();k++) cerr<<" "<<vec[k];
	    cerr<<endl;
	}
    }
    load_ref(argv[2],Chr_vecExon_map,Chr_single_map);
}
*/
