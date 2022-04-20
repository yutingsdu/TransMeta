#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>

using namespace std;
void load_transref(char* file,  map<string, map<vector<int>, bool> >& intron_trans_map, map<string, map<pair<int,int>, bool> >& Chr_Junc_map) //map: chrom-> intron,trans
{
    ifstream in(file);
    istringstream istr;
    string s;

    string chr, strand, lable;
    string tranid;
    string temp;
    int exon_l,exon_r;
    vector<int> vecExon;

    getline(in,s);
    istr.str(s);

    istr>>chr>>temp>>lable>>exon_l>>exon_r>>temp>>strand;
    while(istr>>temp) if( temp == "transcript_id") istr>>tranid;
    if(lable == "exon"){ vecExon.push_back(exon_l); vecExon.push_back(exon_r);}
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
	      vecExon.erase(vecExon.begin()); vecExon.pop_back();
	      /*
	      cerr<<tranid<<endl;
	      for(int i=0;i<vecExon.size();i++) cerr<<vecExon[i]<<" ";
	      cerr<<endl;
	      */
	      string cs = chr + strand;
	      for(int i=0;i<int(vecExon.size())-1;)
	      {
	        pair<int,int> J = make_pair(vecExon[i],vecExon[i+1]);
                if(Chr_Junc_map.find(cs) == Chr_Junc_map.end())
                {
                     map<pair<int,int>, bool> m;
                     m[J] = true;
                     Chr_Junc_map[cs] = m;
                }
                else Chr_Junc_map[cs][J] = true;
		i+=2;
	      }
	      
	      chr += strand;
	      if(intron_trans_map.find(chr) == intron_trans_map.end())
	      {
		map<vector<int>,bool> m;
		m[vecExon] = true;
		intron_trans_map[chr] = m;
	      } else {
		intron_trans_map[chr][vecExon] = true;
	      }
	    }
	    vecExon.clear();
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	    chr = curr_chr; strand = curr_strand;
	    tranid = current_id;
	    
	}
    }

    if(vecExon.size() != 2) {
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
     //for(int i=0;i<vecExon.size();i++) cerr<<vecExon[i]<<" ";
     //cerr<<endl;
     chr += strand;
    
     if(intron_trans_map.find(chr) == intron_trans_map.end())
     {
        map<vector<int>,bool> m;
        m[vecExon] = true;
        intron_trans_map[chr] = m;
      } else {
        intron_trans_map[chr][vecExon] = true;
      }
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
