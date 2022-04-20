#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>

using namespace std;
void load_GTF(bool flag, char* file,  
		map<string, map< pair<int,int>,bool> >& Chr_Junc_map,
		map<string, map< pair<int,int>,bool> >& Chr_Range_map, //chrom-> range
		map<string, map< vector<int>,bool> >& Chr_vecExon_map) //map: chrom->trans
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
	    //vecExon.erase(vecExon.begin()); vecExon.pop_back();
	    string cs = chr + strand;
	    //if(cs == "chrY-")
		    //cout<<"vecExon: "<<vecExon.size()<<endl;
	    pair<int,int> range = make_pair(vecExon.front(),vecExon.back());
	    if(Chr_Range_map.find(cs) == Chr_Range_map.end())
	    {
	        map<pair<int,int>, bool >  m;
		m[range] = true;
		Chr_Range_map[cs] = m;
	    }
	    else
	        Chr_Range_map[cs][range] = true; 

	    if(Chr_vecExon_map.find(cs) == Chr_vecExon_map.end())
	    {
	        map<vector<int>,bool> m;
		m[vecExon] = true;
		Chr_vecExon_map[cs] = m;
	    }
	    else Chr_vecExon_map[cs][vecExon] = true;
	    if(flag)
	    {
		for(size_t j=1;j<vecExon.size()-1;)
		{
		    int jl=vecExon[j],jr=vecExon[j+1];
		    pair<int,int> junc = make_pair(jl,jr);
		    if(Chr_Junc_map.find(cs) == Chr_Junc_map.end())
		    {
			map<pair<int,int>, bool >  m;
			m[junc] = true;
			Chr_Junc_map[cs]=m;
		    }
		    else Chr_Junc_map[cs][junc] = true;
		    j=j+2;
		}
	    }
	    vecExon.clear();
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	    chr = curr_chr; strand = curr_strand;
	    tranid = current_id;
	    
	}
    }

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
     //vecExon.erase(vecExon.begin()); vecExon.pop_back();
     string cs = chr + strand;
     pair<int,int> range = make_pair(vecExon.front(),vecExon.back());
     //cout<<"here: "<<range.first<<" "<<range.second<<endl;
     if(Chr_Range_map.find(cs) == Chr_Range_map.end())
     {
         map<pair<int,int>, bool >  m;
	 m[range] = true;
	 Chr_Range_map[cs] = m;
     }
     else Chr_Range_map[cs][range] = true;
     if(Chr_vecExon_map.find(cs) == Chr_vecExon_map.end())
     {
         map<vector<int>,bool> m;
	 m[vecExon] = true;
	 Chr_vecExon_map[cs] = m;
     }
     else Chr_vecExon_map[cs][vecExon] = true;
    return;
}
 
