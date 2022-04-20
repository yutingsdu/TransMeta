
#ifndef UTILITY_H
#define UTILITY_H
/*
   ADD_JUNCTION
   DELETE_ABNORMAL_JUNCTION
   JUNCTION DESCIDE EXON
  */
#include<iostream>
#include<string>
#include<algorithm>
#include<vector>
#include<fstream>
#include<sstream>
#include<map>
#include<cstring>
#include<stdlib.h>
#include "process.h"
#include"graph_division.h"
#include"pair_path.h"

#include <boost/unordered_map.hpp>
//#include"describe_graph.h"
	
using namespace std;
// ** parameters **
extern int MERGE_PARA;
extern float multi_map_frac;
extern int gr_length;
//extern multimap<string, vector< pair<int,int> > > readid_cigar_map;
//extern map<string, pair<cigar_t,cigar_t> > readid_pairInfo_map;
/*
string read_revcomp(string s)
{
 string revstring="";
    if(s == "*") return s;
    for(int i=s.length() - 1;i>=0;i--)
    {
        char c = s[i];
    char revchar;
    switch (c) {
      case 'g':
        revchar = 'C';
        break;
      case 'G':
        revchar = 'C';
        break;
      case 'a':
        revchar = 'T';
        break;
      case 'A':
        revchar = 'T';
        break;
      case 't':
        revchar = 'A';
        break;
      case 'T':
        revchar = 'A';
        break;
      case 'c':
        revchar = 'G';
        break;
      case 'C':
        revchar = 'G';
        break;
      default:
        revchar = 'N';
    }
    revstring += revchar;
    }
    return revstring;
}
*/
void add_junction (int xs_flag,vector<int>& junc_l,vector<int>& junc_r, vector<double>& junc_cov,
		  vector<double>& junc_cov_plus,vector<double>& junc_cov_minus,
		  vector<int>& read_junc_l, vector<int>& read_junc_r,int coverage_used_NH 
		  )
{
    
       for(int i=0;i<read_junc_l.size();i++)
      {
         if(junc_l.size()==0)
         {
	      junc_l.push_back(read_junc_l[i]); junc_r.push_back(read_junc_r[i]);
	      double d=1.0000/(coverage_used_NH); junc_cov.push_back(d); 

	      if(xs_flag==1) {junc_cov_plus.push_back(d);junc_cov_minus.push_back(0);}//plus_cov
	      else if(xs_flag==0) {junc_cov_plus.push_back(0);junc_cov_minus.push_back(d);}
	      else { junc_cov_plus.push_back(d/2);junc_cov_minus.push_back(d/2); }
	      continue;
         }

         if(junc_l.size()>0){
          if(read_junc_l[i]==junc_l.back() && read_junc_r[i]==junc_r.back())
	  {
	    	  double d=1.00000/(coverage_used_NH);  junc_cov.back()=junc_cov.back()+d; 

		  if(xs_flag==1) {junc_cov_plus.back()=junc_cov_plus.back()+d;junc_cov_minus.back()=junc_cov_minus.back()+0;}//plus_cov
		  else if(xs_flag==0) {junc_cov_plus.back()=junc_cov_plus.back()+0;junc_cov_minus.back()=junc_cov_minus.back()+d;}
		  else { junc_cov_plus.back()=junc_cov_plus.back()+d/2; junc_cov_minus.back()=junc_cov_minus.back()+d/2;}
		  continue;
          }

          if(read_junc_l[i]>junc_l.back()||(read_junc_l[i]==junc_l.back()&&read_junc_r[i]>junc_r.back()))
	  {
	          junc_l.push_back(read_junc_l[i]);  junc_r.push_back(read_junc_r[i]);
	          double d=1.0000/(coverage_used_NH); junc_cov.push_back(d); 

	          if(xs_flag==1) {junc_cov_plus.push_back(d);junc_cov_minus.push_back(0);}//plus_cov
	          else if(xs_flag==0) {junc_cov_plus.push_back(0);junc_cov_minus.push_back(d);}
		  else { junc_cov_plus.push_back(d/2); junc_cov_minus.push_back(d/2);}
		  continue;
          }

          if(read_junc_l[i]<junc_l.back()||(read_junc_l[i]==junc_l.back()&&read_junc_r[i]<junc_r.back()))
	  {
	          for(int k=junc_l.size()-1;k>=0;k--){
		    if(read_junc_l[i]==junc_l[k]&&read_junc_r[i]==junc_r[k])
		    {
		      double d=1.00000/(coverage_used_NH);junc_cov[k]=junc_cov[k]+d; 

		      if(xs_flag==1) {junc_cov_plus[k]=junc_cov_plus[k]+d;junc_cov_minus[k]=junc_cov_minus[k]+0;}//plus_cov
		      else if(xs_flag==0) {junc_cov_plus[k]=junc_cov_plus[k]+0;junc_cov_minus[k]=junc_cov_minus[k]+d;}
		      else { junc_cov_plus[k]=junc_cov_plus[k]+d/2; junc_cov_minus[k]=junc_cov_minus[k]+d/2;}
		      break;
		    }
		    if(read_junc_l[i]>junc_l[k]||(read_junc_l[i]==junc_l[k]&&read_junc_r[i]>junc_r[k]))
		    {
		      junc_l.insert(junc_l.begin()+k+1,read_junc_l[i]);junc_r.insert(junc_r.begin()+k+1,read_junc_r[i]);
		      double d=1.00000/(coverage_used_NH);junc_cov.insert(junc_cov.begin()+k+1,d); 

		      if(xs_flag==1) {junc_cov_plus.insert(junc_cov_plus.begin()+k+1,d);junc_cov_minus.insert(junc_cov_minus.begin()+k+1,0);}//plus_cov
		      else if(xs_flag==0){junc_cov_plus.insert(junc_cov_plus.begin()+k+1,0);junc_cov_minus.insert(junc_cov_minus.begin()+k+1,d);}
		      else { junc_cov_plus.insert(junc_cov_plus.begin()+k+1,d/2); junc_cov_minus.insert(junc_cov_minus.begin()+k+1,d/2);}
		      break;
		    }
		    if(k==0)
		    {//
			junc_l.insert(junc_l.begin(),read_junc_l[i]); junc_r.insert(junc_r.begin(),read_junc_r[i]);

			double d=1.00000/(coverage_used_NH); junc_cov.insert(junc_cov.begin(),d);

			if(xs_flag==1) {junc_cov_plus.insert(junc_cov_plus.begin(),d);junc_cov_minus.insert(junc_cov_minus.begin(),0);}//plus_cov
                      	else if(xs_flag==0){junc_cov_plus.insert(junc_cov_plus.begin(),0);junc_cov_minus.insert(junc_cov_minus.begin(),d);}
			else { junc_cov_plus.insert(junc_cov_plus.begin(),d/2); junc_cov_minus.insert(junc_cov_minus.begin(),d/2);}

		    }//if(k==0)
	          }//for(int k=junc_l.size()-1;k>=0;k--)	
	    }//if(read_junc_l[i]<junc_l.back()||(read_junc_l[i]==junc_l.back()&&read_junc_r[i]<junc_r.back()))

         }//add new junction while exists some juntcion finish

    }
    return;
}

void delete_abnormal_junction(int strand,
			      vector<int>& junc_l,vector<int>& junc_r,
			      vector<double>& junc_cov, 
			      vector<double>& junc_cov_plus,vector<double>& junc_cov_minus,
			      bool& abnormal_flag)
{

    if(strand == 1) //minus
    {
	for(int i=0;i<junc_cov.size();i++)
	    junc_cov[i] = junc_cov[i] - junc_cov_plus[i];
    }
    else //plus
    {
	for(int i=0;i<junc_cov.size();i++)
	    junc_cov[i] = junc_cov[i] - junc_cov_minus[i];
    }
    double abnormal_count = 0;
    //for(int i=0;i<junc_cov.size();i++)
    for(int i=0;i<junc_cov.size();)//(EW1
    {
	if(junc_cov[i] <= 0)
	{
	    abnormal_count ++;
	    junc_l.erase(junc_l.begin() + i); junc_r.erase(junc_r.begin() + i);
	    junc_cov.erase(junc_cov.begin() + i);
	    junc_cov_plus.erase(junc_cov_plus.begin() + i);
	    junc_cov_minus.erase(junc_cov_minus.begin() + i);
	}
	else i++;//NEW1
    }

    if(double(abnormal_count/junc_cov.size()) >0.01) abnormal_flag = true;
    return;
}

void merge_and_get_exon(vector<int>& exon_l, vector<int>& exon_r, vector<double>& v_Gene,int gene_l)
{
    vector<int> v_nu;
    for(int a=0;a<v_Gene.size();a++)
    {
        if(v_Gene[a]==0) v_nu.push_back(1);
        else if(v_Gene[a]!=0)
        {
            if(v_nu.size()>0&&v_nu.size()<=MERGE_PARA) 
	    {
                for(int b=a-v_nu.size();b<a;b++){  v_Gene[b]++; }
             }
            v_nu.clear();
         }
    }
    for(int c=0;c<v_Gene.size();c++)
    {
	 if(c==0|| v_Gene[c-1]==0&&v_Gene[c]!=0) { exon_l.push_back(gene_l+c);}
	if(c==v_Gene.size()-1|| v_Gene[c]!=0&&v_Gene[c+1]==0){ exon_r.push_back(gene_l+c);}
    }

    return;
}
void junction_decide_exon(vector<int>& junc_l,vector<int>& junc_r,
			  vector<int>& exon_l,vector<int>& exon_r)
{
    
    if(junc_l.size()>0){//junction exists
	vector<int> junc_site_l,junc_site_r,junc_site;

	vector<int> junc_r_sort = junc_r;//NEW
	sort(junc_r_sort.begin(),junc_r_sort.end());//NEW
 	for(int k1=0;k1<junc_l.size();k1++)
 	{
	    for(int k2=0;k2<exon_l.size();k2++)
	    {
	        if(junc_l[k1]>=exon_l[k2]&&junc_l[k1]<=exon_r[k2]){
		    if(junc_site_l.size()==0) {junc_site_l.push_back(junc_l[k1]);}
		    if(junc_site_l.size()>0 && junc_site_l.back() != junc_l[k1]) {junc_site_l.push_back(junc_l[k1]);}
		}
		if(junc_r_sort[k1]>exon_l[k2]&&junc_r_sort[k1]<=exon_r[k2]){
		    if(junc_site_r.size()==0) junc_site_r.push_back(junc_r_sort[k1]);
		    if(junc_site_r.size()>0 && junc_site_r.back() != junc_r_sort[k1]) {junc_site_r.push_back(junc_r_sort[k1]);}
		}
	    }
	}
	if(junc_site_l.size()>0||junc_site_r.size()>0){//partial junction exists
	    for(int k1=0;k1<junc_site_r.size();k1++)
	    {
		junc_site_l.push_back(junc_site_r[k1]);
	    }
	    //junc_site=Sort(junc_site_l);
	    junc_site = junc_site_l;//NEW
	    sort(junc_site.begin(),junc_site.end());//NEW

	    vector<int> exon_lt,exon_rt;
	    exon_lt=exon_l; exon_rt=exon_r;
	    exon_l.clear();exon_r.clear();
	    vector<int> v_nu_exon;
	    v_nu_exon.push_back(0);

	    for(int i1=0;i1<junc_site.size();i1++)
	    {
		vector<int> v_in_site;
		for(int i2=v_nu_exon.back();i2<exon_lt.size();i2++)
		{
                    if(junc_site[i1]>exon_rt[i2])
		    {
       	                v_nu_exon.push_back(i2+1);
               	        exon_l.push_back(exon_lt[i2]);
                       	exon_r.push_back(exon_rt[i2]);
               	    }
		    if(junc_site[i1]>exon_lt[i2]&&junc_site[i1]<=exon_rt[i2])
		    {
		  	v_nu_exon.push_back(i2+1);
			for(int i3=i1;i3<junc_site.size();i3++){
			    if(junc_site[i3]>exon_lt[i2]&& junc_site[i3]<=exon_rt[i2]){ v_in_site.push_back(junc_site[i3]);}
			    if(junc_site[i3]>exon_rt[i2]){break;}
			}
			for(int m=0;m<v_in_site.size();m++){
			    if(m==0){ exon_l.push_back(exon_lt[i2]);exon_r.push_back(v_in_site[m]-1);}
			    if(m>0&&v_in_site[m]-1-v_in_site[m-1]>=0){ exon_l.push_back(v_in_site[m-1]); exon_r.push_back(v_in_site[m]-1);}
			    if(m==v_in_site.size()-1) {exon_l.push_back(v_in_site[m]); exon_r.push_back(exon_rt[i2]);}
			}
			break;
		    }
		}
	    }
	    for(int j1=v_nu_exon.back();j1<exon_lt.size();j1++){
	 	exon_l.push_back(exon_lt[j1]);
		exon_r.push_back(exon_rt[j1]);
	    }
	}//partial junction finish
    }//junction exist finish
}


void delete_exon_by_read_number(vector<int>& junc_l,vector<int>& junc_r, vector<double>& junc_cov,
				vector<double>& junc_cov_plus,vector<double>& junc_cov_minus,
				vector<int>& exon_l,vector<int>& exon_r,
				vector<int>& read_beg,vector<int>& read_fin,vector<int>& Read_nh
				)
{
 
    if(multi_map_frac == 1) return;   
    //delete exon depend on read number 
    vector<int> dele_exon_l,dele_exon_r;
    for(int i=0;i<exon_l.size();i++)
    {
        int read_num=0;
        int nh_num=0;
        for(int j=0;j<read_beg.size();j++)
	{
           if( (read_beg[j]>=exon_l[i] && read_beg[j]<=exon_r[i]) || (read_fin[j]>=exon_l[i] && read_fin[j]<=exon_r[i]) ){
                  if(Read_nh[j]>1) {
                       read_num++;nh_num++;
                   }
                   else read_num++;
           }
           if(read_beg[j]>exon_r[i]) break;
         }
         if((double)nh_num/read_num>multi_map_frac)
	 {
                    dele_exon_l.push_back(exon_l[i]);dele_exon_r.push_back(exon_r[i]);
                    exon_l[i]=0;  exon_r[i]=0;
         }
     }
	   
     for(int i=0;i<exon_l.size();)
     {
          if(exon_l[i]==0){exon_l.erase(exon_l.begin()+i);exon_r.erase(exon_r.begin()+i);}
          else i++;
      }

	    //since delete some exon, some junction need to delete
      for(int i=0;i<junc_l.size();)
      {
        if(In_exon(junc_l[i]-1,dele_exon_l,dele_exon_r) || In_exon(junc_r[i],dele_exon_l,dele_exon_r))
	{
             junc_l.erase(junc_l.begin()+i); junc_r.erase(junc_r.begin()+i);
             junc_cov.erase(junc_cov.begin()+i); 
	     junc_cov_plus.erase(junc_cov_plus.begin()+i); junc_cov_minus.erase(junc_cov_minus.begin()+i);
         }
         else i++;
      }
      return;
}

void process_graph_pair_strand( int max_map,int gene_l,
		  		vector<int>& junc_l,vector<int>& junc_r,vector<double>& junc_cov,
			  	vector<double>& junc_cov_plus, vector<double>& junc_cov_minus,
			  	vector<int>& exon_l,vector<int>& exon_r,vector<double>& exon_cov,vector<double>& v_Gene,
		  		vector<int>& seg_l,vector<int>& seg_r, vector<int>& seg_NH,
			  	vector<int>& read_beg, vector<int>& read_fin, vector<int>& Read_nh,
			  	string Chr,int strand,
				boost::unordered_map<string, pair<cigar_t,cigar_t> >& tu_readid_pairInfo_map
	 	   	      )
{
   if(max_map-gene_l<=gr_length ) return;
   bool abnormal_flag = false;
   delete_abnormal_junction(strand,junc_l,junc_r,junc_cov,junc_cov_plus,junc_cov_minus,abnormal_flag); 
   for(int i1=0;i1<junc_l.size();)
   {
		
                if(junc_cov[i1] <= 1 )//GR
                //if(junc_cov[i1] < 1 )//NEW CHANGE
                {
                    junc_l.erase(junc_l.begin()+i1);junc_r.erase(junc_r.begin()+i1);
                    junc_cov.erase(junc_cov.begin()+i1);
                    junc_cov_plus.erase(junc_cov_plus.begin()+i1);junc_cov_minus.erase(junc_cov_minus.begin()+i1);
                }

                else if(abnormal_flag && junc_cov[i1] == 1)
                {
                    junc_l.erase(junc_l.begin()+i1);junc_r.erase(junc_r.begin()+i1);
                    junc_cov.erase(junc_cov.begin()+i1);
                    junc_cov_plus.erase(junc_cov_plus.begin()+i1);junc_cov_minus.erase(junc_cov_minus.begin()+i1);
                }
                else i1++;
		
		
/*
		if(junc_cov[i1] <= 1 )//NEW star
		{
		    junc_l.erase(junc_l.begin()+i1);junc_r.erase(junc_r.begin()+i1);
		    junc_cov.erase(junc_cov.begin()+i1);
		    junc_cov_plus.erase(junc_cov_plus.begin()+i1);junc_cov_minus.erase(junc_cov_minus.begin()+i1);
		}
		else i1++;
*/
    }

    merge_and_get_exon(exon_l,exon_r,v_Gene,gene_l);

    delete_exon_by_read_number(junc_l,junc_r,junc_cov,junc_cov_plus,junc_cov_minus,exon_l,exon_r,read_beg,read_fin,Read_nh);

    junction_decide_exon(junc_l,junc_r,exon_l,exon_r);


    int line = 2;//strand 
    int XS_plus = 0,XS_minus = 0;
    process_gene(exon_l,exon_r,junc_l,junc_r,seg_l,seg_r,seg_NH,exon_cov,junc_cov,v_Gene,
		 strand,line,XS_plus,XS_minus,Chr,gene_l,tu_readid_pairInfo_map,abnormal_flag);

    return;
}

void process_graph_pair_unstrand(int max_map,int gene_l,
		  		vector<int>& junc_l,vector<int>& junc_r,vector<double>& junc_cov,
			  	vector<double>& junc_cov_plus, vector<double>& junc_cov_minus,
			  	vector<int>& exon_l,vector<int>& exon_r,vector<double>& exon_cov,vector<double>& v_Gene,
		  		vector<int>& seg_l,vector<int>& seg_r, vector<int>& seg_NH,
			  	vector<int>& read_beg, vector<int>& read_fin, vector<int>& Read_nh,
			  	string Chr,int strand,int XS_plus,int XS_minus,
				boost::unordered_map<string, pair<cigar_t,cigar_t> >& tu_readid_pairInfo_map
	 	   	      )
{
   if(max_map-gene_l<=gr_length ) return;
   //delete_abnormal_junction(strand,junc_l,junc_r,junc_cov,junc_cov_plus,junc_cov_minus,abnormal_flag);

   for(int i1=0;i1<junc_l.size();)
   {
                //if(junc_cov[i1] < 1 ) //GR
                if(junc_cov[i1] <= 1 )//NEWTTT //0 && 10.5 noref
		//if(junc_cov[i1] < 1 )//NEWTTT
                {
                    junc_l.erase(junc_l.begin()+i1);junc_r.erase(junc_r.begin()+i1);
                    junc_cov.erase(junc_cov.begin()+i1);
                    junc_cov_plus.erase(junc_cov_plus.begin()+i1);junc_cov_minus.erase(junc_cov_minus.begin()+i1);
                }
                else i1++;
     }

    merge_and_get_exon(exon_l,exon_r,v_Gene,gene_l);

    delete_exon_by_read_number(junc_l,junc_r,junc_cov,junc_cov_plus,junc_cov_minus,exon_l,exon_r,read_beg,read_fin,Read_nh);

    junction_decide_exon(junc_l,junc_r,exon_l,exon_r);

    int line = 1;//strand 

    if(XS_plus>0 && XS_minus>0)
    {
        vector<int>junc_plus_l,junc_plus_r,junc_minus_l,junc_minus_r;
        vector<int> exon_plus_l,exon_plus_r,exon_minus_l,exon_minus_r;
        vector<int> same_false_junc;
        vector<double> junc_cov_p,junc_cov_m;

        Graph_division graph_division;
        graph_division.division(exon_l,exon_r,junc_l,junc_r,junc_cov,junc_cov_plus,junc_cov_minus);

        junc_plus_l=graph_division.get_junc_plus_l();   junc_plus_r=graph_division.get_junc_plus_r();
        junc_minus_l=graph_division.get_junc_minus_l(); junc_minus_r=graph_division.get_junc_minus_r();

        exon_plus_l=graph_division.get_exon_plus_l();   exon_plus_r=graph_division.get_exon_plus_r();
        exon_minus_l=graph_division.get_exon_minus_l(); exon_minus_r=graph_division.get_exon_minus_r();

        same_false_junc=graph_division.get_same_false_junc();

        junc_cov_p=graph_division.get_junc_cov_plus();  junc_cov_m=graph_division.get_junc_cov_minus();

        double rate_1=(double)(XS_plus/(XS_plus+XS_minus));
        double rate_2=(double)(XS_minus/(XS_plus+XS_minus));

        XS_plus=1;XS_minus=0;
        //process_unstrand_gene(exon_plus_l,exon_plus_r,junc_plus_l,junc_plus_r,seg_l,seg_r,seg_NH,exon_cov,junc_cov_p,v_Gene,
        //                                 strand,line,XS_plus,XS_minus,ch_note,gene_l,rate_1,same_false_junc,tu_readid_pairInfo_map);

	//10.5-noref
        process_gene(exon_plus_l,exon_plus_r,junc_plus_l,junc_plus_r,seg_l,seg_r,seg_NH,exon_cov,junc_cov_p,v_Gene,
                                         strand,line,XS_plus,XS_minus,Chr,gene_l,tu_readid_pairInfo_map);
        XS_plus=0;XS_minus=1;
        //process_unstrand_gene(exon_minus_l,exon_minus_r,junc_minus_l,junc_minus_r,seg_l,seg_r,seg_NH,exon_cov,junc_cov_m,v_Gene,
        //                      strand,line,XS_plus,XS_minus,ch_note,gene_l,rate_2,same_false_junc,tu_readid_pairInfo_map);   

	//10.5-noref
        process_gene(exon_minus_l,exon_minus_r,junc_minus_l,junc_minus_r,seg_l,seg_r,seg_NH,exon_cov,junc_cov_m,v_Gene,
                              strand,line,XS_plus,XS_minus,Chr,gene_l,tu_readid_pairInfo_map);   
   }
   else //XS_plus == 0 || XS_minus == 0
	process_gene(exon_l,exon_r,junc_l,junc_r,seg_l,seg_r,seg_NH,exon_cov,junc_cov,v_Gene,
                     strand,line,XS_plus,XS_minus,Chr,gene_l,tu_readid_pairInfo_map);

    return;
}

void initiate_next_gene(vector<int>& read_junc_l, vector<int>& read_junc_r,vector<int>&read_seg_l,vector<int>&read_seg_r,
			string read_id,string read_map,int read_start, int read_last_site,int read_NH,string read_pair,
			int& xs_flag,
			int& max_map,int& gene_l,
                        vector<int>& junc_l,vector<int>& junc_r,vector<double>& junc_cov,
                        vector<double>& junc_cov_plus, vector<double>& junc_cov_minus,
                        vector<int>& exon_l,vector<int>& exon_r,vector<double>& exon_cov,vector<double>& v_Gene,
                         vector<int>& seg_l,vector<int>& seg_r, vector<int>& seg_NH,
                        vector<int>& read_beg, vector<int>& read_fin, vector<int>& Read_nh,
                        string Chr,int& strand,int& XS_plus,int& XS_minus,
			boost::unordered_map<string, pair<cigar_t,cigar_t> >& tu_readid_pairInfo_map
			)
{
    //readid_cigar_map.clear();
    //readid_pairInfo_map.clear();
    tu_readid_pairInfo_map.clear();
    //Gene_sequence.clear();
    XS_plus=0;XS_minus=0;

    read_beg.clear();read_fin.clear();Read_nh.clear();
    read_beg.push_back(read_start);read_fin.push_back(read_last_site);Read_nh.push_back(read_NH);

    exon_l.clear();exon_r.clear();junc_l.clear();junc_r.clear();junc_cov.clear(),seg_l.clear();seg_r.clear();exon_cov.clear();
    junc_cov_plus.clear();junc_cov_minus.clear(); seg_NH.clear();  v_Gene.clear();

    gene_l=read_start;
    max_map=read_last_site;

    for(int i=0;i<=max_map-gene_l;i++){
	 v_Gene.push_back(0);
	//Gene_sequence.push_back('a');
    }

    int coverage_used_NH = read_NH;
    //initiate seg
    for(int j=0;j<read_seg_l.size();j++)
    {
                for(int k=read_seg_l[j]-gene_l;k<=read_seg_r[j]-gene_l;k++)
                {
                    double d=1.00000/(coverage_used_NH); v_Gene[k]=v_Gene[k]+d;
                }
                seg_l.push_back(read_seg_l[j]); seg_r.push_back(read_seg_r[j]);
                seg_NH.push_back(coverage_used_NH);
   }

   if(read_pair == "=")
  {
       cigar_t vec_p;
       for(int j = 0;j<read_seg_l.size();j++)
       {
            pair<int,int> p = make_pair(read_seg_l[j],read_seg_r[j]);
            vec_p.push_back(p);
       }

      cigar_t temp;
      pair<cigar_t,cigar_t> p = make_pair(vec_p,temp);
      tu_readid_pairInfo_map[read_id] = p;
   }
   //initiate junction
   for(int i=0;i<read_junc_l.size();i++)
   {
                junc_l.push_back(read_junc_l[i]); junc_r.push_back(read_junc_r[i]);
                double d=1.00000/(coverage_used_NH); junc_cov.push_back(d);

                if(xs_flag == 1){junc_cov_plus.push_back(d);junc_cov_minus.push_back(0);}
                else if(xs_flag == 0) {junc_cov_plus.push_back(0);junc_cov_minus.push_back(d);}
                else { junc_cov_plus.push_back(d/2); junc_cov_minus.push_back(d/2);}
   }

}
#endif
