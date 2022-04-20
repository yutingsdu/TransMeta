#ifndef GTF2FASTA_H
#define GTF2FASTA_H

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include <algorithm>

using namespace std;
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
class GTF2FASTA
{
private:
  const char * path_gtf;
  const char * path_genome;
  const char * path_out;
public:
  vector<string> RefName;
  vector<string> RefSequence;
  ofstream out;

public:
  GTF2FASTA(char* f1,char* f2,char*f3)
  {
	path_gtf = f1; path_genome = f2; path_out = f3;
  }

  void load_genome()
  {

    ifstream in(path_genome);
    istringstream istr;//NEW


    string genome = "" ,s_temp;

    getline(in,s_temp);

    istr.str(s_temp);//NEW
    istr>>s_temp;//NEW
    istr.clear();//NEW
    RefName.push_back(s_temp.substr(1));

    while(getline(in,s_temp))
    {
	if(s_temp[0] == '>') 
	{
	    RefSequence.push_back(genome);
	    istr.str(s_temp);//NEW
	    istr>>s_temp;//NEW
	    istr.clear();//NEW
	    RefName.push_back(s_temp.substr(1));
	    genome.clear();
	    genome = "";
	}
	else genome.append(s_temp);
    }
    RefSequence.push_back(genome);

  }

  void get_genome_sequence(vector<int>& exon, string strand, string chr)
  {
	int refid = -1;
	for(int i=0;i<RefName.size();i++) 
	{
	   if(RefName[i] == chr) 
	   {
		refid = i;
		break;
	   }

	}
	if(strand == "+")
	{
	    for(int i=0;i<exon.size() - 1; )
	    {
		out<<RefSequence[refid].substr(exon[i] - 1,exon[i+1] - exon[i] + 1);
		i += 2;
	    }
	    out<<endl;
	}
	else
	{
	    string rev="";
	    for(int i=0;i<exon.size() - 1;)
	    {
		rev .append(RefSequence[refid].substr(exon[i] - 1,exon[i+1] - exon[i] + 1));
		i += 2;
	    }
	    out<<read_revcomp(rev);
	    out<<endl;
	}
  }
  void process()
  {
	load_genome();
	//cerr<<RefName.size()<<" "<<RefSequence.size()<<endl;

	out.open(path_out);

	ifstream in(path_gtf);
	string s;
	istringstream istr;

	vector< int > exon;
	
	string current_strand,current_chr,current_trid;
	getline(in,s);
	istr.str(s);
	if(s.empty())return;
	istr>>current_chr>>s>>s>>s>>s>>s>>current_strand>>s>>s>>s>>s>>current_trid;
	istr.clear();

	out<<">"<<current_trid.substr(1,current_trid.length() - 3)<<endl;
    	while(getline(in,s))
    	{
	    istr.str(s);
	    string lable,chr,strand,trid;
	    int n1,n2;
	    istr>>chr>>s>>lable>>n1>>n2>>s>>strand>>s>>s>>s>>s>>trid;
	    istr.clear();
	    if(lable == "transcript"){

		get_genome_sequence(exon,current_strand,current_chr);
		exon.clear();
		current_strand = strand;
		current_chr = chr;
	 	current_trid = trid;
		out<<">"<<current_trid.substr(1,current_trid.length() - 3)<<endl;
	    }

            else{
	        exon.push_back(n1); exon.push_back(n2);
            }
	}
	get_genome_sequence(exon,current_strand,current_chr);

	return;
  }

};

#endif
