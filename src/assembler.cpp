// Some of this program was borrowed from Bayesembler and modified by Juntao Liu.
//
#include "Utility.h"
#include <assembler.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <list>
#include <vector>
#include <time.h>
#include <boost/random.hpp>
#include <set>
#include <map>
#include <utils.h>
#include <boost/program_options.hpp>
#include <api/BamWriter.h>
#include <limits>
using namespace std;
extern string suffix;
extern ofstream out_used_read;
extern ofstream out_junction_read;
extern ofstream out_graph;
extern map<pair<int,int>, vector<string> > plus_junction_readid_map;
extern map<pair<int,int>, vector<string> > minus_junction_readid_map;
extern map<pair<int,int>, vector<string> > unstranded_junction_readid_map;
extern string bam_file;
extern string strand_specific;
extern string out_dir;

extern string out_name;
extern int MERGE_PARA;
extern int anchor_length;

//used in single end data;
extern int min_delete_graph_cov;
extern int Path_length;
extern int single_strand_specific;

extern bool SingleEndFlag;

bool correcr_name_flag;
void getDir(string dir)
{
        bool retcode = boost::filesystem::create_directories(dir);

	string jr,ur;
	jr = dir + "/junction.reads." + suffix;

	ur = ur = dir + "/used.reads." + suffix;

	string gr;
	gr = dir + "/MyGraph.simplified." + suffix + ".graph";
	//out_junction_read.open(const_cast<char*>(jr.c_str()));
	//out_used_read.open(const_cast<char*>(ur.c_str()));
	out_graph.open(const_cast<char*>(gr.c_str()));
	return;
}
string reverse_(string s)
{
    std::reverse(s.begin(),s.end());
    return s;
}
string Assembler::read_revcomp(string s)
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
Assembler::Assembler(OptionsContainer options_variables_in) {

	options_variables = options_variables_in;
}
Assembler::Assembler()
{

}
bool Assembler::graphComparePos(GraphInfo first, GraphInfo second) {
    
    return first.left < second.left;
}


bool Assembler::graphCompareReadCount(pair<uint,uint> first, pair<uint,uint> second) {
    
    return first.second > second.second;
}

string Assembler::bamTosamString(BamTools::BamAlignment alignment)
{
    stringstream samString;
    samString << alignment.Name<<" ";
    samString << 0<<" ";
    samString << alignment.RefID<<" ";
    samString << alignment.Position+1<<" ";
    samString << alignment.MapQuality <<" ";
    samString << generateCigarString(alignment.CigarData)<<" ";
    samString << "= " ; //Pair
    samString << alignment.MatePosition<<" ";
    uint NH;
    alignment.GetTag("NH",NH);
    stringstream str; str<<NH;
    string s_NH = "NH:i:" + str.str();
    //cerr<<s_NH<<endl;

    string xs;
    if(alignment.GetTag("XS",xs))
    {
      	string s_xs = "XS:A:" + xs.substr(0,1);
	samString << s_xs<<" " <<s_NH;
    }
    else samString <<s_NH;

    vector<string> vec = alignment.GetTagNames();
    //for(int i=0;i<vec.size();i++) cerr<<vec[i]<<" ";cerr<<endl;
    
  
    return samString.str();
    
}
void Assembler::FindJunc(vector<BamTools::CigarOp> & cigar_data,
			 int read_start,
			 vector<int>& read_junc_l,vector<int>& read_junc_r,vector<int>& read_seg_l,vector<int>& read_seg_r,
			 int& last_site
			) 
{
    int num1 = 0, num2 = 0, sum = 0;
    for (vector<BamTools::CigarOp>::iterator it = cigar_data.begin(); it != cigar_data.end(); it++) {
        if ((it->Type == 'M') or (it->Type == 'D')) {
            //read_cigar << it->Length << "M";
	    num1 += it->Length;
	    num2 += it->Length;

        } else if (it->Type == 'N') {
	    num2 += it->Length;
	    read_junc_l.push_back(read_start + num1);
	    read_junc_r.push_back(read_start + num2);
	    num1 += it->Length;
            //read_cigar << it->Length << "N";

        } else continue;
/*
	else if (it->Type == 'S') {
            //read_cigar << it->Length << "S";
        }
        else if (it->Type == 'H') {
            //read_cigar << it->Length << "H";
        }
         else if (it->Type == 'I') {

        }
*/
   }
   last_site = num1;

   read_seg_l.push_back(read_start);
   if(read_junc_l.size() > 0)
   {
	read_seg_r.push_back(read_junc_l[0] - 1);
	for(int k = 1;k < read_junc_l.size();k++)
	{
	    read_seg_l.push_back(read_junc_r[k-1]); read_seg_r.push_back(read_junc_l[k]-1);
	}

	read_seg_l.push_back(read_junc_r.back()); read_seg_r.push_back(read_start + last_site-1);
   }
   else read_seg_r.push_back(read_start + last_site-1);

   last_site += read_start - 1;
}
void Assembler::add_junction_readid_map(string Strand, string id, vector<int> read_junc_l,vector<int> read_junc_r)
{
    if(Strand == "+" || Strand == ".")
    {
 	for(size_t i=0;i<read_junc_l.size();i++)
	{
	    pair<int,int> junc = make_pair(read_junc_l[i] - 1, read_junc_r[i]);
	    if(plus_junction_readid_map.find(junc) == plus_junction_readid_map.end()){
		vector<string>v(1,id);
	        plus_junction_readid_map[junc] = v;
	    }
	    else plus_junction_readid_map[junc].push_back(id);
	}
    }
    if(Strand == "-" || Strand == ".")
    {
 	for(size_t i=0;i<read_junc_l.size();i++)
	{
	    pair<int,int> junc = make_pair(read_junc_l[i] - 1, read_junc_r[i]);
	    if(minus_junction_readid_map.find(junc) == minus_junction_readid_map.end()){
		vector<string>v(1,id);
	        minus_junction_readid_map[junc] = v;
	    }
	    else minus_junction_readid_map[junc].push_back(id);
	}
    }
}
//void Assembler::run_TransRef(Transcomb_utility& TU, BamTools::BamAlignment& alignment, string Strand)
void Assembler::run_TransRef(Transcomb_utility& TU, BamYu& alignment, string Strand)
{

//   int read_start = alignment.Position + 1;
    int read_start = alignment.start_pos + 1;   
//   string read_id = alignment.Name, read_ch = RefVec[alignment.RefID].RefName;
    string read_id = alignment.read_id, read_ch = RefVec[alignment.ref_id].RefName;

   string read_map;// = generateCigarString(alignment.CigarData);
   string read_pair = "";
   //if(alignment.IsPaired()) read_pair = "=";
   if(alignment.paired) read_pair = "=";
   

   //uint NH;alignment.GetTag("NH",NH);
   uint NH = alignment.NH;
   int read_NH = (int)NH;
   int coverage_used_NH = read_NH;
   //if(FrameFlag) 
	   //out_used_read<<read_id<<'\n';//ipac2
   //if(read_NH > 1) return;

   string ID=read_id;
   //int coverage_used_NH = 1;

   int xs_flag = -1;
   /*
   string xs; 
   if( alignment.GetTag("XS",xs) )
   {
	if(xs.substr(0,1) == "+") { TU.XS_plus++; xs_flag = 1; }
	else { TU.XS_minus++; xs_flag = 0; }
   }
   */
   if(alignment.XS == "+") {TU.XS_plus++; xs_flag = 1;}
   else if(alignment.XS == "-") {TU.XS_minus++; xs_flag = 0;}

   //if(read_id != "end"){read_id.erase(read_id.end()-1);read_id.erase(read_id.end()-1);}


   vector<int> read_junc_l,read_junc_r,read_seg_l,read_seg_r;
   int read_last_site ;
   //FindJunc(alignment.CigarData ,read_start,read_junc_l,read_junc_r,read_seg_l,read_seg_r,read_last_site);
   FindJunc(alignment.cigar ,read_start,read_junc_l,read_junc_r,read_seg_l,read_seg_r,read_last_site);

   //if( !( read_id=="end" || (read_seg_r[0]-read_seg_l[0]+1>anchor_length && read_seg_r.back()-read_seg_l.back()+1>anchor_length) ) ) return;
   if( ! (read_seg_r[0]-read_seg_l[0]+1>anchor_length && read_seg_r.back()-read_seg_l.back()+1>anchor_length) )  return;


   if(TU.ch_note.size() == 0){
	TU.ch_note.push_back(read_ch);
	TU.Chr = read_ch; //9.9
	TU.gene_l = read_start;
	TU.max_map = read_last_site;
   }
   //cout<<read_ch<<" "<<TU.Chr<<" -- "<<read_start<<" "<< TU.max_map + MERGE_PARA <<'\n';
   //if(read_ch == TU.ch_note.back() && read_start <= TU.max_map + MERGE_PARA && read_id != "end")//current read in this gene!
   if(read_ch == TU.Chr && read_start <= TU.max_map + MERGE_PARA )//current read in this gene!
   {
	add_junction_readid_map(Strand,ID,read_junc_l,read_junc_r);//ipac2
	add_junction(xs_flag,TU.junc_l,TU.junc_r,TU.junc_cov,TU.junc_cov_plus,TU.junc_cov_minus,
		     read_junc_l,read_junc_r,coverage_used_NH);

	//TU.read_beg.push_back(read_start); TU.read_fin.push_back(read_last_site); TU.Read_nh.push_back(read_NH);

	if(read_last_site > TU.max_map) TU.max_map = read_last_site;
	
	for(int i = TU.v_Gene.size();i <= TU.max_map - TU.gene_l;i++){
	  TU.v_Gene.push_back(0);
	  //TU.Gene_sequence.push_back('a');
	}
	double d=1.00000/(coverage_used_NH);
	for(int j = 0;j<read_seg_l.size();j++){
	    for(int k = read_seg_l[j] - TU.gene_l;k <= read_seg_r[j] - TU.gene_l; k++) {
		//double d=1.00000/(coverage_used_NH); 
		TU.v_Gene[k] = TU.v_Gene[k] + d;
	    }
	    TU.seg_l.push_back(read_seg_l[j]); TU.seg_r.push_back(read_seg_r[j]);TU.seg_NH.push_back(coverage_used_NH);
	}

	//if(alignment.IsPaired()) 
	if(alignment.paired)
	{
	    cigar_t vec_p;
	    for(int j = 0;j<read_seg_l.size();j++){
		pair<int,int> p = make_pair(read_seg_l[j],read_seg_r[j]);
		vec_p.push_back(p);
	    }
	    boost::unordered_map<string, pair<cigar_t,cigar_t> >::iterator it = TU.tu_readid_pairInfo_map.find(read_id);
	    if(it == TU.tu_readid_pairInfo_map.end())
	    {
		cigar_t temp;
		pair<cigar_t,cigar_t> p = make_pair(vec_p,temp); //DIFFERENT !!!!!!!!
		TU.tu_readid_pairInfo_map[read_id] = p;
	    }
	    else it->second.second = vec_p;
	}
	else
	{
	    cigar_t vec_p;
	    for(int j = 0;j<read_seg_l.size();j++){
		pair<int,int> p = make_pair(read_seg_l[j],read_seg_r[j]);
		vec_p.push_back(p);
	    }
	    boost::unordered_map<string, pair<cigar_t,cigar_t> >::iterator it = TU.tu_readid_pairInfo_map.find(read_id);
	    if(it == TU.tu_readid_pairInfo_map.end())
	    {
		cigar_t temp;
		pair<cigar_t,cigar_t> p = make_pair(vec_p,vec_p); //DIFFERENT !!!!!!!!
		TU.tu_readid_pairInfo_map[read_id] = p;
	    }
	}

   }//current read in curent gene!

   else 
   {
	/*
	cerr<<"***"<<endl;
        for(int i=0;i<TU.junc_l.size();i++)
                cerr<<strand<<" "<<TU.ch_note.back()<<" "<<TU.junc_l[i]-1<<" "<<TU.junc_r[i]<<"  "
                <<TU.junc_cov[i]<<"  + "<<TU.junc_cov_plus[i]<<"  - "<<TU.junc_cov_minus[i]<<endl;
	*/
	

	if(line == 2)//strand 
	    process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
				      TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene, TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
				      TU.Chr, strand,
				      TU.tu_readid_pairInfo_map);
	else if(line == 1){
	    process_graph_pair_unstrand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
					TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
					TU.Chr, strand, TU.XS_plus,TU.XS_minus,
					TU.tu_readid_pairInfo_map);
	}
	initiate_next_gene(read_junc_l, read_junc_r,read_seg_l, read_seg_r,read_id, read_map, read_start,read_last_site,read_NH, read_pair,
			   xs_flag,
			   TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
			   TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,
			   TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
			   TU.Chr, strand, TU.XS_plus,TU.XS_minus,
			   TU.tu_readid_pairInfo_map);
	if(Strand == "+") plus_junction_readid_map.clear();
	else if(Strand == "-") minus_junction_readid_map.clear();
	add_junction_readid_map(Strand,ID,read_junc_l,read_junc_r);
   }
   TU.Chr = read_ch;
   //TU.ch_note.erase(TU.ch_note.begin());
   //TU.ch_note.push_back(read_ch);

   return;
}


double Assembler::writeUniqueReads(BamTools::BamWriter * writer, UniqueReadsYu * unique_reads, FirstReadsYu * first_reads) {
    double nd_nm_pe_reads = 0;

    uint cur_nh_tag;
    uint cur_hi_tag;

    UniqueReadsYu::iterator uit = unique_reads->begin();

    uint left_most_position;

    if (first_reads->empty()) {

        left_most_position = numeric_limits<uint>::max();
    
    } else {

        left_most_position = first_reads->begin()->first;
    }
         
    while (uit->first < left_most_position and !unique_reads->empty()) {
        cur_hi_tag = 0;
    
        //assert(uit->second->GetTag("NH", cur_nh_tag));  //NEW// 
        //assert(uit->second->GetTag("HI", cur_hi_tag) or (cur_nh_tag == 1)); //NEW//

        if (cur_nh_tag == 1) {

            nd_nm_pe_reads += 1;
        }
        //assert(writer->SaveAlignment(*(uit->second)));
	//cerr<<bamTosamString(*(uit->second))<<endl;
	if( line == 2)
	{
	    if(strand == 1)  run_TransRef(minus_info,*(uit->second),"-");
	    else  run_TransRef(plus_info,*(uit->second),"+");
	} else if(line == 1)
	{
	    run_TransRef(unstrand_info,*(uit->second),"."); 
	}
        delete uit->second;
        unique_reads->erase(uit++);

    }
         
    return nd_nm_pe_reads; 
}


//void Assembler::addUniqueReads(UniqueReads * unique_reads, ReadPairs * read_pairs) {
void Assembler::addUniqueReads(UniqueReadsYu * unique_reads, ReadPairsYu * read_pairs) {
    
    for (ReadPairsYu::iterator it = read_pairs->begin(); it != read_pairs->end(); it++) {

        //unique_reads->insert(unique_reads->end(), pair<uint, BamTools::BamAlignment*>(it->second.first->Position, it->second.first));
        //unique_reads->insert(unique_reads->end(), pair<uint, BamTools::BamAlignment*>(it->second.second->Position, it->second.second));

        unique_reads->insert(unique_reads->end(), pair<uint, BamYu*>(it->second.first->start_pos, it->second.first));
        unique_reads->insert(unique_reads->end(), pair<uint, BamYu*>(it->second.second->start_pos, it->second.second));
	/*
	BamTools::BamAlignment alignment1 = *(it->second.first);
	BamTools::BamAlignment alignment2 = *(it->second.second);
	bool junc_flag = false;
	vector<BamTools::CigarOp> cigar_data = alignment1.CigarData;
	for (vector<BamTools::CigarOp>::iterator it = cigar_data.begin(); it != cigar_data.end(); it++) {
	    if (it->Type == 'N'){
	 	junc_flag = true;
		break;
	    }
	}
	if(!junc_flag)
	{
	    vector<BamTools::CigarOp> cigar_data2 = alignment2.CigarData;
	    for (vector<BamTools::CigarOp>::iterator it = cigar_data2.begin(); it != cigar_data2.end(); it++) {
		if (it->Type == 'N'){
		    junc_flag = true;
		    break;
		}
	    }
	}
	//if(!junc_flag) continue;
	if(alignment1.IsFirstMate())
	{
	
	   //out1<<">"<<alignment1.Name<<'\n';
	   //out1<<alignment1.QueryBases<<'\n';
	   //if(alignment1.IsReverseStrand()) out1<<read_revcomp(alignment1.QueryBases)<<'\n';
	   //else out1<<alignment1.QueryBases<<'\n';

	   //out2<<">"<<alignment2.Name<<'\n';
	   //out2<<alignment2.QueryBases<<'\n';
	   //if(alignment2.IsReverseStrand()) out2<<read_revcomp(alignment2.QueryBases)<<'\n';
	   //else out2<<alignment2.QueryBases<<'\n';
	    
	}
	else 
	{
	   //out2<<">"<<alignment1.Name<<'\n';
	   //out2<<alignment1.QueryBases<<'\n';
	   //if(alignment1.IsReverseStrand()) out2<<read_revcomp(alignment1.QueryBases)<<'\n';
           //else out2<<alignment1.QueryBases<<'\n';

           //out1<<">"<<alignment2.Name<<'\n';
	   //out1<<alignment2.QueryBases<<'\n';
	   //if(alignment2.IsReverseStrand()) out1<<read_revcomp(alignment2.QueryBases)<<'\n';
           //else out1<<alignment2.QueryBases<<'\n';
	    
	}
	*/
    }

    read_pairs->clear();
}

string Assembler::generateCigarString(vector<BamTools::CigarOp> & cigar_data) {

    stringstream read_cigar; 
    for (vector<BamTools::CigarOp>::iterator it = cigar_data.begin(); it != cigar_data.end(); it++) {

        if ((it->Type == 'M') or (it->Type == 'D')) {

            read_cigar << it->Length << "M";

        } else if (it->Type == 'N') {

            read_cigar << it->Length << "N";
        
        } else if (it->Type == 'S') {
	    read_cigar << it->Length << "S";
	}
	else if (it->Type == 'H') {
	    read_cigar << it->Length << "H";
	}
	 else if (!(it->Type == 'I')) {

            cerr << "ERROR: Unhandled cigar string symbol '" << it->Type << "'!" << endl;
            exit(-1);
        }
   }
   return read_cigar.str();
}


//void Assembler::markDuplicates(BamTools::BamAlignment & current_alignment, FirstReads * first_reads, ReadPairs * read_pairs) {
void Assembler::markDuplicates(BamYu & current_alignment, FirstReadsYu * first_reads, ReadPairsYu * read_pairs) 
{

    uint cur_hi_tag = 0;
    //current_alignment.GetTag("HI", cur_hi_tag); //NEW//

    ReadId ri;
    //ri.name = current_alignment.Name;
    ri.name = current_alignment.read_id;
    if(correcr_name_flag) ri.name = ri.name.substr(0,ri.name.length() - 2); //NEW !!! IMPORTANT !!!

    ri.hi_tag = cur_hi_tag;

//    FirstReads::iterator first_reads_it = first_reads->find(current_alignment.MatePosition);
    FirstReadsYu::iterator first_reads_it = first_reads->find(current_alignment.mate_pos);

    if (first_reads_it == first_reads->end()) 
    {

//        FirstReads::iterator cur_pos_first_reads_it = first_reads->find(current_alignment.Position);
        FirstReadsYu::iterator cur_pos_first_reads_it = first_reads->find(current_alignment.start_pos);

        if (cur_pos_first_reads_it == first_reads->end()) 
	{

//          ReadIDs temp_read_ids;
            ReadIDsYu temp_read_ids;
//          assert(temp_read_ids.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);
            assert(temp_read_ids.insert(pair<ReadId,BamYu*>(ri, new BamYu(current_alignment))).second);
//          assert(first_reads->insert(pair<uint, ReadIDs>(current_alignment.Position, temp_read_ids)).second);
	    assert(first_reads->insert(pair<uint, ReadIDsYu>(current_alignment.start_pos, temp_read_ids)).second);

        } 
	else 
	{

//          assert(cur_pos_first_reads_it->second.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);
            assert(cur_pos_first_reads_it->second.insert(pair<ReadId, BamYu*>(ri, new BamYu(current_alignment))).second);
        }

    } 
    else 
    { 
        
//        ReadIDs::iterator pos_first_reads_it = first_reads_it->second.find(ri);
          ReadIDsYu::iterator pos_first_reads_it = first_reads_it->second.find(ri);

        if (pos_first_reads_it == first_reads_it->second.end()) 
	{

//            assert(current_alignment.Position == current_alignment.MatePosition);
            assert(current_alignment.start_pos == current_alignment.mate_pos);
//            assert(first_reads_it->second.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);            
            assert(first_reads_it->second.insert(pair<ReadId, BamYu*>(ri, new BamYu(current_alignment))).second);            
        
        } 
	else 
	{

            PairInfo pi;
            //pi.pos = pos_first_reads_it->second->Position;
            pi.pos = pos_first_reads_it->second->start_pos;
            //pi.insert = abs(pos_first_reads_it->second->InsertSize);
            pi.insert = abs(pos_first_reads_it->second->insert);

//            pi.cigar = generateCigarString(pos_first_reads_it->second->CigarData);
            pi.cigar = generateCigarString(pos_first_reads_it->second->cigar);
//            pi.m_cigar = generateCigarString(current_alignment.CigarData);
            pi.m_cigar = generateCigarString(current_alignment.cigar);

//            assert(current_alignment.MatePosition == pi.pos);
            assert(current_alignment.mate_pos == pi.pos);

            ReadPairsYu::iterator read_pairs_it = read_pairs->find(pi);

            if (read_pairs_it == read_pairs->end()) 
	    {

//                assert(read_pairs->insert(pair<PairInfo, ReadPair>(pi, ReadPair(pos_first_reads_it->second, new BamTools::BamAlignment(current_alignment)))).second);
                assert(read_pairs->insert(pair<PairInfo, ReadPairYu>(pi, ReadPairYu(pos_first_reads_it->second, new BamYu(current_alignment)))).second);
            
            } 
	    else 
	    {

                uint cur_first_nh_tag;
//               assert(pos_first_reads_it->second->GetTag("NH", cur_first_nh_tag));
		cur_first_nh_tag=pos_first_reads_it->second->NH;               

                uint cur_second_nh_tag;
//                assert(current_alignment.GetTag("NH", cur_second_nh_tag));
		cur_second_nh_tag=current_alignment.NH;

		uint cur_nh_tag = (cur_first_nh_tag < cur_second_nh_tag) ? cur_first_nh_tag:cur_second_nh_tag;//3 YU NEW in tophat.sam they are equal


                uint prev_nh_tag_1;
                //assert(read_pairs_it->second.first->GetTag("NH", prev_nh_tag_1));
		prev_nh_tag_1=read_pairs_it->second.first->NH;

		uint prev_nh_tag_2;
		//assert(read_pairs_it->second.second->GetTag("NH",prev_nh_tag_2));
		prev_nh_tag_2=read_pairs_it->second.second->NH;

		uint prev_nh_tag = (prev_nh_tag_1 < prev_nh_tag_2) ? prev_nh_tag_1:prev_nh_tag_2;//3 YU NEW in tophat.sam they are equal


		if(cur_nh_tag == 1 and prev_nh_tag > 1) 
		{//NEW 

                    delete read_pairs_it->second.first;
                    delete read_pairs_it->second.second;
                    
//                    read_pairs_it->second = ReadPair(pos_first_reads_it->second, new BamTools::BamAlignment(current_alignment));
                    read_pairs_it->second = ReadPairYu(pos_first_reads_it->second, new BamYu(current_alignment));
                
                } 
		else 
		{

                    delete pos_first_reads_it->second;
                }
            }
            
            first_reads_it->second.erase(pos_first_reads_it);

            if (first_reads_it->second.empty()) 
	    {

                first_reads->erase(first_reads_it);
            }
        }
	
    }
}
void  Assembler::ReadNameModify()
{
    BamTools::BamReader reader;
    reader.Open(bam_file);
    assert(reader.IsOpen());
    int i = 0;
    BamTools::BamAlignment current_alignment;
    while (reader.GetNextAlignment(current_alignment)) {
        i++;
        string name = current_alignment.Name;
        string s1 = name.substr(name.length() - 2);
        if(s1 == ".1" || s1 == ".2")
        {
            if(i>100)
            {
                correcr_name_flag = true;
                break;
            }
        }
        else{
         correcr_name_flag = false;
         break;
        }
    }
    reader.Close();
}


void Assembler::IfSingleEnd()
{
    BamTools::BamReader reader;
    reader.Open(bam_file);
    assert(reader.IsOpen());
    int i = 0;
    BamTools::BamAlignment current_alignment;
    while (reader.GetNextAlignment(current_alignment)) {
	if(current_alignment.IsPrimaryAlignment() ) 
	{
	    i++;
	    if(!(current_alignment.IsPaired())){ //not paired
	       SingleEndFlag = true;
	       break;
	    }
	    else {
	        if(i>100){
		  SingleEndFlag = false;
		  break;
		}
	   } 
	}

    }
    reader.Close();
}


void Assembler::generateSpliceGraphs()
{ 
    ReadNameModify();
    IfSingleEnd();
    //string ru = out_dir + "/reads.used";
    //if(FrameFlag) out_reads.open(const_cast<char*>(ru.c_str()));
    // Segment bamfile into positive and negative strands if stranded
    BamTools::BamReader reader;
    assert(reader.Open(bam_file));
    assert(reader.IsOpen());

    unsigned long total_mapped_pair_reads = 0;
    unsigned long nd_nm_mapped_paired_end_reads = 0;

    unsigned long num_reads_paired = 0;

    vector<pair<list<GraphInfo>, string> > all_graphs; 
    int graph_count = 0;

    const BamTools::SamHeader sam_header = reader.GetHeader();
    const BamTools::RefVector references = reader.GetReferenceData();

    RefVec = references;//NEW
    //cout<<"RefVec.size():"<<RefVec.size()<<endl;
    //for(size_t i=0;i<RefVec.size();i++) cout<<RefVec[i].RefName<<" ";
    //cout<<endl;
    /*
    for(size_t i=0;i<RefVec.size();i++){
	string dir = out_dir + "/" + RefVec[i].RefName;
        bool retcode = boost::filesystem::create_directories(dir);
    }
    */
    //cout << "[" << getLocalTime() << "] Removing duplicate reads ..." << endl; 

    if( SingleEndFlag )//single end 
    {
	string output = "MyReads.fa";
	//out3.open(output.c_str(),ios::app); //8.14NEW

	BamTools::BamWriter unstranded_writer;
	BamTools::BamAlignment current_alignment;
	int cur_reference = -1;
	if( single_strand_specific == 0) //unstranded
	{
	  line = 1 ; strand = 0;
	  while(reader.GetNextAlignment(current_alignment)) 
	  {
	      if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment())
	      {
	    	string Chr = RefVec[current_alignment.RefID].RefName;
		if(cur_reference == -1 )
		{
		    string dir = out_dir + "/" + Chr;
		    getDir(dir);
		    cur_reference = current_alignment.RefID;
		}
		else if(cur_reference != current_alignment.RefID){
		    cur_reference = current_alignment.RefID;
		    Transcomb_utility TU = unstrand_info;
		    process_graph_pair_unstrand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
                                              TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
                                              TU.Chr, strand, TU.XS_plus,TU.XS_minus,TU.tu_readid_pairInfo_map);
		    ClearTU(unstrand_info);
                    plus_junction_readid_map.clear();
                    minus_junction_readid_map.clear();
                    out_graph.flush();
                    out_graph.close();
		    string dir = out_dir + "/" + Chr;
		    getDir(dir);
		}
		BamYu current_alignment_yu(current_alignment);
		run_TransRef(unstrand_info,current_alignment_yu,".");
	      }
	  }
	} else if( single_strand_specific == 1) //forward
	{
	    line = 2;
	    while(reader.GetNextAlignment(current_alignment))
	    {
		if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment())
		{
		    BamYu current_alignment_yu(current_alignment);
                    string Chr = RefVec[current_alignment.RefID].RefName;
                    if(cur_reference == -1 )
                    {
                         string dir = out_dir + "/" + Chr;
                         getDir(dir);
                         cur_reference = current_alignment.RefID;
                    }
		    else if(cur_reference != current_alignment.RefID)
	            {
			    cur_reference = current_alignment.RefID;
			    Transcomb_utility TU = plus_info;
			    strand = 2;
			    process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
			                                TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
			                                TU.Chr, strand,TU.tu_readid_pairInfo_map);
			    TU = minus_info;
			    strand = 1;
			    process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
			                                TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
			                                TU.Chr, strand,TU.tu_readid_pairInfo_map);
			    ClearTU(plus_info);
			    ClearTU(minus_info);
			    plus_junction_readid_map.clear();
			    minus_junction_readid_map.clear();
			    out_graph.flush();
			    out_graph.close();
			    string dir = out_dir + "/" + Chr;
			    getDir(dir);
		    }
		    if( !current_alignment.IsReverseStrand() ) 
		    {
			strand = 2;
			run_TransRef(plus_info,current_alignment_yu,"+"); //plus 
		    }
		    else{
			strand = 1;
			run_TransRef(minus_info,current_alignment_yu,"-"); //minus
		    }
		}
	    }
	} else if(single_strand_specific == 2)//reverse
	{
	    line = 2;
	    while(reader.GetNextAlignment(current_alignment))
	    {
		if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment())
		{
		    BamYu current_alignment_yu(current_alignment);
                    string Chr = RefVec[current_alignment.RefID].RefName;
                    if(cur_reference == -1 )
                    {
                         string dir = out_dir + "/" + Chr;
                         getDir(dir);
                         cur_reference = current_alignment.RefID;
                    }
		    else if(cur_reference != current_alignment.RefID)
	            {
			    cur_reference = current_alignment.RefID;
			    Transcomb_utility TU = plus_info;
			    strand = 2;
			    process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
			                                TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
			                                TU.Chr, strand,TU.tu_readid_pairInfo_map);
			    TU = minus_info;
			    strand = 1;
			    process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
			                                TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
			                                TU.Chr, strand,TU.tu_readid_pairInfo_map);
			    ClearTU(plus_info);
			    ClearTU(minus_info);
			    plus_junction_readid_map.clear();
			    minus_junction_readid_map.clear();
			    out_graph.flush();
			    out_graph.close();
			    string dir = out_dir + "/" + Chr;
			    getDir(dir);
		    }
		    if( !current_alignment.IsReverseStrand() ) 
		    {
			strand = 1;
			run_TransRef(minus_info,current_alignment_yu,"-"); //minus
		    }
		    else
		    {
			strand = 2;
			run_TransRef(plus_info,current_alignment_yu,"+"); //plus
		    }
		}
	    }
	}
	else 
	{
	    cerr<<"[Error] : reads are single-end!"<<endl 
		<<"          -s argument need to  either \"single_unstranded\", \"single_forward\" or \"single_reverse\" for single-end reads"<<endl<<endl;
	    exit(1);
	}

	return;
    }
     
    string output1 = "MyReads_1.fa", output2 = "MyReads_2.fa";
    //out1.open(output1.c_str(),ios::app);//NEW //8.14NEW
    //out2.open(output2.c_str(),ios::app);//NEW//8.14NEW

    if(single_strand_specific != -1)
    {
	cerr<<"[Error] : reads are paired-end!"<<endl
	    <<"          -s argument need to be either \"unstranded\", \"first\" or \"second\" for parired-end reads"<<endl<<endl;
	exit(1);
    } 
    

    if (strand_specific == "unstranded") {        
	line = 1; strand = 0;
        
        BamTools::BamWriter unstranded_writer;
        //assert(unstranded_writer.Open(bam_nd_pe_unstranded_file_name, sam_header, references));

        BamTools::BamAlignment current_alignment;
        
        FirstReads first_reads;
        ReadPairs read_pairs;
        UniqueReads unique_reads;

        FirstReadsYu first_reads_Yu;
        ReadPairsYu read_pairs_Yu;
        UniqueReadsYu unique_reads_Yu;

        int current_position = -1;
        int cur_reference = -1;

        while (reader.GetNextAlignment(current_alignment)) {

            assert(current_alignment.IsPaired());
            assert(!current_alignment.IsFailedQC());

            if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment() and current_alignment.IsMateMapped()) {

                total_mapped_pair_reads += 1;

                // Check that read pair is OK
                if ((current_alignment.IsReverseStrand() != current_alignment.IsMateReverseStrand()) and (current_alignment.RefID == current_alignment.MateRefID)) {
		    BamYu current_alignment_yu(current_alignment);
                    // Same threshold as CEM
                    if (abs(current_alignment.InsertSize) > 700000) {

                        continue;
                    }
		    string Chr = RefVec[current_alignment.RefID].RefName;
		    if(cur_reference == -1 )
		    {
		        string dir = out_dir + "/" + Chr;
			getDir(dir);
			cur_reference = current_alignment.RefID;
		    }

		    if(cur_reference != current_alignment.RefID)
		    {
			    
		        cur_reference = current_alignment.RefID;
			current_position = current_alignment.Position;
			addUniqueReads(&unique_reads_Yu, &read_pairs_Yu);
			nd_nm_mapped_paired_end_reads += writeUniqueReads(&unstranded_writer, &unique_reads_Yu, &first_reads_Yu);
			Transcomb_utility TU = unstrand_info;

			process_graph_pair_unstrand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
                                        TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
                                        TU.Chr, strand, TU.XS_plus,TU.XS_minus,TU.tu_readid_pairInfo_map);
			ClearTU(unstrand_info);
			plus_junction_readid_map.clear();
			minus_junction_readid_map.clear();
			assert(unique_reads_Yu.empty());
			assert(first_reads_Yu.empty());
			out_graph.flush();
                        out_graph.close();
			string dir = out_dir + "/" + Chr;
			getDir(dir);
			
		    }
                    if (cur_reference == current_alignment.RefID && (current_position != current_alignment.Position) ){ //or (cur_reference != current_alignment.RefID)) 

                        num_reads_paired += read_pairs.size();         

                        //addUniqueReads(&unique_reads, &read_pairs);
			addUniqueReads(&unique_reads_Yu, &read_pairs_Yu);
                        //nd_nm_mapped_paired_end_reads += writeUniqueReads(&unstranded_writer, &unique_reads, &first_reads);
                        nd_nm_mapped_paired_end_reads += writeUniqueReads(&unstranded_writer, &unique_reads_Yu, &first_reads_Yu);

                        if (cur_reference != current_alignment.RefID) {
        
                            cur_reference = current_alignment.RefID;

                            assert(unique_reads_Yu.empty());
                            assert(first_reads_Yu.empty());

                        } else {

                            assert(current_position < current_alignment.Position);
                        }
                    }

                    //markDuplicates(current_alignment, &first_reads, &read_pairs);
		    markDuplicates(current_alignment_yu, &first_reads_Yu, &read_pairs_Yu);
                    current_position = current_alignment.Position;
                }
            }
        }

        num_reads_paired += read_pairs.size();

        addUniqueReads(&unique_reads_Yu, &read_pairs_Yu);
        nd_nm_mapped_paired_end_reads += writeUniqueReads(&unstranded_writer, &unique_reads_Yu, &first_reads_Yu);

	Transcomb_utility TU = unstrand_info;
	process_graph_pair_unstrand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
                                        TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene,TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
                                        TU.Chr, strand, TU.XS_plus,TU.XS_minus,TU.tu_readid_pairInfo_map);

        assert(unique_reads_Yu.empty());
        assert(first_reads_Yu.empty());

        //unstranded_writer.Close();

    }
  else 
  {//stranded

	line = 2;
            
        BamTools::BamWriter plus_writer;
        //assert(plus_writer.Open(bam_nd_pe_plus_file_name, sam_header, references));
        BamTools::BamWriter minus_writer;
        //assert(minus_writer.Open(bam_nd_pe_minus_file_name, sam_header, references));

        pair<BamTools::BamWriter *, BamTools::BamWriter *> writer_pair;
	bool strand_specific_first = false;

        if (strand_specific == "first") {

	    strand_specific_first = true;

            writer_pair.first = &plus_writer;
            writer_pair.second = &minus_writer;
        
        } else {
            
            writer_pair.first = &minus_writer;
            writer_pair.second = &plus_writer;   
        }

        BamTools::BamAlignment current_alignment;

        FirstReads first_reads_1;
        FirstReads first_reads_2;
	FirstReadsYu first_reads_1_Yu;
	FirstReadsYu first_reads_2_Yu;

        ReadPairs read_pairs_1;
        ReadPairs read_pairs_2;
	ReadPairsYu read_pairs_1_Yu;
	ReadPairsYu read_pairs_2_Yu;

        UniqueReads unique_reads_1;
        UniqueReads unique_reads_2;
        
	UniqueReadsYu unique_reads_1_Yu;
	UniqueReadsYu unique_reads_2_Yu;

        int current_position_1 = -1;
        int current_position_2 = -1;
        
        int cur_reference_1 = -1;
        int cur_reference_2 = -1;
	int cur_reference = -1;
        while (reader.GetNextAlignment(current_alignment)) {

            assert(current_alignment.IsPaired());
            assert(!current_alignment.IsFailedQC());

            if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment() and current_alignment.IsMateMapped()) {
                total_mapped_pair_reads += 1;

                // Check that read pair is OK
                if ((current_alignment.IsReverseStrand() != current_alignment.IsMateReverseStrand()) and (current_alignment.RefID == current_alignment.MateRefID)) {
                    // Same threshold as CEM
                    if (abs(current_alignment.InsertSize) > 700000) {

                        continue;
                    }
		    string Chr = RefVec[current_alignment.RefID].RefName;
		    BamYu current_alignment_yu(current_alignment);
		    //cout<<"YU: "<<current_alignment_yu.ref_id<<" "<<current_alignment_yu.start_pos<<'\n';
		    //if(cur_reference_1 != current_alignment.RefID)
		    if(cur_reference == -1 )
		    {
			string dir = out_dir + "/" + Chr;
			getDir(dir);

			cur_reference = current_alignment.RefID;
		    }
		    else if(cur_reference != current_alignment.RefID)
		    {
		        cur_reference = current_alignment.RefID;
			current_position_1 = -1;
			current_position_2 = -1;
			Transcomb_utility TU;

			if( strand_specific_first ) strand = 2;//plus
			else strand = 1;//minus
			addUniqueReads(&unique_reads_1_Yu, &read_pairs_1_Yu);
			nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.first, &unique_reads_1_Yu, &first_reads_1_Yu);
			if(strand == 2) TU = plus_info;
			else TU = minus_info;
			process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
                                      TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene, TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
                                      TU.Chr, strand,TU.tu_readid_pairInfo_map);
			if( strand_specific_first ) strand = 1;//minus
			else strand = 2;//plus
			addUniqueReads(&unique_reads_2_Yu, &read_pairs_2_Yu);
			nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.second, &unique_reads_2_Yu, &first_reads_2_Yu);
			if(strand == 2) TU = plus_info;
			else TU = minus_info;
			process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
                                      TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene, TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
                                      TU.Chr, strand,TU.tu_readid_pairInfo_map);

			assert(unique_reads_1_Yu.empty());
			assert(first_reads_1_Yu.empty());

			assert(unique_reads_2_Yu.empty());
			assert(first_reads_2_Yu.empty());

			plus_junction_readid_map.clear();
			minus_junction_readid_map.clear();
			ClearTU(plus_info);
			ClearTU(minus_info);

			out_graph.flush();
			out_graph.close();

			string dir = out_dir + "/" + Chr;
			getDir(dir);
		    }
                    // Plus strand
		    //cout<<"** "<<Chr<<" "<<current_alignment_yu.start_pos<<'\n';
                    if ( ( current_alignment.IsFirstMate() and current_alignment.IsReverseStrand() 
			   and (current_alignment.Position >= current_alignment.MatePosition) )  
			or ( !current_alignment.IsFirstMate() and !current_alignment.IsReverseStrand() 
			     and (current_alignment.Position <= current_alignment.MatePosition))
		       ) 
		    {
                        if (cur_reference == current_alignment.RefID and (current_position_1 != current_alignment.Position) ){//or (cur_reference_1 != current_alignment.RefID)) 

                            num_reads_paired += read_pairs_1.size(); 

                            //addUniqueReads(&unique_reads_1, &read_pairs_1);
			    addUniqueReads(&unique_reads_1_Yu, &read_pairs_1_Yu);

			    if( strand_specific_first ) strand = 2;//plus
			    else strand = 1;//minus

//                            nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.first, &unique_reads_1, &first_reads_1);
                            nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.first, &unique_reads_1_Yu, &first_reads_1_Yu);

                            assert(current_position_1 < current_alignment.Position);
                        }

                        //markDuplicates(current_alignment, &first_reads_1, &read_pairs_1); 
                        markDuplicates(current_alignment_yu, &first_reads_1_Yu, &read_pairs_1_Yu); 
			current_position_1 = current_alignment.Position;
                    }// !!! plus strand finish !!!
		    // Minus strand 
		    else if (  ( current_alignment.IsFirstMate() and !current_alignment.IsReverseStrand() 
			        and (current_alignment.Position <= current_alignment.MatePosition) ) 
			     or ( !current_alignment.IsFirstMate() and current_alignment.IsReverseStrand() 
				  and (current_alignment.Position >= current_alignment.MatePosition) )
			    ) 
		    {
                        if (cur_reference == current_alignment.RefID and (current_position_2 != current_alignment.Position) ){ //or (cur_reference_2 != current_alignment.RefID)) 
                            num_reads_paired += read_pairs_2.size(); 

                            //addUniqueReads(&unique_reads_2, &read_pairs_2);
			    addUniqueReads(&unique_reads_2_Yu, &read_pairs_2_Yu);

			    if( strand_specific_first ) strand = 1;//minus
			    else strand = 2;//plus
//                            nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.second, &unique_reads_2, &first_reads_2);
                            nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.second, &unique_reads_2_Yu, &first_reads_2_Yu);


			    assert(current_position_2 < current_alignment.Position);
                        }

                        //markDuplicates(current_alignment, &first_reads_2, &read_pairs_2); 
			markDuplicates(current_alignment_yu, &first_reads_2_Yu, &read_pairs_2_Yu);
                        current_position_2 = current_alignment.Position;           
                    }// !!! minus strand finish !!!
		    
                }
            }
        }

        num_reads_paired += read_pairs_1.size(); 

	//plus when strand specific first
        //addUniqueReads(&unique_reads_1, &read_pairs_1);
	addUniqueReads(&unique_reads_1_Yu, &read_pairs_1_Yu);

	if( strand_specific_first ) strand = 2;//plus
        else strand = 1;//minus
        nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.first, &unique_reads_1_Yu, &first_reads_1_Yu);//writer.first plus(specific first)  

	Transcomb_utility TU;
	if(strand == 1) TU = minus_info;
	else TU = plus_info;
	process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
                                      TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene, TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
                                      TU.Chr, strand,TU.tu_readid_pairInfo_map);

        assert(unique_reads_1_Yu.empty());
        assert(first_reads_1_Yu.empty()); 

        num_reads_paired += read_pairs_2.size(); 

	//minus when sstrand specific second
        //addUniqueReads(&unique_reads_2, &read_pairs_2);
	addUniqueReads(&unique_reads_2_Yu, &read_pairs_2_Yu);

	if( strand_specific_first ) strand = 1;//minus
	else strand = 2;//plus

        nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.second, &unique_reads_2_Yu, &first_reads_2_Yu);
	
	if(strand == 1) TU = minus_info;
	else TU = plus_info;

	process_graph_pair_strand(TU.max_map,TU.gene_l, TU.junc_l, TU.junc_r,TU.junc_cov,TU.junc_cov_plus, TU.junc_cov_minus,
                                      TU.exon_l,TU.exon_r, TU.exon_cov, TU.v_Gene, TU.seg_l, TU.seg_r, TU.seg_NH,TU.read_beg, TU.read_fin, TU.Read_nh,
                                      TU.Chr, strand,TU.tu_readid_pairInfo_map);


        assert(unique_reads_2_Yu.empty());
        assert(first_reads_2_Yu.empty());
	out_graph<<endl;
        out_graph.close();
        //plus_writer.Close();
        //minus_writer.Close();
  }
//  if(FrameFlag) out_reads.close();
}
