// Some of this program was borrowed from Bayesembler and modified by Juntao Liu.
#ifndef ASSEMBLER_H
#define ASSEMBLER_H

//#include "transcomb_main.h"
//#include "Utility.h"
#include <utils.h> 
//#include "./new/transcomb_main.h"
//#include <fragmentLengthModel.h>
#include <queue>
#include <boost/thread.hpp> 
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>
#include <api/BamWriter.h>
#include <map>
#include <tr1/unordered_set>


#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

using namespace std;

class Assembler {
       
    public:
        
        Assembler(OptionsContainer);
	Assembler();
//    	void runBayesemblerThreaded();
//	vector<pair<list<GraphInfo>, string> >  generateSpliceGraphs(double *);
	void generateSpliceGraphs();
	BamTools::RefVector RefVec;
	string generateCigarString(vector<BamTools::CigarOp> &);//keep
    private:
	class BamYu
	{
	    public:
		int start_pos,mate_pos,insert;
		string read_id;
		//string chr;
		uint ref_id;
		//string Cigar;
		vector<BamTools::CigarOp> cigar;
		string pair_flag;
		uint NH;
		string XS;
		bool paired;
	        BamYu(BamTools::BamAlignment& alignment)
		{
		    start_pos = alignment.Position;
		    mate_pos = alignment.MatePosition;
		    insert = alignment.InsertSize;
		    read_id = alignment.Name;
		    ref_id = alignment.RefID;
		    //chr = RefVec[alignment.RefID].RefName;
		    cigar= alignment.CigarData;
		    alignment.GetTag("NH",NH);

		    if(alignment.GetTag("XS",XS) )
			    XS = XS.substr(0,1);
		    else XS = "";
		    paired = alignment.IsPaired();
		    return;
		}
	};
        struct ReadId {
            
            string name;
            uint hi_tag;

            bool operator == (const ReadId & ri) const { 

                return (name == ri.name and hi_tag == ri.hi_tag);
            }
        };

        struct ReadIdHasher {
            
            size_t operator()(const ReadId & ri) const {

                return ((boost::hash<string>()(ri.name) ^ (boost::hash<uint>()(ri.hi_tag) << 1)) >> 1);
            }
        };

        struct PairInfo {
            
            uint pos;
            uint insert; // e_mpos
            string cigar;
            string m_cigar;

            bool operator == (const PairInfo & pi) const { 

                return (pos == pi.pos and insert == pi.insert and cigar == pi.cigar and m_cigar == pi.m_cigar);
            }
        };

        struct PairInfoHasher {
            
            size_t operator()(const PairInfo & pi) const {

                return (((((boost::hash<uint>()(pi.pos) ^ (boost::hash<uint>()(pi.insert) << 1)) >> 1) ^ (boost::hash<string>()(pi.cigar) >> 1)) << 1) ^ (boost::hash<string>()(pi.m_cigar) >> 1));
            }
        };

        struct SingleInfo {
            
            uint pos;
            string cigar;

            bool operator == (const SingleInfo & si) const { 

                return (pos == si.pos and cigar == si.cigar);
            }
        };

        struct SingleInfoHasher {
            
            size_t operator()(const SingleInfo & si) const {

                return ((boost::hash<uint>()(si.pos) ^ (boost::hash<string>()(si.cigar) << 1)) >> 1);
            }
        };

        typedef pair<BamTools::BamAlignment*, BamTools::BamAlignment*> ReadPair;
	typedef pair<BamYu*,BamYu*> ReadPairYu;
        typedef tr1::unordered_map <ReadId, BamTools::BamAlignment*, ReadIdHasher> ReadIDs;
	typedef tr1::unordered_map <ReadId, BamYu*, ReadIdHasher> ReadIDsYu;

        typedef map <uint, tr1::unordered_map <ReadId, BamTools::BamAlignment*, ReadIdHasher> > FirstReads;
	typedef map <uint, tr1::unordered_map <ReadId, BamYu*, ReadIdHasher> > FirstReadsYu;

        typedef tr1::unordered_map <PairInfo, ReadPair, PairInfoHasher > ReadPairs;     
	typedef tr1::unordered_map <PairInfo, ReadPairYu, PairInfoHasher > ReadPairsYu;

        typedef multimap <uint, BamTools::BamAlignment*> UniqueReads;
	typedef multimap <uint, BamYu*> UniqueReadsYu;

        static bool graphComparePos(GraphInfo first, GraphInfo second);//keep
        static bool graphCompareReadCount(pair<uint,uint> first, pair<uint,uint> second);//keep

        //void generateBamIndex(string);

        //double writeUniqueReads(BamTools::BamWriter *, UniqueReads *, FirstReads *);//keep
        double writeUniqueReads(BamTools::BamWriter *, UniqueReadsYu *, FirstReadsYu *);//keep

//        void addUniqueReads(UniqueReads *, ReadPairs *);//keep
	void addUniqueReads(UniqueReadsYu *, ReadPairsYu *);//keep
	void add_current_alinment_as_UniqueReads(UniqueReads * ,BamTools::BamAlignment & );//YUTING NEW
	
//    	string generateCigarString(vector<BamTools::CigarOp> &);//keep
//        void markDuplicates(BamTools::BamAlignment &, FirstReads *, ReadPairs *);//keep
	void markDuplicates(BamYu &, FirstReadsYu *, ReadPairsYu *);        
	

	ofstream out1,out2,out3;

    	OptionsContainer options_variables;

        struct ReadInfo {
            
            uint position;
            int fragment_length;
            vector<BamTools::CigarOp> cigar_string;
            double probability;
        };

	//TransRef
	
	void ReadNameModify();
	void IfSingleEnd();
	typedef  vector<pair<int,int> > cigar_t;
	
	void FindJunc(vector<BamTools::CigarOp> & cigar_data, int read_start,
                        vector<int>& read_junc_l,vector<int>& read_junc_r,
			vector<int>& read_seg_l,vector<int>& read_seg_r,
	                int& last_site);


	string bamTosamString(BamTools::BamAlignment alignment);

	struct Transcomb_utility
	{
	  vector<int> exon_l,exon_r,junc_l,junc_r,seg_l,seg_r,seg_NH;// **seg_l, seg_r used to compute partial_exon_junction coverage;**
	  vector<double> exon_cov,junc_cov,v_Gene,junc_cov_plus,junc_cov_minus;
	  vector<int> read_beg,read_fin,Read_nh;// ** used in delete exon depend on read number **
	  int XS_plus,XS_minus;
	  int max_map;
	  int gene_l;
	  vector<string> ch_note;
	  string Chr;
	  //vector<int> strand_note;
	  //typedef multimap <string,vector<string> >::iterator iter;
	  //istringstream istr;
	  boost::unordered_map<string, pair<cigar_t,cigar_t> > tu_readid_pairInfo_map;
	  vector<char> Gene_sequence;
	  
	  

	  Transcomb_utility& operator = (Transcomb_utility& TU) {
	    exon_l= TU.exon_l; exon_r = TU.exon_r; junc_l = TU.junc_l; junc_r = TU.junc_r;
	    seg_l = TU.seg_l; seg_r = TU.seg_r; seg_NH = TU.seg_NH;

	    exon_cov = TU.exon_cov; junc_cov = TU.junc_cov; v_Gene = TU.v_Gene;
	    junc_cov_plus = TU.junc_cov_plus; junc_cov_minus = TU.junc_cov_minus;

	    read_beg = TU.read_beg; read_fin = TU.read_fin; Read_nh = TU.Read_nh;

	    XS_plus = TU.XS_plus; XS_minus = TU.XS_minus;
	    max_map = TU.max_map; gene_l = TU.gene_l; ch_note = TU.ch_note; 
	    Chr = TU.Chr;
	    //strand_note = TU.strand_note;

	    tu_readid_pairInfo_map = TU.tu_readid_pairInfo_map;
	    //Gene_sequence = TU.Gene_sequence;
	  }

	};

	int line,strand;

	Transcomb_utility plus_info,minus_info,unstrand_info;
	//void run_TransRef(Transcomb_utility& , BamTools::BamAlignment&,string);
	void run_TransRef(Transcomb_utility& , BamYu&,string);
	string read_revcomp(string s);

	void ClearTU(Transcomb_utility& TU){
	    TU.exon_l.clear(); TU.exon_r.clear();TU.junc_l.clear();TU.junc_r.clear();
	    TU.exon_cov.clear(); TU.junc_cov.clear(); TU.v_Gene.clear(); 
	    TU.junc_cov_plus.clear(); TU.junc_cov_minus.clear();

	    TU.seg_l.clear(); TU.seg_r.clear(); TU.seg_NH.clear();
	    TU.read_beg.clear();TU.read_fin.clear(); TU.Read_nh.clear();
	    TU.ch_note.clear();
	    TU.tu_readid_pairInfo_map.clear();
	    TU.Gene_sequence.clear();
	    return;
	}
	void add_junction_readid_map(string Strand, string id, vector<int> read_junc_l,vector<int> read_junc_r);//ipac2
/*
	vector<int> exon_l,exon_r,junc_l,junc_r,seg_l,seg_r,seg_NH;// **seg_l, seg_r used to compute partial_exon_junction coverage;**
        vector<double> exon_cov,junc_cov,v_Gene,junc_cov_plus,junc_cov_minus;
        vector<int> read_beg,read_fin,Read_nh;// ** used in delete exon depend on read number **

        int XS_plus,XS_minus;
        int max_map;
        int gene_l;
        vector<string> ch_note;
        vector<int> strand_note;
        int line,strand;

        multimap <string,vector<string> > m_id_map;//multimap <read_id, start+cigar>

        typedef multimap <string,vector<string> >::iterator iter;

        istringstream istr;
*/
};

#endif
    
    

