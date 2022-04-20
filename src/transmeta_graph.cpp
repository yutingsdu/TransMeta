// Some of this program was borrowed from Bayesembler and modified by Juntao Liu.
#include <string>
#include <time.h>
#include <boost/random.hpp>
#include <utils.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> 
#include <assembler.h>
#include <unistd.h>
#include <getopt.h>
#include "pathsearch/gtf2fa.h"
#include "pathsearch/abundance.h"
using namespace std;

string suffix="";
map<pair<int,int>, vector<string> > plus_junction_readid_map;
map<pair<int,int>, vector<string> > minus_junction_readid_map;
map<pair<int,int>, vector<string> > unstranded_junction_readid_map;

ofstream out_junction_read;
ofstream out_used_read;
ofstream out_graph;

int pack_graph_num=0;
int unpack_graph_num=0;
bool PackingFlag = false;
bool MyFlag = false;
string bam_file ="";
string strand_specific = "unstranded";
string genome_file = "";

int single_strand_specific = -1;

bool SingleEndFlag = false;
bool SFlag = false;
bool SSFlag = false;
bool Help = false;

// *** parameter ***
int Path_length = 500;
double AVERAGE_REMOVE_RATE = 0.03;
double UNBALANCE_RATE = 0.03;
int gr_length = 200;
int trim = 0; //revise Borrow
int MERGE_PARA = 0;//200; //GR
float SEED_Filter = 1.01;
int anchor_length = 0;
float multi_map_frac = 1;

int min_delete_graph_cov = 2;// in single end = 3

string min_cov = "0";

bool version_flag = false;


int rg_index = -1; //then graph from 0
int trans_id;

string out_name = "";
string out_dir = "iPAC_outdir";

string usage();
static const char* short_options = "b:o:hs:g:l:L:d:D:e:f:F:a:m:vt:";
static struct option long_options[] = {
    {"bam-file",		1,  0,  'b'},
    {"out-dir", 		1,  0,  'o'},
    {"help",			0,  0,  'h'},
    {"strand-specific",		1,  0,  's'},
    {"genome-file",     	1,  0,  'g'},
    {"min-trans-length",	1,  0,  'l'},
    {"min-gene-length",		1,  0,  'L'},
    {"min-average-junc-ratio",  1,  0,  'd'},
    {"min-unbalance-ratio",	1,  0,  'D'},
    {"min-gap-length",		1,  0,  'e'},
    {"min-trans-cov",		1,  0,  'f'},
    {"min-seed-coverage",	1,  0,  'F'},
    {"min-anchor-length",	1,  0,  'a'},
    {"min-multiHit-reads-ratio",1,  0,  'm'},
    {"Single-Cell-Seq", 	0,  0,  310},
    {"trim",			1,  0,  't'},
    {"version",			0,  0,  'v'},
    {"single_end_ss",           1,  0,  320},
    {"suffix",			1,  0,  340},
    {0,0,0,0}
  };
int parse_options(int argc, char* argv[]) {
    int option_index = 0;
    int next_option;

    do {
	next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
	switch (next_option) {
	    case -1:     /* Done with options. */
		break;
	    case 'b':
		bam_file = optarg;
		break;
	    case 'o': //"out-dir"
		out_dir = optarg;
		break;
	    case 's':  //"strand-specific"
		strand_specific = optarg;
		break;
	    case 'g': //"genome-file"
		genome_file = optarg;
		break;
	    case 'l': // "min-trans-length" 
		Path_length = atoi(optarg); //300
		break;
	    case 'L': //"min-gene-length"
		gr_length = atoi(optarg); //200 
		break;
	    case 'd'://"min-average-junc-ratio"
		AVERAGE_REMOVE_RATE = atof(optarg);// 0.1 
		break;
	    case 'D': //"min-unbalance-ratio"
		UNBALANCE_RATE = atof(optarg); //0.1
		break;
	    case 'e': //"min-gap-length"
		MERGE_PARA = atoi(optarg); //200
		break;
	    case 'f': //"min-trans-cov"
		min_cov = optarg; //0
		break;
	    case 'F': //"min-seed-coverage"
		SEED_Filter = atof(optarg);//0
		break;
	    case 'a': //"min-anchor-length"
		anchor_length = atoi(optarg);//0
		break;
	    case 'm': //"min-multiHit-reads-ratio"
		multi_map_frac = atof(optarg);//1
		break;
	    case 't':// trim 
		trim = atoi(optarg);//1
		break;
	    case 'h':
		Help = true;
		break;
	    case 'v':
		version_flag = true;
		break;
	    case 320://
		single_strand_specific = atoi(optarg); //single_end_ss
		break;
	    case 340:
		suffix = optarg;;
		break;
	    default:
		cout<<usage();
		exit(1);
	}
    } while(next_option != -1);

    if(Help){
	//cout << usage() ;
	exit(1);
    }
    if(version_flag)
    {
        cout<<"The current version is: v.2.0"<<endl;
	exit(1);
    }
    if( bam_file == ""){
	cerr<<"Error : --bam-file option needs an argument!! "<<endl;
	//cout<<usage();
	exit(1);
    }
    /*
    if( genome_file == ""){
	cerr<<"Error : --genome-file option needs an argument!! "<<endl;
	//cout<<usage();
	exit(1);
    }
    */
    if(strand_specific != "unstranded" && strand_specific != "first" && strand_specific != "second"
	&& strand_specific != "single_unstranded" && strand_specific != "single_forward" && strand_specific != "single_reverse")
    {
	cerr<<"[Error] :  -s argument need to be either \"unstranded\", \"first\" or \"second\" for parired-end reads"<<endl
	    <<"            OR \"single_unstranded\", \"single_forward\" or \"single_reverse\" for single-end reads" << endl<<endl;
	exit(1);
    }
    if(strand_specific == "single_unstranded") single_strand_specific = 0;
    else if(strand_specific == "single_forward") single_strand_specific = 1;
    else if(strand_specific == "single_reverse" ) single_strand_specific = 2;
    return 0;
}
string usage()
{
    std::stringstream use_info;
    use_info<<endl
	    <<"================================================================"<<endl;
    use_info<<"ipac v.2.0 usage:"<<"\n"<<endl<<endl;
    use_info<<"** Required **"<<endl;

    use_info<<"--bam-file/-b <string>: BAM file produced by Tophat or Tophat2 or Hisat"<<endl<<endl;
    use_info<<"--genome-file/-g <string>: genome file"<<endl<<endl;

    use_info<<"--strand-specific/-s <string>: Strand-specific RNA-Seq reads orientation "<<endl
            <<"	      1) Use <unstranded> to indicate RNA-seq reads are non-strand-specific" <<endl
	    <<"	      2) Use <first> to indicate fr-first-stranded RNA-seq reads."<<endl
	    <<"	      3) Use <second> to indicate fr-second-stranded RNA-seq reads."<<endl<<endl;

    use_info<<"--Single-Cell-Seq: Requied if reads are from Single Cell Sequence."<<endl<<endl;
    use_info<<"----------------------------------------------------------------"<<endl;

    use_info<<endl<<"** Options **"<<endl<<endl
    <<"--help/-h				: Output ipac Help Information."<<endl<<endl
    <<"--out-dir/-o <string>			: Name of directory for output, default: iPAC_outdir/"<<endl<<endl
    <<"--min-trans-cov/-f <float>		: Minimum expression level estimated by abundance analysis for output, default: >0."<<endl<<endl
    <<"--min-seed-coverage/-F <float>		: Minimum seed coverage used for generate a new transcript."<<endl<<endl
    <<"--min-trans-length/-l <int> 		: Minimum assembled transcript length, default: 300."<<endl<<endl
    <<"--min-gene-length/-L <int> 		: Minimum  gene length, default: 200."<<endl<<endl
    <<"--min-average-junc-ratio/-d <float>	: Minimum junction coverage fraction by average junction coverage, default: 0.1."<<endl<<endl
    <<"--min-unbalance-ratio/-D <float>	: Minimum farction of unbalanced junction, default: 0.1."<<endl<<endl
    <<"--min-gap-length/-e <int>		: Minimum gap length between two exons, default: 200."<<endl<<endl
    <<"--min-anchor-length/-a <int>		: Minimum anchor length for junctions, default: 1."<<endl<<endl
    <<"--min-multiHit-reads-ratio/-m <float>	: Fraction of exon allowed to be covered by multi-hit reads, default: 1."<<endl<<endl
    <<"--version/-v  				: Report the current version of ipac and exit."<<endl<<endl;

    use_info<<"----------------------------------------------------------------"<<endl;
   use_info<<endl<<"** Note **"<<endl<<endl
   <<"A typical command of ipac might be:"<<endl
   <<"    ipac -b File.bam -g genome.fa -s unstranded"<<endl
   <<"    (If your data are strand-specific, it is recommended to set -s option.)"<<endl<<endl
   <<"If your data are from Single Cell Sequence,the command might be:"<<endl
   <<"    ipac -b File.bam -g genome.fa --Single-Cell-Seq"<<endl<<endl;

   use_info<<endl
            <<"================================================================"<<endl;
   return (use_info.str());
}

typedef  vector<pair<int,int> > cigar_t;



//MAIN 
int main (int argc, char * argv[]) {

    OptionsContainer options_variables;

    int parse_ret = parse_options(argc, argv);
    if (parse_ret)
      return parse_ret;
    if(trim == 1) SFlag = true;//spk
    else if(trim == 2){//SC
	SFlag = true;
	MERGE_PARA = 50;
    }
    if(strand_specific == "unstranded"){ //GR

    }
    //SFlag = false;//3.1 
    //if(SSFlag)
    //	MERGE_PARA = 50;
/*
    if( (access( out_dir.c_str(), 0 )) != -1 )
    {
	cout << "[Warning] "<<out_dir<<" exists already. It will be overwritten.\n";

	FILE* stream;

	string removeDir = "rm -rf " + out_dir;
	stream = popen(removeDir.c_str(),"r");
	pclose(stream);
    }

    bool retcode = boost::filesystem::create_directories(out_dir);
*/

    out_name ="./" + out_dir + "/raw." + suffix + ".gtf";
//cout<<out_name<<endl;
    Assembler assembler(options_variables);
    assembler.generateSpliceGraphs();
    out_junction_read.close();
    out_used_read.close();
    return 0;
}
