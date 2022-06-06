#include "./path_search.h"
#include "boost/icl/interval_map.hpp"
#include "boost/icl/split_interval_map.hpp"
using namespace boost;
using namespace std;

extern ofstream out_gtf;
int MERGE_PARA = 200;
string out_name;
ofstream out_graph;
string mode;
extern bool unstranded;
int pack_graph_num=0;
int unpack_graph_num=0;
bool SFlag = false;
bool PackingFlag = false;
bool MyFlag = false;
int Path_length = 500;
double AVERAGE_REMOVE_RATE = 0.01;//0.03
double UNBALANCE_RATE = 0.03;
double SEED_Filter = 1.01;
int rg_index = 0;
int trans_id;
typedef pair<int,int> junc_t;
typedef pair<int,int> exon_t;
double SampleSize = 1;
struct junc_info{
    string sample;
    double coverage;
};
class RawGraph
{
  //private:
  public:
    string Mode;//R|I|U 
    string strand, chr;
    int start_pos, end_pos;
    vector<double> VecGene;
    vector<pair<int,int> > junction;
    vector<double> junc_cov;

    vector<int> exon_l,exon_r,junc_l,junc_r;
    vector<pair<int,int> > exon,raw_ref_exon;
    vector<double> exon_cov, raw_ref_exon_cov;

    vector<path_t> pair_path;
    vector<double> pair_path_cov;
    vector< vector<int> > pair_path_exon;
  public:
    RawGraph(){}
    RawGraph(vector<pair<int,int> >junction_, vector<double> junc_cov_, 
	     int sp, vector<double> vg,
	     string str, string ch, string mode,
	     vector<path_t> pair_path_, vector<double> pair_path_cov_,
	     vector<pair<int,int> > raw_exon_, vector<double> raw_exon_cov_)
    {
	junction = junction_;
	junc_cov = junc_cov_;
  	start_pos = sp;
	VecGene = vg;
	end_pos = sp + vg.size() - 1;
	strand = str;
	chr = ch;
	Mode = mode;
  	if(1 || Mode == "R") 
	{
	    raw_ref_exon = raw_exon_;
	    raw_ref_exon_cov = raw_exon_cov_;
	}
  	exon = raw_exon_;
	exon_cov = raw_exon_cov_;

	pair_path = pair_path_;
	pair_path_cov = pair_path_cov_;
	for(size_t i=0;i<pair_path.size();i++)
	{
	    vector<int> path_exon;
	    path_t p = pair_path[i];
	    for(size_t j=0;j<p.size();j++)
	    {
	         path_exon.push_back(exon[p[j]].first);
		 path_exon.push_back(exon[p[j]].second);
	    }
	    pair_path_exon.push_back(path_exon);
	}
			
    }
  void show()
  {
    cout<<"pos: "<<start_pos<<" "<<end_pos<<endl;
    cout<<"Nodes"<<endl;
    for(int i=0;i<exon_l.size();i++) cout<<exon_l[i]<<" "<<exon_r[i]<<": "<<exon_cov[i]<<endl;
  }
  void get_exon_and_junction()
  {
	exon_l.clear(); exon_r.clear(); 
	junc_l.clear(); junc_r.clear();
	for(size_t i=0;i<exon.size();i++)
  	{
	    exon_l.push_back(exon[i].first);
            exon_r.push_back(exon[i].second);
	}
	for(size_t i=0;i<junction.size();i++)
        {
            junc_l.push_back(junction[i].first + 1);
            junc_r.push_back(junction[i].second);
        }
	return;
  }
  void merge_and_get_exon()
  {
    exon_l.clear(); exon_r.clear(); exon.clear();
    exon_cov.clear();
    junc_l.clear(); junc_r.clear();
    vector<double> v_Gene = VecGene;
    int gene_l = start_pos;
    vector<int> v_nu;
    for(int a=0;a<v_Gene.size();a++)
    {
        if(v_Gene[a]==0) v_nu.push_back(1);
        else if(v_Gene[a]!=0)
        {
            if(v_nu.size()>0&&v_nu.size()<= 50) //MERGE_PARA) 
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
    for(size_t i=0;i<exon_l.size();i++) exon.push_back(make_pair(exon_l[i], exon_r[i]));
    return;
  }

  void junction_decide_exon()
  {
    for(size_t i=0;i<junction.size();i++)
    {
         junc_l.push_back(junction[i].first + 1);
         junc_r.push_back(junction[i].second);
    }
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
    for(int i1=0;i1<exon_l.size();i1++)
    {
          double  k=0.00000;
          for(int i2=exon_l[i1]-start_pos;i2<=exon_r[i1]-start_pos;i2++){
              k=k+VecGene[i2];
           }
           double cov=k/(exon_r[i1]-exon_l[i1]+1);
           exon_cov.push_back(cov);
    }
}
    void Clear()
    {
	start_pos = -1;
	end_pos = -1;
	VecGene.clear();
	junction.clear(); junc_cov.clear();
	exon.clear(); exon_cov.clear();
	raw_ref_exon.clear();
	raw_ref_exon_cov.clear();
	pair_path.clear();
	return;
	
    }
};

map<string, vector<RawGraph> > plus_ref, minus_ref, unstranded_ref;
map<string,bool> skipped_reads;
void add_refinfo(string strd, string chr, RawGraph rg, map<string, vector<RawGraph> >& ref)
{
    string s = strd+chr;
    if(ref.find(s) == ref.end())
    {
	vector<RawGraph> V(1,rg);
	ref[s] = V;
    }
    else {
	ref[s].push_back(rg);
    }
    return;
}
void load_unmapped_reads(char* file)
{
    ifstream in(file);
    string s;
    int i=0;
    while(getline(in,s))
    {
	if(i % 1000000 == 0) cerr<<"loading "<<i<<" lines..."<<endl;
	i++;
	if(s[0] != '@') continue;
	s = s.substr(1,s.length() - 3);
	skipped_reads[s] = true;
    }
    cout<<"skipped_reads size: "<<skipped_reads.size()<<endl;
    return;
}
bool sorter(const RawGraph R1, const RawGraph R2)
{
    return (R1.start_pos<R2.start_pos
	    || (R1.start_pos==R2.start_pos && R1.end_pos<R2.end_pos)
	    || (R1.start_pos==R2.start_pos && R1.end_pos==R2.end_pos && R1.Mode < R2.Mode));
}

void combine(RawGraph& G1, RawGraph& G2)
{
    int shift = G2.start_pos - G1.start_pos;
    if(G2.end_pos > G1.end_pos)
    {
	for(size_t i=0;i<G2.end_pos - G1.end_pos;i++)
	    G1.VecGene.push_back(0);
	G1.end_pos = G2.end_pos;
    }
    for(size_t i=0;i<G2.VecGene.size();i++)
    {
	G1.VecGene[shift + i] += G2.VecGene[i];
    }

    if(1 || G2.Mode == "R")
    {
	for(size_t i=0;i<G2.raw_ref_exon.size();i++)
	{
	    G1.raw_ref_exon.push_back(G2.raw_ref_exon[i]);
	    G1.raw_ref_exon_cov.push_back(G2.raw_ref_exon_cov[i]);
	}
    }

    vector<junc_t> Juncs = G1.junction;
    for(size_t i=0;i<G2.junction.size();i++) Juncs.push_back(G2.junction[i]);

    sort(Juncs.begin(),Juncs.end());
    Juncs.erase(unique(Juncs.begin(),Juncs.end()), Juncs.end());

    vector<double> Covs;
    for(size_t i=0;i<Juncs.size();i++)
    {
	double cov = 0;
	vector<junc_t>::iterator it = find(G1.junction.begin(),G1.junction.end(),Juncs[i]);
	if( it != G1.junction.end()) cov += G1.junc_cov[ it - G1.junction.begin() ];
	it = find(G2.junction.begin(),G2.junction.end(),Juncs[i]);
	if( it != G2.junction.end()) cov += G2.junc_cov[ it - G2.junction.begin() ];

	Covs.push_back(cov);
    }
    vector<junc_t> Exons = G1.exon;
    for(size_t i=0;i<G2.exon.size();i++) Exons.push_back(G2.exon[i]);
    sort(Exons.begin(),Exons.end());
    Exons.erase(unique(Exons.begin(),Exons.end()), Exons.end());
    vector<double> ECovs;
    for(size_t i=0;i<Exons.size();i++)
    {
	double cov = 0;
	vector<junc_t>::iterator it1 = find(G1.exon.begin(),G1.exon.end(),Exons[i]);
	if( it1 != G1.exon.end()) cov += G1.exon_cov[ it1 - G1.exon.begin() ];
	vector<junc_t>::iterator it2 = find(G2.exon.begin(),G2.exon.end(),Exons[i]);
	if( it2 != G2.exon.end()) cov += G2.exon_cov[ it2 - G2.exon.begin() ];
	if(it1 != G1.exon.end() && it2 != G2.exon.end()) cov = double(cov/2.0);
	ECovs.push_back(cov);
    }


    vector<vector<int> > PairPathExon = G1.pair_path_exon;
    for(size_t i=0;i<G2.pair_path_exon.size();i++) PairPathExon.push_back(G2.pair_path_exon[i]);
    sort(PairPathExon.begin(),PairPathExon.end());
    PairPathExon.erase(unique(PairPathExon.begin(),PairPathExon.end()),PairPathExon.end());

    G1.junction = Juncs;
    G1.junc_cov = Covs;
    G1.exon = Exons;
    G1.exon_cov = ECovs;
    G1.pair_path_exon = PairPathExon;

    //if(G1.Mode != "U" && (G1.Mode == "R" || G2.Mode == "R"))
    if(G1.Mode != "U" && (G1.Mode  != G2.Mode))
        G1.Mode = "U";
	
    G2.Clear();
    return;
    
}
void get_pair( vector<vector<int> > PairPathExon,vector<path_t>& Unused_pair_paths,vector<double> & Unused_pair_paths_cov,
		vector<int> exon_l, vector<int> exon_r,
		vector<int> junc_l, vector<int> junc_r, vector<double> Vec_weights)
{
    //cout<<"PairPathExon size: "<<PairPathExon.size()<<endl;
    vector<pair<int,int> > exon,junc;
    for(int i=0;i<exon_l.size();i++)
    {   
        pair<int,int> p = make_pair(exon_l[i],exon_r[i]);
        exon.push_back(p);
    }   
    for(int i = 0;i<junc_l.size();i++)
    {   
        pair<int,int> p = make_pair(junc_l[i],junc_r[i]);
        junc.push_back(p);
    }
    for(size_t i=0;i<PairPathExon.size();i++)
    {
	vector<int> exons = PairPathExon[i];
	path_t p;
	if(exons.size() < 6 ) continue;
	for(size_t i1=0;i1<exons.size();)
	{
	    int el = exons[i1], er = exons[i1+1];
	    //for each exon in PAIR, Do: 

	    for(size_t j=0;j<exon_l.size();j++)
	    {
	    	if(exon_l[j] <= el && exon_r[j] > el)
	    	{
		  for(size_t k=j;k<exon_l.size();k++)
		  {
		     if(exon_r[k]<=er) p.push_back(k); 
		     else break;
		  }
		  break;
	    	}
	    }
	    i1 += 2;
	}
	bool nojunc_flag = false;
	double cov = 10000000;
	//double cov = 0;
	for(size_t j=0;j<p.size()-1;j++)
	{
	    junc_t junc_ = make_pair(exon_r[p[j]] + 1, exon_l[p[j+1]]);
	    vector<pair<int,int> >::iterator it = find(junc.begin(),junc.end(),junc_);
	    if(it == junc.end())
	    {
	        nojunc_flag = true;
		break;
	    }
	    int J = it - junc.begin();
	 //   cout<<p[j]<<" "<<p[j+1]<<" "<<junc_.first<<" "<<junc_.second<<" J: "<<J<<" cov: "<<cov<<endl;
	    if(cov > Vec_weights[J]) cov = Vec_weights[J];
	}
	if(!nojunc_flag){
	    Unused_pair_paths.push_back(p);
	    Unused_pair_paths_cov.push_back(cov);
	}
    }
    

}
void search(RawGraph& G)
{
    rg_index++;
    trans_id = 1;
    if(G.Mode == "U")
    {
        G.merge_and_get_exon();
        G.junction_decide_exon();
    }
    else
	G.get_exon_and_junction();
    

    vector<vector<int> > Vec_edges;
    vector<double> Vec_weights;

    vector<path_t> Unused_pair_paths;
    vector<double> Unused_pair_paths_cov;
    vector<int> Unused_junctions;
    vector<string> Node_seq_1;
    vector<int>Node_seq_2,Node_seq_3;
    vector<double>Node_seq_4;
    vector<int> junc_l = G.junc_l, junc_r = G.junc_r;
    vector<int> exon_l = G.exon_l, exon_r = G.exon_r;
    vector<vector<int> > PairPathExon = G.pair_path_exon;
    for(int i=0;i<junc_l.size();i++)
    {
        string edge;
        vector<int> edge_;
        stringstream ss;
        for(int j=0;j<exon_l.size();j++)
        {
            if(junc_l[i]-1==exon_r[j]) {
                edge_.push_back(j);
            }
            if(junc_r[i]==exon_l[j]) {
                edge_.push_back(j);
                Vec_edges.push_back(edge_);
                Vec_weights.push_back(G.junc_cov[i]);
            }
        }
//        vec_edges.push_back(edge);
    }
    for(int i=0;i<exon_l.size();i++)
    {
        Node_seq_1.push_back(G.chr);
        Node_seq_2.push_back(exon_l[i]);
        Node_seq_3.push_back(exon_r[i]);
        Node_seq_4.push_back(G.exon_cov[i]);
    }


    if(G.Mode == "U")
    {
      for(size_t i=0;i<exon_l.size()-1;i++)
      {
	if(exon_l[i+1] - exon_r[i] == 1 && G.exon_cov[i] != 0 && G.exon_cov[i+1] != 0) //get partial cov
	{
	    pair<int,int> junc = make_pair(exon_r[i], exon_l[i+1]);
	    if(find(G.junction.begin(),G.junction.end(),junc) != G.junction.end()) continue;
	    bool flag = false;
	    double C = 0;
	    for(size_t j=0;j<G.raw_ref_exon.size();j++)
	    {
		if(junc.first>G.raw_ref_exon[j].first && junc.first < G.raw_ref_exon[j].second)
		{
		    flag = true;
		    C = G.raw_ref_exon_cov[j];
		    break;
		}
	    }
	    if( !flag ) continue;
	    continue;
	    double cov = (G.exon_cov[i] + G.exon_cov[i+1])/2.0;
	    //double cov = C;
	    vector<int> edge;
	    edge.push_back(i); edge.push_back(i+1);
	    Vec_edges.push_back(edge);
	    Vec_weights.push_back(cov);

	}
      }
    }
    get_pair( PairPathExon,Unused_pair_paths,Unused_pair_paths_cov,
	      exon_l, exon_r,
	      junc_l, junc_r, Vec_weights);

    cout<<"Graph: "<<rg_index<<" mode:"<<G.Mode<<"-----"<<endl;
    cout<<"Edges: "<<endl;
    for(int i=0;i<Vec_edges.size();i++) cout<<Vec_edges[i][0]<<"->"<<Vec_edges[i][1]<<": "<<Vec_weights[i]<<endl;
    cout<<"Nodes"<<endl;
    for(int i=0;i<Node_seq_1.size();i++) cout<<i<<": "<<G.strand<<" "<<Node_seq_1[i]<<" "<<Node_seq_2[i]<<" "<<Node_seq_3[i]<<" "<<Node_seq_4[i]<<endl;
    cout<<"Pair"<<endl;
    for(size_t i=0;i<Unused_pair_paths.size();i++) {
        vector<int> pp = Unused_pair_paths[i];
        for(size_t j=0;j<pp.size();j++) cout<<pp[j]<<" ";
        cout<<": "<<Unused_pair_paths_cov[i]<<endl;
    }

    //if(rg_index == 19) exit(0);


    vector<int> edge_out, edge_in;
    vector<double> edge_weight;
    for(int i=0;i<Vec_edges.size();i++)
    {
 	if(Vec_weights[i] < 2) continue; //empirical
	edge_out.push_back(Vec_edges[i][0]);
	edge_in.push_back(Vec_edges[i][1]);
	edge_weight.push_back(Vec_weights[i]);
    }
    mode = G.Mode;

    SEED_Filter = 2;
    mode = "I";
    PathSearch ps( edge_out, edge_in, edge_weight,
                   Unused_pair_paths, Unused_pair_paths_cov,
                   Node_seq_1, Node_seq_2, Node_seq_3,  Node_seq_4,
		   G.strand,G.chr);
    ps.path_search();
    
}

void combine_and_search(vector<RawGraph>& VecGraph)
{
    sort(VecGraph.begin(),VecGraph.end(),sorter);
    for(size_t i=0;i<VecGraph.size();i++)
    {
	if(VecGraph[i].start_pos == -1) continue;
	for(size_t j=i+1;j<VecGraph.size();j++)
	{
	    if(VecGraph[i].end_pos + 0 > VecGraph[j].start_pos )
		combine(VecGraph[i],VecGraph[j]);
	    else
	    {
		search(VecGraph[i]);
		break;
	    }
	    	 
	}
    }
}
class MergeSample{
  public: 
	  string Strand, Chr;
	  vector<vector<string> > GraphOfSamples;
	  vector< vector<pair<int,int> > > CoordinatesOfGraphs;
	  vector< pair<int,int> > GeneRanges;
	  map<junc_t, vector< pair<int,double> > > JunctionCoverage_map;
	  typedef map<junc_t, vector< pair<int,double> > >::iterator junction_iter;

	  map<exon_t, vector< pair<int,double> > > RawExonCoverage_map;
	  typedef map<exon_t, vector<pair<int,double> > >::iterator exon_iter;
	  //map<exon_t, vector< pair<int,double> > > ExonCoverage_map;
	  typedef icl::right_open_interval<int32_t> ROI;
//	  typedef icl::split_interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> split_interval_map;

	  typedef icl::split_interval_map<int32_t, double, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> split_interval_map;

	  typedef split_interval_map::const_iterator split_exon_iter;

	  split_interval_map ExonCoverage_map;

	  string sample_list, Dir;
	  MergeSample(){}
	  MergeSample(string sample_list_,string Dir_)
	  {
		sample_list = sample_list_;
		Dir = Dir_;
		ifstream in(const_cast<char*>(sample_list.c_str()));
		string sample;
	      	while(getline(in,sample)) //each sample
		{
		    sample = Dir + "/" + sample;
		    char* file = const_cast<char*>(sample.c_str());
		    ifstream insample(file);
		    string temp;
		    vector<string> graph_of_one_sample;
		    while(getline(insample,temp)) {
			    graph_of_one_sample.push_back(temp);
		    }
		    //cout<<"graph_of_one_sample size: "<<graph_of_one_sample.size()<<endl;
		    GraphOfSamples.push_back(graph_of_one_sample);

		}
		//cout<<"GraphOfSamples size: "<<GraphOfSamples.size()<<endl;
		SampleSize = double(1.0*GraphOfSamples.size());
	  }
	  void get_graph_range(vector<string>& file,string strand, int sample_id,
			       vector<pair<int,int> >& coordinates_of_one_sample);

	  void load_graph(vector<string>& file,string strand, vector<pair<int,int> >& coordinates_of_one_sample);

  	  void process_junction();
	  void sample_cluster();
	  void process(string STRD);
	  void process_each_gene();
	  double get_coverage(vector< pair<int,double> > juncInfo);

	  void get_partial_junction( 	vector<pair<int,int> >& junction,
					vector< vector< pair<int,double> > >& junction_cov,
					vector<pair<int,int> >& raw_exon, 
					vector< vector< pair<int,double> > >& raw_exon_cov,
					vector<pair<int,int> >& exon);
	  bool IsMetaGraph(vector<junc_t> junctions,vector< vector< pair<int,double> > >junction_cov);
	  void get_meta_graph(vector<pair<int,int> > junction,
			      vector< vector< pair<int,double> > >junction_cov,
			      vector<pair<int,int> > raw_exon,
			      vector< vector< pair<int,double> > > raw_exon_cov,
			      vector<pair<int,int> > exon,
			      vector<double> exon_cov);

	  int create_split(split_interval_map &imap, int32_t p);

	  void test_split_interval_map();

};

void MergeSample::get_graph_range(vector<string>& file,string strand, int sample_id,
				  vector<pair<int,int> >& coordinates_of_one_sample)//graph corordinates
{
    istringstream istr;
    string s,temp;
    bool edge_flag = false, node_flag = false, pair_flag = false;
    string strd = "", chr = "";
    int start_pos=0;
    vector<int> exon_l,exon_r,junc_l,junc_r,edge_out,edge_in;
    vector<double> exon_cov, junc_cov, edge_weight;
    vector<pair<int,int> > junction, exon;
    vector<path_t> ppaths;
    vector<double> ppaths_cov;
    for(size_t i=0;i<file.size();i++)
    {
	string s = file[i];
        if( s == "Edges"){ edge_flag = true; node_flag = false;pair_flag = false;continue;}
        if( s == "Nodes") { edge_flag = false; node_flag = true;pair_flag = false;continue;}
        if( s == "Pair") {edge_flag = false; node_flag = false;pair_flag = true;continue;}
        if( s[0] == '#') {
	    if(edge_out.size() > 0 && strd == strand)
	    {
	       	pair<int,int> codnt=make_pair(exon_l.front(),exon_r.back());
		coordinates_of_one_sample.push_back(codnt);
	      	for(int i=0;i<edge_out.size();i++)
		{
		    junc_t J = make_pair(exon_r[edge_out[i]], exon_l[edge_in[i]]);
		    pair<int,double> sc = make_pair(sample_id,edge_weight[i]);
		    if(JunctionCoverage_map.find(J) == JunctionCoverage_map.end())
		    {
		        vector<pair<int,double> > v (1,sc);
			JunctionCoverage_map[J] = v;
		    }
		    else JunctionCoverage_map[J].push_back(sc);

		}
		for(size_t i=0;i<exon_l.size();i++)
		{
		    ExonCoverage_map += make_pair(ROI(exon_l[i], exon_r[i]), exon_cov[i]);
		    exon_t E = make_pair(exon_l[i], exon_r[i]);
		    double cov = exon_cov[i];
		    pair<int,double> sc = make_pair(sample_id,exon_cov[i]);
		    if(RawExonCoverage_map.find(E) == RawExonCoverage_map.end())
		    {
			 vector<pair<int,double> > v (1,sc);
		        RawExonCoverage_map[E] = v;
		    }
		    else RawExonCoverage_map[E].push_back(sc);
		    
		}

	    }
	    junc_l.clear();junc_r.clear();junc_cov.clear(); 
	    junction.clear();
	    exon_l.clear();exon_r.clear();exon_cov.clear();
	    exon.clear();
	    edge_out.clear();edge_in.clear();edge_weight.clear();
	    ppaths.clear(); ppaths_cov.clear();
            edge_flag = false; node_flag = false;pair_flag = false;continue;
        }

        if(edge_flag) {

            istr.str(s);
            int out,in;
            vector<int> v;
            double cov;
            istr>>out>>temp>>in>>temp>>cov;
            istr.clear();
            edge_out.push_back(out); edge_in.push_back(in);
            edge_weight.push_back(cov);
	}
        if(node_flag)
        {
 	    istr.str(s);
            int left, right;
            double cov;
            istr>>strd>>chr>>left>>right>>cov;
	    Chr = chr;
	    if( cov == 0) cov = 0.01;
            istr.clear();
            exon_l.push_back(left); exon_r.push_back(right);exon_cov.push_back(cov);
	    exon.push_back(make_pair(left,right));

        }  
	if(pair_flag) continue;
    }
    file.clear();
    return;
}
void MergeSample::load_graph(vector<string>& file,string strand, vector<pair<int,int> >& coordinates_of_one_sample) //graph corordinates
{
    istringstream istr;
    string s,temp;
    bool edge_flag = false, node_flag = false, pair_flag = false;
    string strd = "", chr = "";
    vector<double> gene_vec;
    int start_pos=0;
    vector<int> exon_l,exon_r,junc_l,junc_r,edge_out,edge_in;
    vector<double> exon_cov, junc_cov, edge_weight;
    vector<pair<int,int> > junction, exon;
    vector<path_t> ppaths;
    vector<double> ppaths_cov;
    for(size_t i=0;i<file.size();i++)
    {
	string s = file[i];
        if( s == "Edges"){ edge_flag = true; node_flag = false;pair_flag = false;continue;}
        if( s == "Nodes") { edge_flag = false; node_flag = true;pair_flag = false;continue;}
        if( s == "Pair") {edge_flag = false; node_flag = false;pair_flag = true;continue;}
        if( s[0] == '#') {
	    if(edge_out.size() > 0 && strd == strand)
	    {
	       	pair<int,int> codnt=make_pair(exon_l.front(),exon_r.back());
		coordinates_of_one_sample.push_back(codnt);
	      	for(int i=0;i<edge_out.size();i++)
		    junction.push_back(make_pair(exon_r[edge_out[i]], exon_l[edge_in[i]]));

	      	RawGraph rg(junction,edge_weight,start_pos,gene_vec,strd,chr,mode,ppaths,ppaths_cov, exon, exon_cov);
	        //if(strd == "+") plus_graphs.push_back(rg);
	        //else if(strd == "-") minus_graphs.push_back(rg);
		//else if(strd == ".") unstranded_graphs.push_back(rg);
	    }
	    junc_l.clear();junc_r.clear();junc_cov.clear(); gene_vec.clear();
	    junction.clear();
	    exon_l.clear();exon_r.clear();exon_cov.clear();
	    exon.clear();
	    edge_out.clear();edge_in.clear();edge_weight.clear();
	    ppaths.clear(); ppaths_cov.clear();
            edge_flag = false; node_flag = false;pair_flag = false;continue;
        }

        if(edge_flag)
        {
            istr.str(s);
            int out,in;
            vector<int> v;
            double cov;
            istr>>out>>temp>>in>>temp>>cov;
            istr.clear();
            edge_out.push_back(out); edge_in.push_back(in);
            edge_weight.push_back(cov);
        }
        if(node_flag)
        {
 	    istr.str(s);
            int left, right;
            double cov;
            istr>>strd>>chr>>left>>right>>cov;
	    if( cov == 0) cov = 0.01;
            istr.clear();
            exon_l.push_back(left); exon_r.push_back(right);exon_cov.push_back(cov);
	    exon.push_back(make_pair(left,right));
	    /*
	    if(gene_vec.empty()) start_pos = left;
	    int S = gene_vec.size();
	    for(int i=gene_vec.size(); i < left-start_pos;i++) gene_vec.push_back(0); //intron 
	    for(int i=gene_vec.size(); i <= right-start_pos;i++) gene_vec.push_back(cov); //exon 
	    */

        }  
	if(pair_flag)
        {
            istr.str(s);
            double cov;
            vector<int> v;
            while(istr>>temp)
            {
                if(temp == ":"){
                    istr>>cov;
                    break;
                }
                else
                    v.push_back(atoi(temp.c_str()) );
            }
            istr.clear();
            ppaths.push_back(v);
            ppaths_cov.push_back(cov);
        }
    }
    return;
}

void MergeSample::sample_cluster()
{
    vector<size_t> vec_index(CoordinatesOfGraphs.size(),0);
    bool Stop = true;
    pair<int,int> range = make_pair(-1,-1);
    while(1) //for each gene(index)
    {
        //pair<int,int> range = CoordinatesOfGraphs[vec_index.front()].front();
	/*
	cout<<"index: ";
	for(int j = 0;j<vec_index.size();j++) cout<<vec_index[j]<<" ";
	cout<<endl;
	*/
	bool get_one_range = true;
	for(size_t i=0;i<CoordinatesOfGraphs.size();i++) //each sample
	{
	   if(vec_index[i] >= CoordinatesOfGraphs[i].size()) continue; //each gene
	   else
	   {
		Stop = false;
	   	if(range.second == -1 || CoordinatesOfGraphs[i][vec_index[i]].first < range.second)  //each gene
	   	{

		    if(range.second == -1 || CoordinatesOfGraphs[i][vec_index[i]].first < range.first) 
			    range.first = CoordinatesOfGraphs[i][vec_index[i]].first;

		    if(range.second == -1 || CoordinatesOfGraphs[i][vec_index[i]].second > range.second) 
			     range.second = CoordinatesOfGraphs[i][vec_index[i]].second;
		    get_one_range = false;
		    vec_index[i]++;
	   	}
		/*
		if(i == CoordinatesOfGraphs.size() - 1 && get_one_range){
			cout<<"range: "<<range.first<<" "<<range.second<<endl;
			cout<<"index_new: ";
			for(int j = 0;j<vec_index.size();j++) cout<<vec_index[j]<<" ";
			cout<<endl;
			break;
			range = make_pair(-1,-1);
		}
		else get_one_range = true;
		*/
	   } 
	}
	if(get_one_range && range.first != -1 && range.second != -1){
	    /*
	    cout<<"range: "<<range.first<<" "<<range.second<<endl;
	    cout<<"index_new: ";
	    for(int j = 0;j<vec_index.size();j++) cout<<"("<<j<<")"<<vec_index[j]<<" ";
	    cout<<endl;
	    */
	    //break;
	    GeneRanges.push_back(range);
	    range = make_pair(-1,-1);
	}

	if(Stop) break;
	Stop = true;

    }
}
/*
int MergeSample::create_split(split_interval_map &imap, int32_t p)
{
     split_exon_iter it = imap.find(p);
     if(it == imap.end()) return 0;
     int32_t l = lower(it->first);
     int32_t r = upper(it->first);
     int32_t w = it->second;
     assert(l <= p);
     assert(r >= p);
     if(l == p || r == p) return 0;
     imap -= make_pair(ROI(l, r), w);
     imap += make_pair(ROI(l, p), w);
     imap += make_pair(ROI(p, r), w);
     return 0;
}
*/
void MergeSample::test_split_interval_map()
{

	split_interval_map imap;

        imap += make_pair(ROI(6, 7), 3.0);
        imap += make_pair(ROI(1, 3), 3.0);
        imap += make_pair(ROI(1, 2), 1.0);
        imap += make_pair(ROI(2, 5), 2.0);

 //       create_split(imap, 4);

	split_exon_iter it;

	for(it = imap.begin(); it != imap.end(); it++)
	{
	         //printf("interval: [%d,%d) -> %d\n", lower(it->first), upper(it->first), it->second);
		 cout<<"interval: "<<lower(it->first)<<" "<<upper(it->first)<<": "<<it->second<<endl;
	}
	it = imap.find(2);
	if(it  != imap.end())
		cout<<"interval: "<<lower(it->first)<<" "<<upper(it->first)<<": "<<it->second<<endl;
	else cout<<"not find"<<endl;

	it = imap.upper_bound(ROI(4, 6));



}
void MergeSample::process(string STRD)
{
	//test_split_interval_map();
	//return;
	ifstream in(const_cast<char*>(sample_list.c_str()));
	string sample;
	for(size_t i=0;i<GraphOfSamples.size();i++)
	{
	
	    vector<string> one_sample = GraphOfSamples[i];
	    vector<pair<int,int> > coordinates,exons;
	    //cout<<"sample: "<<i<<" ";
	    get_graph_range(one_sample,STRD,i,coordinates);//get graph info
	    Strand = STRD;

	    CoordinatesOfGraphs.push_back(coordinates);

	}
	//cerr<<Chr<<" "<<Strand<<endl;
	//cout<<"CoordinatesOfGraphs size: "<<CoordinatesOfGraphs.size()<<endl;
	//cout<<"JunctionCoverage_map size(): "<<JunctionCoverage_map.size()<<endl;
	//cout<<"ExonCoverage_map size(): "<<ExonCoverage_map.size()<<endl;
	//cout<<"RawExonCoverage_map size(): "<<RawExonCoverage_map.size()<<endl;
	/*
	for(size_t i=0;i<CoordinatesOfGraphs.size();i++)
	{
	    vector<pair<int,int> > coords = CoordinatesOfGraphs[i];
	    cout<<"sample_"<<i<<": ";
	    for(size_t j=0;j<10;j++) cout<<coords[j].first<<" "<<coords[j].second<<"; ";
	    cout<<endl;

	    
	}
	*/
	/*
	split_exon_iter it;
	for(it = ExonCoverage_map.begin(); it != ExonCoverage_map.end();it++)
		cout<<"exon: "<<lower(it->first)<<"-"<<upper(it->first)<<": "<<it->second<<endl;
	*/

	sample_cluster();

	process_each_gene();
	//cout<<"GeneRanges: "<<GeneRanges.size()<<endl;

	int A=0,B=0,C=0,D=0,E=0;
	for(junction_iter i = JunctionCoverage_map.begin();i!= JunctionCoverage_map.end();i++)
	{
	    /*
	    cout<<i->first.first<<"->"<<i->first.second<<" ";
	    for(size_t j=0;j<i->second.size();j++) cout<<i->second[j].first<<","<<i->second[j].second<<"; ";
	    cout<<endl;
	    */
	    if(i->second.size()  == 1 ) E++;
	    else if(i->second.size() < 10) A++;
	    else if(i->second.size() < 30) B++;
	    else if(i->second.size() < 50) C++;
	    else if(i->second.size() <= 73) D++;

	}
	//cout<<"==1,<10,<30,<50,<73: "<<E<<","<<A<<","<<B<<","<<C<<","<<D<<endl;
	//for(int i=0;i<1;){}
}
double MergeSample::get_coverage(vector< pair<int,double> > juncInfo)
{
    //do Kalman Filter here
    double cov = 0.0;
    double cov2 = 0.0;
    for(size_t i=0;i<juncInfo.size();i++)
    {
        cov += juncInfo[i].second;
	cov2 += juncInfo[i].second*juncInfo[i].second;
    }

    //cov = cov/(double(1.0*juncInfo.size()));
    cov = double(cov/(1.0*SampleSize));
  //  cov = double(sqrt(cov2)/(1.0*SampleSize)); //revision
    return cov;
}
void MergeSample::get_partial_junction( vector<pair<int,int> >& junction,
					vector< vector< pair<int,double> > >& junction_cov,
					vector<pair<int,int> >& raw_exon, 
					vector< vector< pair<int,double> > >& raw_exon_cov,
					vector<pair<int,int> >& exon)
{
    for(size_t i=0;i<exon.size()-1;i++)
    {
        if(exon[i+1].first - exon[i].second == 1)//partial
	{
	    junc_t J = make_pair(exon[i].second, exon[i+1].first);
	    if(find(junction.begin(),junction.end(),J) == junction.end())//not record
	    {
	       vector< pair<int,double> > sample_cov;
	       bool flag = false;;
	       for(size_t j=0;j<raw_exon.size();j++)
	       {
	           if(J.first > raw_exon[j].first && J.second <  raw_exon[j].second)
		   {
		       flag = true;
		       for(size_t k =0;k<raw_exon_cov[j].size();k++) 
			       sample_cov.push_back(raw_exon_cov[j][k]);
		   }
		   if(raw_exon[j].first >= J.second) break;
	       }
		if(flag)
		{
	       	   junction.push_back(J);
	           junction_cov.push_back(sample_cov);
		}
	    }
	}
    }
    return;
}

bool MergeSample::IsMetaGraph(vector<junc_t> junctions, vector< vector< pair<int,double> > >junction_cov) //pair-><sapleid,cov>
{
    vector<int> edgesNumber_of_eachSample(SampleSize,0);
    for(size_t i=0;i<junction_cov.size();i++)
    {
	if(junctions[i].second - junctions[i].first == 1) continue;
        for(size_t j=0;j<junction_cov[i].size();j++)
	{
	    int sampleid = junction_cov[i][j].first;
	    edgesNumber_of_eachSample[sampleid]++;
	}
    }
    int biggest = *max_element(edgesNumber_of_eachSample.begin(),edgesNumber_of_eachSample.end());
    int smallest = *min_element(edgesNumber_of_eachSample.begin(),edgesNumber_of_eachSample.end());

    //if(double(1.0*smallest)/(1.0*biggest) <= 0.0) return false;
    int Number = 0; //number of sample that contians no edge int this graph;
    for(size_t i=0;i<edgesNumber_of_eachSample.size();i++)
    {
        if(edgesNumber_of_eachSample[i] == 0) Number++;
    }
    if(Number > 10) return false;
    return true;
}
void MergeSample::get_meta_graph(vector<pair<int,int> > junction,
		                 vector< vector< pair<int,double> > >junction_cov,
		                 vector<pair<int,int> > raw_exon,
		                 vector< vector< pair<int,double> > > raw_exon_cov,
		                 vector<pair<int,int> > exon,
		                 vector<double> exon_cov)
{
    /*
    cout<<"Edges: "<<endl;
    for(int a=0;a<junction.size();a++)
    {
	cout<<junction[a].first<<" "<<junction[a].second<<": ";
	for(int b=0;b<junction_cov[a].size();b++) cout<<junction_cov[a][b].first<<"-"<<junction_cov[a][b].second<<"; ";
	cout<<endl;
    }
    cout<<"Nodes: "<<exon.size()<<endl;
    
    for(int a=0;a<exon.size();a++)
	cout<<exon[a].first<<" "<<exon[a].second<<" "<<exon_cov[a]<<endl;
    
    cout<<"RawNodes: "<<endl;
    for(int a=0;a<raw_exon.size();a++)
    {
	cout<<"raw: "<<raw_exon[a].first<<" "<<raw_exon[a].second<<": ";
	for(int b=0;b<raw_exon_cov[a].size();b++) cout<<raw_exon_cov[a][b].first<<"-"<<raw_exon_cov[a][b].second<<"; ";
	cout<<endl;
    }
    */
    vector<int>junc_l,junc_r;
    for(size_t i=0;i<junction.size();i++)
    {
        junc_l.push_back(junction[i].first);
	junc_r.push_back(junction[i].second);
    }
    for(size_t i=0;i<exon.size()-1;)
    {
        if(exon[i+1].first == exon[i].second) //partial [10-20),[20,30)
	{
	    if( find(junc_l.begin(),junc_l.end(),exon[i].second) == junc_l.end() 
		&& find(junc_r.begin(),junc_r.end(),exon[i+1].first) == junc_r.end()) //not splice position
		{
		    exon[i].second = exon[i+1].second;
		    exon_cov[i] += exon_cov[i+1];
		    exon.erase(exon.begin() + i + 1);
		    exon_cov.erase(exon_cov.begin() + i + 1);
		}
	    else if( find(junc_l.begin(),junc_l.end(),exon[i].second) != junc_l.end() ){ 
	        exon[i+1].first++;// [10-20) [21,30)
		i++;
	    }
	    else if(find(junc_r.begin(),junc_r.end(),exon[i+1].first) != junc_r.end())
	    {
	        exon[i].second--; // [10-19) ,[20,30)
		i++;
	    }
	    else i++;
	}
	else i++;
    }
    for(size_t i=0;i<exon.size();)
    {
        if(exon[i].first > exon[i].second){
	    exon.erase(exon.begin() + i);
	    exon_cov.erase(exon_cov.begin() + i);
	}
	else i++;
    }
    /*
     cout<<"Nodes-new: "<<exon.size()<<endl;
     for(int a=0;a<exon.size();a++) 
         cout<<a<<" "<<exon[a].first<<" "<<exon[a].second<<" "<<exon_cov[a]<<endl;
	 */
    get_partial_junction(junction,junction_cov, raw_exon,raw_exon_cov,exon);

     vector<int> exon_l,exon_r;
     vector<string> exon_chr;
     for(size_t i=0;i<exon.size();i++){
         exon_l.push_back(exon[i].first);
	 exon_r.push_back(exon[i].second);
	 exon_chr.push_back(Chr);
	 exon_cov[i] = exon_cov[i]/double(SampleSize);
	 exon_cov[i] = exon_cov[i]/double(SampleSize);

     }
     //get_edges
     vector<int> edge_out, edge_in;
     vector<double> edge_weight;
     //if( !IsMetaGraph(junction,junction_cov) ) return;
     for(size_t i=0;i<junction.size();i++){
         vector<int>::iterator a,b;
	 a = find(exon_r.begin(),exon_r.end(),junction[i].first);
	 b = find(exon_l.begin(),exon_l.end(),junction[i].second);
	 if( a != exon_r.end() && b != exon_l.end())
	 {
	     edge_out.push_back(a - exon_r.begin());
	     edge_in.push_back(b - exon_l.begin());
	     double cov = get_coverage(junction_cov[i]);
	     edge_weight.push_back(cov);
	 }
     }
     /*
     cout<<"Edges: "<<endl;
     for(size_t i=0;i<edge_out.size();i++) cout<<edge_out[i]<<" -> "<<edge_in[i]<<": "<<edge_weight[i]<<endl;
     cout<<"Graph"<<endl<<endl;
     */
     vector< vector<int> > Unused_pair_paths;
     vector<double> Unused_pair_paths_cov;
     PathSearch ps( edge_out, edge_in, edge_weight,
                   Unused_pair_paths, Unused_pair_paths_cov,
                   exon_chr,exon_l,exon_r,exon_cov,
                   Strand,Chr);
     ps.junctions = junction;
     ps.junctions_MappingInfo = junction_cov;
     ps.path_search();

}
void MergeSample::process_each_gene()
{
	//out_graph.open("Graph.txt",ios::app);
	junction_iter j=JunctionCoverage_map.begin();
	split_exon_iter k = ExonCoverage_map.begin();
	exon_iter m = RawExonCoverage_map.begin();
	for(size_t i=0;i<GeneRanges.size();i++)
	{
	    //rg_index++;
	    trans_id = 1;
	    //cout<<"Gene Range: "<<GeneRanges[i].first<<" "<<GeneRanges[i].second<<endl;
	    int l = GeneRanges[i].first, r = GeneRanges[i].second;
	    //cout<<"Gene-"<<i<<" range: "<<l<<" "<<r<<endl;
	    vector<pair<int,int> > junction;
	    vector<pair<int,int> > exon,raw_exon;
	    vector<double> exon_cov;
	    vector< vector< pair<int,double> > >junction_cov, raw_exon_cov;

	    for(;j != JunctionCoverage_map.end();j++)
	    {
	    	if((j->first.first >= l && j->first.first <= r)
		   && (j->first.second >= l && j->first.second <= r))
		{
		    //junction.push_back(make_pair(j->first.first,j->first.second));
		    junction.push_back(j->first);
		    junction_cov.push_back(j->second);
		}
		else{
		    break;
		}
	    }
	    for(;k != ExonCoverage_map.end();k++)
	    {
		if(lower(k->first) >= l && upper(k->first) <= r)
		{
		    exon.push_back(make_pair(lower(k->first),upper(k->first)));
		    exon_cov.push_back(k->second);
		}
		else{
		    break;
		}
	    }
	    for(;m != RawExonCoverage_map.end();m++)
	    {
	        if((m->first.first >= l && m->first.first <= r)
		   && (m->first.second >= l && m->first.second <= r))
		{
		    raw_exon.push_back(m->first);
		    raw_exon_cov.push_back(m->second);
		}
		else{
		    break;
		}    
	    }
	    get_meta_graph(junction,junction_cov, raw_exon, raw_exon_cov, exon, exon_cov);
	}
	//out_graph.close();
	return;
}
