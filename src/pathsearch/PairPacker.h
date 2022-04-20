#ifndef IPACACKER_H
#define IPACACKER_H

// PairPacker.h

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include<algorithm>
#include<cstring>
using namespace std;
extern bool SFlag;
class PairPacker
{
	private:
	vector<vector<int> > Path_cover;
	vector<string> Node_seq_1;
	vector<int> Node_seq_2,Node_seq_3;
        vector<double> Node_seq_4;
	int Graph_num;
	int Trans_id;
	vector<int> Single_nodes;// from junction_graph.Single_nodes
	char * Output_path;
	int Path_length;
	int Strand;
	public:
	PairPacker(vector<vector<int> > path_cover,vector<string> node_seq_1,vector<int> node_seq_2,vector<int> node_seq_3,vector<double> node_seq_4,int graph_num,int trans_id,vector<int> single_nodes,char * output_path,int path_length,int strand);
	~PairPacker();
	void PairPacker_Output();
	int Curr_id;
};
PairPacker::PairPacker(vector<vector<int> > path_cover,vector<string> node_seq_1,vector<int> node_seq_2,vector<int> node_seq_3,vector<double> node_seq_4,int graph_num,int trans_id,vector<int> single_nodes,char * output_path,int path_length,int strand)
{
	Path_cover=path_cover;
	Node_seq_1=node_seq_1;
	Node_seq_2=node_seq_2;
	Node_seq_3=node_seq_3;
	Node_seq_4=node_seq_4;
	Graph_num=graph_num;
	Trans_id=trans_id;
	Single_nodes=single_nodes;
	Output_path=output_path;
	Path_length=path_length;
	Strand=strand;
	return;
}
PairPacker::~PairPacker()
{
	return;
}
void PairPacker::PairPacker_Output()
{
	Curr_id=Trans_id;
	size_t i,j,k;
	int strand_id;
	string Symble=".+-";
	if (Strand==0)
	{
		strand_id=0;
	}
	else if (Strand==1)
	{
		strand_id=1;
	}
	else
	{
		strand_id=2;
	}
	ofstream file;
	//file.open(Output_path);
        file.open(Output_path,ios_base::app);


	//ofstream fa("TransComb.fa",ios::app);//NEW
	/*
	for(int i = 0;i<Path_cover.size();i++){
	    vector<int> vec = Path_cover[i];
	    for(int j =0;j<vec.size();j++) cerr<<vec[j]<<"->";
	    cerr<<endl;
	}
	cerr<<"Graph"<<Graph_num<<endl;
	*/
/*
	sort(Path_cover.begin(),Path_cover.end());
 	Path_cover.erase(unique(Path_cover.begin(),Path_cover.end()),Path_cover.end());
*/
	double L = 500, COV1 = 5, COV2=10;//Yu-Borrow
/*
	if(SFlag){
	    L = 800;
	    COV1 = 6;
	    COV2 = 20;
	}
*/
	for (i=0;i<Path_cover.size();i++)
        {
		//break;//8.1NEW
		k=0;
		double cov = 0;
                for (j=0;j<Path_cover[i].size();j++)
                {
                        k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
			cov += Node_seq_4[Path_cover[i][j]]*(Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1);
                }
		cov = double(cov/k);
		if (k>=size_t(Path_length) || (Path_cover.size() == 1 && cov > COV2))
	   	//if (!(k>=size_t(Path_length)) && Path_cover.size() == 1 && cov > 6)//NEW-borrow
		{
			file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"
			    <<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"
			    <<1000<<"	"<<Symble[strand_id]
			    <<"	.	gene_id "<<"\""<<"PAIRP."<<Graph_num<<"\""<<"; "
			    <<"transcript_id "<<"\""<<"PAIRP."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
			//NEW
			//fa<<">IPAC."<<Graph_num<<"."<<Curr_id<<" gene="<<"IPAC."<<Graph_num<<" "<<Symble[strand_id]<<'\n';

			for (j=0;j<Path_cover[i].size()-1;j++)
               		{
				file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"
				    <<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"
				    <<1000<<"	"<<Symble[strand_id]
				    <<"	.	gene_id "<<"\""<<"PAIRP."<<Graph_num<<"\""<<"; "
				    <<"transcript_id "<<"\""<<"PAIRP."<<Graph_num<<"."<<Curr_id<<"\""<<"; "
				    <<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;

				//fa<<exon_sequence[Path_cover[i][j]];

			}
			file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"
			    <<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"
			    <<1000<<"	"<<Symble[strand_id]
			    <<"	.	gene_id "<<"\""<<"PAIRP."<<Graph_num<<"\""<<"; "
			    <<"transcript_id "<<"\""<<"PAIRP."<<Graph_num<<"."<<Curr_id<<"\""<<"; "
			    <<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;

			//NEW
			//fa<<exon_sequence[Path_cover[i][Path_cover[i].size()-1]]<<'\n';
			
			Curr_id++;
		}
	}
	for (i=0;i<Single_nodes.size();i++)
        {
	//	break;
		int length = Node_seq_3[Single_nodes[i]]-Node_seq_2[Single_nodes[i]]+1;
		double cov = Node_seq_4[Single_nodes[i]];
//		if (Node_seq_3[Single_nodes[i]]-Node_seq_2[Single_nodes[i]]+1>=Path_length && Node_seq_4[Single_nodes[i]]>20)
		if((Node_seq_3[Single_nodes[i]]-Node_seq_2[Single_nodes[i]]+1>=L
		   && Node_seq_4[Single_nodes[i]]>=COV1) || cov > COV2)
		//if( !(Node_seq_3[Single_nodes[i]]-Node_seq_2[Single_nodes[i]]+1>=500
		   //&& Node_seq_4[Single_nodes[i]]>=5)  && cov >10)
		{
			file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"
			    <<"transcript"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"
			    <<1000<<"	"<<Symble[strand_id]
			    <<"	.	gene_id "<<"\""<<"PAIRP."<<Graph_num<<"\""<<"; "
			    //<<"transcript_id "<<"\""<<"PAIRP."<<Graph_num<<"."<<Curr_id<<"-"<<length<<"-"<<cov<<"\";"<<endl;
			    <<"transcript_id "<<"\""<<"PAIRP."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;


			//NEW
			//fa<<"IPAC."<<Graph_num<<"."<<Curr_id<<" gene="<<"IPAC."<<Graph_num<<" "<<Symble[strand_id]<<'\n';


                        file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"
			    <<"exon"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"
			    <<1000<<"	"<<Symble[strand_id]
			    <<"	.	gene_id "<<"\""<<"PAIRP."<<Graph_num<<"\""<<"; "
			    //<<"transcript_id "<<"\""<<"PAIRP."<<Graph_num<<"."<<Curr_id<<"-"<<length<<"-"<<cov<<"\""<<"; "
			    <<"transcript_id "<<"\""<<"PAIRP."<<Graph_num<<"."<<Curr_id<<"\""<<"; "
			    <<"exon_number "<<"\""<<1<<"\""<<";"<<endl;

			//NEW
			//fa<<exon_sequence[Single_nodes[i]]<<'\n';
			Curr_id++;
		}
	}
	file.close();
	return;
}// end of PairPacker_Output

#endif


