#ifndef PAIRPACKER_H
#define PAIRPACKER_H

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
	vector<vector<double> > Edges_cov;
	string Trans_name;
	public:
	PairPacker(vector<vector<int> > path_cover,vector<string> node_seq_1,vector<int> node_seq_2,vector<int> node_seq_3,vector<double> node_seq_4,int graph_num,int trans_id,vector<int> single_nodes,char * output_path,int path_length,int strand,vector<vector<double> > edges_cov,string trans_name);
	~PairPacker();
	void PairPacker_Output();
	void PairPacker_Output_combine();
	void PairPacker_Output_short();
	int Curr_id;
};
PairPacker::PairPacker(vector<vector<int> > path_cover,vector<string> node_seq_1,vector<int> node_seq_2,vector<int> node_seq_3,vector<double> node_seq_4,int graph_num,int trans_id,vector<int> single_nodes,char * output_path,int path_length,int strand,vector<vector<double> > edges_cov,string trans_name)
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
	Edges_cov=edges_cov;
	Trans_name=trans_name;
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
        file.open(Output_path,ios_base::app);
    if (Path_cover.size() == 1)
    {
	    for (i=0;i<Path_cover.size();i++)
            {
		k=0;
                for (j=0;j<Path_cover[i].size();j++)
                {
                        k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
                }
		if (k>=size_t(Path_length) && k>=500)
		{
			file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
			for (j=0;j<Path_cover[i].size()-1;j++)
               		{
				file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			}
			file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			Curr_id++;
		}
	    }
    }
    if (Path_cover.size() >=2)
    {
	for (i=0;i<Path_cover.size();i++)
        {
		k=0;
                for (j=0;j<Path_cover[i].size();j++)
                {
                        k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
                }
		if (k>=size_t(Path_length) && k>=500)
		{
			file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
			for (j=0;j<Path_cover[i].size()-1;j++)
               		{
				file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			}
			file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			Curr_id++;
		}
	}
    }
/*
	for (i=0;i<Single_nodes.size();i++)
        {
		if (Node_seq_3[Single_nodes[i]]-Node_seq_2[Single_nodes[i]]+1>=500 && Node_seq_4[Single_nodes[i]]>=5)
		{
			file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
                        file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<1<<"\""<<";"<<endl;
			Curr_id++;
		}
	}
*/
	file.close();
	return;
}// end of PairPacker_Output

void PairPacker::PairPacker_Output_combine()
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
        file.open(Output_path,ios_base::app);
    if (Path_cover.size() == 1)
    {
	    for (i=0;i<Path_cover.size();i++)
            {
		k=0;
                for (j=0;j<Path_cover[i].size();j++)
                {
                        k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
                }
		if (k>=0)
		{
			file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
			for (j=0;j<Path_cover[i].size()-1;j++)
               		{
				file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			}
			file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			Curr_id++;
		}
	    }
    }
    if (Path_cover.size() >=2)
    {
	for (i=0;i<Path_cover.size();i++)
        {
		k=0;
                for (j=0;j<Path_cover[i].size();j++)
                {
                        k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
                }
		if (k>=0)
		{
			file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
			for (j=0;j<Path_cover[i].size()-1;j++)
               		{
				file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			}
			file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			Curr_id++;
		}
	}
    }
/*
	for (i=0;i<Single_nodes.size();i++)
        {
		if (Node_seq_3[Single_nodes[i]]-Node_seq_2[Single_nodes[i]]+1>=500 && Node_seq_4[Single_nodes[i]]>=5)
		{
			file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
                        file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<1<<"\""<<";"<<endl;
			Curr_id++;
		}
	}
*/
	file.close();
	return;
}// end of PairPacker_Output
void PairPacker::PairPacker_Output_short()
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
        file.open(Output_path,ios_base::app);
    if (Path_cover.size() == 1)
    {
	if (Path_cover[0].size() > 2) return;
	map<vector<double>,double> hash_edge_cov;
	for (i=0;i<Edges_cov.size();i++)
	{
		vector<double> edge;
		edge.push_back(Edges_cov[i][0]);
		edge.push_back(Edges_cov[i][1]);
		double cov = Edges_cov[i][2];
		hash_edge_cov[edge] = cov;
	}
	i=0;
	int is_single=0;
	vector<double> edge;
	edge.push_back(Edges_cov[0][0]);
	edge.push_back(Edges_cov[0][1]);
	double max_edge_cov = hash_edge_cov[edge];
	edge.clear();
	k=0;
        for (j=0;j<Path_cover[i].size();j++)
        {
                k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
        }
	for (j=0;j<Path_cover[i].size()-1;j++)
	{
		if (Node_seq_3[Path_cover[i][j]]+1 == Node_seq_2[Path_cover[i][j+1]])
		{
			is_single++;
		}
		edge.push_back(Path_cover[i][j]);
		edge.push_back(Path_cover[i][j+1]);
		if (hash_edge_cov[edge] > max_edge_cov)
		{
			max_edge_cov = hash_edge_cov[edge];
		}
	}
	if (is_single == Path_cover[i].size()-1)
	{
	    if (max_edge_cov >= 10 && k>=500)
	    {
		file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
		for (j=0;j<Path_cover[i].size()-1;j++)
             	{
			file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
		}
		file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
		Curr_id++;
	    }
	}
	else
	{
	    for (i=0;i<Path_cover.size();i++)
            {
		k=0;
                for (j=0;j<Path_cover[i].size();j++)
                {
                        k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
                }
		double min_node_cov=Node_seq_4[Path_cover[i][0]];
		for (j=0;j<Path_cover[i].size();j++)
		{
			if (min_node_cov > Node_seq_4[Path_cover[i][j]])
			{
				min_node_cov=Node_seq_4[Path_cover[i][j]];
			}
		}
		if (k>=size_t(Path_length) && k>=500)
		{
			file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
			for (j=0;j<Path_cover[i].size()-1;j++)
               		{
				file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			}
			file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			Curr_id++;
		}
	    }
	}
    }
    if (Path_cover.size() >=2)
    {
	for (i=0;i<Path_cover.size();i++)
        {
		if (Path_cover[i].size() > 2) continue;
		k=0;
                for (j=0;j<Path_cover[i].size();j++)
                {
                        k=k+Node_seq_3[Path_cover[i][j]]-Node_seq_2[Path_cover[i][j]]+1;
                }
		double min_node_cov=Node_seq_4[Path_cover[i][0]];
		for (j=0;j<Path_cover[i].size();j++)
		{
			if (min_node_cov > Node_seq_4[Path_cover[i][j]])
			{
				min_node_cov=Node_seq_4[Path_cover[i][j]];
			}
		}
		if (k>=size_t(Path_length) && k>=500)
		{
			file<<Node_seq_1[Path_cover[i][0]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Path_cover[i][0]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
			for (j=0;j<Path_cover[i].size()-1;j++)
               		{
				file<<Node_seq_1[Path_cover[i][j]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][j]]<<"	"<<Node_seq_3[Path_cover[i][j]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			}
			file<<Node_seq_1[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<Node_seq_3[Path_cover[i][Path_cover[i].size()-1]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<j+1<<"\""<<";"<<endl;
			Curr_id++;
		}
	}
    }

	for (i=0;i<Single_nodes.size();i++)
        {
		if (Node_seq_3[Single_nodes[i]]-Node_seq_2[Single_nodes[i]]+1>=500 && Node_seq_4[Single_nodes[i]]>=20)
		{
			file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"<<"transcript"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\";"<<endl;
                        file<<Node_seq_1[Single_nodes[i]]<<"	"<<"PairPacker"<<"	"<<"exon"<<"	"<<Node_seq_2[Single_nodes[i]]<<"	"<<Node_seq_3[Single_nodes[i]]<<"	"<<1000<<"	"<<Symble[strand_id]<<"	.	gene_id "<<"\""<<Trans_name<<"."<<Graph_num<<"\""<<"; "<<"transcript_id "<<"\""<<Trans_name<<"."<<Graph_num<<"."<<Curr_id<<"\""<<"; "<<"exon_number "<<"\""<<1<<"\""<<";"<<endl;
			Curr_id++;
		}
	}

	file.close();
	return;
}// end of PairPacker_Output

#endif


