#ifndef GET_BOUND
#define GET_BOUND

//get_genes_bound.h

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include <algorithm>

using namespace std;
class Get_gene_bound
{
	private:
	const char * path_gtf; 
	public:
	void Process();
	Get_gene_bound(const char * p);
	vector<vector<int> > Gene_bound;
	vector<int> Gene_chr_num;
	vector<string> Genome_chr;
	void variant_clear();
};
Get_gene_bound::Get_gene_bound(const char * p)
{
	path_gtf=p;
	return;
}
void Get_gene_bound::variant_clear()
{
		Gene_bound.clear();
		Gene_chr_num.clear();
		Genome_chr.clear();
}
void Get_gene_bound::Process()
{
istringstream istr;
ifstream ifs(path_gtf);
int left_pos, right_pos;
int gene_min=1000000000;
int gene_max=1;
string gene_id_old="no";
string chr_old="no";
string strand_old="no";
vector<int> gene_left,gene_right;
string temp, chr_genome, gene_id, chr, strand;
getline(ifs,temp);
istr.str(temp);
istr>>chr>>temp>>temp>>left_pos>>right_pos>>temp>>strand>>temp>>temp>>gene_id;
istr.clear();
vector<string> gene_chr, gene_strand;
gene_min=left_pos;
gene_max=right_pos;
gene_id_old=gene_id;
chr_old=chr;
strand_old=strand;
while (getline(ifs,temp))
{
        istr.str(temp);
        istr>>chr>>temp>>temp>>left_pos>>right_pos>>temp>>strand>>temp>>temp>>gene_id;
        istr.clear();
        if (gene_id == gene_id_old)
        {
                if (left_pos<gene_min)
                {
                        gene_min=left_pos;
                }
                if (right_pos>gene_max)
                {
                        gene_max=right_pos;
                }
        }
        else
        {
                gene_chr.push_back(chr_old);
                gene_left.push_back(gene_min);
                gene_right.push_back(gene_max);
                gene_strand.push_back(strand_old);
                gene_min=left_pos;
                gene_max=right_pos;
                gene_id_old=gene_id;
                chr_old=chr;
                strand_old=strand;
        }
}
ifs.close();
gene_chr.push_back(chr);
gene_left.push_back(gene_min);
gene_right.push_back(gene_max);
gene_strand.push_back(strand);
string curr_gene_chr=gene_chr[0];
string curr_gene_strand=gene_strand[0];
int curr_gene_left=gene_left[0];
int curr_gene_right=gene_right[0];
vector<int> gene_bound_temp;
int t=0;
Gene_chr_num.push_back(t);
for (size_t i=1;i<gene_left.size();i++)
{
        if (gene_chr[i]==curr_gene_chr && gene_strand[i]==curr_gene_strand && gene_left[i]<=curr_gene_right)
        {
                if (gene_right[i]>curr_gene_right) curr_gene_right=gene_right[i];
        }
        else //(gene_left[i]>curr_gene_right)
        {
                t++;
                if (gene_chr[i]!=curr_gene_chr)
                {
                        Gene_chr_num.push_back(t);
			Genome_chr.push_back(curr_gene_chr);
                }
                gene_bound_temp.push_back(curr_gene_left);
                gene_bound_temp.push_back(curr_gene_right);
                Gene_bound.push_back(gene_bound_temp);
                gene_bound_temp.clear();

                curr_gene_chr=gene_chr[i];
                curr_gene_strand=gene_strand[i];
                curr_gene_left=gene_left[i];
                curr_gene_right=gene_right[i];
        }
}
Gene_chr_num.push_back(t+1);
Genome_chr.push_back(curr_gene_chr);
gene_bound_temp.push_back(curr_gene_left);
gene_bound_temp.push_back(curr_gene_right);
Gene_bound.push_back(gene_bound_temp);
gene_bound_temp.clear();
}// end of Process() for gene bound

class Get_exon_bound
{
    private:
    const char * path_gtf;
	const char * path_genome;
	const char * path_gene;
	const char * path_chr;
	vector<vector<int> > Gene_bound;
	vector<int> Gene_chr_num;
	vector<string> Genome_chr;
    public:
	ofstream file_gene;
	ofstream file_chr;
	void variant_clear();
	vector<vector<int> > node_id_on_pesudo(int curr_pos, vector<vector<int> >& exons);
	void get_chr_from_genome();
	void get_genome_seq(string chr, string strand, vector<vector<int> >& exons);
	vector<vector<int> > exon_split(vector<vector<int> >& exons);
    void Process();
    Get_exon_bound(const char * p1,const char *p2,const char *p3,const char *p4, vector<vector<int> >& gene_bound, vector<int>& gene_chr_num, vector<string>& genome_chr);
    vector<vector<int> > Exon_bound;
    vector<int> Exon_chr_num;
	vector<string> chrs;
	vector<vector<string> > chr_genome_vec;
    int length;
};
Get_exon_bound::Get_exon_bound(const char * p1,const char *p2, const char *p3,const char *p4,vector<vector<int> >& gene_bound, vector<int>& gene_chr_num, vector<string>& genome_chr)
{
    path_gtf=p1;
	path_genome=p2;
	path_gene=p3;
	path_chr=p4;
	file_gene.open(path_gene);
	file_chr.open(path_chr);
	Gene_bound=gene_bound;
	Gene_chr_num=gene_chr_num;
	Genome_chr=genome_chr; 
       return;
}
void Get_exon_bound::variant_clear()
{
		Gene_bound.clear();
		Gene_chr_num.clear();
		Genome_chr.clear();
		Exon_bound.clear();
		Exon_chr_num.clear();
		chrs.clear();
		chr_genome_vec.clear();
}
vector<vector<int> > Get_exon_bound::node_id_on_pesudo(int curr_pos, vector<vector<int> >& exons)
{
// The computation is dependent on the number of "NNN...".
	vector<vector<int> > output_vec;
	vector<int> vec_temp;
	int pos=curr_pos;
	int pos_temp;
	for (size_t i=0;i<exons.size();i++)
	{
		vec_temp.push_back(pos);
		pos_temp=exons[i][1]-exons[i][0]+pos;
		vec_temp.push_back(pos_temp);
		output_vec.push_back(vec_temp);
		vec_temp.clear();
		//pos=pos_temp+27;
		if (i<exons.size()-1 && exons[i][1]+1!=exons[i+1][0])
		{
			pos=pos_temp+27;
		}
		else
		{
			pos=pos_temp+1;
		}
	}
	return output_vec;
}
void Get_exon_bound::get_chr_from_genome()
{
	ifstream ifs_1(path_genome);
	string temp;
        vector<string> chr_genome_vec_tmp;
  //      int t=0;
        getline(ifs_1,temp);
        chrs.push_back(temp.substr(1,temp.size()-1));
        getline(ifs_1,temp);
        chr_genome_vec_tmp.push_back(temp);
        length=temp.length();
        while (getline(ifs_1,temp))
        {
                if (temp[0]=='>')
                {
                        chr_genome_vec.push_back(chr_genome_vec_tmp);
                        chr_genome_vec_tmp.clear();
                        chrs.push_back(temp.substr(1,temp.size()-1));
                }
                else
                {
                        chr_genome_vec_tmp.push_back(temp);
                }
        }
        chr_genome_vec.push_back(chr_genome_vec_tmp);
        chr_genome_vec_tmp.clear();
}
void Get_exon_bound::get_genome_seq(string chr, string strand, vector<vector<int> >& exons)
{
        int chr_id=-1;
	for (size_t ii=0;ii<chrs.size();ii++)
        {
        	if (chrs[ii] == chr)
                {
			chr_id=ii;
			break;
		}
	}
	for (size_t i=0;i<exons.size();i++)
	{
		int num_1=exons[i][0]/length;
		int num_2=exons[i][0]%length;
		vector<int> vec_1,vec_2;
		if (num_2 == 0)
                {
                        vec_1.push_back(num_1-1);
                        vec_1.push_back(length-1);
                }
                else
                {
                        vec_1.push_back(num_1);
                        vec_1.push_back(num_2-1);
                }
		num_1=exons[i][1]/length;
                num_2=exons[i][1]%length;
                if (num_2 == 0)
                {
                        vec_2.push_back(num_1-1);
                        vec_2.push_back(length-1);
                }
                else
                {
                        vec_2.push_back(num_1);
                        vec_2.push_back(num_2-1);
                }
		if (vec_2[0]>vec_1[0])
		{
                	file_chr<<chr_genome_vec[chr_id][vec_1[0]].substr(vec_1[1],length-vec_1[1]);
                	for (int j=vec_1[0]+1;j<vec_2[0];j++)
                	{
                        	file_chr<<chr_genome_vec[chr_id][j];
                	}
			file_chr<<chr_genome_vec[chr_id][vec_2[0]].substr(0,1+vec_2[1]);
		}
		else
		{
			file_chr<<chr_genome_vec[chr_id][vec_1[0]].substr(vec_1[1],vec_2[1]-vec_1[1]+1);
		}
		if (i<exons.size()-1 && exons[i][1]+1==exons[i+1][0])
			continue;
		if (strand == "+")
		{
			file_chr<<"GTNNNNNNNNNNNNNNNNNNNNNNAG";
		}
		else
		{
			file_chr<<"CTNNNNNNNNNNNNNNNNNNNNNNAC";
		}
	}//for (int i=0;i<exons.size();i++)
return;
}//get_genome_seq()
vector<vector<int> > Get_exon_bound::exon_split(vector<vector<int> >& exons)
{

	vector<vector<int> > output_vec;
	vector<int> output_vec_temp;
	size_t i,j;
	vector<int> curr_exon=exons[0];
	vector<int> split_pos_L;
	vector<int> split_pos_R;
	vector<int> split_pos;
	vector<int>::iterator iter,pos_L,pos_R;
	split_pos_L.push_back(curr_exon[0]);
	split_pos_R.push_back(curr_exon[1]);
	for (i=1;i<exons.size();i++)
	{
		if (exons[i][0]>curr_exon[1])
		{
			sort(split_pos_L.begin(),split_pos_L.end());
			iter=unique(split_pos_L.begin(),split_pos_L.end());
			split_pos_L.erase(iter,split_pos_L.end());
			sort(split_pos_R.begin(),split_pos_R.end());
			iter=unique(split_pos_R.begin(),split_pos_R.end());
			split_pos_R.erase(iter,split_pos_R.end());
			split_pos=split_pos_L;
			for (j=0;j<split_pos_R.size();j++)
			{
				split_pos.push_back(split_pos_R[j]);
			}
			sort(split_pos.begin(),split_pos.end());
	/*
			for (j=0;j<split_pos.size();j++)
                        {
				cout<<split_pos[j]<<"  ";
			}
			cout<<endl;
*/
			for (j=0;j<split_pos.size()-1;j++)
			{
			    
				pos_L=find(split_pos_L.begin(),split_pos_L.end(),split_pos[j]);
				pos_R=find(split_pos_R.begin(),split_pos_R.end(),split_pos[j]);
				if (pos_L != split_pos_L.end() && pos_R != split_pos_R.end())
				{
					output_vec_temp.push_back(split_pos[j]);
					output_vec_temp.push_back(split_pos[j]);
					output_vec.push_back(output_vec_temp);
					output_vec_temp.clear();
					output_vec_temp.push_back(split_pos[j]+1);
				}
				else if (pos_L != split_pos_L.end())
				{
					output_vec_temp.push_back(split_pos[j]);
				}
				else
				{
					output_vec_temp.push_back(split_pos[j]+1);
				}
				pos_L=find(split_pos_L.begin(),split_pos_L.end(),split_pos[j+1]);
				if (pos_L != split_pos_L.end())
				{
					output_vec_temp.push_back(split_pos[j+1]-1);
				}
				else
				{
					output_vec_temp.push_back(split_pos[j+1]);
				}
				if (output_vec_temp[1]>=output_vec_temp[0])
				{
					output_vec.push_back(output_vec_temp);
				}
				output_vec_temp.clear();
			}
			curr_exon=exons[i];
			split_pos_L.clear();
			split_pos_R.clear();
			split_pos_L.push_back(curr_exon[0]);
			split_pos_R.push_back(curr_exon[1]);
		}
		else
		{
			split_pos_L.push_back(exons[i][0]);
			split_pos_R.push_back(exons[i][1]);
			if (exons[i][1]>curr_exon[1])
			{
				curr_exon[1]=exons[i][1];
			}
		}
	}
			sort(split_pos_L.begin(),split_pos_L.end());
			iter=unique(split_pos_L.begin(),split_pos_L.end());
			split_pos_L.erase(iter,split_pos_L.end());
			sort(split_pos_R.begin(),split_pos_R.end());
			iter=unique(split_pos_R.begin(),split_pos_R.end());
			split_pos_R.erase(iter,split_pos_R.end());
			split_pos=split_pos_L;
			for (j=0;j<split_pos_R.size();j++)
			{
				split_pos.push_back(split_pos_R[j]);
			}
			sort(split_pos.begin(),split_pos.end());
			for (j=0;j<split_pos.size()-1;j++)
			{
			    
				pos_L=find(split_pos_L.begin(),split_pos_L.end(),split_pos[j]);
				pos_R=find(split_pos_R.begin(),split_pos_R.end(),split_pos[j]);
				if (pos_L != split_pos_L.end() && pos_R != split_pos_R.end())
				{
					output_vec_temp.push_back(split_pos[j]);
					output_vec_temp.push_back(split_pos[j]);
					output_vec.push_back(output_vec_temp);
					output_vec_temp.clear();
					output_vec_temp.push_back(split_pos[j]+1);
				}
				else if (pos_L != split_pos_L.end())
				{
					output_vec_temp.push_back(split_pos[j]);
				}
				else
				{
					output_vec_temp.push_back(split_pos[j]+1);
				}
				pos_L=find(split_pos_L.begin(),split_pos_L.end(),split_pos[j+1]);
				if (pos_L != split_pos_L.end())
				{
					output_vec_temp.push_back(split_pos[j+1]-1);
				}
				else
				{
					output_vec_temp.push_back(split_pos[j+1]);
				}
				if (output_vec_temp[1]>=output_vec_temp[0])
				{
					output_vec.push_back(output_vec_temp);
				}
				output_vec_temp.clear();
			}
	vector<vector<int> >::iterator iter_output;
	sort(output_vec.begin(),output_vec.end());
	iter_output=unique(output_vec.begin(),output_vec.end());
	output_vec.erase(iter_output,output_vec.end());
	return output_vec;
}//exon_split(vector<vector<int> >& exons)
 
void Get_exon_bound::Process()
{
size_t i;
int j,chr_num;
ifstream ifs(path_gtf);
istringstream istr;
string temp, chr_genome, gene_id, chr, strand, exon, trans_id;
string curr_chr,curr_strand;
int left_pos, right_pos;
chr_num=1;
getline(ifs,temp);
istr.str(temp);
istr>>chr>>temp>>exon>>left_pos>>right_pos>>temp>>strand>>temp>>temp>>gene_id>>temp>>trans_id;
istr.clear();
for (i=0;i<Gene_chr_num.size()-1;i++)
{
    curr_chr=Genome_chr[i];
    vector<string>::iterator pos_chr;
    pos_chr=find(chrs.begin(),chrs.end(),curr_chr);
    if (pos_chr == chrs.end())
    {
	//cout<<"Error: cannot detect \""<<curr_chr<<"\" in reference genome."<<endl;
	return;
    }
    curr_strand=strand;
    file_chr<<">"<<curr_chr<<curr_strand<<endl;
    int curr_pos=1;
    for (j=Gene_chr_num[i];j<Gene_chr_num[i+1];j++)
    {
	vector<vector<int> > exon_vec_raw, exon_vec_split;
	vector<string> exon_trans_id;
	vector<vector<int> >::iterator iter;
	vector<int> exon_vec_temp;
	if (exon == "exon" && right_pos<=Gene_bound[j][1])
	{
		exon_vec_temp.push_back(left_pos);
		exon_vec_temp.push_back(right_pos);
		exon_vec_raw.push_back(exon_vec_temp);
		exon_trans_id.push_back(trans_id);
		exon_vec_temp.clear();
	}
	while (getline(ifs,temp))
	{
		istr.str(temp);
        	istr>>chr>>temp>>exon>>left_pos>>right_pos>>temp>>strand>>temp>>temp>>gene_id>>temp>>trans_id;
        	istr.clear();
		if (curr_chr == chr && exon == "exon" && right_pos<=Gene_bound[j][1])
		{
			exon_vec_temp.push_back(left_pos);
			exon_vec_temp.push_back(right_pos);
			exon_vec_raw.push_back(exon_vec_temp);
			exon_trans_id.push_back(trans_id);
			exon_vec_temp.clear();
		}
		else if (left_pos>Gene_bound[j][1] || curr_chr != chr)
		{
			break;
		}
	} //while (getline(ifs,temp))
	vector<vector<int> > exon_vec_search_edge=exon_vec_raw;
	sort(exon_vec_raw.begin(),exon_vec_raw.end());
	iter=unique(exon_vec_raw.begin(),exon_vec_raw.end());
	exon_vec_raw.erase(iter,exon_vec_raw.end());
	exon_vec_split=exon_split(exon_vec_raw);
// Searching for edges...
	vector<int> Edges_L,Edges_R;
	vector<vector<int> > Paths;
	vector<int> path_temp;
	int add_id=1;
	for (size_t k=0;k<exon_vec_search_edge.size()-1;k++)
	{
		if (exon_trans_id[k]==exon_trans_id[k+1])
		{
			int temp_edge_search_L=exon_vec_search_edge[k][1];
			int temp_edge_search_R=exon_vec_search_edge[k+1][0];
			int tmp_left=0;
			int tmp_right=0;
			for (size_t m=0;m<exon_vec_split.size();m++)
			{
				if (tmp_left == 1 && tmp_right ==1)
				{
					break;
				}
				if (exon_vec_split[m][1]==temp_edge_search_L)
				{
					tmp_left=1;
					Edges_L.push_back(m);
				    if(add_id == 1)
				    {
					add_id=0;
					for (size_t t=m;t>=0;t--)
					{
						if (exon_vec_split[t][0]==exon_vec_search_edge[k][0])
						{
							for (size_t n=t;n<=m;n++)
							{
								path_temp.push_back(n);
							}
							break;
						}
					}
				    }
				}
				else if (exon_vec_split[m][0]==temp_edge_search_R)
				{
					tmp_right=1;
					Edges_R.push_back(m);
					for (size_t t=m;t<exon_vec_split.size();t++)
					{
						if (exon_vec_split[t][1]==exon_vec_search_edge[k+1][1])
						{
							for (size_t n=m;n<=t;n++)
							{
								path_temp.push_back(n);
							}
							break;
						}
					}
				}
			}
		} //if (exon_trans_id[k]==exon_trans_id[k+1])
		else
		{
			if (path_temp.size()>0)
			{
				Paths.push_back(path_temp);
				path_temp.clear();
				add_id=1;
			}
			else
			{
				int temp_edge_search_L=exon_vec_search_edge[k][0];
				int temp_edge_search_R=exon_vec_search_edge[k][1];
				int tmp_left=-1;
				int tmp_right=-1;
				for (size_t m=0;m<exon_vec_split.size();m++)
				{
					if (tmp_left != -1 && tmp_right != -1)
					{
						break;
					}
					if (exon_vec_split[m][0]==temp_edge_search_L)
					{
						tmp_left=m;
					}
					if (exon_vec_split[m][1]==temp_edge_search_R)
					{
						tmp_right=m;
					}
				}
				for (size_t m=tmp_left;m<=tmp_right;m++)
				{
					path_temp.push_back(m);
				}
				Paths.push_back(path_temp);
				path_temp.clear();
			}
		}
	} //for (int k=0;k<exon_vec_search_edge.size()-1;k++)
	if (path_temp.size()>0)
	{
		Paths.push_back(path_temp);
		path_temp.clear();
	}
	else
	{
		int k=exon_vec_search_edge.size()-1;
		int temp_edge_search_L=exon_vec_search_edge[k][0];
		int temp_edge_search_R=exon_vec_search_edge[k][1];
		int tmp_left=-1;
		int tmp_right=-1;
		for (size_t m=0;m<exon_vec_split.size();m++)
		{
			if (tmp_left != -1 && tmp_right != -1)
			{
				break;
			}
			if (exon_vec_split[m][0]==temp_edge_search_L)
			{
				tmp_left=m;
			}
			if (exon_vec_split[m][1]==temp_edge_search_R)
			{
				tmp_right=m;
			}
		}
		for (size_t m=tmp_left;m<=tmp_right;m++)
		{
			path_temp.push_back(m);
		}
		Paths.push_back(path_temp);
		path_temp.clear();
	}
	for (size_t k=0;k<exon_vec_split.size()-1;k++)
	{
		if (exon_vec_split[k][1]+1 == exon_vec_split[k+1][0])
		{
			Edges_L.push_back(k);
			Edges_R.push_back(k+1);
		}
	}
	vector<vector<int> > Edge_L_R;
	vector<vector<int> >::iterator iter_edge;
	vector<int> temp_edge;
	for (size_t k=0;k<Edges_R.size();k++)
	{
		temp_edge.push_back(Edges_L[k]);
		temp_edge.push_back(Edges_R[k]);
		Edge_L_R.push_back(temp_edge);
		temp_edge.clear();
	}
	sort(Edge_L_R.begin(),Edge_L_R.end());
	iter_edge=unique(Edge_L_R.begin(),Edge_L_R.end());
	Edge_L_R.erase(iter_edge,Edge_L_R.end());
	vector<vector<int> > Node_id=node_id_on_pesudo(curr_pos,exon_vec_split);
	curr_pos=Node_id[Node_id.size()-1][1]+27;
	file_gene<<"** Chr **"<<endl;
	file_gene<<curr_chr<<curr_strand<<endl;
	file_gene<<"** Nodes **"<<endl;
	for (size_t k=0;k<exon_vec_split.size();k++)
	{
		file_gene<<exon_vec_split[k][0]<<" "<<exon_vec_split[k][1]<<" "<<Node_id[k][0]<<" "<<Node_id[k][1]<<endl;
	}
	file_gene<<"** Edges **"<<endl;
	for (size_t k=0;k<Edge_L_R.size();k++)
	{
		file_gene<<Edge_L_R[k][0]<<" "<<Edge_L_R[k][1]<<endl;
	}
	file_gene<<"** Paths **"<<endl;
	for (size_t k=0;k<Paths.size();k++)
	{
		for (size_t t=0;t<Paths[k].size();t++)
		{
			file_gene<<Paths[k][t]<<" ";
		}
		file_gene<<endl;
	}
	file_gene<<"** End **"<<endl;
	get_genome_seq(curr_chr,curr_strand,exon_vec_split);

    } //for (j=Gene_chr_num[i];j<Gene_chr_num[i+1];j++)
    file_chr<<endl;
} //for (i=0;i<Gene_chr_num.size();i++)
ifs.clear();
}//Process()
#endif




//get_genome_seq(string chr, string strand, vector<vector<int> >& exons)













