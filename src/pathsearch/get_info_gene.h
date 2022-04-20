#ifndef GET_INFO_GENE_H
#define GET_INFO_GENE_H
// get_info_gene.h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>

using namespace std;
typedef struct struct_1
{
	string chr;
	int pos_left;
	int pos_right;
} Hash_key;
class Gene_info
{
private:
  const char * path_of_gene_info;
  const char * path_of_bam;
public:
  Gene_info(const char * p1, const char *p2, double cov_uniq);
  vector<string> Gene_chrs;
  vector<vector<vector<double> > > Gene_nodes;
  vector<vector<vector<double> > > Gene_edges;
  vector<vector<vector<int> > > Gene_edges_partial;
  vector<vector<vector<int> > > Gene_paths;
  map<vector<double>, vector<double> > Hash_edges;
  map<vector<double>, vector<double> > Hash_left, Hash_right;
  vector<vector<vector<int> > > Gene_pair_paths;
  vector<int> Chr_int;
  vector<string> Chr_raw;
  vector<vector<vector<int> > > Pair_vec_multi; //[chr_id; left_start; left_num; left_multi; right_start; right_num; right_multi]
  vector<vector<vector<int> > > Pair_vec_junc; //[gene_id; * * * -1 * * *], where * means the node numbers
  vector<vector<vector<int> > > Map_need_to_do; // these reads may belong to partial exons
  double Cov_uniq;
  void variant_clear_1();
  void variant_clear_2();
public:
  void Generate_genes();
  int Identify_cigar(string cigar);
  vector<vector<int> > Identify_cigar(int start_pos, string cigar);
  int Check_flag(int flag);
  void Update_genes_from_sam();
  void process_map_need_to_do();
  void process_pair_vec_junc();
  void process_pair_vec_multi();
  vector<vector<int> > search_paths(int node_start, int node_end, int gene_id);
  void process_gene_paths();
  vector<vector<int> > enumerate_pair_paths(vector<int> Pair_edge_left,vector<int> Pair_edge_right,vector<int> Pair_edge_pos,int max);
  void show_genes();
  void show_hash();
  void show_hash_add();
  void show_cigar();
  void show_map_need_to_do();
  void show_pair_vec_multi();
  void show_pair_vec_junc();
  void show_gene_pair_paths();
  void test_pair_path(char * path);
  void test();
  void multi_clear();
  void junc_clear();
  void map_clear();
};
void Gene_info::multi_clear()
{
		Pair_vec_multi.clear();
}
void Gene_info::junc_clear()
{
		Pair_vec_junc.clear();
}
void Gene_info::map_clear()
{
		Map_need_to_do.clear();
}
void Gene_info::variant_clear_1()
{
		Hash_left.clear();
		Hash_right.clear();
}
void Gene_info::variant_clear_2()
{
		Hash_edges.clear();
		Gene_edges_partial.clear();
		Gene_nodes.clear();
		Gene_edges.clear();
		Gene_paths.clear();
		Gene_pair_paths.clear();
}

void Gene_info::test()
{
		vector<vector<int> > a=search_paths(8,10,11);
		for (size_t i=0;i<a.size();i++)
		{
				for (size_t j=0;j<a[i].size();j++)
				{
						cout<<a[i][j]<<" ";
				}
				cout<<endl;
		}
}
Gene_info::Gene_info(const char * p1, const char *p2, double cov_uniq)
{
	path_of_gene_info=p1;
	path_of_bam=p2;
	Cov_uniq=cov_uniq;
	return;
}
void Gene_info::Generate_genes()
{
	ifstream ifs(path_of_gene_info);
        string temp,curr_chr;
        int val_1,val_2,val_3,val_4;
        vector<double> curr_node_temp;
	vector<double> curr_edge_temp;
	vector<int> curr_edge_partial_temp;
	vector<int> curr_path_temp;
        vector<vector<double> > curr_node;
	vector<vector<double> > curr_edge;
	vector<vector<int> > curr_edge_partial;
	vector<vector<int> > curr_path;
        istringstream istr;
	double gene_id=0;
	double chr_id=0;
	vector<vector<int> > pair_temp;
        while (getline(ifs,temp))
        {
                if (temp == "** Chr **")
                {
			Pair_vec_junc.push_back(pair_temp);
			Gene_pair_paths.push_back(pair_temp);
                        getline(ifs,temp);
                        curr_chr=temp;
			Gene_chrs.push_back(curr_chr);
			if (Chr_raw.size() == 0)
			{
				Chr_raw.push_back(curr_chr);
				Chr_int.push_back(chr_id);
				Pair_vec_multi.push_back(pair_temp);
				Map_need_to_do.push_back(pair_temp);
			}
			else if (curr_chr != Chr_raw[Chr_raw.size()-1])
			{
				Chr_raw.push_back(curr_chr);
				chr_id++;
				Chr_int.push_back(chr_id);
				Pair_vec_multi.push_back(pair_temp);
				Map_need_to_do.push_back(pair_temp);
			}
                }
                if (temp == "** Nodes **")
                {
                        curr_node.clear();
			vector<double> key_left,key_right,value_left,value_right;
                        while (1)
                        {
                                getline(ifs,temp);
                                if (temp == "** Edges **")
                                        break;
                                istr.str(temp);
                                istr>>val_1>>val_2>>val_3>>val_4;
                                istr.clear();
                                curr_node_temp.push_back(val_1);
                                curr_node_temp.push_back(val_2);
                                curr_node_temp.push_back(val_3);
                                curr_node_temp.push_back(val_4);
                                curr_node_temp.push_back(0);
                                curr_node.push_back(curr_node_temp);
                                curr_node_temp.clear();
                        }
			Gene_nodes.push_back(curr_node);
			for (size_t i=0;i<curr_node.size()-1;i++)
			{
				key_left.push_back(chr_id);
				key_left.push_back(curr_node[i][3]);
				value_left.push_back(gene_id);
				value_left.push_back(i);
				Hash_left[key_left]=value_left;
				key_left.clear();
				value_left.clear();
			}
			for (size_t i=1;i<curr_node.size();i++)
			{
				key_right.push_back(chr_id);
                                key_right.push_back(curr_node[i][2]);
                                value_right.push_back(gene_id);
                                value_right.push_back(i);
                                Hash_right[key_right]=value_right;
				key_right.clear();
				value_right.clear();
			}
                }
                if (temp == "** Edges **")
                {
			double VAL_1, VAL_2;
			int edge_id=0;
			curr_edge.clear();
			curr_edge_partial.clear();
			vector<double> key,value;
			while (1)
			{
				value.clear();
				getline(ifs,temp);
				if (temp == "** Paths **")
					break;
				istr.str(temp);
				istr>>VAL_1>>VAL_2;
				istr.clear();
				curr_edge_temp.push_back(VAL_1);
				curr_edge_temp.push_back(VAL_2);
				curr_edge_temp.push_back(0);
				curr_edge.push_back(curr_edge_temp);
				curr_edge_temp.clear();
				// adding to hash...
				if (curr_node[VAL_2][2]-curr_node[VAL_1][3] == 1)
				{
					curr_edge_partial_temp.push_back(int(VAL_1));
					curr_edge_partial_temp.push_back(int(VAL_2));
					curr_edge_partial.push_back(curr_edge_partial_temp);
					curr_edge_partial_temp.clear();
				}
				key.push_back(chr_id);
				key.push_back(curr_node[VAL_1][3]);
				key.push_back(curr_node[VAL_2][2]);
				value.push_back(gene_id);
				value.push_back(edge_id);
				if (Hash_edges.find(key) != Hash_edges.end())
				{
					cout<<"Error: "<<curr_chr<<" "<<value[0]<<" "<<value[1]<<" "<<value[2]<<endl;
                                        return;
				}
				Hash_edges[key]=value;
				key.clear();
				value.clear();
				edge_id++;
			}
			if (curr_node.size() == 1)
			{
				curr_edge_temp.push_back(0);
				curr_edge.push_back(curr_edge_temp);
				curr_edge_temp.clear();
			}
			Gene_edges.push_back(curr_edge);
			Gene_edges_partial.push_back(curr_edge_partial);
		}
		if (temp == "** Paths **")
		{
			curr_path.clear();
			while (1)
                        {
                                getline(ifs,temp);
                                if (temp == "** End **")
                                        break;
				int path_start_pos=0;
				for (size_t i=0;i<temp.size();i++)
				{
					if (temp[i]==' ')
					{
						curr_path_temp.push_back(atoi(temp.substr(path_start_pos,i-path_start_pos).c_str()));
						path_start_pos=i+1;
					}
				}
				curr_path.push_back(curr_path_temp);
				curr_path_temp.clear();
			}
			Gene_paths.push_back(curr_path);
			gene_id++;
		}
        }//while (getline(ifs,temp))
        ifs.close();
	return;
}//Generate_genes()
void Gene_info::show_cigar()
{
	string x="62M19801N8M5S";
	vector<vector<int> > y = Identify_cigar(12,x);
	for (size_t i=0;i<y.size();i++)
	{
		cout<<y[i][0]<<"--"<<y[i][1]<<endl;
	}
	x="8M2I14M1D3M";
	cout<<Identify_cigar(x)<<endl;
}
int Gene_info::Identify_cigar(string cigar)
{
	int sum_M=0;
	int curr_num=0;
	for (size_t i=0;i<cigar.size();i++)
	{
		if (cigar[i] == 'M' || cigar[i] == 'D')
		{
			sum_M=sum_M+atoi(cigar.substr(curr_num,i-curr_num).c_str());
			curr_num=i+1;
		}
		else if (cigar[i] == 'I' || cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P')
		{
			curr_num=i+1;
		}
	}
	return sum_M;
}

vector<vector<int> > Gene_info::Identify_cigar(int start_pos, string cigar)
{
	int curr_pos_left=start_pos;
	int curr_pos_right=0;
	int curr_num=0;
	int module_id=0;
	vector<int> Pos_left,Pos_right,Module_N;
	for (size_t i=0;i<cigar.size();i++)
	{
		if (cigar[i] == 'M' || cigar[i] == 'D' || cigar[i] == 'N')
		{
			curr_pos_right=atoi(cigar.substr(curr_num,i-curr_num).c_str())+curr_pos_left-1;
			Pos_left.push_back(curr_pos_left);
			Pos_right.push_back(curr_pos_right);
//			cout<<curr_pos_left<<"  "<<curr_pos_right<<endl;
			if (cigar[i] == 'N') Module_N.push_back(module_id);
			module_id++;
			curr_pos_left=curr_pos_right+1;
			curr_num=i+1;
		}
		else if (cigar[i] == 'I' || cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P')
		{
			curr_num=i+1;
		}
	}
	vector<vector<int> > output;
	vector<int> output_temp;
	for (size_t i=0;i<Module_N.size();i++)
	{
		output_temp.push_back(Pos_right[Module_N[i]-1]);
		output_temp.push_back(Pos_left[Module_N[i]+1]);
		output.push_back(output_temp);
		output_temp.clear();
	}
	return output;
}
int Gene_info::Check_flag(int flag)
{
	int is_pair=0;
	int temp=flag;
    	int num;
    	vector<int> binary;
    	while(temp!=0)
	{
        	num=temp%2;
        	binary.push_back(num);
        	temp=temp/2;
    	}
    	for(int i=binary.size();i<11;i++) binary.push_back(0);
	if (binary[1] == 1) is_pair=1;
	return is_pair;
}
void Gene_info::Update_genes_from_sam()
{
	string temp,read_id;
	string read_id_1,chr_id_1,map_qual_1,map_cigar_1,equal_star_1;
	string read_id_2,chr_id_2,map_qual_2,map_cigar_2,equal_star_2;
	int flag_1, map_start_1, pair_start_1, map_length_1;
	int flag_2, map_start_2, pair_start_2, map_length_2;
	ifstream ifs(path_of_bam);
	istringstream istr;
	vector<int> pair_vec, pair_vec_L, pair_vec_R;
	long long line=0;
	int N_L=0;
	int N_R=0;
	while (getline(ifs,temp))
	{
		line++;
		if (temp.substr(0,3) == "@PG")
			break;
	}
	vector<string> sam_line;
	getline(ifs,temp);
	line++;
	string curr_line=temp;
	sam_line.push_back(curr_line);
	istr.str(temp);
	istr>>read_id;
	istr.clear();
	string curr_read=read_id;
	while (getline(ifs,temp))
	{
		line++;
		if (line%100000==0) cerr<<line<<endl;
		sam_line.push_back(temp);
		istr.str(temp);
		istr>>read_id;
		istr.clear();
		curr_read=read_id;
		while (getline(ifs,temp))
		{
			line++;
			if (line%100000==0) cerr<<line<<endl;
			istr.str(temp);
                        istr>>read_id;
                        istr.clear();
			if (curr_read != read_id)
			{
				curr_line=temp;
				break;
			}
			sam_line.push_back(temp);
		}
		int pair_num=int(sam_line.size())/2;
		vector<int> read_flag_1, read_flag_2; //[x y], where x=0: 1st read, x=1: 2nd read; y=0: forward, y=2: reverse.
		double multi_1, multi_2;
		for (int i=0;i<pair_num;i++)
		{
			istr.str(sam_line[2*i]);
			istr>>read_id_1>>flag_1>>chr_id_1>>map_start_1>>map_qual_1>>map_cigar_1>>equal_star_1>>pair_start_1>>map_length_1;
			multi_1=0;
			while (istr>>temp)
                	{
                        	if (map_cigar_1 == "*" || multi_1>0)
                        	{
                                	break;
                        	}
                        	if (temp.substr(0,5)=="NH:i:")
                        	{
                                	multi_1=atof(temp.substr(5,temp.size()-5).c_str());
                        	}
                	}
			istr.clear();
			istr.str(sam_line[2*i+1]);
                        istr>>read_id_2>>flag_2>>chr_id_2>>map_start_2>>map_qual_2>>map_cigar_2>>equal_star_2>>pair_start_2>>map_length_2;
                        multi_2=0;
                        while (istr>>temp)
                        {
                                if (map_cigar_2 == "*" || multi_2>0)
                                {
                                        break;
                                }
                                if (temp.substr(0,5)=="NH:i:")
                                {
                                        multi_2=atof(temp.substr(5,temp.size()-5).c_str());
                                }
                        }
                        istr.clear();
			int N_count_1=count(map_cigar_1.begin(),map_cigar_1.end(),'N');
			int N_count_2=count(map_cigar_2.begin(),map_cigar_2.end(),'N');
			pair_vec_L.clear();
			pair_vec_R.clear();
			if (N_count_1 == 0 && multi_1 >=1)
                	{
                        	N_L=0;
                        	for (size_t j=0;j<Chr_raw.size();j++)
                       	 	{
                                	if (Chr_raw[j]==chr_id_1)
                                	{
                                        	pair_vec_L.push_back(j);
                                        	break;
                                	}
                        	}
                        	pair_vec_L.push_back(map_start_1);
                        	pair_vec_L.push_back(Identify_cigar(map_cigar_1));
                        	pair_vec_L.push_back(multi_1);
                        }
			else if (N_count_1>0)
			{
				N_L=1;
                        	vector<vector<int> > Junc_pos = Identify_cigar(map_start_1,map_cigar_1);
                        	for (size_t i=0;i<Junc_pos.size();i++)
                        	{
                                	vector<double> key,key_left,key_right;
                                	for (size_t j=0;j<Chr_raw.size();j++)
                               		{
                                        	if (Chr_raw[j]==chr_id_1)
                                        	{
                                                	key.push_back(j);
							key_left.push_back(j);
							key_right.push_back(j);
                                                	break;
                                       		}
                                	}
                                	key.push_back(Junc_pos[i][0]);
                                	key.push_back(Junc_pos[i][1]);
                                	if (Hash_edges.find(key) != Hash_edges.end())
					{
                                		vector<double> value = Hash_edges[key];
						if (pair_vec_L.size() == 0) 
						{
							pair_vec_L.push_back(int(value[0])); // record gene_id;
						}
						Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi_1;
						if (pair_vec_L.size() == 1)
						{
							pair_vec_L.push_back(int(Gene_edges[value[0]][value[1]][0]));
						}
						pair_vec_L.push_back(int(Gene_edges[value[0]][value[1]][1]));
						if (value[0] != pair_vec_L[0])
						{
							pair_vec_L.clear();
						}
					}
					else
					{
						key_left.push_back(Junc_pos[i][0]);
						key_right.push_back(Junc_pos[i][1]);
						if (Hash_left.find(key_left) != Hash_left.end() && Hash_right.find(key_right) != Hash_right.end())
						{
							vector<double> value_left=Hash_left[key_left];
							vector<double> value_right=Hash_right[key_right];
							vector<double> value;
							if (value_left[0]==value_right[0])
							{
								vector<double> edge_temp;
								edge_temp.push_back(value_left[1]);
								edge_temp.push_back(value_right[1]);
								edge_temp.push_back(1/multi_1);
								Gene_edges[value_left[0]].push_back(edge_temp);
								value.push_back(value_left[0]);
								value.push_back(Gene_edges[value_left[0]].size()-1);
								Hash_edges[key]=value;
								if (pair_vec_L.size() == 0)
								{
									pair_vec_L.push_back(int(value[0])); // record gene_id;
								}
								if (pair_vec_L.size() == 1)
                                                		{
                                                        		pair_vec_L.push_back(int(Gene_edges[value[0]][value[1]][0]));
                                                		}
								pair_vec_L.push_back(int(Gene_edges[value[0]][value[1]][1]));
								if (value[0] != pair_vec_L[0] || pair_vec_L[pair_vec_L.size()-2]!=Gene_edges[value[0]][value[1]][0])
                                                		{
                                                        		pair_vec_L.clear();
                                                		}
							}
						}
					}
                        	} //for (int i=0;i<Junc_pos.size();i++)
			}//else if (N_count_1>0)
			if (N_count_2 == 0 && multi_2 >=1)
                	{
                        	N_R=0;
                        	for (size_t j=0;j<Chr_raw.size();j++)
                       	 	{
                                	if (Chr_raw[j]==chr_id_2)
                                	{
                                        	pair_vec_R.push_back(j);
                                        	break;
                                	}
                        	}
                        	pair_vec_R.push_back(map_start_2);
                        	pair_vec_R.push_back(Identify_cigar(map_cigar_2));
                        	pair_vec_R.push_back(multi_2);
                        }
			else if (N_count_2>0)
			{
				N_R=1;
                        	vector<vector<int> > Junc_pos = Identify_cigar(map_start_2,map_cigar_2);
                        	for (size_t i=0;i<Junc_pos.size();i++)
                        	{
                                	vector<double> key,key_left,key_right;
                                	for (size_t j=0;j<Chr_raw.size();j++)
                               		{
                                        	if (Chr_raw[j]==chr_id_2)
                                        	{
                                                	key.push_back(j);
							key_left.push_back(j);
							key_right.push_back(j);
                                                	break;
                                       		}
                                	}
                                	key.push_back(Junc_pos[i][0]);
                                	key.push_back(Junc_pos[i][1]);
                                	if (Hash_edges.find(key) != Hash_edges.end())
					{
	                                	vector<double> value = Hash_edges[key];
						if (pair_vec_R.size() == 0)
						{
                        	                        pair_vec_R.push_back(int(value[0])); // record gene_id;
						}
						Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi_2;
						if (pair_vec_R.size() == 1)
						{
							pair_vec_R.push_back(int(Gene_edges[value[0]][value[1]][0]));
						}
						pair_vec_R.push_back(int(Gene_edges[value[0]][value[1]][1]));
						if (value[0] != pair_vec_R[0])
						{
							pair_vec_R.clear();
						}
					}
					else
					{
						key_left.push_back(Junc_pos[i][0]);
						key_right.push_back(Junc_pos[i][1]);
						if (Hash_left.find(key_left) != Hash_left.end() && Hash_right.find(key_right) != Hash_right.end())
						{
							vector<double> value_left=Hash_left[key_left];
							vector<double> value_right=Hash_right[key_right];
							vector<double> value;
							if (value_left[0]==value_right[0])
							{
								vector<double> edge_temp;
								edge_temp.push_back(value_left[1]);
								edge_temp.push_back(value_right[1]);
								edge_temp.push_back(1/multi_1);
								Gene_edges[value_left[0]].push_back(edge_temp);
								value.push_back(value_left[0]);
								value.push_back(Gene_edges[value_left[0]].size()-1);
								Hash_edges[key]=value;
								if (pair_vec_R.size() == 0)
								{
									pair_vec_R.push_back(int(value[0])); // record gene_id;
								}
								if (pair_vec_R.size() == 1)
                                                		{
                                                        		pair_vec_R.push_back(int(Gene_edges[value[0]][value[1]][0]));
                                                		}
								pair_vec_R.push_back(int(Gene_edges[value[0]][value[1]][1]));
								if (value[0] != pair_vec_R[0] || pair_vec_R[pair_vec_R.size()-2]!=Gene_edges[value[0]][value[1]][0])
                                                		{
                                                        		pair_vec_R.clear();
                                                		}
							}
						}
					}
                        	}
			}
			if (N_L == 0 && N_R == 0 && Check_flag(flag_1) == 1 && pair_vec_L.size()>0 && pair_vec_R.size()>0)
                	{
                        	pair_vec.clear();
                                string chr_temp = Chr_raw[pair_vec_L[0]];
                                if (chr_temp[chr_temp.size()-1] == '+' && pair_vec_L[1]<=pair_vec_R[1])
                                {
                                       	pair_vec.push_back(pair_vec_L[0]);
                                       	pair_vec.push_back(pair_vec_L[1]);
                                       	pair_vec.push_back(pair_vec_R[1]);
                                       	pair_vec.push_back(pair_vec_L[2]);
                                       	pair_vec.push_back(pair_vec_R[2]);
                                       	pair_vec.push_back(pair_vec_L[3]);
                                       	pair_vec.push_back(pair_vec_R[3]);
                                       	Pair_vec_multi[pair_vec_L[0]].push_back(pair_vec);
                                }
                                else if (chr_temp[chr_temp.size()-1] == '-' && pair_vec_L[1]>=pair_vec_R[1])
                                {
                                       	pair_vec.push_back(pair_vec_L[0]);
                                       	pair_vec.push_back(pair_vec_R[1]);
                                        pair_vec.push_back(pair_vec_L[1]);
                                    	pair_vec.push_back(pair_vec_R[2]);
                                        pair_vec.push_back(pair_vec_L[2]);
                                        pair_vec.push_back(pair_vec_R[3]);
                                        pair_vec.push_back(pair_vec_L[3]);
                                        Pair_vec_multi[pair_vec_L[0]].push_back(pair_vec);
                                }
				pair_vec_L.clear();
				pair_vec_R.clear();
                	}
			else if (N_L == 1 && N_R == 0 && Check_flag(flag_1) == 1 && pair_vec_L.size()>0 && pair_vec_R.size()>0)
	                {
				pair_vec.clear();
				int gene_id_tmp=pair_vec_L[0];
				int node_id_tmp=-1;
				for (size_t j=0;j<Gene_nodes[gene_id_tmp].size();j++)
				{
					if (pair_vec_R[1]+2>=Gene_nodes[gene_id_tmp][j][2] && pair_vec_R[1]<=Gene_nodes[gene_id_tmp][j][3])
					{
						node_id_tmp=j;
						break;
					}
				}
				if (node_id_tmp != -1) //still -1 means the start_pos is not in this gene, but in other gene;
				{
					pair_vec_R.clear();
					pair_vec_R.push_back(node_id_tmp);
					int end_pos=pair_vec_R[1]+pair_vec_R[2]-1;
					double exon_cov=double(pair_vec_R[2])/double(pair_vec_R[3])/double(Gene_nodes[gene_id_tmp][node_id_tmp][3]-Gene_nodes[gene_id_tmp][node_id_tmp][2]+1);
					if (Gene_edges_partial[gene_id_tmp].size()==0)
					{
						Gene_nodes[gene_id_tmp][node_id_tmp][4]+=exon_cov;
					}
					else
					{
						int num_id=0;
						vector<double> key;
						vector<double> value;
						key.push_back(pair_vec_R[0]);
						key.push_back(0);
						key.push_back(0);
						for (int j=node_id_tmp;j<int(Gene_nodes[gene_id_tmp].size())-1;j++)
						{
							if (end_pos>Gene_nodes[gene_id_tmp][j][3])
							{
								num_id++;
								key[1]=Gene_nodes[gene_id_tmp][j][3];
								key[2]=Gene_nodes[gene_id_tmp][j+1][2];
								if (Hash_edges.find(key) == Hash_edges.end())
									continue;
								value = Hash_edges[key];
								Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi_2;
								pair_vec_R.push_back(Gene_edges[value[0]][value[1]][1]);
							}
							else
							{
								break;
							}
						}
						if (num_id==0)
						{
							Gene_nodes[gene_id_tmp][node_id_tmp][4]+=exon_cov;
						}
					}
					if (chr_id_1[chr_id_1.size()-1] == '+' && map_start_1 <= map_start_2 && multi_1==1 && multi_2==1)
					{
						pair_vec=pair_vec_L;
						pair_vec.push_back(-1);
						for (size_t j=0;j<pair_vec_R.size();j++)
						{
							pair_vec.push_back(pair_vec_R[j]);
						}
						Pair_vec_junc[pair_vec[0]].push_back(pair_vec);
					}
					if (chr_id_1[chr_id_1.size()-1] == '-' && map_start_1 >= map_start_2 && multi_1==1 && multi_2==1)
					{
						pair_vec.push_back(pair_vec_L[0]); //gene_id
						for (size_t j=0;j<pair_vec_R.size();j++)
                                                {
                                                        pair_vec.push_back(pair_vec_R[j]);
                                                }
						pair_vec.push_back(-1);
						for (size_t j=1;j<pair_vec_L.size();j++)
						{
							pair_vec.push_back(pair_vec_L[j]);
						}
						Pair_vec_junc[pair_vec[0]].push_back(pair_vec);
					}
					pair_vec_L.clear();
					pair_vec_R.clear();
				}
        	        }
                	else if (N_L == 0 && N_R == 1 && Check_flag(flag_1) == 1 && pair_vec_L.size()>0 && pair_vec_R.size()>0)
                	{
				pair_vec.clear();
				int gene_id_tmp=pair_vec_R[0];
				int node_id_tmp=-1;
				for (size_t j=0;j<Gene_nodes[gene_id_tmp].size();j++)
				{
					if (pair_vec_L[1]+2>=Gene_nodes[gene_id_tmp][j][2] && pair_vec_L[1]<=Gene_nodes[gene_id_tmp][j][3])
					{
						node_id_tmp=j;
						break;
					}
				}
				if (node_id_tmp != -1) //still -1 means the start_pos is not in this gene, but in other gene;
				{
					pair_vec_L.clear();
					pair_vec_L.push_back(node_id_tmp);
					int end_pos=pair_vec_L[1]+pair_vec_L[2]-1;
					double exon_cov=double(pair_vec_L[2])/double(pair_vec_L[3])/double(Gene_nodes[gene_id_tmp][node_id_tmp][3]-Gene_nodes[gene_id_tmp][node_id_tmp][2]+1);
					if (Gene_edges_partial[gene_id_tmp].size()==0)
					{
						Gene_nodes[gene_id_tmp][node_id_tmp][4]+=exon_cov;
					}
					else
					{
						int num_id=0;
						vector<double> key;
						vector<double> value;
						key.push_back(pair_vec_L[0]);
						key.push_back(0);
						key.push_back(0);
						for (int j=node_id_tmp;j<int(Gene_nodes[gene_id_tmp].size())-1;j++)
						{
							if (end_pos>Gene_nodes[gene_id_tmp][j][3])
							{
								num_id++;
								key[1]=Gene_nodes[gene_id_tmp][j][3];
								key[2]=Gene_nodes[gene_id_tmp][j+1][2];
								if (Hash_edges.find(key) == Hash_edges.end())
									continue;
								value = Hash_edges[key];
								Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi_2;
								pair_vec_L.push_back(Gene_edges[value[0]][value[1]][1]);
							}
							else
							{
								break;
							}
						}
						if (num_id == 0)
						{
							Gene_nodes[gene_id_tmp][node_id_tmp][4]+=exon_cov;
						}
					}
					if (chr_id_1[chr_id_1.size()-1] == '+' && map_start_1 <= map_start_2 && multi_1==1 && multi_2==1)
                                        {
                                                pair_vec.push_back(pair_vec_R[0]); //gene_id
						for (size_t j=0;j<pair_vec_L.size();j++)
                                                {
                                                        pair_vec.push_back(pair_vec_L[j]);
                                                }
                                                pair_vec.push_back(-1);
                                                for (size_t j=1;j<pair_vec_R.size();j++)
                                                {
                                                        pair_vec.push_back(pair_vec_R[j]);
                                                }
                                                Pair_vec_junc[pair_vec[0]].push_back(pair_vec);
                                        }
                                        if (chr_id_1[chr_id_1.size()-1] == '-' && map_start_1 >= map_start_2 && multi_1==1 && multi_2==1)
                                        {
                                                pair_vec=pair_vec_R;
                                                pair_vec.push_back(-1);
						for (size_t j=0;j<pair_vec_L.size();j++)
                                                {
                                                        pair_vec.push_back(pair_vec_L[j]);
                                                }
                                                Pair_vec_junc[pair_vec[0]].push_back(pair_vec);
                                        }
					pair_vec_L.clear();
                                        pair_vec_R.clear();
				}
                	}
                	else if (N_L == 1 && N_R == 1 && Check_flag(flag_1) == 1 && pair_vec_L.size()>0 && pair_vec_R.size()>0 && pair_vec_L[0] == pair_vec_R[0])
                	{
				// pair_vec_L.size()=0 means that the junction read does not span a junction in Hash_edges;
				pair_vec.clear();
				if (chr_id_1[chr_id_1.size()-1] == '+' && map_start_1 <= map_start_2 && multi_1==1 && multi_2==1)
				{
					pair_vec=pair_vec_L;
					pair_vec.push_back(-1);
					for (size_t j=1;j<pair_vec_R.size();j++)
					{
						pair_vec.push_back(pair_vec_R[j]);
					}
					Pair_vec_junc[pair_vec[0]].push_back(pair_vec);
				}
				if (chr_id_1[chr_id_1.size()-1] == '-' && map_start_1 >= map_start_2 && multi_1==1 && multi_2==1)
                                {
                                        pair_vec=pair_vec_R;
                                        pair_vec.push_back(-1);
                                        for (size_t j=1;j<pair_vec_L.size();j++)
                                        {
                                                pair_vec.push_back(pair_vec_L[j]);
                                        }
                                        Pair_vec_junc[pair_vec[0]].push_back(pair_vec);
                                }
				pair_vec_L.clear();
                                pair_vec_R.clear();
                	}
			if (N_L == 0 && pair_vec_L.size()>0)
			{
				Map_need_to_do[pair_vec_L[0]].push_back(pair_vec_L);
			}
			if (N_R == 0 && pair_vec_R.size()>0)
			{
				Map_need_to_do[pair_vec_R[0]].push_back(pair_vec_R);
			}
		} //for (int i=0;i<pair_num;i++) : for each pair
		int pair_num_last=int(sam_line.size())%2;
		if (pair_num_last == 1)
		{
			istr.str(sam_line[sam_line.size()-1]);
			istr>>read_id_1>>flag_1>>chr_id_1>>map_start_1>>map_qual_1>>map_cigar_1>>equal_star_1>>pair_start_1>>map_length_1;
			multi_1=0;
			while (istr>>temp)
                	{
                        	if (map_cigar_1 == "*" || multi_1>0)
                        	{
                                	break;
                        	}
                        	if (temp.substr(0,5)=="NH:i:")
                        	{
                                	multi_1=atoi(temp.substr(5,temp.size()-5).c_str());
                        	}
                	}
			istr.clear();
			int N_count_1=count(map_cigar_1.begin(),map_cigar_1.end(),'N');
			pair_vec_L.clear();
			pair_vec_R.clear();
			if (N_count_1 == 0 && multi_1 >=1)
                	{
                        	N_L=0;
                        	for (size_t j=0;j<Chr_raw.size();j++)
                       	 	{
                                	if (Chr_raw[j]==chr_id_1)
                                	{
                                        	pair_vec_L.push_back(j);
                                        	break;
                                	}
                        	}
                        	pair_vec_L.push_back(map_start_1);
                        	pair_vec_L.push_back(Identify_cigar(map_cigar_1));
                        	pair_vec_L.push_back(multi_1);
                        }
			if (N_count_1>0)
			{
                        	vector<vector<int> > Junc_pos = Identify_cigar(map_start_1,map_cigar_1);
                        	for (size_t i=0;i<Junc_pos.size();i++)
                        	{
                                	vector<double> key;
                                	vector<double> key_left,key_right;
                                	for (size_t j=0;j<Chr_raw.size();j++)
                               		{
                                        	if (Chr_raw[j]==chr_id_1)
                                        	{
                                                	key.push_back(j);
                                                	break;
                                       		}
                                	}
                                	key.push_back(Junc_pos[i][0]);
                                	key.push_back(Junc_pos[i][1]);
                                	if (Hash_edges.find(key) != Hash_edges.end())
					{
                                		vector<double> value = Hash_edges[key];
						Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi_1;
					}
					else
					{
						key_left.push_back(Junc_pos[i][0]);
						key_right.push_back(Junc_pos[i][1]);
						if (Hash_left.find(key_left) != Hash_left.end() && Hash_right.find(key_right) != Hash_right.end())
						{
							vector<double> value_left=Hash_left[key_left];
							vector<double> value_right=Hash_right[key_right];
							vector<double> value;
							if (value_left[0]==value_right[0])
							{
								vector<double> edge_temp;
								edge_temp.push_back(value_left[1]);
								edge_temp.push_back(value_right[1]);
								edge_temp.push_back(1/multi_1);
								Gene_edges[value_left[0]].push_back(edge_temp);
								value.push_back(value_left[0]);
								value.push_back(Gene_edges[value_left[0]].size()-1);
								Hash_edges[key]=value;
							}
						}
					}
                        	} //for (int i=0;i<Junc_pos.size();i++)
			}
			if (N_L == 0 && pair_vec_L.size()>0)
			{
				Map_need_to_do[pair_vec_L[0]].push_back(pair_vec_L);
			}
		} //if (pair_num_last == 1)
		sam_line.clear();
		sam_line.push_back(curr_line);

	} //while (getline(ifs,temp));
	ifs.close();
}

vector<vector<int> > Gene_info::search_paths(int node_start, int node_end, int gene_id)
{
	vector<vector<int> > Enu_paths_temp,Enu_paths;
        vector<int> enu_temp,add_temp;
	enu_temp.push_back(node_start);
	Enu_paths.push_back(enu_temp);
        enu_temp.clear();
		while (1)
        {
                for (size_t i=0;i<Enu_paths.size();i++)
                {
                        add_temp.clear();
                        enu_temp=Enu_paths[i];
                        for (size_t j=0;j<Gene_edges[gene_id].size();j++)
                        {
                                if (Gene_edges[gene_id][j][0]==enu_temp[enu_temp.size()-1] && Gene_edges[gene_id][j][1]<=node_end)
                                {
                                        add_temp.push_back(Gene_edges[gene_id][j][1]);
                                }
                        }
                        if (add_temp.size()>0)
                        {
                                for (size_t j=0;j<add_temp.size();j++)
                                {
                                        enu_temp.push_back(add_temp[j]);
                                        Enu_paths_temp.push_back(enu_temp);
                                        enu_temp=Enu_paths[i];
                                }
                        }
                        else
                        {
                                Enu_paths_temp.push_back(Enu_paths[i]);
                        }
                }
				if (Enu_paths == Enu_paths_temp)
				{
					break;
				}
				else
				{
                	Enu_paths.clear();
					for (size_t i=0;i<Enu_paths_temp.size();i++)
					{
						int sum_i=0;
						for (size_t j=1;j<Enu_paths_temp[i].size();j++)
						{
							sum_i=sum_i+Gene_nodes[gene_id][Enu_paths_temp[i][j]][3]-Gene_nodes[gene_id][Enu_paths_temp[i][j]][2]+1;
						}
						if (Enu_paths_temp[i][Enu_paths_temp[i].size()-1]==node_end || sum_i<=300)
						{
							Enu_paths.push_back(Enu_paths_temp[i]);
						}
					}
					Enu_paths_temp.clear();
				}
        }//while (1)
		Enu_paths_temp.clear();
		for (size_t i=0;i<Enu_paths.size();i++)
		{
				if (Enu_paths[i][Enu_paths[i].size()-1] == node_end)
				{
						Enu_paths_temp.push_back(Enu_paths[i]);
				}
		}
	return Enu_paths_temp;
}

void Gene_info::process_gene_paths() //do this after process_pair_vec_junc();
{
	for (size_t i=0;i<Gene_paths.size();i++)
	{
		if (Gene_paths[i].size()>0) //there exist paths;
		{
			//cerr<<"Processing Gene Path: "<<i<<" ..."<<endl;
			vector<vector<vector<int> > > paths;
			vector<vector<int> > sub_path,Edges, Paths_reserved;
			vector<int> sub_path_temp, edge, uniq_paths, Edges_left, Edges_right, Edges_temp;
			int uniq_id,uniq_path;
			vector<vector<int> >::iterator pos;
			vector<int>::iterator iter;
			vector<int> Single_exon_paths_reserved;
			for (size_t j=0;j<Gene_paths[i].size();j++)
			{
				if (Gene_paths[i][j].size() == 1)
				{
					if (Gene_nodes[i][Gene_paths[i][j][0]][4]>0)
					{
						Single_exon_paths_reserved.push_back(Gene_paths[i][j][0]);
					}
					continue;
				}
				for (size_t k=0;k<Gene_paths[i][j].size()-1;k++)
				{
					sub_path_temp.push_back(Gene_paths[i][j][k]);
					sub_path_temp.push_back(Gene_paths[i][j][k+1]);
					sub_path.push_back(sub_path_temp);
					sub_path_temp.clear();
				}
				paths.push_back(sub_path);
				sub_path.clear();
			}
			edge.push_back(0);
			edge.push_back(0);
			//vector<double> path_edge_R;
			for (size_t j=0;j<Gene_edges[i].size();j++)
			{
				if (Gene_edges[i][j][2]==0)
				{
					continue;
				}
				edge[0]=Gene_edges[i][j][0];
				edge[1]=Gene_edges[i][j][1];
				uniq_id=0;
				for (size_t k=0;k<paths.size();k++)
				{
					pos=find(paths[k].begin(),paths[k].end(),edge);
					if (pos != paths[k].end())
					{
						uniq_path=k;
						uniq_id++;
					}
					if (uniq_id>=2)
						break;
				}
				if (uniq_id == 1 && Gene_edges[i][j][2]<=Cov_uniq)
				{
					uniq_paths.push_back(uniq_path);
					//path_edge_R.push_back(Gene_edges[i][j][2]);
				}
			}
			sort(uniq_paths.begin(),uniq_paths.end());
			iter=unique(uniq_paths.begin(),uniq_paths.end());
			uniq_paths.erase(iter,uniq_paths.end());
			for (size_t j=0;j<uniq_paths.size();j++)
			{
				for (size_t k=0;k<paths[uniq_paths[j]].size();k++)
				{
					for (size_t t=0;t<Gene_edges[i].size();t++)
					{
						if (paths[uniq_paths[j]][k][0]==Gene_edges[i][t][0] && paths[uniq_paths[j]][k][1]==Gene_edges[i][t][1])
						{
							if (Gene_edges[i][t][2]==0)
							{
								//cerr<<"Gene_"<<i<<"_Path_"<<uniq_paths[j]<<": "<<path_edge_R[j]<<endl;
								Gene_edges[i][t][2]=Gene_edges[i][t][2]+1;
							}
							break;
						}
					}
				}
			}
			for (size_t j=0;j<Gene_edges[i].size();j++)
			{
				if (Gene_edges[i][j][2]>0)
				{
					Edges_temp.push_back(Gene_edges[i][j][0]);
					Edges_temp.push_back(Gene_edges[i][j][1]);
					Edges.push_back(Edges_temp);
					Edges_temp.clear();
				}
			}
			for (size_t j=0;j<paths.size();j++)
			{
				int path_keep=1;
				for (size_t k=0;k<paths[j].size();k++)
				{
					edge=paths[j][k];
					pos=find(Edges.begin(),Edges.end(),edge);
					if (pos == Edges.end())
					{
						path_keep=0;
						break;
					}
				}
				if (path_keep==1)
				{
					Paths_reserved.push_back(Gene_paths[i][j]);
				}
			}
			vector<int> single_temp;
			single_temp.push_back(0);
			for (size_t j=0;j<Single_exon_paths_reserved.size();j++)
			{
				single_temp[0]=Single_exon_paths_reserved[j];
				Paths_reserved.push_back(single_temp);
			}
			Gene_paths[i]=Paths_reserved;
		} //if (Gene_paths[i].size()>0)
	} //for (int i=0;i<Gene_paths.size();i++)
} //end of process_gene_paths()

void Gene_info::process_pair_vec_junc()
{
	vector<vector<int> >::iterator iter;
	for (size_t i=0;i<Pair_vec_junc.size();i++) //for each gene
	{
	    if (Pair_vec_junc[i].size()>0) // if it is not single-exon
	    {
	    	cerr<<"Processing gene: "<<i<<" ..."<<endl;
		sort(Pair_vec_junc[i].begin(),Pair_vec_junc[i].end());
		vector<vector<int> > Pair_temp = Pair_vec_junc[i];
		vector<int> Pair_vec_junc_num;
		int pair_num=1;
		iter=unique(Pair_vec_junc[i].begin(),Pair_vec_junc[i].end());
		Pair_vec_junc[i].erase(iter,Pair_vec_junc[i].end()); // unique of pairs for each gene
		for (size_t j=1;j<Pair_temp.size();j++)
		{
			if (Pair_temp[j] != Pair_temp[j-1])
			{
				Pair_vec_junc_num.push_back(pair_num); //number of each unique pair
				pair_num=1;
			}
			else
			{
				pair_num++;
			}
		}
		Pair_vec_junc_num.push_back(pair_num); // last one num
		for (size_t j=0;j<Pair_vec_junc[i].size();j++)
		{
			vector<int> Pair_L,Pair_R,Pair_path;
			int k_id=0;
			for (size_t k=1;k<Pair_vec_junc[i][j].size();k++)
			{
				if (Pair_vec_junc[i][j][k] == -1)
				{
					k_id=k+1;
					break;
				}
				else
					Pair_L.push_back(Pair_vec_junc[i][j][k]);
			}
			for (size_t k=k_id;k<Pair_vec_junc[i][j].size();k++)
			{
				Pair_R.push_back(Pair_vec_junc[i][j][k]);
			}
			k_id=-1;
			for (size_t k=0;k<Pair_L.size();k++)
			{
				if (Pair_R[0] == Pair_L[k])
				{
					k_id=k;
					break;
				}
			}
			int stop_signal=0;
			if (k_id != -1 && Pair_R.size()>=Pair_L.size()-k_id) //k_id = -1 means L and R have no overlap;
			{
				Pair_path=Pair_L;
				for (size_t k=k_id+1;k<Pair_L.size();k++)
				{
					if (Pair_L[k] != Pair_R[k-k_id])
					{
						stop_signal=1;
						break;
					}
				}
				if (stop_signal==0)
				{
					for (size_t k=Pair_L.size()-k_id;k<Pair_R.size();k++)
					{
						Pair_path.push_back(Pair_R[k]);
					}
					Gene_pair_paths[i].push_back(Pair_path);
				}
			}
			else
			{
				Pair_path=Pair_L;
				vector<vector<int> > Inter_paths=search_paths(Pair_L[Pair_L.size()-1],Pair_R[0],i);
				if (Inter_paths.size() == 1)
				{
					vector<int> Inter_path=Inter_paths[0];
					int t_id=0;
					for (size_t k=1;k<Inter_path.size();k++)
					{
						Pair_path.push_back(Inter_path[k]);
						for (size_t t=t_id;t<Gene_edges[i].size();t++)
						{
							if (Gene_edges[i][t][0]==Inter_path[k-1] && Gene_edges[i][t][1]==Inter_path[k])
							{
								Gene_edges[i][t][2]=Gene_edges[i][t][2]+Pair_vec_junc_num[j];
								t_id=t+1;
								break;
							}
						}
					}
					for (size_t k=1;k<Pair_R.size();k++)
					{
						Pair_path.push_back(Pair_R[k]);
					}
					Gene_pair_paths[i].push_back(Pair_path);
				}
			}
		} //for (int j=0;j<Pair_vec_junc[i].size();j++)
		sort(Gene_pair_paths[i].begin(),Gene_pair_paths[i].end());
		iter=unique(Gene_pair_paths[i].begin(),Gene_pair_paths[i].end());
                Gene_pair_paths[i].erase(iter,Gene_pair_paths[i].end());
		vector<int> vec_rm;
		for (size_t j=0;j<Gene_pair_paths[i].size();j++)
		{
			for (size_t k=0;k<Gene_pair_paths[i].size();k++)
			{
				int j_size=Gene_pair_paths[i][j].size();
				int k_size=Gene_pair_paths[i][k].size();
				vector<int> vec_j,vec_k;
				if (j<k && j_size>k_size)
				{
					vec_k=Gene_pair_paths[i][k];
					for (int t=0;t<j_size-k_size+1;t++)
					{
						for (int m=t;m<t+k_size;m++)
						{
							vec_j.push_back(Gene_pair_paths[i][j][m]);
						}
						if (vec_k == vec_j)
						{
							vec_rm.push_back(k);
							break;
						}
						vec_j.clear();
					}
				}
				if (j<k && j_size<k_size)
				{
					vec_j=Gene_pair_paths[i][j];
					for (int t=0;t<k_size-j_size+1;t++)
					{
						for (int m=t;m<t+j_size;m++)
						{
							vec_k.push_back(Gene_pair_paths[i][k][m]);
						}
						if (vec_k == vec_j)
						{
							vec_rm.push_back(j);
							break;
						}
						vec_k.clear();
					}
				}
			}
		}
		vector<vector<int> > Paths;
		vector<int>::iterator pos_rm;
		for (size_t j=0;j<Gene_pair_paths[i].size();j++)
		{
			pos_rm=find(vec_rm.begin(),vec_rm.end(),j);
			if (pos_rm == vec_rm.end() && Gene_pair_paths[i][j].size()>=3)
			{
				Paths.push_back(Gene_pair_paths[i][j]);
			}
		}
		Gene_pair_paths[i]=Paths; // there is no abundance any more from here;
		// Begin to connect pair_paths...
	/*
		vector<int> Pair_edge_left,Pair_edge_right,Pair_edge_pos;
		for (int j=0;j<Gene_pair_paths[i].size();j++)
		{
			for (int k=j+1;k<Gene_pair_paths[i].size();k++)
			{
				int conn_start=-1;
				vector<int> vec_j,vec_k;
				if (Gene_pair_paths[i][j][0]<Gene_pair_paths[i][k][0])
				{
					for (int t=0;t<Gene_pair_paths[i][j].size();t++)
					{
						if (Gene_pair_paths[i][j][t]==Gene_pair_paths[i][k][0])
						{
							conn_start=t;
							break;
						}
					}
					if (conn_start == -1 ||conn_start == Gene_pair_paths[i][j].size()-1 || Gene_pair_paths[i][j].size()-conn_start >= Gene_pair_paths[i][k].size())
						continue;
					for (int t=conn_start;t<Gene_pair_paths[i][j].size();t++)
					{
						vec_j.push_back(Gene_pair_paths[i][j][t]);
						vec_k.push_back(Gene_pair_paths[i][k][t-conn_start]);
					}
					if (vec_j == vec_k)
					{
						Pair_edge_left.push_back(j);
						Pair_edge_right.push_back(k);
						Pair_edge_pos.push_back(Gene_pair_paths[i][j].size()-conn_start);
					}
				}
				if (Gene_pair_paths[i][j][0]>Gene_pair_paths[i][k][0])
				{
					for (int t=0;t<Gene_pair_paths[i][k].size();t++)
                                        {
                                                if (Gene_pair_paths[i][k][t]==Gene_pair_paths[i][j][0])
                                                {
                                                        conn_start=t;
                                                        break;
                                                }
                                        }
                                        if (conn_start == -1 || conn_start == Gene_pair_paths[i][k].size()-1 || Gene_pair_paths[i][k].size()-conn_start >= Gene_pair_paths[i][j].size())
                                                continue;
                                        for (int t=conn_start;t<Gene_pair_paths[i][k].size();t++)
                                        {
						vec_k.push_back(Gene_pair_paths[i][k][t]);
						vec_j.push_back(Gene_pair_paths[i][j][t-conn_start]);
                                        }
                                        if (vec_j == vec_k)
                                        {
                                                Pair_edge_left.push_back(k);
                                                Pair_edge_right.push_back(j);
                                                Pair_edge_pos.push_back(Gene_pair_paths[i][k].size()-conn_start);
                                        }
				}
			}
		}
		// Begin to enumerate all pair connecting paths
		vector<vector<int> > Enu_paths=enumerate_pair_paths(Pair_edge_left,Pair_edge_right,Pair_edge_pos,Gene_pair_paths[i].size());
		vector<int> path_temp;
		int pos_id;
		Paths.clear();
		for (int j=0;j<Enu_paths.size();j++)
		{
			path_temp=Gene_pair_paths[i][Enu_paths[j][0]];
			for (int k=1;k<Enu_paths[j].size();k++)
			{
				for (int t=0;t<Pair_edge_pos.size();t++)
				{
					if (Pair_edge_left[t]==Enu_paths[j][k-1] && Pair_edge_right[t]==Enu_paths[j][k])
					{
						pos_id=t;
						break;
					}
				}
				for (int t=Pair_edge_pos[pos_id];t<Gene_pair_paths[i][Enu_paths[j][k]].size();t++)
				{
					path_temp.push_back(Gene_pair_paths[i][Enu_paths[j][k]][t]);
				}
			}
			Paths.push_back(path_temp);
		}
		Gene_pair_paths[i]=Paths;
	*/
	    } //if (Pair_vec_junc[i].size()>0)
	} //for (int i=0;i<Pair_vec_junc.size();i++)
}
void Gene_info::test_pair_path(char * path)
{
		vector<vector<int> > Gene_pair_paths_i;
		ifstream ifs(path);
		string temp;
		vector<int> curr_path;
		while (getline(ifs,temp))
		{
			int path_start_pos=0;
                        for (size_t j=0;j<temp.size();j++)
			{
				if (temp[j]==' ')
				{
					curr_path.push_back(atoi(temp.substr(path_start_pos,j-path_start_pos).c_str()));
					path_start_pos=j+1;
				}
			}
			Gene_pair_paths_i.push_back(curr_path);
			curr_path.clear();
		}
		vector<int> Pair_edge_left,Pair_edge_right,Pair_edge_pos;
		for (size_t j=0;j<Gene_pair_paths_i.size();j++)
		{
			for (size_t k=j+1;k<Gene_pair_paths_i.size();k++)
			{
				int conn_start=-1;
				vector<int> vec_j,vec_k;
				if (Gene_pair_paths_i[j][0]<Gene_pair_paths_i[k][0])
				{
					for (size_t t=0;t<Gene_pair_paths_i[j].size();t++)
					{
						if (Gene_pair_paths_i[j][t]==Gene_pair_paths_i[k][0])
						{
							conn_start=t;
							break;
						}
					}
					if (conn_start == -1 ||conn_start == int(Gene_pair_paths_i[j].size())-1 || int(Gene_pair_paths_i[j].size())-conn_start >= int(Gene_pair_paths_i[k].size()))
						continue;
					for (size_t t=conn_start;t<Gene_pair_paths_i[j].size();t++)
					{
						vec_j.push_back(Gene_pair_paths_i[j][t]);
						vec_k.push_back(Gene_pair_paths_i[k][t-conn_start]);
					}
					if (vec_j == vec_k)
					{
						Pair_edge_left.push_back(j);
						Pair_edge_right.push_back(k);
						Pair_edge_pos.push_back(Gene_pair_paths_i[j].size()-conn_start);
					}
				}
				if (Gene_pair_paths_i[j][0]>Gene_pair_paths_i[k][0])
				{
					for (size_t t=0;t<Gene_pair_paths_i[k].size();t++)
                                        {
                                                if (Gene_pair_paths_i[k][t]==Gene_pair_paths_i[j][0])
                                                {
                                                        conn_start=t;
                                                        break;
                                                }
                                        }
                                        if (conn_start == -1 || conn_start == int(Gene_pair_paths_i[k].size())-1 || int(Gene_pair_paths_i[k].size())-conn_start >= int(Gene_pair_paths_i[j].size()))
                                                continue;
                                        for (size_t t=conn_start;t<Gene_pair_paths_i[k].size();t++)
                                        {
						vec_k.push_back(Gene_pair_paths_i[k][t]);
						vec_j.push_back(Gene_pair_paths_i[j][t-conn_start]);
                                        }
                                        if (vec_j == vec_k)
                                        {
                                                Pair_edge_left.push_back(k);
                                                Pair_edge_right.push_back(j);
                                                Pair_edge_pos.push_back(Gene_pair_paths_i[k].size()-conn_start);
                                        }
				}
			}
		}
		// Begin to enumerate all pair connecting paths
		vector<vector<int> > Enu_paths=enumerate_pair_paths(Pair_edge_left,Pair_edge_right,Pair_edge_pos,Gene_pair_paths_i.size());
		vector<int> path_temp;
		int pos_id=0;
		vector<vector<int> > Paths;
		for (size_t j=0;j<Enu_paths.size();j++)
		{
			path_temp=Gene_pair_paths_i[Enu_paths[j][0]];
			for (size_t k=1;k<Enu_paths[j].size();k++)
			{
				for (size_t t=0;t<Pair_edge_pos.size();t++)
				{
					if (Pair_edge_left[t]==Enu_paths[j][k-1] && Pair_edge_right[t]==Enu_paths[j][k])
					{
						pos_id=t;
						break;
					}
				}
				for (size_t t=Pair_edge_pos[pos_id];t<Gene_pair_paths_i[Enu_paths[j][k]].size();t++)
				{
					path_temp.push_back(Gene_pair_paths_i[Enu_paths[j][k]][t]);
				}
			}
			Paths.push_back(path_temp);
		}
		cout<<"size = "<<Paths.size()<<endl;
}

vector<vector<int> > Gene_info::enumerate_pair_paths(vector<int> Pair_edge_left,vector<int> Pair_edge_right,vector<int> Pair_edge_pos,int max)
{
	vector<int> temp_in,temp_out,S,T,Single_pair;
	int max_num=max;
	int i,j,t;
	for (i=0;i<max_num;i++)
        {
                temp_in.clear();
                temp_out.clear();
		for (j=0;j<int(Pair_edge_left.size());j++)
                {
                        if (Pair_edge_right[j]==i)
                        {
                                temp_in.push_back(j);
                        }
                        else if (Pair_edge_left[j]==i)
                        {
                                temp_out.push_back(j);
                        }
                }
                if (temp_in.size()==0 && temp_out.size()>0)
                {
                        S.push_back(i);
                }
                else if (temp_in.size()>0 && temp_out.size()==0)
                {
                        T.push_back(i);
                }
                else if (temp_in.size()==0 && temp_out.size()==0)
                {
                        Single_pair.push_back(i);
                }
        }
	vector<vector<int> > Enu_paths,Enu_paths_temp;
        vector<int> enu_temp,add_temp;
        vector<int> update_index;
	vector<int>::iterator pos;
        for (i=0;i<int(S.size());i++)
        {
                enu_temp.push_back(S[i]);
                Enu_paths.push_back(enu_temp);
                enu_temp.clear();
        }
        for (i=0;i<int(Enu_paths.size());i++)
        {
                pos=find(T.begin(),T.end(),Enu_paths[i][Enu_paths[i].size()-1]);
                if (pos==T.end())
                {
                        update_index.push_back(i);
                }
        }
	t=0;
	while (update_index.size()>0)
        {
		cout<<"t = "<<t<<endl;
		t++;
                for (i=0;i<int(Enu_paths.size());i++)
                {
                        add_temp.clear();
                        enu_temp=Enu_paths[i];
                        for (j=0;j<int(Pair_edge_left.size());j++)
                        {
                                if (Pair_edge_left[j]==enu_temp[enu_temp.size()-1])
                                {
                                        add_temp.push_back(Pair_edge_right[j]);
                                }
                        }
                        if (add_temp.size()>0)
                        {
                                for (j=0;j<int(add_temp.size());j++)
                                {
                                        enu_temp.push_back(add_temp[j]);
                                        Enu_paths_temp.push_back(enu_temp);
                                        enu_temp=Enu_paths[i];
                                }
                        }
                        else
                        {
                                Enu_paths_temp.push_back(Enu_paths[i]);
                        }
                }
                Enu_paths=Enu_paths_temp;
		cout<<t<<": "<<Enu_paths.size()<<endl;
                Enu_paths_temp.clear();
                update_index.clear();
                for (i=0;i<int(Enu_paths.size());i++)
                {
                        pos=find(T.begin(),T.end(),Enu_paths[i][Enu_paths[i].size()-1]);
                        if (pos==T.end())
                        {
                                update_index.push_back(i);
                        }
                }
	}//while (update_index.size()>0)
	vector<int> single_temp;
	for (i=0;i<int(Single_pair.size());i++)
	{
		single_temp.push_back(Single_pair[i]);
		Enu_paths.push_back(single_temp);
		single_temp.clear();
	}
	return Enu_paths;
}
/*
void Gene_info::check_pair_vec_junc()
{
	for (int i=0;i<Pair_vec_junc.size();i++)
	{
		for (int j=0;j<Pair_vec_junc[i].size();j++)
		{
			
		}
	}
}
*/
void Gene_info::process_map_need_to_do()
{
	for (size_t i=0;i<Map_need_to_do.size();i++)
        {
		sort(Map_need_to_do[i].begin(),Map_need_to_do[i].end());
	}
	int gene_index=0;
	vector<int> pair_vec;
	for (size_t i=0;i<Map_need_to_do.size();i++)
        {
		for (size_t j=0;j<Map_need_to_do[i].size();j++)
                {
			int node_index=-1;
			for (size_t k=gene_index;k<Gene_nodes.size();k++)
			{
				if (Chr_raw[Map_need_to_do[i][j][0]] == Gene_chrs[k] && Map_need_to_do[i][j][1]+2>=Gene_nodes[k][0][2] && Map_need_to_do[i][j][1]<=Gene_nodes[k][Gene_nodes[k].size()-1][3])
				{
					gene_index=k;
					for (int t=0;t<int(Gene_nodes[k].size());t++)
					{
						if (Map_need_to_do[i][j][1]+2>=Gene_nodes[k][t][2] && Map_need_to_do[i][j][1]<=Gene_nodes[k][t][3])
						{
							node_index=t;
							break;
						}
					}
					break;
				}
			}
			int end_pos=Map_need_to_do[i][j][1]+Map_need_to_do[i][j][2]-1;
			if (Gene_edges_partial[gene_index].size()==0 && node_index != -1)
			{
				double exon_cov=double(Map_need_to_do[i][j][2])/double(Map_need_to_do[i][j][3])/double(Gene_nodes[gene_index][node_index][3]-Gene_nodes[gene_index][node_index][2]+1);
				Gene_nodes[gene_index][node_index][4]+=exon_cov;
			}
			else if (Gene_edges_partial[gene_index].size()>0 && node_index != -1)
			{
				int num_id=0;
				vector<double> key;
                                vector<double> value;
				key.push_back(Map_need_to_do[i][j][0]);
				key.push_back(0);
				key.push_back(0);
				for (int t=node_index;t<int(Gene_nodes[gene_index].size())-1;t++)
				{
					if (end_pos>Gene_nodes[gene_index][t][3])
					{
						num_id++;
						key[1]=Gene_nodes[gene_index][t][3];
						key[2]=Gene_nodes[gene_index][t+1][2];
						if (Hash_edges.find(key) == Hash_edges.end())
							continue;
						value = Hash_edges[key];
						double multi=double(Map_need_to_do[i][j][3]);
						Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi;
					}
					else
						break;
				}
				if (num_id == 0)
				{
					double exon_cov=double(Map_need_to_do[i][j][2])/double(Map_need_to_do[i][j][3])/double(Gene_nodes[gene_index][node_index][3]-Gene_nodes[gene_index][node_index][2]+1);
					Gene_nodes[gene_index][node_index][4]+=exon_cov;
				}
			}
		}
	}
}
void Gene_info::process_pair_vec_multi()
{
	// sort Pair_vec_multi based on each chr.
	int gene_index_1=0;
	int gene_index_2=0;
	vector<int> pair_vec;
	for (size_t i=0;i<Pair_vec_multi.size();i++)
        {
		sort(Pair_vec_multi[i].begin(),Pair_vec_multi[i].end());
		cerr<<"Pair_vec_multi "<<i<<": "<<Pair_vec_multi[i].size()<<endl;
		for (size_t j=0;j<Pair_vec_multi[i].size();j++)
                {
			double multi_1=0;
			double multi_2=0;
			int node_index_1=-1;
			int node_index_2=-1;
			pair_vec.clear();
			for (size_t k=gene_index_1;k<Gene_nodes.size();k++)
			{
				string curr_chr=Chr_raw[Pair_vec_multi[i][j][0]];
				int start_pos=Pair_vec_multi[i][j][1];
				int node_last_pos=Gene_nodes[k][Gene_nodes[k].size()-1][3];
				if (curr_chr == Gene_chrs[k] && start_pos+2>=Gene_nodes[k][0][2] && start_pos<=node_last_pos)
				{
					gene_index_1=k;
					gene_index_2=k;
					for (int t=0;t<int(Gene_nodes[k].size());t++)
					{
						if (start_pos+2>=Gene_nodes[k][t][2] && start_pos<=Gene_nodes[k][t][3])
						{
							node_index_1=t;
							break;
						}
					}
					break;
				}
			}
			for (size_t k=gene_index_2;k<Gene_nodes.size();k++)
                        {
				string curr_chr=Chr_raw[Pair_vec_multi[i][j][0]];
				int start_pos=Pair_vec_multi[i][j][2];
				int node_last_pos=Gene_nodes[k][Gene_nodes[k].size()-1][3];
				if (curr_chr == Gene_chrs[k] && start_pos+2>=Gene_nodes[k][0][2] && start_pos<=node_last_pos)
				{
					gene_index_2=k;
					for (int t=0;t<int(Gene_nodes[k].size());t++)
                                        {
                                                if (start_pos+2>=Gene_nodes[k][t][2] && start_pos<=Gene_nodes[k][t][3])
                                                {
                                                        node_index_2=t;
                                                        break;
                                                }
                                        }
					break;
				}
			}
			pair_vec.push_back(gene_index_1);
			pair_vec.push_back(node_index_1);
			int end_pos_1=Pair_vec_multi[i][j][1]+Pair_vec_multi[i][j][3]-1;
			int end_pos_2=Pair_vec_multi[i][j][2]+Pair_vec_multi[i][j][4]-1;
			double exon_cov_1=double(Pair_vec_multi[i][j][3])/double(Pair_vec_multi[i][j][5])/double(Gene_nodes[gene_index_1][node_index_1][3]-Gene_nodes[gene_index_1][node_index_1][2]+1);
			double exon_cov_2=double(Pair_vec_multi[i][j][4])/double(Pair_vec_multi[i][j][6])/double(Gene_nodes[gene_index_2][node_index_2][3]-Gene_nodes[gene_index_2][node_index_2][2]+1);
			if (Gene_edges_partial[gene_index_1].size()==0)
			{
				Gene_nodes[gene_index_1][node_index_1][4]+=exon_cov_1;
			}
			else
			{
				int num_id=0;
				vector<double> key;
                                vector<double> value;
				key.push_back(Pair_vec_multi[i][j][0]);
                                key.push_back(0);
                                key.push_back(0);
                                for (int t=node_index_1;t<int(Gene_nodes[gene_index_1].size())-1;t++)
                                {
                                        if (end_pos_1>Gene_nodes[gene_index_1][t][3])
                                        {
						num_id++;
                                                key[1]=Gene_nodes[gene_index_1][t][3];
                                                key[2]=Gene_nodes[gene_index_1][t+1][2];
                                                if (Hash_edges.find(key) == Hash_edges.end())
                                                        continue;
                                                value = Hash_edges[key];
                                                double multi=double(Pair_vec_multi[i][j][5]);
						multi_1=double(Pair_vec_multi[i][j][5]);
                                                Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi;
						pair_vec.push_back(Gene_edges[value[0]][value[1]][1]);
                                        }
                                        else
                                                break;
                                }
				if (num_id == 0)
				{
					Gene_nodes[gene_index_1][node_index_1][4]+=exon_cov_1;
				}
			}//if (Gene_edges_partial[gene_index_1].size()>0)
			pair_vec.push_back(-1);
			pair_vec.push_back(node_index_2);
			if (Gene_edges_partial[gene_index_2].size()==0)
			{
				Gene_nodes[gene_index_2][node_index_2][4]+=exon_cov_2;
			}
			else
			{
				int num_id=0;
				vector<double> key;
                                vector<double> value;
				key.push_back(Pair_vec_multi[i][j][0]);
                                key.push_back(0);
                                key.push_back(0);
				for (int t=node_index_2;t<int(Gene_nodes[gene_index_2].size())-1;t++)
                                {
                                        if (end_pos_2>Gene_nodes[gene_index_2][t][3])
                                        {
						num_id++;
                                                key[1]=Gene_nodes[gene_index_2][t][3];
                                                key[2]=Gene_nodes[gene_index_2][t+1][2];
                                                if (Hash_edges.find(key) == Hash_edges.end())
                                                        continue;
                                                value = Hash_edges[key];
                                                double multi=double(Pair_vec_multi[i][j][6]);
						multi_2=double(Pair_vec_multi[i][j][6]);
                                                Gene_edges[value[0]][value[1]][2]=Gene_edges[value[0]][value[1]][2]+1/multi;
						pair_vec.push_back(Gene_edges[value[0]][value[1]][1]);
                                        }
                                        else
                                                break;
                                }
				if (num_id == 0)
				{
					Gene_nodes[gene_index_2][node_index_2][4]+=exon_cov_2;
				}
			}//if (Gene_edges_partial[gene_index_2].size()>0)
			if (gene_index_1 == gene_index_2)
			{
				if (pair_vec.size() == 4 && pair_vec[1] != pair_vec[3] && multi_1==1 && multi_2==1)
				{
					Pair_vec_junc[gene_index_1].push_back(pair_vec);
				}
				if (pair_vec.size() > 4 && multi_1==1 && multi_2==1)
				{
					Pair_vec_junc[gene_index_1].push_back(pair_vec);
				}
				if (Gene_nodes[gene_index_1].size() == 1)
				{
					Gene_edges[gene_index_1][0][0]=Gene_edges[gene_index_1][0][0]+(double(Pair_vec_multi[i][j][3])/multi_1+double(Pair_vec_multi[i][j][4])/multi_2)/double(Gene_nodes[gene_index_1][0][3]-Gene_nodes[gene_index_1][0][2]+1);
				}
			}
		}//for (int j=0;j<Pair_vec_multi[i].size();j++)
	}//for (int i=0;i<Pair_vec_multi.size();i++)
}
void Gene_info::show_pair_vec_multi()
{
	for (size_t i=0;i<Pair_vec_multi.size();i++)
	{
		for (size_t j=0;j<Pair_vec_multi[i].size();j++)
		{
			for (size_t k=0;k<Pair_vec_multi[i][j].size();k++)
			{
				cout<<Pair_vec_multi[i][j][k]<<" ";
			}
			cout<<endl;
		}
		cout<<"************************"<<endl;
	}
}
void Gene_info::show_pair_vec_junc()
{
	for (size_t i=0;i<Pair_vec_junc.size();i++)
        {
                for (size_t j=0;j<Pair_vec_junc[i].size();j++)
                {
                        for (size_t k=0;k<Pair_vec_junc[i][j].size();k++)
                        {
                                cout<<Pair_vec_junc[i][j][k]<<" ";
                        }
                        cout<<endl;
                }
                cout<<"************************"<<endl;
        }
}
void Gene_info::show_gene_pair_paths()
{
	for (size_t i=0;i<Gene_pair_paths.size();i++)
	{
		for (size_t j=0;j<Gene_pair_paths[i].size();j++)
		{
			cout<<"gene_"<<i<<": ";
			for (size_t k=0;k<Gene_pair_paths[i][j].size();k++)
			{
				cout<<Gene_pair_paths[i][j][k]<<" ";
			}
			cout<<endl;
		}
	}
}
void Gene_info::show_hash_add()
{
	map<vector<double>, vector<double> >::iterator i=Hash_left.begin();
	for(;i!=Hash_left.end();i++)
        {
		vector<double> first=i->first;
		cout<<Chr_raw[first[0]]<<": "<<first[1]<<endl;
                vector<double> second=i->second;
                cout<<second[0]<<": "<<second[1]<<endl;
        }
	i=Hash_right.begin();
	for(;i!=Hash_right.end();i++)
        {
                vector<double> first=i->first;
                cout<<Chr_raw[first[0]]<<": "<<first[1]<<endl;
                vector<double> second=i->second;
                cout<<second[0]<<": "<<second[1]<<endl;
        }
}
void Gene_info::show_hash()
{
	map<vector<double>, vector<double> >::iterator i=Hash_edges.begin();
	for(;i!=Hash_edges.end();i++)
        {
		vector<double> first=i->first;
		cout<<Chr_raw[first[0]]<<": "<<first[1]<<"->"<<first[2]<<endl;
                vector<double> second=i->second;
                cout<<second[0]<<": "<<second[1]<<endl;
        }
}
void Gene_info::show_map_need_to_do()
{
	for (size_t i=0;i<Map_need_to_do.size();i++)
	{
		for (size_t j=0;j<Map_need_to_do[i].size();j++)
		{
			for (size_t k=0;k<Map_need_to_do[i][j].size();k++)
			{
				cout<<Map_need_to_do[i][j][k]<<" ";
			}
			cout<<endl;
		}
	}
}
void Gene_info::show_genes()
{
	size_t i,j,k;
	for (i=0;i<Gene_chrs.size();i++)
	{
		cout<<"Gene_"<<i<<endl;
		cout<<"** Chr **"<<endl;
		cout<<Gene_chrs[i]<<endl;
		cout<<"** Nodes **"<<endl;
		for (j=0;j<Gene_nodes[i].size();j++)
		{
			for (k=0;k<Gene_nodes[i][j].size();k++)
			{
				cout<<Gene_nodes[i][j][k]<<" ";
			}
			cout<<endl;
		}
		cout<<"** Edges **"<<endl;
		if (Gene_nodes[i].size()==1) 
		{
			cout<<Gene_edges[i][0][0]<<endl;
		}
		else
		{
			for (j=0;j<Gene_edges[i].size();j++)
			{
				cout<<Gene_edges[i][j][0]<<" "<<Gene_edges[i][j][1]<<": "<<Gene_edges[i][j][2]<<endl;
			}
		}
		cout<<"** Partial Edges **"<<endl;
		for (j=0;j<Gene_edges_partial[i].size();j++)
		{
			cout<<Gene_edges_partial[i][j][0]<<" "<<Gene_edges_partial[i][j][1]<<endl;
		}
		cout<<"** Paths **"<<endl;
		for (j=0;j<Gene_paths[i].size();j++)
		{
			for (k=0;k<Gene_paths[i][j].size();k++)
			{
				cout<<Gene_paths[i][j][k]<<" ";
			}
			cout<<endl;
		}
		cout<<"** End **"<<endl;
	}
}

#endif

