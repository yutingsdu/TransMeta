#ifndef GET_INFO_GENE_SINGLE_H
#define GET_INFO_GENE_SINGLE_H

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
class Gene_info_single
{
private:
  const char * path_of_gene_info;
  const char * path_of_bam;
public:
  Gene_info_single(const char * p1, const char *p2, double cov_uniq);
  vector<string> Gene_chrs;
  vector<vector<vector<double> > > Gene_nodes;
  vector<vector<vector<double> > > Gene_edges;
  vector<vector<vector<int> > > Gene_edges_partial;
  vector<vector<vector<int> > > Gene_paths;
  map<vector<double>, vector<double> > Hash_edges;
  map<vector<double>, vector<double> > Hash_left, Hash_right;
  vector<int> Chr_int;
  vector<string> Chr_raw;
  vector<vector<vector<int> > > Map_need_to_do; // these reads may belong to partial exons
  double Cov_uniq;
public:
  void Generate_genes();
  int Identify_cigar(string cigar);
  vector<vector<int> > Identify_cigar(int start_pos, string cigar);
  void Update_genes_from_sam();
  void process_map_need_to_do();
  vector<vector<int> > search_paths(int node_start, int node_end, int gene_id);
  void process_gene_paths();
  void show_genes();
  void show_hash();
  void show_hash_add();
  void show_cigar();
  void show_map_need_to_do();
  void variant_clear_1();
  void variant_clear_2();
  void map_clear();
};
void Gene_info_single::map_clear()
{
		Map_need_to_do.clear();
}
void Gene_info_single::variant_clear_1()
{
		Hash_left.clear();
		Hash_right.clear();
}
void Gene_info_single::variant_clear_2()
{
		Hash_edges.clear();
		Gene_edges_partial.clear();
		Gene_nodes.clear();
		Gene_edges.clear();
		Gene_paths.clear();
}
Gene_info_single::Gene_info_single(const char * p1, const char *p2, double cov_uniq)
{
	path_of_gene_info=p1;
	path_of_bam=p2;
	Cov_uniq=cov_uniq;
	return;
}
void Gene_info_single::Generate_genes()
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
                        getline(ifs,temp);
                        curr_chr=temp;
			Gene_chrs.push_back(curr_chr);
			if (Chr_raw.size() == 0)
			{
				Chr_raw.push_back(curr_chr);
				Chr_int.push_back(chr_id);
				Map_need_to_do.push_back(pair_temp);
			}
			else if (curr_chr != Chr_raw[Chr_raw.size()-1])
			{
				Chr_raw.push_back(curr_chr);
				chr_id++;
				Chr_int.push_back(chr_id);
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
int Gene_info_single::Identify_cigar(string cigar)
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

vector<vector<int> > Gene_info_single::Identify_cigar(int start_pos, string cigar)
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
void Gene_info_single::Update_genes_from_sam()
{
	string temp;
	string read_id_1,chr_id_1,map_qual_1,map_cigar_1,equal_star_1;
	int flag_1, map_start_1, pair_start_1, map_length_1;
	ifstream ifs(path_of_bam);
	istringstream istr;
	vector<int> pair_vec_L;
	long long line=0;
	while (getline(ifs,temp))
	{
		line++;
		if (temp.substr(0,3) == "@PG")
			break;
	}
	while (getline(ifs,temp))
	{
		line++;
		if (line%100000==0) cerr<<line<<endl;
		double multi_1;
		istr.str(temp);
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
		int N_count_1=count(map_cigar_1.begin(),map_cigar_1.end(),'N');
		pair_vec_L.clear();
		if (N_count_1 == 0 && multi_1 >=1)
                {
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
			Map_need_to_do[pair_vec_L[0]].push_back(pair_vec_L);
		}
		else if (N_count_1>0)
		{
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
		}//else if (N_count_1>0)
	} //while (getline(ifs,temp));
	ifs.close();
} //end of function

vector<vector<int> > Gene_info_single::search_paths(int node_start, int node_end, int gene_id)
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

void Gene_info_single::process_gene_paths() 
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

void Gene_info_single::process_map_need_to_do()
{
	int gene_index=0;
	vector<int> pair_vec;
	for (size_t i=0;i<Map_need_to_do.size();i++)
        {
		sort(Map_need_to_do[i].begin(),Map_need_to_do[i].end());
		cerr<<"Map_need_to_do "<<i<<": "<<Map_need_to_do[i].size()<<endl;
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
void Gene_info_single::show_hash_add()
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
void Gene_info_single::show_hash()
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
void Gene_info_single::show_map_need_to_do()
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
void Gene_info_single::show_genes()
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

