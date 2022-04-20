#ifndef RECOVER_PATHS
#define RECOVER_PATHS

//recover_paths.h
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include<algorithm>
#include<cmath>
using namespace std;

class Recover_Junction_Paths
{
	private:
	vector<int> DAG_left;
	vector<int> DAG_right;
	vector<double> DAG_weights;
	vector<vector<int> > Junction_Path_cover;
	int max_node_num;
	vector<int> Node_left, Node_right;
	public:
	Recover_Junction_Paths(vector<int> dag_left,vector<int> dag_right,vector<double> dag_weights,vector<vector<int> > junction_path_cover,int max_node,vector<int> node_left,vector<int> node_right);
	~Recover_Junction_Paths();
	void Get_Path_cover();
	void Get_seeds();
	vector<vector<int> > Path_Cover;
	vector<int> Seeds; //index of egdes
	int raw_size;
};
Recover_Junction_Paths::Recover_Junction_Paths(vector<int> dag_left,vector<int> dag_right,vector<double> dag_weights,vector<vector<int> > junction_path_cover,int max_node,vector<int> node_left,vector<int> node_right)
{
	DAG_left=dag_left;
	DAG_right=dag_right;
	DAG_weights=dag_weights;
	Junction_Path_cover=junction_path_cover;
	max_node_num=max_node;
	Node_left=node_left;
	Node_right=node_right;
	return;
}
Recover_Junction_Paths::~Recover_Junction_Paths()
{
	return;
}
void Recover_Junction_Paths::Get_Path_cover()
{
	size_t i,j;
	vector<int> Path_temp;
	for (i=0;i<Junction_Path_cover.size();i++)
	{
		Path_temp.push_back(DAG_left[Junction_Path_cover[i][0]]);
		for (j=0;j<Junction_Path_cover[i].size();j++)
		{
			Path_temp.push_back(DAG_right[Junction_Path_cover[i][j]]);
		}
		Path_Cover.push_back(Path_temp);
		Path_temp.clear();
	}
	return;
}

#endif


