#ifndef GET_JUNCTION_GRAPH
#define GET_JUNCTION_GRAPH

//get_junction_graph.h

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include<algorithm>
#include<cmath>
#include"rand_in_out.h"
using namespace std;

class Get_Junction_Graph
{
	private:
	vector<vector<int> > Vec_edges;
	vector<double> Vec_weights;
	vector<vector<int> > Pair_paths;
	int Node_num;
	public:
	Get_Junction_Graph(vector<vector<int> > vec_edges,vector<double> vec_weights,vector<vector<int> > pair_paths,int node_num);
	~Get_Junction_Graph();
	void Construct_Junction_Graph();// attention: an edge in junction graph represents two conseccutive edges in Raw graph.
	vector<int> Edges_left;// left node of an edge in junction graph
	vector<int> Edges_right;// right node of an edge in junction graph
	vector<double> Weights;//weights of an edge in junction graph
	vector<int> DAG_left;// denote the two nodes( left one) of a junction edge (a node in junction graph).
	vector<int> DAG_right;// denote the two nodes (right one) of a junction edge (a node in junction graph).
	vector<double> DAG_weights;
	vector<int> Cons_left;
	vector<int> Cons_right;
	vector<int> S,T;
	vector<int> Single_nodes;
	int junction_num,max_node_num;
};
Get_Junction_Graph::Get_Junction_Graph(vector<vector<int> > vec_edges,vector<double> vec_weights,vector<vector<int> > pair_paths,int node_num)
{
	Vec_edges=vec_edges;
	Vec_weights=vec_weights;
	Pair_paths=pair_paths;
	Node_num=node_num;
	return;
}
Get_Junction_Graph::~Get_Junction_Graph()
{
	return;
}
void Get_Junction_Graph::Construct_Junction_Graph()
{
	string temp;
	int i,j,k,t;
	int max_num,temp_max,power;
	double min_weights,sum_weights,temp_sum;
	vector<string> edge_left,edge_right,weight;
	vector<double> temp_weights;
	vector<int> Pair_edges,Pair_nodes,Pair_num;
	vector<int>::iterator pos_left;
	vector<int>::iterator pos_right;
	vector<int>::iterator pos;
	vector<int> temp_in,temp_out,temp_cons;
	vector<vector<int> > Cons_edges;
	vector<vector<int> >::iterator iter;
	vector<double> TEMP_weights;
	int index_in,index_out;
	double index_weights=100000000000000000;
	junction_num=-1;
	max_node_num=-1;
    if (Vec_edges.size()==0)
    {
	for (i=0;i<Node_num;i++)
	{
		Single_nodes.push_back(i);
	}
    }
    if (Vec_edges.size()>=1)
    {
	// Processing pairs
        for (i=0;i<int(Vec_edges.size());i++)
        {
                DAG_left.push_back(Vec_edges[i][0]);
                DAG_right.push_back(Vec_edges[i][1]);
                DAG_weights.push_back(Vec_weights[i]);
        }
	vector<vector<int> > Cons_edges;
	for (i=0;i<int(Pair_paths.size());i++)
	{
		for (j=0;j<int(Pair_paths[i].size())-2;j++)
		{
			int left_id=-1;
			int right_id=-1;
			vector<int> cons_temp;
			for (k=0;k<int(DAG_left.size());k++)
			{
				if (DAG_left[k]==Pair_paths[i][j] && DAG_right[k]==Pair_paths[i][j+1])
				{
					left_id=k;
				}
				if (DAG_left[k]==Pair_paths[i][j+1] && DAG_right[k]==Pair_paths[i][j+2])
				{
					right_id=k;
				}
				if (left_id != -1 && right_id != -1)
				{
					break;
				}
			}
			cons_temp.push_back(left_id);
			cons_temp.push_back(right_id);
			Cons_edges.push_back(cons_temp);
		}
	}
	sort(Cons_edges.begin(),Cons_edges.end());
	vector<vector<int> >::iterator iter=unique(Cons_edges.begin(),Cons_edges.end());
	Cons_edges.erase(iter,Cons_edges.end());
	for (i=0;i<int(Cons_edges.size());i++)
        {
                Cons_left.push_back(Cons_edges[i][0]);
                Cons_right.push_back(Cons_edges[i][1]);
        }
        Cons_edges.clear();
	// Processing each node for in-coming and out-going
	for (i=0;i<Node_num;i++)
	{
		temp_in.clear();
		temp_out.clear();
		for (j=0;j<int(DAG_weights.size());j++)
		{
			if (DAG_right[j]==i)
			{
				temp_in.push_back(j);
			}
			else if (DAG_left[j]==i)
			{
				temp_out.push_back(j);
			}
		}
                if (temp_in.size()>1 && temp_out.size()>1)
                {
                    if (temp_in.size()*temp_out.size()>20 && temp_in.size()>=temp_out.size())
                    {
                        vector<double> weights_in,weights_out;
                        for (j=0;j<int(temp_in.size());j++)
                        {
                                weights_in.push_back(DAG_weights[temp_in[j]]);
                        }
                        for (j=0;j<int(temp_out.size());j++)
                        {
                                weights_out.push_back(DAG_weights[temp_out[j]]);
                        }
                        Rand_compute rand(weights_in,weights_out,temp_in,temp_out,3);
                        rand.Search_rand();
                        vector<int> Best_rand=rand.Best_rand;
                        vector<vector<int> > Packing_results;
                        vector<int> packing_temp;
                        vector<vector<int> >::iterator packing_pos;
                        for (j=0;j<int(temp_in.size());j++)
                        {
                                packing_temp.push_back(temp_in[j]);
                                packing_temp.push_back(Best_rand[j]);
                                Packing_results.push_back(packing_temp);
                                packing_temp.clear();
                        }
                        for (j=0;j<int(temp_in.size());j++)
                        {
                                for (k=0;k<int(temp_out.size());k++)
                                {
                                        Edges_left.push_back(temp_in[j]);
                                        Edges_right.push_back(temp_out[k]);
                                        packing_temp.push_back(temp_in[j]);
                                        packing_temp.push_back(temp_out[k]);
                                        packing_pos=find(Packing_results.begin(),Packing_results.end(),packing_temp);
                                        packing_temp.clear();
                                        if (packing_pos == Packing_results.end())
                                        {
                                                Weights.push_back(1);
                                        }
                                        else
                                        {
                                                Weights.push_back(0);
                                        }
                                }
                        }
                    }
                    if (temp_in.size()*temp_out.size()>20 && temp_in.size()<temp_out.size())
                    {
                        vector<double> weights_in,weights_out;
                        for (j=0;j<int(temp_in.size());j++)
                        {
                                weights_in.push_back(DAG_weights[temp_in[j]]);
                        }
                        for (j=0;j<int(temp_out.size());j++)
                        {
                                weights_out.push_back(DAG_weights[temp_out[j]]);
                        }
                        Rand_compute rand(weights_out,weights_in,temp_out,temp_in,3);
                        rand.Search_rand();
                        vector<int> Best_rand=rand.Best_rand;
                        vector<vector<int> > Packing_results;
                        vector<int> packing_temp;
                        vector<vector<int> >::iterator packing_pos;
                        for (j=0;j<int(temp_out.size());j++)
                        {
                                packing_temp.push_back(Best_rand[j]);
                                packing_temp.push_back(temp_out[j]);
                                Packing_results.push_back(packing_temp);
                                packing_temp.clear();
                        }
                        for (j=0;j<int(temp_in.size());j++)
                        {
                                for (k=0;k<int(temp_out.size());k++)
                                {
                                        Edges_left.push_back(temp_in[j]);
                                        Edges_right.push_back(temp_out[k]);
                                        packing_temp.push_back(temp_in[j]);
                                        packing_temp.push_back(temp_out[k]);
                                        packing_pos=find(Packing_results.begin(),Packing_results.end(),packing_temp);
                                        packing_temp.clear();
                                        if (packing_pos == Packing_results.end())
                                        {
                                                Weights.push_back(1);
                                        }
                                        else
                                        {
                                                Weights.push_back(0);
                                        }
                                }
                        }
                    }
                    if (temp_in.size()*temp_out.size()<=20)
                    {
                        if (temp_in.size()>temp_out.size())
                        {
                                temp_max=temp_in.size();
                        }
                        else
                        {
                                temp_max=temp_out.size();
                        }
                        power=1;
                        for (j=0;j<int(temp_in.size()*temp_out.size());j++)
                        {
                                power=2*power;
                        }
                        for (j=0;j<power;j++)
                        {
                                temp_weights.clear();
                                int a=j;
                                for (k=0;a>0;k++)
                                {
                                        temp_weights.push_back(a%2);
                                        a/=2;
                                }
                                if (k!=int(temp_in.size()*temp_out.size()))
                                {
                                        for (t=0;t<int(temp_in.size()*temp_out.size())-k;t++)
                                        {
                                                temp_weights.push_back(0);
                                        }
                                }
                                temp_sum=0;
                                for (k=0;k<int(temp_weights.size());k++)
                                {
                                        temp_sum=temp_sum+temp_weights[k];
                                }
                                if (temp_sum==temp_max)
                                {
                                        index_in=0;
                                        index_out=0;
                                        for (k=0;k<int(temp_in.size());k++)
                                        {
                                                temp_sum=0;
                                                for (t=0;t<int(temp_out.size());t++)
                                                {
                                                        temp_sum=temp_sum+temp_weights[temp_out.size()*k+t];
                                                }
                                                if (temp_sum==0)
                                                {
                                                        break;
                                                }
                                                else
                                                {
                                                        index_in++;
                                                }
                                        }
                                        for (k=0;k<int(temp_out.size());k++)
                                        {
                                                temp_sum=0;
                                                for (t=0;t<int(temp_in.size());t++)
                                                {
                                                        temp_sum=temp_sum+temp_weights[k+temp_out.size()*t];
                                                }
                                                if (temp_sum==0)
                                                {
                                                        break;
                                                }
                                                else
                                                {
                                                        index_out++;
                                                }
                                        }
                                        if (index_in==int(temp_in.size()) && index_out==int(temp_out.size()))
                                        {
                                                min_weights=0;
                                                if (temp_in.size()>=temp_out.size())
                                                {
                                                        for (k=0;k<int(temp_out.size());k++)
                                                        {
                                                                sum_weights=0;
                                                                for (t=0;t<int(temp_in.size());t++)
                                                                {
                                                                        sum_weights=sum_weights+temp_weights[k+temp_out.size()*t]*DAG_weights[temp_in[t]];
                                                                }
                                                                min_weights=min_weights+(1-sum_weights/DAG_weights[temp_out[k]])*(1-sum_weights/DAG_weights[temp_out[k]]);
                                                        }
                                                        if (min_weights<index_weights)
                                                        {
                                                                index_weights=min_weights;
                                                                TEMP_weights.clear();
                                                                TEMP_weights=temp_weights;
                                                        }
                                                }
                                                else
                                                {
                                                        for (k=0;k<int(temp_in.size());k++)
                                                        {
                                                                sum_weights=0;
                                                                for (t=0;t<int(temp_out.size());t++)
                                                                {
                                                                        sum_weights=sum_weights+temp_weights[temp_out.size()*k+t]*DAG_weights[temp_out[t]];
                                                                }
                                                                min_weights=min_weights+(1-sum_weights/DAG_weights[temp_in[k]])*(1-sum_weights/DAG_weights[temp_in[k]]);
                                                        }
                                                        if (min_weights<index_weights)
                                                        {
                                                                index_weights=min_weights;
                                                                TEMP_weights.clear();
                                                                TEMP_weights=temp_weights;
                                                        }
                                                }
                                        }
                                }//end of if (temp_sum==temp_max)
                        }//end of for (j=0;j<pow((double)2,temp_in.size()*temp_out.size());j++)
                        for (j=0;j<int(temp_in.size());j++)
                        {
                                for (k=0;k<int(temp_out.size());k++)
                                {
                                        Edges_left.push_back(temp_in[j]);
                                        Edges_right.push_back(temp_out[k]);
                                }
                        }
                        for (j=0;j<int(TEMP_weights.size());j++)
                        {
                                Weights.push_back(1-TEMP_weights[j]);
                        }
                        index_weights=100000000000000000;
                    }
                }//end of if (temp_in.size()>1 && temp_out.size()>1)

		else if ((temp_in.size()==1 || temp_out.size()==1) && (temp_in.size()>0 && temp_out.size()>0))
		{
			for (j=0;j<int(temp_in.size());j++)
			{
				for (k=0;k<int(temp_out.size());k++)
				{
					Edges_left.push_back(temp_in[j]);
					Edges_right.push_back(temp_out[k]);
					Weights.push_back(0);
				}
			}
		}
		else if (temp_in.size()==0 && temp_out.size()>0)
		{
			for (j=0;j<int(temp_out.size());j++)
			{
				S.push_back(temp_out[j]);
				pos=find(DAG_left.begin(),DAG_left.end(),DAG_right[temp_out[j]]);
				if (pos==DAG_left.end() && temp_out[j]>junction_num)
				{
					junction_num=temp_out[j];
				}
			}
		}
		else if (temp_in.size()>0 && temp_out.size()==0)
		{
			for (j=0;j<int(temp_in.size());j++)
			{
				T.push_back(temp_in[j]);
				pos=find(DAG_right.begin(),DAG_right.end(),DAG_left[temp_in[j]]);
				if (pos==DAG_right.end() && temp_in[j]>junction_num)
				{
					junction_num=temp_in[j];
				}
			}
		}
		else if (temp_in.size()==0 && temp_out.size()==0)
		{
			Single_nodes.push_back(i);
		}
	}
    }//end of if (vec_edges.size()>=1)
/*
//////////////////////////////////////////////////////////////////////////////////////
        cout<<"**** Output ****"<<endl;
	cout<<"DAG Edges:"<<endl;
	for (i=0;i<DAG_weights.size();i++)
	{
		cout<<DAG_left[i]<<"->"<<DAG_right[i]<<": "<<DAG_weights[i]<<endl;
	}
	cout<<"***************************************"<<endl;
	cout<<"Junction Edges:"<<endl;
	for (i=0;i<Edges_left.size();i++)
	{
		cout<<Edges_left[i]<<"->"<<Edges_right[i]<<": "<<Weights[i]<<endl;
	}
	cout<<"S:"<<endl;
	for (i=0;i<S.size();i++)
	{
		cout<<S[i]<<"; ";
	}
	cout<<endl;
	cout<<"T:"<<endl;
	for (i=0;i<T.size();i++)
	{
		cout<<T[i]<<"; ";
	}
	cout<<endl;
	cout<<"Single Nodes:"<<endl;
	for (i=0;i<Single_nodes.size();i++)
	{
		cout<<Single_nodes[i]<<"; ";
	}
	cout<<endl;
*/
} //end of "Construct_Junction_Graph()"

#endif























