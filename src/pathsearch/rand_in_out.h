#ifndef RAND
#define RAND

//emcompute.h
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<stdlib.h>
#include<algorithm>
#include<cmath>
#include <ctime>
#include <cstdlib>

using namespace std;
class Rand_compute
{
	private:
	vector<double> Weights_left;
	vector<double> Weights_right;
	vector<int> Edges_left;
	vector<int> Edges_right;
	int Num_of_running;
	public:
	vector<int> best_rand;
	vector<int> Best_rand;
	double curr_num;
	long long curr_iter;
	Rand_compute(vector<double> weights_left,vector<double> weights_right,vector<int> nodes_left,vector<int> nodes_right,int num);
	~Rand_compute();
	void Search_rand();
};
Rand_compute::Rand_compute(vector<double> weights_left,vector<double> weights_right,vector<int> nodes_left,vector<int> nodes_right,int num)
{
	Weights_left=weights_left;
	Weights_right=weights_right;
	Edges_left=nodes_left;
	Edges_right=nodes_right;
	Num_of_running=num;
	return;
}
Rand_compute::~Rand_compute()
{
	return;
}
double random(double start, double end)
{
	return start+(end-start)*rand()/(RAND_MAX + 1.0);
}
void Rand_compute::Search_rand()
{
// Begin here ...
srand(unsigned(time(0)));
long long iter,max_iter;
int i,j,k,t,r,m,n;
int xx=0;
vector<int> rand_num;
vector<int> temp,Temp;
vector<int> best_index_left,best_index_right;
vector<vector<int> > Edges_rand;
vector<int>:: iterator pos;
curr_num=100000000;
//max_iter=pow(10,7);
max_iter=10000000;
vector<double> Weights_rand;
double eps=0;
n=Weights_left.size();
m=Weights_right.size();
for (i=0;i<m;i++)
{
	Weights_rand.push_back(0);
}
    vector<double> Weights_right_sort=Weights_right;
    vector<int> Weights_right_sort_index;
    sort(Weights_right_sort.begin(),Weights_right_sort.end());
    for (i=0;i<m;i++)
    {
	for (j=0;j<m;j++)
	{
		if (Weights_right_sort[i]==Weights_right[j])
		{
			Weights_right_sort_index.push_back(j);
			Weights_right[j]=0;
			break;
		}
	}
    }
    //if (pow(n,5)>max_iter) 
    if (n*n*n*n*n>max_iter)
    {
	curr_iter=max_iter;
    }
    else
    {
	//curr_iter=pow(n,5);
	curr_iter=n*n*n*n*n;
    }
    //cout<<"curr_iter="<<curr_iter<<endl;

    for (iter=0;iter<curr_iter*Num_of_running;iter++)// 20 simulations
    {
	rand_num.clear();
	for (j=0;j<n;j++)
	{
		r=int(random(0,m));
		rand_num.push_back(r);
	}
	t=0;
	for (j=0;j<m;j++)
	{
		pos=find(rand_num.begin(),rand_num.end(),j);
		if (pos!=rand_num.end())
		{
			t++;
		}
	}
	if (t==m)
	{
		xx++;
		for (j=0;j<m;j++)
		{
			double sum_temp=0;
			for (k=0;k<n;k++)
			{
				if (rand_num[k]==j)
				{
					sum_temp=sum_temp+Weights_left[k];
				}
			}
			Weights_rand[j]=sum_temp;
		}
	//sorting	
	vector<double> Weights_rand_raw=Weights_rand;
	vector<int> Weights_rand_sort_index;
	sort(Weights_rand.begin(),Weights_rand.end());
	for (j=0;j<m;j++)
	{
		for (k=0;k<m;k++)
		{
			if (Weights_rand[j]==Weights_rand_raw[k])
			{
				Weights_rand_sort_index.push_back(k);
				Weights_rand_raw[k]=0;
				break;
			}
		}
	}
	eps=0;
	for (j=0;j<m;j++)
	{
		//eps=eps+(Weights_rand[j]-Weights_right_sort[j])*(Weights_rand[j]-Weights_right_sort[j]);
		eps=eps+(1-Weights_rand[j]/Weights_right_sort[j])*(1-Weights_rand[j]/Weights_right_sort[j]);
	}
	if (eps<curr_num)
	{
		best_rand=rand_num;
		curr_num=eps;
		best_index_left=Weights_rand_sort_index;
		best_index_right=Weights_right_sort_index;
	}
	rand_num.clear();
	}//if (t==m)
    }// simulation num 20
Best_rand=best_rand;
for (j=0;j<n;j++)
{
	for (k=0;k<m;k++)
	{
		if (best_index_left[k]==best_rand[j])
		{
			Best_rand[j]=Edges_right[best_index_right[k]];
		}
	}
}
}//end of final function
#endif







