#ifndef SORT_GTF_H
#define SORT_GTF_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iterator>

using namespace std;
class Sort_GTF
{
	private:
	const char * path_in;
	const char * path_out;
	public:
	Sort_GTF(const char * p1,const char * p2);
	void Process();
};
Sort_GTF::Sort_GTF(const char * p1,const char *p2)
{
	path_in=p1;
	path_out=p2;
}
void Sort_GTF::Process()
{
string temp;
vector<string> Chrs;
int t=0;
vector<vector<string> > Trans_chr;
vector<string> Trans_chr_temp;
/*
ifstream ifs_chr(path_chr);
while (getline(ifs_chr,temp))
{
	Chrs[temp]=t;
	t++;
	Trans_chr.push_back(Trans_chr_temp);
}
ifs_chr.close();
*/
ifstream ifs_chr(path_in);
istringstream istr;
string a1,a2,a3,a6,a7,a8,a9,a10,a11,a12;
int a4,a5;
getline(ifs_chr,temp);
istr.str(temp);
istr>>a1;
istr.clear();
string curr_chr=a1;
Chrs.push_back(a1);
while (getline(ifs_chr,temp))
{
	istr.str(temp);
	istr>>a1;
	istr.clear();
	if (a1 != curr_chr)
	{
		Chrs.push_back(a1);
		curr_chr=a1;
	}
}
ifs_chr.close();
sort(Chrs.begin(),Chrs.end());
vector<string>::iterator pos;
pos=unique(Chrs.begin(),Chrs.end());
Chrs.erase(pos,Chrs.end());
map<string,int > Hash_Chrs;
for (size_t i=0;i<Chrs.size();i++)
{
	Hash_Chrs[Chrs[i]]=i;
	Trans_chr.push_back(Trans_chr_temp);
}

ofstream file(path_out);
ifstream ifs(path_in);
while (getline(ifs,temp))
{
	istr.str(temp);
	istr>>a1;
	istr.clear();
	Trans_chr[Hash_Chrs[a1]].push_back(temp);
}
ifs.close();

for (size_t i=0;i<Trans_chr.size();i++)
{
	map<string,int > hash_trans_ids;
	vector<vector<int> > trans_start_pos;
	vector<int> start_temp;
	vector<string> vec_trans_ids, vec_trans_ids_sorted;
	map<string,vector<string> > Hash_id_all;
	int k=0;
	int mode=1;
	for (size_t j=0;j<Trans_chr[i].size();j++)
	{
		istr.str(Trans_chr[i][j]);
		istr>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12;
		istr.clear();
		vector<string> hash_temp;
		if (Hash_id_all.find(a12) == Hash_id_all.end())
		{
			hash_temp.push_back(Trans_chr[i][j]);
			Hash_id_all[a12]=hash_temp;
		}
		else
		{
			hash_temp=Hash_id_all[a12];
			hash_temp.push_back(Trans_chr[i][j]);
			Hash_id_all[a12]=hash_temp;
		}
		if (hash_trans_ids.find(a12) == hash_trans_ids.end())
		{
			hash_trans_ids[a12]=a4;
		}
		else if (a4 < hash_trans_ids[a12])
		{
			hash_trans_ids[a12]=a4;
		}
/*
		if (hash_trans_ids.find(a12) == hash_trans_ids.end())
		{
			hash_trans_ids[a12]=a4;
			vec_trans_ids.push_back(a12);
			start_temp.push_back(a4);
			start_temp.push_back(k);
			trans_start_pos.push_back(start_temp);
			start_temp.clear();
			k++;
		}
*/
	}
	map<string,int >::iterator iter=hash_trans_ids.begin();
	for(;iter!=hash_trans_ids.end();iter++)
        {
		vec_trans_ids.push_back(iter->first);
		start_temp.push_back(iter->second);
		start_temp.push_back(k);
		trans_start_pos.push_back(start_temp);
		start_temp.clear();
		k++;
        }
	sort(trans_start_pos.begin(),trans_start_pos.end());
	for (size_t j=0;j<trans_start_pos.size();j++)
	{
		vec_trans_ids_sorted.push_back(vec_trans_ids[trans_start_pos[j][1]]);
	}
	for (size_t j=0;j<vec_trans_ids_sorted.size();j++)
	{
		vector<string> hash_temp;
                hash_temp=Hash_id_all[vec_trans_ids_sorted[j]];
		if (hash_temp.size()>=2)
		{
			istr.str(hash_temp[0]);
                        istr>>a1>>a2>>a3>>a4;
                        istr.clear();
                        istr.str(hash_temp[1]);
                        istr>>a1>>a2>>a3>>a5;
                        istr.clear();
			if (a4<a5)
			{
				mode=1;
			}
			else
			{
				mode=2;
			}
			break;
		}
	}
	for (size_t j=0;j<vec_trans_ids_sorted.size();j++)
	{
		vector<string> hash_temp;
		hash_temp=Hash_id_all[vec_trans_ids_sorted[j]];
		if (hash_temp.size() == 1)
		{
			file<<hash_temp[0]<<endl;
		}
		else 
		{
			if (mode==1)
			{
				for (size_t k=0;k<hash_temp.size();k++)
				{
					file<<hash_temp[k]<<endl;
				}
			}
			else
			{
				for (size_t k=0;k<hash_temp.size();k++)
                                {
                                        file<<hash_temp[hash_temp.size()-k-1]<<endl;
                                }
			}
		}
	}
}

file.close();
return;
}

#endif


