#ifndef COMBINE_GTF
#define COMBINE_GTF

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<map>
#include<set>
#include<iterator>
#include<stdlib.h>

using namespace std;
class Combine_gtf
{
	private:
	const char * path_1;
	const char * path_2;
	const char * path_3;
	public:
	ofstream file;
	Combine_gtf(const char * p1,const char * p2,const char * p3);
	void Process();

};
Combine_gtf::Combine_gtf(const char * p1,const char * p2,const char * p3)
{
	path_1=p1;
	path_2=p2;
	path_3=p3;
	file.open(path_3);
}
void Combine_gtf::Process()
{
string temp;
ifstream ifs(path_1);
while (getline(ifs,temp))
{
	file<<temp<<endl;
}
ifs.close();
ifs.open(path_2);
while (getline(ifs,temp))
{
	file<<temp<<endl;
}
ifs.close();
return;
}

#endif

















