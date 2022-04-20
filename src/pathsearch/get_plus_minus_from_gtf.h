#ifndef GET_PLUS_MINUS_FROM_GTF_H
#define GET_PLUS_MINUS_FROM_GTF_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <algorithm>
#include <sstream>

using namespace std;
class Get_plus_minus
{
	private:
	const char * path_in;
	const char * path_out;
	string strand;
	public:
	Get_plus_minus(const char * p1,const char * p2,string s);
	void Process();
};
Get_plus_minus::Get_plus_minus(const char * p1,const char *p2,string s)
{
	path_in=p1;
	path_out=p2;
	strand=s;
}
void Get_plus_minus::Process()
{
ifstream ifs(path_in);
ofstream file(path_out);
string temp;
istringstream istr;
string a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12;
while (getline(ifs,temp))
{
	istr.str(temp);
	istr>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12;
	istr.clear();
	if (a7==strand && a3=="exon")
	file<<a1<<"	"<<a2<<"	"<<a3<<"	"<<a4<<"	"<<a5<<"	"<<a6<<"	"<<a7<<"	"<<a8<<"	"<<a9<<" "<<a10<<" "<<a11<<" "<<a12<<endl;
}
file.close();
ifs.close();
return;
}

#endif


