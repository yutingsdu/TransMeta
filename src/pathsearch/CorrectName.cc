#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
using namespace std;
int main(int argc,char*argv[]){
	if (argc!=6)
	{
		cout<<"====================================================================="<<endl;
		cout<<"Usage of CorrectName:"<<"\n"<<endl;
		cout<<"  CorrectName <fa|fq> <reads_1> <reads_2> <processed_reads_1> <processed_reads_2>"<<"\n"<<endl;
		cout<<"** Required **"<<"\n"<<endl;
		cout<<"  <fa|fq>: RNA-seq data type."<<"\n"<<endl;
		cout<<"         If it is stored by fasta format, please type <fa>"<<endl;
		cout<<"         If it is stored by fastq format, please type <fq>"<<"\n"<<endl;
		cout<<"  <reads_1>: Path to left RNA-seq reads."<<"\n"<<endl;
		cout<<"  <reads_2>: Path to right RNA-seq reads."<<"\n"<<endl;
		cout<<"  <processed_reads_1>: Path to write modified left reads."<<"\n"<<endl;
		cout<<"  <processed_reads_2>: Path to write modified right reads."<<"\n"<<endl;
		cout<<"** Note **"<<endl;
		cout<<"If fastq RNA-seq data: reads_1.fq and reads_2.fq need to process"<<endl;
		cout<<"The command would be:"<<"\n"<<endl;
		cout<<"./CorrectName fq reads_1.fq reads_2.fq New_reads_1.fq New_reads_2.fq"<<endl;
		cout<<"====================================================================="<<endl;
		return 0;
	}
	string data_type=argv[1];
	if (data_type=="fa")
	{
		ifstream in_1(argv[2]);
		ifstream in_2(argv[3]);
		ofstream out_1(argv[4]);
		ofstream out_2(argv[5]);
		string temp_1,temp_2,name_1,name_2;
		istringstream istr_1,istr_2;
		while (getline(in_1,temp_1))
		{
			getline(in_2,temp_2);
			istr_1.str(temp_1);
			istr_2.str(temp_2);
			istr_1>>name_1;
			istr_2>>name_2;
			istr_1.clear();
			istr_2.clear();
			name_1[name_1.size()-2]='/';
			name_2[name_2.size()-2]='/';
			out_1<<name_1<<endl;
			out_2<<name_2<<endl;
			getline(in_1,temp_1);
			getline(in_2,temp_2);
			out_1<<temp_1<<endl;
                        out_2<<temp_2<<endl;
		}
	}
	else if (data_type=="fq")
	{
		ifstream in_1(argv[2]);
		ifstream in_2(argv[3]);
		ofstream out_1(argv[4]);
		ofstream out_2(argv[5]);
		string temp_1,temp_2,name_1,name_2;
		istringstream istr_1,istr_2;
		while (getline(in_1,temp_1))
		{
			getline(in_2,temp_2);
			istr_1.str(temp_1);
			istr_2.str(temp_2);
			istr_1>>name_1;
			istr_2>>name_2;
			istr_1.clear();
			istr_2.clear();
			name_1[name_1.size()-2]='/';
			name_2[name_2.size()-2]='/';
			out_1<<name_1<<endl;
			out_2<<name_2<<endl;//@name
			getline(in_1,temp_1);
			getline(in_2,temp_2);
			out_1<<temp_1<<endl;
                        out_2<<temp_2<<endl;//ATCG sequence
			getline(in_1,temp_1);
                        getline(in_2,temp_2);
                        out_1<<"+"<<endl;
                        out_2<<"+"<<endl;//@
			getline(in_1,temp_1);
                        getline(in_2,temp_2);
                        out_1<<temp_1<<endl;
                        out_2<<temp_2<<endl;//sequencing quality
		}
	}//else if (data_type=="fq")
}









