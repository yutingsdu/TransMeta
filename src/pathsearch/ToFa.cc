#include "gtf2fa.h"

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;
int main(int argc,char* argv[])
{
    GTF2FASTA gtf2fasta(argv[1],argv[2],argv[3]);
    gtf2fasta.process();

    //time ../src/newh/ToFa TransComb.gtf /storage/juntaosdu/yuting/New_transcomb/Human_single/genome.fa TransComb_sctie.fa
}
