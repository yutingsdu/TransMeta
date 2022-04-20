#include "abundance.h"
#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

int main(int argc,char* argv[])
{
    Abundance abundance(argv[1],argv[2],argv[3],argv[4]);
    abundance.Process();

    //../src/newh/Estimate ./k_results_Genome/abundance.tsv TransComb.gtf TransComb_final.gtf 1 

}
