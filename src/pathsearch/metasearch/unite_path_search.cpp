#include "./path_search.h"
#include "./unite_path_search.h"
ofstream out_gtf;
ofstream out_info;
string strand;
bool unstranded = false;

int main (int argc, char* argv[])
{
    string Strand = argv[4];
    if(Strand == "unstranded") unstranded = true;
    string file = argv[3];
    out_gtf.open(argv[3],ios::app);
    file = file.substr(0,file.length() - 4) + ".info";
    out_info.open(file.c_str(),ios::app);

    MergeSample ms(argv[1],argv[2]);//graph-list and graph-Dir
    ms.process("+");
    MergeSample ms2(argv[1],argv[2]);//graph-list and graph-Dir
    ms2.process("-");

    out_gtf.close();
    out_info.close();
    return 0;
}
