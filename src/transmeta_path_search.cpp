#include "./metapathsearch/meta_path_search.h"
#include "./metapathsearch/meta_merge_sample.h"
ofstream out_gtf;
ofstream out_info;
string strand;
bool unstranded = false;

int main (int argc, char* argv[])
{
    string Strand = argv[4];
    if(Strand == "unstranded") unstranded = true;
    string dir = argv[3];

    string gtf_gile = dir + "/transmeta-temp.gtf";
    out_gtf.open(gtf_gile.c_str(),ios::app);

    string info_file = dir + "/transmeta.info";
    out_info.open(info_file.c_str(),ios::app);

    MergeSample ms(argv[1],argv[2]);//graph-list and graph-Dir
    ms.process("+");
    MergeSample ms2(argv[1],argv[2]);//graph-list and graph-Dir
    ms2.process("-");

    out_gtf.close();
    out_info.close();
    return 0;
}
