//#include"expression_level.h"

//#include"get_junction_graph_new.h"
//#include"junction_paths_new.h"
//#include"recover_paths.h"
//#include"PairPacker.h"
#include "path_search.h"
using namespace std;

string out_name;
string mode;
int pack_graph_num=0;
int unpack_graph_num=0;
bool SFlag = false;
bool PackingFlag = false;
bool MyFlag = false;
int Path_length = 500;
double AVERAGE_REMOVE_RATE = 0.01;//0.03
double UNBALANCE_RATE = 0.03;
float SEED_Filter = 1.01;
int rg_index = 0;
int trans_id;
vector<pair<int,int> > good_junction;
int main (int argc, char* argv[])
{
    //SFlag = true;
    mode = "I";
    out_name = argv[2];
    load_graph(argv[1]);
    //cout<<pack_graph_num<<" "<<unpack_graph_num<<endl;
    return 0;
}
