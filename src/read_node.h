#include<vector>
using namespace std;
int read_in_node(int read_start,vector<int> exon_l,vector<int> exon_r){
        for(int i=0;i<exon_l.size();i++){
                if(read_start>=exon_l[i] && read_start<=exon_r[i]){ return i;}
        }
	return -1;
}

