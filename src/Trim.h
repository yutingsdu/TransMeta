#include<vector>
#include<algorithm>
#define CHI_WIN 100
#define CHI_THR 100
using namespace std;
vector<double> Sorted(vector<double>v){
        for(int i=0;i<v.size()-1;i++){
            double index=i;
            for(int j=i+1;j<v.size();j++){
                if(v[j]<v[index]) {index=j;}
            }
            if(index!=i){
                double a=v[i];v[i]=v[index];v[index]=a;
            }
        }
        return v;
}
vector<double> Inserted(vector<double> v,double n){
    int index;
    if(n>=v.back()) v.push_back(n);
    else{
	if(n<=v[0]){
	    vector<double> v_t=v;
	    v.clear();
	    v.push_back(n);
	    for(int i=0;i<v_t.size();i++) v.push_back(v_t[i]);
	}
	else{	
            for(int i=v.size()-1;i>0;i--){
	        if(v[i]>n&&v[i-1]<=n){
		    index=i;
		    break;
		}
            }
	    v.push_back(0);
	    for(int i=v.size()-1;i>index;i--)
		v[i]=v[i-1];
	    v[index]=n;
	    
	}
    }
    return v;
}
double compute_chi(vector<double> winleft,vector<double> winright,double sumleft,double sumright){
    double chi=0;
    for(int i=0;i<CHI_WIN;i++){
	double mul=(winleft[i]+winright[i])/(sumleft+sumright);
	double mur=mul;
	mul*=sumleft;
	mur*=sumright;
	chi+= (winleft[i]-mul)/mul+(winright[i]-mur)/mur;
    }
    return chi;
   
}

vector<int> find_trims(int gene_l,int exon_l,int exon_r,vector<double> v_Gene){
    int c=0;
    double sumleft=0;
    double sumright=0;
    double maxsinkabundance=0;
    double maxsourceabundance=0;
    double sinkabundance=0;
    double sourceabundance=0;
    int sinkend=0;
    int sourcestart=0;
    vector<int> seg_node;
    vector<double> winleft,winright;
    if(exon_r-exon_l<2*CHI_WIN-1) return seg_node;
    for(int i=exon_l;i<=exon_r;i++){
	if(i-exon_l<2*CHI_WIN-1){
	    if(i-exon_l<CHI_WIN){
		sumleft=sumleft+v_Gene[i-gene_l];
		winleft.push_back(v_Gene[i-gene_l]);
	    }
	    else{
		sumright=sumright+v_Gene[i-gene_l];
		winright.push_back(v_Gene[i-gene_l]);
		if(i-exon_l==2*CHI_WIN-2) {
		    winleft=Sorted(winleft);
		    winright=Sorted(winright);
		}
	    }
	}
	else {
	    sumright=sumright+v_Gene[i-gene_l];
	    winright=Inserted(winright,v_Gene[i-gene_l]);//insert?
	    //winright.push_back(v_Gene[i-gene_l]);
	    double chi=0;
	    if(sumleft!=sumright) chi=compute_chi(winleft,winright,sumleft,sumright);
	    if(chi>CHI_THR) {
		if(sumleft>sumright) {
		    sinkabundance=(sumleft-sumright)/CHI_WIN;
		    if(maxsinkabundance<sinkabundance){
			sinkend=i-CHI_WIN;
			maxsinkabundance=sinkabundance;
		    }
		}
		else if(sumright>sumleft) {
		    sourceabundance=(sumright-sumleft)/CHI_WIN;
		    if(maxsourceabundance<sourceabundance) {
			sourcestart=i-CHI_WIN+1;
			maxsourceabundance=sourceabundance;
		    }
		}
	    }
	    sumleft=sumleft-v_Gene[i-gene_l-2*CHI_WIN+1];
	    vector<double>::iterator it1=find(winleft.begin(),winleft.end(),v_Gene[i-gene_l-2*CHI_WIN+1]);
	    winleft.erase(it1);
	    sumleft=sumleft+v_Gene[i-gene_l-CHI_WIN+1];
	    winleft=Inserted(winleft,v_Gene[i-gene_l-CHI_WIN+1]);//insert?
	    //winleft.push_back(v_Gene[i-gene_l-CHI_WIN+1]);
	    sumright=sumright-v_Gene[i-gene_l-CHI_WIN+1];
	    vector<double>::iterator it2=find(winright.begin(),winright.end(),v_Gene[i-gene_l-CHI_WIN+1]);
	    winright.erase(it2);
	    int d=(int)chi;
	    if(d>c)  c=d;
	    //c=(int)chi;
	}	
    }
    
    seg_node.push_back(sourcestart);
    seg_node.push_back(sinkend);
    seg_node.push_back(c);
    return seg_node;
}
