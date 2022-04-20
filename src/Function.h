#ifndef FUNCTION_H
#define FUNCTION_H
using namespace std;
bool In_exon(int junc,vector<int>exon_l,vector<int>exon_r){
        for(int i=0;i<exon_l.size();i++){
                if(junc<=exon_r[i] && junc>=exon_l[i]) return true;
        }
        return false;
}
int Read_last_site(int Read_start_site,string Read_map){
    int num=0;
    int sum=0;
    for(int i=0;i<Read_map.length();i++){
        if(Read_map[i]>='0'&&Read_map[i]<='9'){
            int k=(int)(Read_map[i]-'0');
            sum=sum*10+k;
        }
        else{
            if(Read_map[i]=='M'||Read_map[i]=='I'||Read_map[i]=='N'||Read_map[i]=='D'){
                num=num+sum;sum=0;
            }
            else sum=0;
        }
    }
    return num+Read_start_site-1;
}

int find(int n,vector<int> v){
    for(int i=0;i<v.size();i++){
        if(v[i]==n) return i;
    }
    return -1;
}
bool is_in(int n,vector<int>v){
    for(int i=0;i<v.size();i++){
        if(n==v[i]) return true;
        else if(i==v.size()) return false;
    }
    return false;
}
bool Has_only_M(string s){
        for(int i=0;i<s.length();i++){
                if(s[i]<'9' && s[i]>'0')continue;
                else if(s[i]!='M') return false;
        }
        return true;
}
bool Has_XS_A(string s){
string::size_type idx=s.find("XS:A");
    if ( idx != string::npos ) return true;
    else return false;
}
bool Has_nh(string s){
    string::size_type idx=s.find("NH:i:");
    if ( idx != string::npos ) return true;
    else return false;
}
bool Has_chr(string s){
    string::size_type idx=s.find("chr");
    if(idx != string::npos ) return true;
    else return false;
}
vector<int> Binary(int n){
     int temp=n;
    int num;
    vector<int> binary;
    while(temp!=0){
        num=temp%2;
        binary.push_back(num);
        temp=temp/2;
    }
    for(int i=binary.size();i<11;i++) binary.push_back(0);
    return binary;
}
vector<int> seg_path(vector<int> exon_l,vector<int>exon_r,int seg_l,int seg_r){
        vector<int> node;
        int start=-1;int end=-1;
        for(int i=0;i<exon_l.size();i++){
                if(seg_l<=exon_r[i] && seg_l>=exon_l[i]) start=i;
                if(seg_r<=exon_r[i] && seg_r>=exon_l[i]) end=i;
        }
        if(start!=-1 && end!=-1){node.push_back(start);node.push_back(end);}
        return node;
}
bool cross_intron(int seg_l,int seg_r,vector<int>intron_l,vector<int>intron_r){
        for(int i=0;i<intron_l.size();i++){
                if((intron_r[i]-intron_l[i]>=0) && ((intron_l[i]>=seg_l && intron_l[i]<=seg_r) || (intron_r[i]>=seg_l && intron_r[i]<=seg_r))) {return true;}
		if((intron_r[i]-intron_l[i]<0) && ((intron_l[i]>seg_l && intron_l[i]<seg_r) || (intron_r[i]>seg_l && intron_r[i]<seg_r))) {return true;}
        }
        return false;
}
bool Has_the_junc(int juncl,int juncr,vector<int> junc_l,vector<int> junc_r){
        for(int i=0;i<junc_l.size();i++){
                if(juncl==junc_l[i]&&juncr==junc_r[i]) return true;
        }
        return false;
}

vector<int> subpath(vector<int> seg,vector<int> exon_l,vector<int>exon_r){
	vector<int> subpath;
	for(int i=0;i<seg.size()-1;){
		for(int j=0;j<exon_l.size();j++){
			if((seg[i]>=exon_l[j]&& seg[i]<=exon_r[j]) || (seg[i+1]>=exon_l[j] && seg[i+1]<=exon_r[j])){
				subpath.push_back(j);
			}
		}
	
		i=i+2;
	}
	sort(subpath.begin(),subpath.end());
	vector<int>::iterator i=unique(subpath.begin(),subpath.end());
	subpath.erase(i,subpath.end());	
	return subpath;
}

vector<int> Subpath(vector<int> v,vector<int> exon_l,vector<int>exon_r,vector<int>intron_l,vector<int>intron_r,vector<int> junc_l,vector<int>junc_r){
	vector<int> read_path;
	if(v.size()==2){
		int segl=v[0];int segr=v[1];
		if(cross_intron(segl,segr,intron_l,intron_r)){read_path.clear();}
		else{
			vector<int> read_path_temp=seg_path(exon_l,exon_r,segl,segr);
			if(read_path_temp.size()!=0 && read_path_temp[0]!=read_path_temp.back()){for(int k=read_path_temp[0];k<=read_path_temp.back();k++)read_path.push_back(k);}
			if(read_path_temp.size()!=0 && read_path_temp[0]==read_path_temp.back()){read_path.push_back(read_path_temp[0]);}
		}
	}
	else {
		for(int j=0;j<v.size()-2;){
			int segl=v[j];int segr=v[j+1];int juncl=v[j+2];int juncr=v[j+3];
//if(v[0]==33087535) cout<<segl<<" "<<segr<<" "<<juncl<<" "<<juncr<<endl;
//if(v[0]==33087535){
//	for(int b=0;b<intron_l.size();b++) cout<<intron_l[b]<<" "<<intron_r[b]<<endl;
//	if(cross_intron(segl,segr,intron_l,intron_r)) cout<<"h"<<endl;
//}
			if(cross_intron(segl,segr,intron_l,intron_r)){read_path.clear();break;}
			else{
				vector<int> read_path_temp=seg_path(exon_l,exon_r,segl,segr);
				if(read_path_temp.size()!=0){for(int k=read_path_temp[0];k<=read_path_temp.back();k++) {if(read_path.size()==0 || k!=read_path.back()) read_path.push_back(k);}}
			}
			if(!Has_the_junc(juncl,juncr,junc_l,junc_r)){read_path.clear();break;}
			else{
				for(int k=0;k<exon_l.size();k++){if(exon_l[k]==juncr) {read_path.push_back(k);}}
			}
			j=j+4;
		}
//if(v[0]==33087535) cout<<"hh "<<read_path.size()<<endl;
		int segl=v[v.size()-2];int segr=v.back();
		if(cross_intron(segl,segr,intron_l,intron_r)){read_path.clear();}
		else{
			vector<int> read_path_temp=seg_path(exon_l,exon_r,segl,segr);
			if(read_path_temp.size()!=0 && read_path.size()!=0){for(int k=read_path_temp[0];k<=read_path_temp.back();k++) {if(read_path.size()==0 || k!=read_path.back()) read_path.push_back(k);}}
		}

	}
	return read_path;
}
bool If_subpath_in_path(vector<int> subpath,vector<int> path){
	if(subpath.size()==0) return false;
	if(subpath.size()==1){
		for(int i=0;i<path.size();i++){
			if(subpath[0]==path[i]) return true;
		}
		return false;
	}
	else{
		int j=find(subpath[0],path);
		if(j==-1) return false;
		else{
			if(path.size()-j<subpath.size()) return false;
			else{
				for(int t=0;t<subpath.size();t++){
					if(subpath[t]!=path[j+t]) return false;
				}
				return true;
			}
		}
	}
}
bool If_junc_small(vector<vector<int> > junc,vector<int> junc_l,vector<int> junc_r,vector<double> junc_cov){
	vector<vector<int> > raw_junc;
	vector<int> junc_temp;
	for(int i=0;i<junc_l.size();i++){
		junc_temp.push_back(junc_l[i]);
		junc_temp.push_back(junc_r[i]);
		raw_junc.push_back(junc_temp);
		junc_temp.clear();
		
	}
	for(int i=0;i<junc.size();i++){
		for(int j=0;j<raw_junc.size();j++){
			if(junc[i]==raw_junc[j]){
				if(junc_cov[j]>1.5) return false;
			}
		}
	}
	return true;
}
#endif
