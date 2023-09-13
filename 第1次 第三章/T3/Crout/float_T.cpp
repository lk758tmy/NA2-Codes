#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory.h>
#include <chrono>
#include <cblas.h>
using namespace std;
using namespace std::chrono;
void matrix_make(float *h,int n){
	float *p=h;
	memset(h,0,n*n*sizeof(float));
	*h=4;
	for(int i=1;i<n;i++){
		p+=(n+1); *p=4;
		*(p-1)=*(p-n)=-1;
	}
	return ;
}
void matrix_print(float *h,int n){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			printf("%.4f ",*h);
			h++;
		}
		printf("\n");
	} printf("\n");
	return ;
}
void GJ_rev(float *m,int n,bool colPivot){
	int maxRowId,rowId;
	float *mp=m,*mp2=m,maxElement,tmp;
	for(int i=0;i<n;i++){
		rowId=i;
		if(colPivot==true){//启用列主元法
            mp=m+n*rowId+i;
            maxElement=abs(*mp); maxRowId=rowId;
            for(int j=i+1;j<n;j++){
                mp+=n;
                if(abs(*mp)>maxElement){
                    maxElement=abs(*mp); maxRowId=j;
                }
            }
            cblas_sswap(n,m+n*rowId,1,m+n*maxRowId,1);
		}
		mp=m+rowId; tmp=-1/(*(m+rowId*n+rowId));
		for(int j=0;j<n;j++,mp+=n){
            if(j==rowId) *mp=-tmp;
            else (*mp)*=tmp;
		}
		mp=m; mp2=m+n*rowId;
        for(int j=0;j<n;j++,mp+=n){
            if(j==rowId) continue;
            tmp=*(mp+rowId); *(mp+rowId)=0;
            cblas_saxpy(n,tmp,mp2,1,mp,1);
            *(mp+rowId)=tmp;
        }
        mp=m+n*rowId; tmp=*(mp+rowId);
        for(int j=0;j<n;j++,mp++){
            if(j==rowId) continue;
            (*mp)*=tmp;
        }
	}
	return ;
}
int main() {
	int n,n2;
	scanf("%d",&n); n2=n*n; //n should be at least 3!
	float *m=(float *)malloc(n2*sizeof(float));
	matrix_make(m,n);
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	GJ_rev(m,n,false);
	e=steady_clock::now();
	//matrix_print(m,n);
	printf("Time(ms): %.3f",(e-s).count()/1000000.0);
	free(m);
	return 0;
}
