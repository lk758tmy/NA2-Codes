#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory.h>
#include <chrono>
#include <cblas.h>
using namespace std;
using namespace std::chrono;
void matrix_make(float *h,int n){
	int n2=n*n; float *p=h;
	memset(h,0,n2*n2*sizeof(float));
	*h=4;
	for(int i=1;i<n2;i++){
		p+=(n2+1); *p=4;
		if(i%n!=0) *(p-1)=*(p-n2)=-1;
	}
	for(int i=n,j=0;i<n2;i++,j++)
		*(h+i*n2+j)=*(h+j*n2+i)=-1;
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
        /*tmp=1/(*(m+n*rowId+rowId));
        for(int j=0;j<n;j++,mp++) (*mp)*=tmp;
        tmp=-tmp; mp=m; mp2=m+n*rowId;
        for(int j=0;j<n;j++,mp+=n){
            if(j==rowId) continue;
            cblas_saxpy(n,-*(mp+rowId),mp2,1,mp,1);
        }*/ //GJ消元化为单位矩阵，但此种方法需要另申请内存将I变换为A逆
	}
	return ;
}
int _main(int n) {
	int n2=n*n; int n4=n2*n2;
	float *m=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	GJ_rev(m,n2,false);
	e=steady_clock::now();
	//matrix_print(m,n2);
	printf("Size: %d\nTime(ms): %.3f\n\n",n,(e-s).count()/1000000.0);
	free(m);
	return 0;
}
int main(){
    for(int n=5;n<=50;n+=5) _main(n);
    return 0;
}
