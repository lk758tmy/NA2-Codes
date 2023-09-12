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
void vector_make(float *h,int n){
	*h=2;
	for(int i=1;i<n-1;i++){h++; *h=1;}
	h++; *h=2;
	for(int i=1;i<n-1;i++){
		h++; *h=1;
		for(int j=1;j<n-1;j++){h++; *h=0;}
		h++; *h=1;
	}
	h++; *h=2;
	for(int i=1;i<n-1;i++){h++; *h=1;}
	h++; *h=2;
	return ;
}
void vector_print(float *h,int n){
	for(int i=0;i<n;i++){
		printf("%.5f ",*h); h++;
	} printf("\n\n");
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
void gauss_colPivot(float *m,float *x,float *c,int n){
	//int *p=(int *)malloc(n*sizeof(int));
	int /**pp=p,*/maxRowId,rowId;
	float *mp=m,*cp=c,maxElement,tmp;
	//for(int i=0;i<n;i++){*pp=i; pp++;}
	for(int i=0;i<n;i++){
		rowId=i;//*(p+i);
		mp=m+n*rowId+i; cp=c+rowId;
		maxElement=abs(*mp); maxRowId=rowId;
		for(int j=i+1;j<n;j++){
			mp+=n;//mp=m+n*(*(p+j))+i;
			if(abs(*mp)>maxElement){
				maxElement=abs(*mp); maxRowId=j;
			}
		}
		//*(p+i)=*(p+maxRowId); *(p+maxRowId)=rowId;
		cblas_sswap(n,m+n*rowId,1,m+n*maxRowId,1);
		tmp=*cp; *cp=*(c+maxRowId); *(c+maxRowId)=tmp;
		mp=m+n*rowId+i;
		tmp=1.0/(*mp); (*cp)*=tmp;
		for(int j=i;j<n;j++){
            (*mp)*=tmp; *mp++;
        }
		mp=m+n*rowId+i; cp=c+rowId; tmp=*cp;
		for(int j=i+1;j<n;j++){//BLAS-2
            cp++; mp+=n;
            *cp-=tmp*(*mp);//*(c+j(*(p+j)))-=(tmp);
			cblas_saxpy(n,-*mp,m+n*rowId,1,m+n*j/*(*(p+j))*/,1);
		}
	}
    //matrix_print(m,n); vector_print(c,n);
	float *xp=x+n-1,*xp2;
	*xp=*(c+n-1);
	for(int i=n-2;i>=0;i--){
        xp--; *xp=*(c+i);
        mp=m+n*(i+1)-1; xp2=x+n-1;
        for(int j=n-1;j>i;j--){
            *xp-=(*xp2)*(*mp);
            mp--; xp2--;
        }
	}
	return ;
}
int main() {
	int n,n2,n4;
	scanf("%d",&n); n2=n*n; n4=n2*n2;
    //n should be at least 3!
	float *m=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	vector_make(c,n);
    //matrix_print(m,n2); vector_print(c,n2);
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	gauss_colPivot(m,x,c,n2);
	e=steady_clock::now();
	//vector_print(x,n2);
	float *xp=x;
	for(int i=0;i<n2;i++){*xp-=1; xp++;}
	printf("Error: %.6f\n",cblas_snrm2(n2,x,1));
	printf("Time(ms): %.3f",(e-s).count()/1000000.0);
	free(m); free(c); free(x);
	return 0;
}
