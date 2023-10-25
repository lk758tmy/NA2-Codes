#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <cblas.h>
using namespace std;
void matrix_make(double *h,int n){
	int n2=n*n; double *p=h;
	for(int i=0;i<n2*n2;i++){*p=0; p++;}
    p=h;
    *h=4;
	for(int i=1;i<n2;i++){
		p+=(n2+1); *p=4;
		if(i%n!=0) *(p-1)=*(p-n2)=-1;
	}
	for(int i=n,j=0;i<n2;i++,j++)
		*(h+i*n2+j)=*(h+j*n2+i)=-1;
	return ;
}
void vector_make(double *h,int n){
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
int CG(double *m,double *x,double *c,int n,double *r,double epsilon){
	double alpha,beta,*p=(double *)malloc(n*sizeof(double));
	double tmp_nrm2r2,*tmp_Amp=(double *)malloc(n*sizeof(double));
	int cnt=0;
	
	for(int i=0;i<n;i++) *(r+i)=*(c+i);
	cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,x,1,-1,r,1);
	for(int i=0;i<n;i++) *(p+i)=-*(r+i);
	cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
	tmp_nrm2r2=cblas_ddot(n,r,1,r,1);
	alpha=tmp_nrm2r2/cblas_ddot(n,p,1,tmp_Amp,1);
	
	for(cnt=1;cnt<=n;cnt++){
		cblas_daxpy(n,alpha,p,1,x,1);
		cblas_daxpy(n,alpha,tmp_Amp,1,r,1);
		beta=cblas_ddot(n,r,1,r,1)/tmp_nrm2r2;
		tmp_nrm2r2*=beta;
		if(tmp_nrm2r2<epsilon) goto CG_END;
		//使用残差2范数的平方做收敛判定，少一些计算量
		cblas_daxpy(n,-1/beta,r,1,p,1);
		for(int j=0;j<n;j++) *(p+j)*=beta;
		cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
		alpha=tmp_nrm2r2/cblas_ddot(n,p,1,tmp_Amp,1);
	}
	
	CG_END: free(p); free(tmp_Amp); return cnt;
}
int _main(int n){
    int n2=n*n, n4=n2*n2;
	double *m=(double *)malloc(n4*sizeof(double));
	matrix_make(m,n);
	double *c=(double *)malloc(n2*sizeof(double));
	double *x=(double *)malloc(n2*sizeof(double));
	double *r=(double *)malloc(n2*sizeof(double));
	for(int i=0;i<n;i++) *(x+i)=0;
	vector_make(c,n);
	printf("%d\t%d\n",n,CG(m,x,c,n2,r,1e-12));
	free(m); free(c); free(x); free(r);
    return 0;
}
int main(){
	printf("Size\tStep\n");
	for(int i=3;i<=60;i+=3) _main(i);
	return 0;
}
