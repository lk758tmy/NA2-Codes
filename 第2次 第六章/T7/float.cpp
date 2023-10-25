#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <cblas.h>
using namespace std;
void matrix_make(float *h,int n){
	int n2=n*n; float *p=h;
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
void precond_make(float *h,int n,float omega){
	int n2=n*n; float *p=h,*q,tmp1,tmp2,tmp3,omega2=omega*omega;
	for(int i=0;i<n2*n2;i++){*p=0; p++;}
	p=h;
	tmp1=16*(1-2*omega+omega2); tmp2=tmp1+omega2; tmp3=4*(omega-omega2);
	for(int i=0;i<n2;i++){
		if(i%n==0) *p=tmp1;
		else{ *p=tmp2; *(p-1)=*(p-n2)=tmp3;}
		p+=(n2+1);
	}
	p=h+n*n2+n; 
	for(int i=n;i<n2;i++){
		*p+=omega2; p+=(n2+1);
	}
	p=h+n*n2; q=h+n;
	for(int i=n,j=0;i<n2;i++,j++){
		*p=*q=-tmp3;
		if(i%n==0) continue;
		*(p-1)=*(q-n2)=-omega2;
		p+=(n2+1); q+=(n2+1);
	}
	return ;
}
int CG(float *m,float *x,float *c,int n,float *r,float epsilon){
	float alpha,beta,*p=(float *)malloc(n*sizeof(float));
	float tmp_nrm2r2,*tmp_Amp=(float *)malloc(n*sizeof(float));
	int cnt;
	
	for(int i=0;i<n;i++) *(r+i)=*(c+i);
	cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,x,1,-1,r,1);
	for(int i=0;i<n;i++) *(p+i)=-*(r+i);
	cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
	tmp_nrm2r2=cblas_sdot(n,r,1,r,1);
	alpha=tmp_nrm2r2/cblas_sdot(n,p,1,tmp_Amp,1);
	
	for(cnt=1;cnt<=n;cnt++){
		cblas_saxpy(n,alpha,p,1,x,1);
		cblas_saxpy(n,alpha,tmp_Amp,1,r,1);
		beta=cblas_sdot(n,r,1,r,1)/tmp_nrm2r2;
		tmp_nrm2r2*=beta;
		if(tmp_nrm2r2<epsilon) goto CG_END;
		cblas_saxpy(n,-1/beta,r,1,p,1);
		for(int j=0;j<n;j++) *(p+j)*=beta;
		cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
		alpha=tmp_nrm2r2/cblas_sdot(n,p,1,tmp_Amp,1);
	}
	
	CG_END: free(p); free(tmp_Amp); return cnt;
}
int main(){
    int n,n2,n4; float omega=1;
	scanf("%d",&n);//n should be at least 3!
	n2=n*n; n4=n2*n2;
	float *m=(float *)malloc(n4*sizeof(float));
	float *q=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	precond_make(q,n,omega);
	
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	float *r=(float *)malloc(n2*sizeof(float));
	for(int i=0;i<n;i++) *(x+i)=0;
	vector_make(c,n);
	
	printf("%d\t%d\n",n,CG(m,x,c,n2,r,1e-12));

	free(m); free(c); free(x); free(r); free(q);
    return 0;
}
