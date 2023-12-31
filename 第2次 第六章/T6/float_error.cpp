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
void CG(float *m,float *x,float *c,int n,float *r,float epsilon){
	float alpha,beta,*p=(float *)malloc(n*sizeof(float));
	float tmp_nrm2r,tmp_nrm2r2,*tmp_Amp=(float *)malloc(n*sizeof(float));
	float *e=(float *)malloc(n*sizeof(float)),error;
	
	for(int i=0;i<n;i++) *(r+i)=*(c+i);
	cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,x,1,-1,r,1);
	for(int i=0;i<n;i++) *(p+i)=-*(r+i);
	cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
	tmp_nrm2r=cblas_snrm2(n,r,1); tmp_nrm2r2=tmp_nrm2r*tmp_nrm2r;
	alpha=tmp_nrm2r2/cblas_sdot(n,p,1,tmp_Amp,1);
	
	printf("Step\tError\tResidual\n");
	for(int i=1;i<=n;i++){
		cblas_saxpy(n,alpha,p,1,x,1);
		cblas_saxpy(n,alpha,tmp_Amp,1,r,1);
		tmp_nrm2r=cblas_snrm2(n,r,1);
		beta=tmp_nrm2r*tmp_nrm2r/tmp_nrm2r2;
		tmp_nrm2r2*=beta;
		for(int j=0;j<n;j++) *(e+j)=*(x+j)-1;
		error=cblas_snrm2(n,e,1);
		printf("%d\t%.6f\t%.6f\n",i,error,tmp_nrm2r);
		if(error<epsilon) goto CG_END;
		cblas_saxpy(n,-1/beta,r,1,p,1);
		for(int j=0;j<n;j++) *(p+j)*=beta;
		cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
		alpha=tmp_nrm2r2/cblas_sdot(n,p,1,tmp_Amp,1);
	}
	
	CG_END: free(p); free(e); free(tmp_Amp); return ;
}
int main(){
    int n,n2,n4;
	scanf("%d",&n);//n should be at least 3!
	n2=n*n; n4=n2*n2;
	float *m=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	float *r=(float *)malloc(n2*sizeof(float));
	for(int i=0;i<n;i++) *(x+i)=0;
	vector_make(c,n);
	
	CG(m,x,c,n2,r,3*1e-5);

	free(m); free(c); free(x); free(r);
    return 0;
}
