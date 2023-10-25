#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <cblas.h>
using namespace std;
int cnt=0; float L1,L2;
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
void Richardson(float *m,float *x,float *c,int n,int step,float *r){
    float tau,residual=0,*e=(float *)malloc(n*sizeof(float));
    for(int k=1;k<=step;k++){
        cnt++; tau=1/(cos((2*k-1.0)*M_PI/2/step)*L1+L2);
        if(residual<1e-8){
            cblas_scopy(n,c,1,r,1);
            cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,-1,m,n,x,1,1,r,1);
        }
        cblas_saxpy(n,tau,r,1,x,1);
        cblas_scopy(n,c,1,r,1);
        cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,-1,m,n,x,1,1,r,1);
        residual=cblas_snrm2(n,r,1);

        for(int i=0;i<n;i++) *(e+i)=*(x+i)-1;
        printf("%d\t%.8f\t%.8f\n",cnt,residual,cblas_snrm2(n,e,1));
    }
    free(e); return ;
}
int main(){
    int n,n2,n4,step,t;
	scanf("%d %d %d",&n,&step,&t); //n should be at least 3!
	n2=n*n; n4=n2*n2; L1=4*sin((n-1.0)*M_PI/2/(n+1)); L2=4;
	float *m=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	float *r=(float *)malloc(n2*sizeof(float));
	vector_make(c,n);
	for(int i=0;i<n2;i++) *(x+i)=*(r+i)=0;

	printf("Steps\tResidual\tError\n");
	for(int i=0;i<t;i++) Richardson(m,x,c,n2,step,r);
	free(m); free(c); free(x); free(r);
    return 0;
}
