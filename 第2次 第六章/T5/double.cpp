#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <cblas.h>
using namespace std;
double niu,kse,rou;
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
void Jacobi_Chebyshev(double *m,double *x,double *c,int n,int step){
	double *y1=(double *)malloc(n*sizeof(double));
	double *y2=(double *)malloc(n*sizeof(double));
	double *r=(double *)malloc(n*sizeof(double));
	double *e=(double *)malloc(n*sizeof(double));
    for(int i=0;i<n;i++){
    	*(y1+i)=0;
        *(y2+i)=*(c+i)-cblas_ddot(n,x,1,m+i*n,1);
        *(y2+i)+=*(x+i)*(*(m+i*n+i));
		*(y2+i)/=*(m+i*n+i);
    }
	rou=2; double tmp=1.0/4/kse/kse;
	
	printf("Step\tResidual\tError\n");
	for(int j=2;j<=step;j++){
		rou=1.0/(1-rou*tmp);
		for(int i=0;i<n;i++){
			*(x+i)=*(c+i)-cblas_ddot(n,y2,1,m+i*n,1);
        	*(x+i)+=*(y2+i)*(*(m+i*n+i));
			*(x+i)*=(rou*niu/(*(m+i*n+i)));
		}
		cblas_daxpy(n,rou*(1-niu),y2,1,x,1);
		cblas_daxpy(n,1-rou,y1,1,x,1);
		cblas_dcopy(n,y2,1,y1,1); cblas_dcopy(n,x,1,y2,1);
		
		cblas_dcopy(n,c,1,r,1);
		cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,x,1,-1,r,1);
		printf("%d\t%.6f",j,cblas_dnrm2(n,r,1));
		for(int i=0;i<n;i++) *(e+i)=*(x+i)-1;
		printf("\t%.6f\n",cblas_dnrm2(n,e,1));
	}
	free(y1); free(y2); free(r); free(e);	
	return ;
}
int main(){
    int n,n2,n4,step;
	scanf("%d %d",&n,&step);//n should be at least 3!
	n2=n*n; n4=n2*n2; niu=1.0/3;
	kse=-0.75/sin((n-1.0)*M_PI/(2*n+2));
	double *m=(double *)malloc(n4*sizeof(double));
	matrix_make(m,n);
	double *c=(double *)malloc(n2*sizeof(double));
	double *x=(double *)malloc(n2*sizeof(double));
	for(int i=0;i<n;i++) *(x+i)=0;
	vector_make(c,n);
	
	Jacobi_Chebyshev(m,x,c,n2,step);

	free(m); free(c); free(x);
    return 0;
}
