#include <cstdio>
#include <random>
#include <cblas.h>
#include <cmath>
#include <time.h>
#define T_S 1000000000.0
using namespace std;
void matrix_build(double *m,int n){
	double a=1;
	for(int i=0;i<n;i++){
		a/=2; *(m+i*n+i)=a;
	}
	for(int i=0;i<80;i++) printf("%.2e\n",*(m+i*n+i)); printf("\n");
	srand(time(0));
	int N=800,j,k; double s,c;
	for(int i=0;i<N;i++){
		j=rand()%n; k=rand()%(n-1);
		if(k>=j) k++;
		while(abs(0.5-(s=(double(rand())/RAND_MAX)))>0.4) ;
		c=sqrt(1-s*s);
		cblas_drot(n,m+j*n,1,m+k*n,1,c,s);
 
		j=rand()%n; k=rand()%(n-1);
		if(k>=j) k++;
		while(abs(0.5-(s=(double(rand())/RAND_MAX)))>0.4) ;
		c=sqrt(1-s*s);
		cblas_drot(n,m+j,n,m+k,n,c,s);
	}
	return ;
}
void CGS(double *m,int n){
	double *q=(double *)malloc(n*n*sizeof(double));
	double *r=(double *)malloc(n*n*sizeof(double));
	for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,q+k,n);

	for(int k=0;k<n;k++){
		for(int j=0;j<k;j++){
			*(r+j*n+k)=cblas_ddot(n,q+j*n,1,q+k*n,1);
			cblas_daxpy(n,-*(r+j*n+k),q+j*n,1,q+k*n,1);
		}	
		*(r+k*n+k)=cblas_dnrm2(n,q+k*n,1);
		for(int i=0;i<n;i++) *(q+k*n+i)/=*(r+k*n+k);
	}
	
	for(int i=0;i<80;i++) printf("%.2e\n",*(r+i*n+i)); printf("\n");
	
	free(q); free(r); return ;
}
void MGS(double *m,int n){
	double *q=(double *)malloc(n*n*sizeof(double));
	double *r=(double *)malloc(n*n*sizeof(double));
	for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,q+k,n);
	
	for(int k=0;k<n;k++){
		*(r+k*n+k)=cblas_dnrm2(n,q+k*n,1);
		for(int i=0;i<n;i++) *(q+k*n+i)/=*(r+k*n+k);
		for(int j=k+1;j<n;j++){
			*(r+k*n+j)=cblas_ddot(n,q+k*n,1,q+j*n,1);
			cblas_daxpy(n,-*(r+k*n+j),q+k*n,1,q+j*n,1);
		}
	}
	
	for(int i=0;i<80;i++) printf("%.2e\n",*(r+i*n+i)); printf("\n");
	
	free(q); free(r); return ;
}
int main(){
	int n=80; double *m=(double *)malloc(6400*sizeof(double));
	matrix_build(m,n);
	CGS(m,n); MGS(m,n);
	free(m); 
	return 0;
}
