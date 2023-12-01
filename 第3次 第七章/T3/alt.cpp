#include <cstdio>
#include <random>
#include <cblas.h>
#include <cmath>
#include <time.h>
#include <lapacke.h>
#define T_S 1000000000.0
using namespace std;
void matrix_build(double *m,int n){
	/*double a=1;
	for(int i=0;i<n;i++){
		a/=2; *(m+i*n+i)=a;
	}*/
	for(int i=0;i<25;i++)
		for(int j=0;j<15;j++) *(m+i*n+j)=pow((i+1.0)/25,j);
	//for(int i=0;i<80;i++) printf("%.2e\n",*(m+i*n+i)); printf("\n");
	
/*	srand(time(0));
	int N=8000,j,k; double s,c;
	for(int i=0;i<N;i++){
		j=rand()%n; k=rand()%(n-1);
		if(k>=j) k++;
		while(abs(0.5-(s=(double(rand())/RAND_MAX)))>0.45) ;
		c=sqrt(1-s*s);
		cblas_drot(n,m+j*n,1,m+k*n,1,c,s);
 
		j=rand()%n; k=rand()%(n-1);
		if(k>=j) k++;
		while(abs(0.5-(s=(double(rand())/RAND_MAX)))>0.45) ;
		c=sqrt(1-s*s);
		cblas_drot(n,m+j,n,m+k,n,c,s);
  }*/
	
	double rcond,anorm,*mm=(double *)malloc(6400*sizeof(double));
	cblas_dcopy(n*n,m,1,mm,1);
	double *w1=(double *)malloc(4*n*sizeof(double));
	int *w2=(int *)malloc(n*sizeof(int)),lda=n,info;
	
	anorm=LAPACK_dlange("1",&n,&n,mm,&lda,w1);
	LAPACK_dgetrf(&n,&n,mm,&lda,w2,&info);
	LAPACK_dgecon("1",&n,mm,&lda,&anorm,&rcond,w1,w2,&info);
	printf("Cond:%.2e\n\n",1.0/rcond);
	free(w1); free(w2);
	
	return ;
}
double Orth(double *q,int n){
	double *e=(double *)malloc(n*n*sizeof(double)),ans;
	for(int i=0;i<n*n;i++) *(e+i)=0;
	for(int i=0;i<n;i++) *(e+i*n+i)=1;
	cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,n,n,n,1,q,n,q,n,-1,e,n);
	ans=cblas_dnrm2(n*n,e,1);
	free(e); return ans;
}
double Error(double *q,double *r,CBLAS_TRANSPOSE t2,double *m,int n){
	double *a=(double *)malloc(n*n*sizeof(double)),ans;
	cblas_dcopy(n*n,m,1,a,1);
	ans=cblas_dnrm2(n*n,a,1);
	cblas_dgemm(CblasRowMajor,CblasTrans,t2,n,n,n,1,q,n,r,n,-1,a,n);
	ans=cblas_dnrm2(n*n,a,1)/ans;
	free(a); return ans;
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
	//for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,q+k,n);
	/*for(int k=0;k<n;k++)
		for(int j=0;j<n;j++) *(q+n*k+j)=*(m+n*j+k);

	double tmp;
	for(int k=0;k<n;k++){
		for(int j=0;j<k;j++){
			//*(r+j*n+k)=cblas_ddot(n,q+j*n,1,q+k*n,1);
			tmp=0;
			for(int t=0;t<n;t++) tmp+=*(q+k*n+t)*(*(q+j*n+t));
			*(r+j*n+k)=tmp;
			//cblas_daxpy(n,-*(r+j*n+k),q+j*n,1,q+k*n,1);
			for(int t=0;t<n;t++)
				*(q+k*n+t)-=*(r+j*n+k)*(*(q+j*n+t));
		}	
		//*(r+k*n+k)=cblas_dnrm2(n,q+k*n,1);
		tmp=0;
		for(int j=0;j<n;j++) tmp+=*(q+k*n+j)*(*(q+k*n+j));
		*(r+k*n+k)=sqrt(tmp);
		for(int i=0;i<n;i++) *(q+k*n+i)/=*(r+k*n+k);
	}*/
	
	//for(int i=0;i<n;i++) printf("%.2e\n",*(r+i*n+i)); printf("\n");
	printf("%.4e %.4e\n",Orth(q,n),Error(q,r,CblasNoTrans,m,n));
	free(q); free(r); return ;
}
void MGS(double *m,int n){
	double *q=(double *)malloc(n*n*sizeof(double));
	double *r=(double *)malloc(n*n*sizeof(double));
	//for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,q+k,n);
	for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,q+k,n);
	for(int k=0;k<n;k++){
		*(r+k*n+k)=cblas_dnrm2(n,q+k*n,1);
		for(int i=0;i<n;i++) *(q+k*n+i)/=*(r+k*n+k);
		for(int j=k+1;j<n;j++){
			*(r+k*n+j)=cblas_ddot(n,q+k*n,1,q+j*n,1);
			cblas_daxpy(n,-*(r+k*n+j),q+k*n,1,q+j*n,1);
		}
	}
	/*for(int k=0;k<n;k++)
		for(int j=0;j<n;j++) *(q+n*k+j)=*(m+n*j+k);
	
	double tmp;
	for(int k=0;k<n;k++){
		//*(r+k*n+k)=cblas_dnrm2(n,q+k*n,1);
		tmp=0;
		for(int j=0;j<n;j++) tmp+=*(q+k*n+j)*(*(q+k*n+j));
		*(r+k*n+k)=sqrt(tmp);
		for(int i=0;i<n;i++) *(q+k*n+i)/=*(r+k*n+k);
		for(int j=k+1;j<n;j++){
			//*(r+k*n+j)=cblas_ddot(n,q+k*n,1,q+j*n,1);
			tmp=0;
			for(int t=0;t<n;t++) tmp+=*(q+k*n+t)*(*(q+j*n+t));
			*(r+k*n+j)=tmp;
			//cblas_daxpy(n,-*(r+k*n+j),q+k*n,1,q+j*n,1);
			for(int t=0;t<n;t++)
				*(q+j*n+t)-=*(r+k*n+j)*(*(q+k*n+t));
		}
	}*/
	
	//for(int i=0;i<n;i++) printf("%.2e\n",*(r+i*n+i)); printf("\n");
	printf("%.4e %.4e\n",Orth(q,n),Error(q,r,CblasNoTrans,m,n));
	free(q); free(r); return ;
}
int main(){
	int n=25;
	//double *m=(double *)malloc(n*n*sizeof(double));
	double *m=(double *)malloc(25*15*sizeof(double));
	matrix_build(m,n);
	CGS(m,n);
	MGS(m,n);
	free(m); 
	return 0;
}
