#include <cstdio>
#include <random>
#include <cblas.h>
#include <cmath>
//#include <lapacke.h>
#include <chrono>
#define T_S 1000000000.0
using namespace std;
using namespace std::chrono;
void matrix_vector_build(double *m,int n){
	double *b=m+n-1;//向量b存储在m的最后一列里
	for(int i=0;i<n*(n-1);i++) *(m+i)=0;
	for(int i=0;i<n-1;i++){
		*(m+i*n+i)=2; *(m+i*n+i+1)=*(m+i*n+i+n)=-1;
		*(b+i*n)=(i+1)/double(n);
	} //b把m多出来的一项覆盖掉了
	*b+=1; *(b+n*(n-2))+=1; *(b+n*(n-1))=0;
	
	srand(system_clock::to_time_t(system_clock::now()));
	int N=rand()%1001+1000,j,k; double s,c;
	for(int i=0;i<N;i++){
		j=rand()%n; k=rand()%(n-1);
		if(k>=j) k++;
		while(abs(0.5-(s=(double(rand())/RAND_MAX)))>0.4) ;
		c=sqrt(1-s*s);
		cblas_drot(n,m+j*n,1,m+k*n,1,c,s);
	}
	return ;
}
double Error(double *x,int n){
	for(int i=0;i<n;i++) *(x+i)-=1;
	return cblas_dnrm2(n,x,1);
}
void NormEqCG(double *m,int n){
	double *c=(double *)malloc(n*sizeof(double));	
	double *x=(double *)malloc(n*sizeof(double));
	for(int i=0;i<n;i++) *(x+i)=0;
	time_point<steady_clock> s,e;
	s=steady_clock::now();

	cblas_dgemv(CblasRowMajor,CblasTrans,n,n-1,1,m,n,m+n-1,n,0,c,1);
	
	double alpha,beta,*p=(double *)malloc(n*sizeof(double));
	double tmp_nrm2r2,*tmp_Amp=(double *)malloc(n*sizeof(double));
	double *tmp_Amp_2=(double *)malloc(n*sizeof(double));
	double *r=(double *)malloc(n*sizeof(double));
	for(int i=0;i<n-1;i++) *(r+i)=-*(c+i);
	for(int i=0;i<n-1;i++) *(p+i)=-*(r+i);

	cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n-1,1,m,n,p,1,0,tmp_Amp_2,1);
	cblas_dgemv(CblasRowMajor,CblasTrans,n,n-1,1,m,n,tmp_Amp_2,1,0,tmp_Amp,1);
	tmp_nrm2r2=cblas_ddot(n-1,r,1,r,1);
	alpha=tmp_nrm2r2/cblas_ddot(n-1,p,1,tmp_Amp,1);
	
	for(int cnt=1;cnt<n;cnt++){
		cblas_daxpy(n-1,alpha,p,1,x,1);
		cblas_daxpy(n-1,alpha,tmp_Amp,1,r,1);
		beta=cblas_ddot(n-1,r,1,r,1)/tmp_nrm2r2;
		tmp_nrm2r2*=beta;
		//if(tmp_nrm2r2<epsilon) goto CG_END;
		cblas_daxpy(n-1,-1/beta,r,1,p,1);
		for(int j=0;j<n-1;j++) *(p+j)*=beta;
		cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n-1,1,m,n,p,1,0,tmp_Amp_2,1);
		cblas_dgemv(CblasRowMajor,CblasTrans,n,n-1,1,m,n,tmp_Amp_2,1,0,tmp_Amp,1);
		alpha=tmp_nrm2r2/cblas_ddot(n-1,p,1,tmp_Amp,1);
	}
	
	e=steady_clock::now();
	printf(",%.3f,%.4e",(e-s).count()/T_S,Error(x,n-1));
	free(r); free(p); free(tmp_Amp); free(x); free(c); free(tmp_Amp_2); return ;
}
void AugNormEq(double *m,int n){
	double *a=(double *)malloc((2*n-1)*(2*n-1)*sizeof(double));	
	double *x=(double *)malloc((2*n-1)*sizeof(double));
	double *c=(double *)malloc((2*n-1)*sizeof(double));
	for(int i=0;i<n;i++){
		*(a+i*(2*n-1)+i)=1; *(c+i)=*(m+n-1+i*n);	
		cblas_dcopy(n-1,m+i*n,1,a+n*(2*n-1)+i,2*n-1);
		cblas_dcopy(n-1,m+i*n,1,a+n+i*(2*n-1),1);
	}
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	
	double alpha,beta,*p=(double *)malloc((2*n-1)*sizeof(double));
	double tmp_nrm2r2,*tmp_Amp=(double *)malloc((2*n-1)*sizeof(double));
	double *r=(double *)malloc((2*n-1)*sizeof(double));
	for(int i=0;i<2*n-1;i++) *(r+i)=-*(c+i);
	for(int i=0;i<n*2-1;i++) *(p+i)=-*(r+i);
	
	cblas_dgemv(CblasRowMajor,CblasNoTrans,2*n-1,2*n-1,1,a,2*n-1,p,1,0,tmp_Amp,1);
	tmp_nrm2r2=cblas_ddot(n*2-1,r,1,r,1);
	alpha=tmp_nrm2r2/cblas_ddot(n*2-1,p,1,tmp_Amp,1);
	
	for(int cnt=1;cnt<2*n;cnt++){
		cblas_daxpy(n*2-1,alpha,p,1,x,1);
		cblas_daxpy(n*2-1,alpha,tmp_Amp,1,r,1);
		beta=cblas_ddot(n*2-1,r,1,r,1)/tmp_nrm2r2;
		tmp_nrm2r2*=beta;
		//if(tmp_nrm2r2<epsilon) goto CG_END;
		cblas_daxpy(n*2-1,-1/beta,r,1,p,1);
		for(int j=0;j<2*n-1;j++) *(p+j)*=beta;
		cblas_dgemv(CblasRowMajor,CblasNoTrans,2*n-1,2*n-1,1,a,2*n-1,p,1,0,tmp_Amp,1);
		alpha=tmp_nrm2r2/cblas_ddot(2*n-1,p,1,tmp_Amp,1);
	}
	
	e=steady_clock::now();
	printf(",%.3f,%.4e",(e-s).count()/T_S,Error(x+n,n-1));
	free(r); free(p); free(tmp_Amp); free(x); free(c); free(a); return ;
}
void MGS(double *m,int n){
	double *q=(double *)malloc(n*n*sizeof(double));
	double *r=(double *)malloc(n*n*sizeof(double));;
	double *x=(double *)malloc((n-1)*sizeof(double));
	for(int k=0;k<n;k++) cblas_dcopy(n,m+k,n,q+n*k,1);
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	
	for(int k=0;k<n;k++){
		*(r+k*n+k)=cblas_dnrm2(n,q+k*n,1);
		for(int i=0;i<n;i++) *(q+k*n+i)/=*(r+k*n+k);
		for(int j=k+1;j<n;j++){
			*(r+k*n+j)=cblas_ddot(n,q+k*n,1,q+j*n,1);
			cblas_daxpy(n,-*(r+k*n+j),q+k*n,1,q+j*n,1);
		}
	}
	for(int k=n-2;k>-1;k--){
		*(x+k)=*(r+(k+1)*n-1);
		for(int j=n-2;j>k;j--) *(x+k)-=*(r+k*n+j)*(*(x+j));
		*(x+k)/=*(r+k*n+k);
	}
	
	e=steady_clock::now();
	printf(",%.3f,%.4e",(e-s).count()/T_S,Error(x,n-1));
	free(q); free(r); free(x); return ;
}
void Householder(double *m,int n){
	double *r=(double *)malloc(n*n*sizeof(double));
	double *x=(double *)malloc((n-1)*sizeof(double));
	double *u=(double *)malloc(n*sizeof(double)),*y=(double *)malloc(n*sizeof(double));
	double a,b; int tmp;
	for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,r+k,n);
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	
	for(int i=0;i<n-1;i++){
		tmp=i*n+i;
		a=cblas_dnrm2(n-i,r+tmp,1);
		if(*(r+tmp)>0) a=-a;
		b=a*(a-*(r+tmp));
		cblas_dcopy(n-i,r+tmp,1,u+i,1);
		*(u+i)-=a;
		cblas_dgemv(CblasColMajor,CblasTrans,n-i,n-i,1,r+tmp,n,u+i,1,0,y+i,1);
		cblas_dger(CblasColMajor,n-i,n-i,-1/b,u+i,1,y+i,1,r+tmp,n);
	}
	for(int k=n-2;k>-1;k--){
		*(x+k)=*(r+(n-1)*n+k);
		for(int j=n-2;j>k;j--) *(x+k)-=*(r+j*n+k)*(*(x+j));
		*(x+k)/=*(r+k*n+k);
	}
	
	e=steady_clock::now();
	printf(",%.3f,%.4e",(e-s).count()/T_S,Error(x,n-1));
	free(x); free(r); free(u); free(y); return ;
}
void Givens(double *m,int n){
	double *r=(double *)malloc(n*n*sizeof(double)),_a,_b,_c,_s;
	double *x=(double *)malloc((n-1)*sizeof(double));
	cblas_dcopy(n*n,m,1,r,1);
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			_b=*(r+j*n+i); //if(abs(_b)<1e-12) continue;
			_a=*(r+i*n+i);
			cblas_drotg(&_a,&_b,&_c,&_s);
			cblas_drot(n-i,r+i*n+i,1,r+j*n+i,1,_c,_s);
		}
	}
	for(int k=n-2;k>-1;k--){
		*(x+k)=*(r+(k+1)*n-1);
		for(int j=n-2;j>k;j--) *(x+k)-=*(r+k*n+j)*(*(x+j));
		*(x+k)/=*(r+k*n+k);
	}
	
	e=steady_clock::now();
	printf(",%.3f,%.4e",(e-s).count()/T_S,Error(x,n-1));
	free(x); free(r); return ;
}
int main(){
	double *m;
	for(int n=500;n<2501;n+=500){
		m=(double *)malloc(n*n*sizeof(double));
		matrix_vector_build(m,n);
		printf("%d",n);
		NormEqCG(m,n); AugNormEq(m,n); MGS(m,n); Householder(m,n); Givens(m,n);
		free(m); 
		printf("\n");
	}
	return 0;
}
