#include <cstdio>
#include <random>
#include <cblas.h>
#include <cmath>
//#include <lapacke.h>
#include <chrono>
#define T_S 1000000000.0
using namespace std;
using namespace std::chrono;
void rand_matrix(double *m,int n){
	for(int i=0;i<n*n;i++) *(m+i)=0;
	for(int i=0;i<n*n;i+=(n+1)) *(m+i)=100;
	double alpha; int nmax=sqrt(n);
	double **cards=(double **)malloc(n*sizeof(double *)),**tmpp,*tmpi;
	for(int j=0;j<n;j++) *(cards+j)=m+n*j;
	for(int i=0;i<nmax;i++){
		for(int j=n-1;j>0;j--){
			tmpp=cards+rand()%j;
			tmpi=*(cards+j); *(cards+j)=*tmpp; *tmpp=tmpi;
		}
		for(int j=0;j<n-1;j+=2){
			alpha=((rand()%2==1)?1:-1);
			alpha*=(rand()%80+21)/100.0;
			cblas_daxpy(n,alpha,*(cards+j),1,*(cards+j+1),1);			
		}
	}
	free(cards); return ;
}
double Orth(double *q,int n){
	double *e=(double *)malloc(n*n*sizeof(double)),ans;
	//double *work=(double *)malloc(n*sizeof(double));
	for(int i=0;i<n*n;i++) *(e+i)=0;
	for(int i=0;i<n;i++) *(e+i*n+i)=1;
	cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,n,n,n,1,q,n,q,n,-1,e,n);
	ans=cblas_dnrm2(n*n,e,1);
	//ans=LAPACK_dlange("M",&n,&n,e,&n,work);
	free(e); /*free(work);*/ return ans;
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
	time_point<steady_clock> s,e;
	for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,q+k,n);
	s=steady_clock::now();

	for(int k=0;k<n;k++){
		for(int j=0;j<k;j++){
			*(r+j*n+k)=cblas_ddot(n,q+j*n,1,q+k*n,1);
			cblas_daxpy(n,-*(r+j*n+k),q+j*n,1,q+k*n,1);
		}	
		*(r+k*n+k)=cblas_dnrm2(n,q+k*n,1);
		for(int i=0;i<n;i++) *(q+k*n+i)/=*(r+k*n+k);
	}
	
	e=steady_clock::now();
	printf(",%.4e,%.3f,%.4e",Orth(q,n),(e-s).count()/T_S,Error(q,r,CblasNoTrans,m,n));
	free(q); free(r); return ;
}
void MGS(double *m,int n){
	double *q=(double *)malloc(n*n*sizeof(double));
	double *r=(double *)malloc(n*n*sizeof(double));
	for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,q+k,n);
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
	
	e=steady_clock::now();
	printf(",%.4e,%.3f,%.4e",Orth(q,n),(e-s).count()/T_S,Error(q,r,CblasNoTrans,m,n));
	free(q); free(r); return ;
}
void Householder(double *m,int n){
	double *q=(double *)malloc(n*n*sizeof(double));
	double *r=(double *)malloc(n*n*sizeof(double));
	double *u=(double *)malloc(n*sizeof(double)),*y=(double *)malloc(n*sizeof(double));
	double a,b; int tmp;
	for(int k=0;k<n;k++) cblas_dcopy(n,m+n*k,1,r+k,n);
	for(int k=0;k<n;k++){
		*(q+k*n+k)=1;
		for(int j=k+1;j<n;j++) *(q+k*n+j)=*(q+j*n+k)=0;
	}
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	
	for(int i=0;i<n;i++){
		tmp=i*n+i;
		a=cblas_dnrm2(n-i,r+tmp,1);
		if(*(r+tmp)>0) a=-a;
		b=a*(a-*(r+tmp));
		cblas_dcopy(n-i,r+tmp,1,u+i,1);
		*(u+i)-=a;
		cblas_dgemv(CblasRowMajor,CblasTrans,n-i,n,1,q+i*n,n,u+i,1,0,y,1);
		cblas_dger(CblasRowMajor,n-i,n,-1/b,u+i,1,y,1,q+i*n,n);
		cblas_dgemv(CblasColMajor,CblasTrans,n-i,n-i,1,r+tmp,n,u+i,1,0,y+i,1);
		cblas_dger(CblasColMajor,n-i,n-i,-1/b,u+i,1,y+i,1,r+tmp,n);
	}
	
	e=steady_clock::now();
	printf(",%.4e,%.3f,%.4e",Orth(q,n),(e-s).count()/T_S,Error(q,r,CblasTrans,m,n));
	free(q); free(r); free(u); free(y); return ;
}
void Givens(double *m,int n){
	double *q=(double *)malloc(n*n*sizeof(double));
	double *r=(double *)malloc(n*n*sizeof(double)),_a,_b,_c,_s;
	cblas_dcopy(n*n,m,1,r,1);
	for(int k=0;k<n;k++){
		*(q+k*n+k)=1;
		for(int j=k+1;j<n;j++) *(q+k*n+j)=*(q+j*n+k)=0;
	}
	time_point<steady_clock> s,e;
	s=steady_clock::now();
	
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			_b=*(r+j*n+i); if(abs(_b)<1e-9) continue;
			_a=*(r+i*n+i);
			cblas_drotg(&_a,&_b,&_c,&_s);
			cblas_drot(n,q+i*n,1,q+j*n,1,_c,_s);
			cblas_drot(n-i,r+i*n+i,1,r+j*n+i,1,_c,_s);
		}
	}
	
	e=steady_clock::now();
	printf(",%.4e,%.3f,%.4e",Orth(q,n),(e-s).count()/T_S,Error(q,r,CblasNoTrans,m,n));
	free(q); free(r); return ;
}
int main(){
	int n; double *m=0;
	freopen("\data.csv","w",stdout);
	printf("No.,Order,CGS,,,MGS,,,H,,,G\n");
	printf(",,Orth.,Time(s),R.Err,Orth.,Time(s),R.Err,Orth.,Time(s),R.Err,Orth.,Time(s),R.Err\n");
	for(int i=1;i<=100;i++){
		srand(system_clock::to_time_t(system_clock::now()));
		n=rand()%1501+500;
		m=(double *)malloc(n*n*sizeof(double));
		rand_matrix(m,n);
		printf("%d,%d",i,n);
		CGS(m,n); MGS(m,n); Householder(m,n); Givens(m,n);
		free(m); 
		printf("\n");
	}
	return 0;
}
