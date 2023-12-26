#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <random>
#include <time.h>
#include <lapacke.h>
using namespace std;
const double epsilon=1e-6;
void func(int n,double *x,double *y){
	double lambda=*(x+n);
	*(y+n)=cblas_ddot(n,x,1,x,1)-1;
	*(y+n-1)=-lambda*(*(x+n-1))+2*(*(x+n-1))-*(x+n-2);
	*y=-lambda*(*x)+2*(*x)-*(x+1);
	for(int i=1;i<n-1;i++){
		*(y+i)=-lambda*(*(x+i));
		*(y+i)+=*(x+i)*2-*(x+i-1)-*(x+i+1);
	}
	return ;
}
void func_d(int n,double *x,double *m){
	double lambda=*(x+n);
	for(int i=0;i<(n+1)*(n+1);i++) *(m+i)=0;
	for(int i=0;i<n;i++){
		*(m+i*(n+1)+i)=2-lambda;
		*(m+i*(n+1)+i+1)=*(m+i*(n+1)+i+n+1)=-1;
		*(m+n*(n+1)+i)=2*(*(x+i));
		*(m+n+i*(n+1))=-(*(x+i));
	}
}
void matrix_make(double *m,int n){
	for(int i=0;i<n-1;i++){
		*(m+i*n+i)=2; *(m+i*n+i+1)=*(m+i*n+i+n)=-1;
		for(int j=2+i;j<n;j++)
			*(m+i*n+j)=*(m+j*n+i)=0;
	}
	*(m+n*n-1)=2;
	return ;
}
int newton(double *x,int n,void(*f)(int,double*,double*),void(*d)(int,double*,double*)){
	int cnt=0;
	double *y=(double *)malloc((n+1)*sizeof(double));
	int *ipiv=(int *)malloc((n+1)*sizeof(int));
	double *m=(double *)malloc((n+1)*(n+1)*sizeof(double)),tmp_f=1;
	printf("%d,,%.10f\n",cnt,*(x+n));
	while(tmp_f>epsilon){
		cnt++;
		f(n,x,y); d(n,x,m); tmp_f=cblas_dnrm2(n+1,y,1);
		//z=x-M^(-1)y => Mz=Mx-y
		cblas_dgemv(CblasRowMajor,CblasNoTrans,n+1,n+1,1,m,n+1,x,1,-1,y,1);
		LAPACKE_dgesv(LAPACK_ROW_MAJOR,n+1,1,m,n+1,ipiv,y,1);
		for(int i=0;i<=n;i++) *(x+ipiv[i]-1)=*(y+i);
		printf("%d,%.10f,%.10f\n",cnt,tmp_f,*(x+n));
		if(cnt==20) break;
	}
	free(y); free(m); free(ipiv); return cnt; 
}
int main(){
	double *x,*tmp_v,tmp,*t;
	srand(time(0));
	for(int n=5;n<9;n+=3){
		printf("%d\n",n);
		x=(double *)malloc((n+1)*sizeof(double));
		tmp_v=(double *)malloc(n*sizeof(double));
		t=(double *)malloc(n*n*sizeof(double));
		matrix_make(t,n);
		for(int i=0;i<n;i++) *(x+i)=double(rand())/RAND_MAX;
		tmp=cblas_dnrm2(n,x,1);
		for(int i=0;i<n;i++) *(x+i)/=tmp;
		//for(int i=0;i<n;i++) *(x+i)=1;
		cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,t,n,x,1,0,tmp_v,1);
		*(x+n)=cblas_ddot(n,x,1,tmp_v,1);
		newton(x,n,func,func_d);
		free(x); free(tmp_v); free(t);
	}
	return 0;
}
