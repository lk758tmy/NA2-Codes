#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
#include <chrono>
#include <lapacke.h>
using namespace std;
using namespace std::chrono;
const double epsilon=1e-10;
void matrix_make(double *m,int n){
	for(int i=0;i<n-1;i++){
		*(m+i*n+i)=2; *(m+i*n+i+1)=*(m+i*n+i+n)=-1;
		for(int j=2+i;j<n;j++)
			*(m+i*n+j)=*(m+j*n+i)=0;
	}
	*(m+n*n-1)=2;
	return ;
}
double maxbar(double *v,int n){
	double ans=0;
	for(int i=0;i<n;i++)
		ans=(abs(*(v+i))>abs(ans))?*(v+i):ans;
	return ans;
}
double sgn(double x){
	if(x>=0) return 1;
	else return -1;
}
void subspace(double *a,int n,double l1r,double l2r,double *u1,double *u2){
	double *v1=(double *)malloc(2*n*sizeof(double)),tmp;
	double *v2=(double *)malloc(2*n*sizeof(double)); //ColMajor
	double *r=(double *)malloc(4*sizeof(double)),l1,l2,l0,e;
	srand(time(0));
	for(int i=0;i<2*n;i++) *(v1+i)=double(rand())/RAND_MAX;
	tmp=1/cblas_dnrm2(n,v1,2);
	for(int i=0;i<n;i++) *(v1+i*2)*=tmp;
	cblas_daxpy(n,-tmp*cblas_ddot(n,v1,2,v1+1,2),v1,2,v1+1,2);
	tmp=1/cblas_dnrm2(n,v1+1,2);
	for(int i=0;i<n;i++) *(v1+1+i*2)*=tmp;
	
	int cnt=0;
	do{
		cnt++; //if(cnt==1000) break;
		l0=(l1>l2)?l1:l2;
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,2,n,1,a,n,v1,2,0,v2,2);
		
		cblas_dcopy(n*2,v2,1,v1,1);//MGS-QR
		*r=cblas_dnrm2(n,v1,2);//Step 1
		for(int i=0;i<n;i++) *(v1+2*i)/=*r;
		*(r+1)=cblas_ddot(n,v1,2,v1+1,2);
		cblas_daxpy(n,-*(r+1),v1,2,v1+1,2);
		*(r+3)=cblas_dnrm2(n,v1+1,2);//Step 2
		for(int i=0;i<n;i++) *(v1+1+2*i)/=*(r+3);
		
		l1=*r; l2=*(r+3);
		//if(cnt%100==0) printf("%d %.10f %.10f\n",cnt,l1,l2);
	}while(abs(l0-((l1>l2)?l1:l2))>epsilon);
	e=cblas_ddot(n,u1,1,v1+1,2)/cblas_dnrm2(n,u1,1)/cblas_dnrm2(n,v1+1,2);
	printf("%.10f,%.10f,%.10f\n",l1,abs(l1-l1r),sqrt(1-e*e));
	e=cblas_ddot(n,u2,1,v1,2)/cblas_dnrm2(n,u2,1)/cblas_dnrm2(n,v1,2);
	printf("%.10f,%.10f,%.10f\n",l2,abs(l2-l2r),sqrt(1-e*e));
	free(v1); free(v2); free(r);return ;
}
int main(){
	double *m,*u1,*u2;
	for(int n=100;n<102;n++){
		printf("%d\n",n);
		m=(double *)malloc(n*n*sizeof(double));
		u1=(double *)malloc(n*sizeof(double));
		u2=(double *)malloc(n*sizeof(double));
		matrix_make(m,n);
		for(int i=0;i<n;i++){
			*(u1+i)=sin((i+1)*(n/(n+1.0))*M_PI);
			*(u2+i)=sin((i+1)*((n-1)/(n+1.0))*M_PI);
		}
		time_point<steady_clock> st,ed;
		st=steady_clock::now();
		subspace(m,n,2+2*cos(2*M_PI/(n+1)),2+2*cos(M_PI/(n+1)),u2,u1);
		ed=steady_clock::now();
		printf("%.3fms\n",(ed-st).count()/1000000.0);		
		free(m); free(u1); free(u2);
	}
	return 0;
}
