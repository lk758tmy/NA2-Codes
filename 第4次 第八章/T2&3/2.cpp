#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
#include <lapacke.h>
using namespace std;
const double epsilon=1e-8;
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
void inverse_power(double *a,int n,double &l,double *u,double *v0,double q){
	int cnt=0,*ipiv=(int *)malloc(n*sizeof(int));
	double *m=(double *)malloc(n*n*sizeof(double));
	double *v1=(double *)malloc(n*sizeof(double));
	double *v2=(double *)malloc(n*sizeof(double)),e,m0,m1=1e9;
	cblas_dcopy(n,v0,1,v1,1);
	for(int i=0;i<n;i++) *(a+i*n+i)-=q;
	
	e=cblas_ddot(n,v1,1,u,1)/cblas_dnrm2(n,v1,1);
	//printf("%d,%.7e,%.7e\n",cnt,abs(l-q-m1),sqrt(1-e*e));
	do{
		cnt++; if(cnt==500) break;
		//gauss_colPivot(a,v2,v1,n);
		cblas_dcopy(n,v1,1,v2,1); cblas_dcopy(n*n,a,1,m,1);
		LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,1,m,n,ipiv,v2,1);
		m0=m1; m1=1/maxbar(v2,n);
		for(int i=0;i<n;i++) *(v1+ipiv[i]-1)=*(v2+i)*m1;
		
		e=1/cblas_dnrm2(n,v1,1);
		e*=cblas_ddot(n,v1,1,u,1);
		//printf("%d,%.7e,%.7e\n",cnt,abs(l-q-m1),sqrt(1-e*e));	
	}while(abs(m1-m0)>epsilon);
	printf("%d,%.7e,%.7e\n",cnt,abs(l-q-m1),sqrt(1-e*e));		
	printf("%d: %.5f\n",cnt,m1+q);
	free(v1); free(v2); free(ipiv); return ;
}
int main(){
	double *a,*v0,*u,l;
	//freopen("\data.csv","w",stdout);
	for(int n=100;n<102;n++){
	for(int q=2;q<4;q++){
		a=(double *)malloc(n*n*sizeof(double));
		v0=(double *)malloc(n*sizeof(double));
		u=(double *)malloc(n*sizeof(double));
		
		if(n==100){
			if(q==2){
				for(int i=0;i<n;i++) *(u+i)=sin((i+1)*51*M_PI/101);
				l=1/cblas_dnrm2(n,u,1);
				for(int i=0;i<n;i++) *(u+i)*=l;
				l=2-2*cos(51*M_PI/101);
				
				/*matrix_make(a,n); 
				for(int i=0;i<n;i++) printf("%.2f ",*(u+i)); printf("\n");
				cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1/l,a,n,u,1,0,v0,1);
				  for(int i=0;i<n;i++) printf("%.2f ",*(v0+i)); printf("\n");*/
			}else{
				for(int i=0;i<n;i++) *(u+i)=sin((i+1)*67/101);
				l=1/cblas_dnrm2(n,u,1);
				for(int i=0;i<n;i++) *(u+i)*=l;
				l=2-2*cos(67*M_PI/101);
			}
		}if(n==101){
			if(q==2){
				for(int i=0;i<n;i++) *(u+i)=sin((i+1)*M_PI/2); // 51/120
				l=1/cblas_dnrm2(n,u,1);
				for(int i=0;i<n;i++) *(u+i)*=l;
				l=2;
			}else{
				for(int i=0;i<n;i++) *(u+i)=sin((i+1)*2*M_PI/3); // 68/102
				l=1/cblas_dnrm2(n,u,1);
				for(int i=0;i<n;i++) *(u+i)*=l;
				l=3;
			}
		}
		
		srand(time(0));
		for(int i=0;i<n;i++) *(v0+i)=double(rand())/RAND_MAX;
		
		matrix_make(a,n); inverse_power(a,n,l,u,v0,q+0.001);
		free(u); free(a);		
	}}

	return 0;
}
