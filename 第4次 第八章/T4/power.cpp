#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
#include <chrono>
using namespace std;
using namespace std::chrono;
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
double sgn(double x){
	if(x>=0) return 1;
	else return -1;
}
void power(double *a,int n,int lda,char s,double &l,double *u){
	//s: 'A'=Atiken 'R'=Rayleigh
	int cnt=0;
	double *v1=(double *)malloc(n*sizeof(double));
	double *v2=(double *)malloc(n*sizeof(double)),m0,m1=1e9,m2,m3,m4;
	srand(time(0));
	for(int i=0;i<n;i++) *(v1+i)=double(rand())/RAND_MAX;
	if(s=='R') cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,a,lda,v1,1,0,v2,1);
	do{
		m0=m1;
		if(s=='R'){
			m1=1/cblas_dnrm2(n,v2,1);
			for(int i=0;i<n;i++) *(v1+i)=*(v2+i)*m1;
			cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,a,lda,v1,1,0,v2,1);
			m1=cblas_ddot(n,v1,1,v2,1);
		}else{
			cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,a,lda,v1,1,0,v2,1);
			m1=maxbar(v2,n);
			for(int i=0;i<n;i++) *(v1+i)=*(v2+i)/m1;
			if(s=='A'){
				if(cnt==1){
					m2=m1;
				}if(cnt==2){
					m3=m2; m2=m1;
				}else{
					m4=m3; m3=m2; m2=m1;
					m1=m4-(m3-m4)*(m3-m4)/(m2-2*m3+m4);
				}
			}
		}
		cnt++;
	}while(abs(m1-m0)>epsilon);
	l=m1; cblas_dcopy(n,v1,1,u,1);
	free(v1); free(v2); return ;
}
void householderTransformSymmetric(double *a,int n,int lda,double *u,double b){
	double *y=(double *)malloc(n*sizeof(double));
	cblas_dgemv(CblasRowMajor,CblasTrans,n,n,1,a,lda,u,1,0,y,1);
	cblas_dger(CblasRowMajor,n,n,-1/b,u,1,y,1,a,lda);
	cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,a,lda,u,1,0,y,1);
	cblas_dger(CblasRowMajor,n,n,-1/b,y,1,u,1,a,lda);	
	free(y);
	return ;
}
void householderTransformLeftVector(double *a,int n,double *u,double b){
	b=-cblas_ddot(n,u,1,a,1)/b;
	cblas_daxpy(n,b,u,1,a,1);
	return ;
}
int main(){
	double *m,*u,*w,*u2,*real,l,l2,alpha,b,e;
	time_point<steady_clock> st,ed;
	//freopen("\data.csv","w",stdout);
	for(int n=100;n<102;n++){
		printf("%d\n",n);
		
		m=(double *)malloc(n*n*sizeof(double));
		u=(double *)malloc(n*sizeof(double));
		u2=(double *)malloc(n*sizeof(double));
		real=(double *)malloc(n*sizeof(double));
		w=(double *)malloc(n*sizeof(double));
		matrix_make(m,n);
		
		st=steady_clock::now();
		for(int i=0;i<n;i++) *(real+i)=sin((i+1)*(n/(n+1.0))*M_PI);
		power(m,n,n,'R',l,u);
		e=cblas_ddot(n,u,1,real,1)/cblas_dnrm2(n,u,1)/cblas_dnrm2(n,real,1);
		e=sqrt(1-e*e);
		printf("%.10f,%.10f,%.10f\n",l,abs(l-2-2*cos(M_PI/(n+1))),e);
		
		alpha=-sgn(*u)*cblas_dnrm2(n,u,1);
		b=alpha*(alpha-*u); *u-=alpha;
		householderTransformSymmetric(m,n,n,u,b);
		cblas_dcopy(n-1,m+1,1,w+1,1);
		/*for(int i=0;i<n;i++){
			for(int j=0;j<n;j++) printf("%.2f ",*(m+i*n+j));
			printf("\n");
		}printf("\n");*/
		
		for(int i=0;i<n;i++) *(real+i)=sin((i+1)*((n-1)/(n+1.0))*M_PI);
		power(m+n+1,n-1,n,'R',l2,u2+1);
		*u2=cblas_ddot(n-1,w+1,1,u2+1,1)/(l2-l);
		householderTransformLeftVector(u2,n,u,b);
		e=cblas_ddot(n,u2,1,real,1)/cblas_dnrm2(n,u2,1)/cblas_dnrm2(n,real,1);
		e=sqrt(1-e*e);
		printf("%.10f,%.10f,%.10f\n",l2,abs(l2-2-2*cos(2*M_PI/(n+1))),e);
		ed=steady_clock::now();
		
		printf("%.3fms\n",(ed-st).count()/1000000.0);
		free(u); free(u2); free(m); free(real); free(w);
	}
	return 0;
}
