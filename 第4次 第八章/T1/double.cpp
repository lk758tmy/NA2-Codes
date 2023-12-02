#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
using namespace std;
const double epsilon=1e-8;
void matrix_make(double *m,int n){
	for(int i=0;i<n-1;i++){
		*(m+i*n+i)=2; *(m+i*n+i+1)=*(m+i*n+i+n)=-1;
		for(int j=2+i;j<n;j++)
			*(m+i*n+j)=*(m+j*n+i)=0;
	}
	return ;
}
double maxbar(double *v,int n){
	double ans=0;
	for(int i=0;i<n;i++)
		ans=(abs(*(v+i))>abs(ans))?*(v+i):ans;
	return ans;
}
void power(double *a,int n,char s,double &l,double *u,double *v0){
	//s: 'A'=Atiken 'R'=Rayleigh(nrm2) 'r'=Rayleigh(maxbar) 
	int cnt=0;
	double *v1=(double *)malloc(n*sizeof(double)),e;
	double *v2=(double *)malloc(n*sizeof(double)),m0,m1=1e9,m2,m3,m4;
	cblas_dcopy(n,v0,1,v1,1);
	//srand(time(0));
	//for(int i=0;i<n;i++) *(v1+i)=double(rand())/RAND_MAX;
	if(s=='R'/*||s=='r'*/) cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,a,n,v1,1,0,v2,1);
	do{
		m0=m1;
		if(s=='R'/*||s=='r'*/){
			/*if(s=='R')*/ m1=1/cblas_dnrm2(n,v2,1);
			//else m1=1/maxbar(v2,n);
			for(int i=0;i<n;i++) *(v1+i)=*(v2+i)*m1;
			cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,a,n,v1,1,0,v2,1);
			m1=cblas_ddot(n,v1,1,v2,1);
			//if(s=='r') m1/=cblas_ddot(n,v1,1,v1,1);
		}else{
			cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,a,n,v1,1,0,v2,1);
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
		if(cnt%10==0){
			if(s=='R') e=1;
			else e=1/cblas_dnrm2(n,v1,1);
			e*=cblas_ddot(n,v1,1,u,1);
			printf("%d,%.7e,%.7e\n",cnt,abs(l-m1),sqrt(1-e*e));			
		}
		cnt++;
	}while(abs(m1-m0)>epsilon);
	printf("%d,%.7e,%.7e\n",cnt,abs(l-m1),sqrt(1-e*e));
	printf("\n");
	free(v1); free(v2); return ;
}
int main(){
	double *m,*v0,*u,l;
	freopen("\data.csv","w",stdout);
	for(int n=100;n<102;n++){
		m=(double *)malloc(n*n*sizeof(double));
		v0=(double *)malloc(n*sizeof(double));
		
		l=sqrt(2.0/(n+1));
		u=(double *)malloc(n*sizeof(double));
		for(int i=0;i<n;i++) *(u+i)=l*sin((i+1)*(1-1.0/(n+1))*M_PI);
		l=1/cblas_dnrm2(n,u,1);
		for(int i=0;i<n;i++) *(u+i)*=l;
		l=2*(1+cos(M_PI/(n+1)));
		
		srand(time(0));
		for(int i=0;i<n;i++) *(v0+i)=double(rand())/RAND_MAX;
		matrix_make(m,n); power(m,n,'O',l,u,v0);
		//实际应用时l、u并不知道，而v0的随机选取应在子程序内部执行
		matrix_make(m,n); power(m,n,'A',l,u,v0);
		//matrix_make(m,n); power(m,n,'r',l,u,v0);
		matrix_make(m,n); power(m,n,'R',l,u,v0);
		free(u); free(m);
	}
	return 0;
}
