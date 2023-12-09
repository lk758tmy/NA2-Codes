#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
using namespace std;
const double epsilon=1e-8;
void matrix_make(double *a,int n){
	double *b=a+n,*c=a+2*n;
	for(int i=0;i<n;i++){
		*(a+i)=2; *(b+i)=*(c+i)=-1;
	}
	return ;
}
double maxbar(double *v,int n){
	double ans=0;
	for(int i=0;i<n;i++)
		ans=(abs(*(v+i))>abs(ans))?*(v+i):ans;
	return ans;
}
void inverse_power_tri(double *a,int n,double *u,double *v0,double q){
	//s: 'A'=Atiken 'R'=Rayleigh(nrm2) 'r'=Rayleigh(maxbar) 
	int cnt=0;
	double *v1=(double *)malloc(n*sizeof(double)),*b=a+n,*c=a+2*n;
	double *v2=(double *)malloc(n*sizeof(double)),e,m0,m1=1e9;
	cblas_dcopy(n,v0,1,v1,1);
	for(int i=0;i<n;i++) *(a+i)-=q;
	for(int i=1;i<n;i++){
		*(b+i)/=*(a+i-1); *(a+i)-=*(b+i)*(*(c+i-1));
	}  //三对角初始化
	for(int i=0;i<n;i++) printf("%.3e ",*(b+i));
	printf("\n");
	
	do{
		cnt++;
		for(int i=1;i<n;i++) *(v1+i)-=*(v1+i-1)*(*(c+i-1));
		*(v2+n-1)=*(v1+n-1)/(*(a+n-1));
		for(int i=n-2;i>-1;i--)
			*(v2+i)=(*(v1+i)-*(c+i)*(*(v2+i+1)))/(*(a+i));
		m0=m1; m1=maxbar(v2,n);
		
		for(int i=0;i<n;i++) *(v2+i)/=m1;
		cblas_dcopy(n,v2,1,v1,1);
		
	}while(abs(m1-m0)>epsilon);
	printf("%d: %.3e\n",cnt,m1);
	for(int i=0;i<n;i++) printf("%.3e ",*(v2+i));
	printf("\n");
	free(v1); free(v2); return ;
}
int main(){
	double *a,*v0,*u,l;
	//freopen("\data.csv","w",stdout);
	//for(int n=100;n<102;n++){
	int n=10;
		a=(double *)malloc(3*n*sizeof(double));
		v0=(double *)malloc(n*sizeof(double));
		
		l=sqrt(2.0/(n+1));
		u=(double *)malloc(n*sizeof(double));
		for(int i=0;i<n;i++) *(u+i)=l*sin((i+1)*(1-1.0/(n+1))*M_PI);
		l=1/cblas_dnrm2(n,u,1);
		for(int i=0;i<n;i++) *(u+i)*=l;
		l=2*(1+cos(M_PI/(n+1)));
		
		srand(time(0));
		for(int i=0;i<n;i++) *(v0+i)=double(rand())/RAND_MAX;
		matrix_make(a,n); inverse_power_tri(a,n,u,v0,2);
		matrix_make(a,n); inverse_power_tri(a,n,u,v0,3);
		free(u); free(a);
	//}
	return 0;
}
