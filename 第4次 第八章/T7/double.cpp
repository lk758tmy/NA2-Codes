#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <algorithm>
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
double sgn(double x){
	if(x>=0) return 1;
	else return -1;
}
int eigen_QR_tri(double *a,int n){
	double x,y,u,alpha,c,s,errf; int cnt=0;
	do{
		alpha=(*(a+n*n-n-2)-*(a+n*n-1))/2;
		u=*(a+n*n-1)+alpha-sgn(alpha)*sqrt(alpha*alpha+*(a+n*n-2)*(*(a+n*n-2)));
		x=*a-u; y=*(a+n);
		for(int i=0;;i++){
			cblas_drotg(&x,&y,&c,&s);
			cblas_drot(n,a+i*n,1,a+(i+1)*n,1,c,s);
			cblas_drot(n,a+i,n,a+i+1,n,c,s);
			if(i==n-2) break;
			x=*(a+n*(i+1)+i); y=*(a+n*(i+2)+i);		
		}
		errf=cblas_dnrm2(n*n,a,1)-cblas_dnrm2(n,a,n+1);
		cnt++;
	}while(errf>n*epsilon);
	printf("%.3e,",errf);
	return cnt;
}
int main(){
	double *m,*eR,*eN;
	for(int n=90;n<110;n++){
		m=(double *)malloc(n*n*sizeof(double));
		eR=(double *)malloc(n*sizeof(double));
		eN=(double *)malloc(n*sizeof(double));
		for(int i=0;i<n;i++) *(eR+i)=2-2*cos((i+1)*M_PI/(n+1));
		
		printf("%d,",n);
		matrix_make(m,n);
		printf("%d,",eigen_QR_tri(m,n));
		cblas_dcopy(n,m,n+1,eN,1);
		sort(eN,eN+n);
		cblas_daxpy(n,-1,eN,1,eR,1);
		printf("%.3e\n",cblas_dnrm2(n,eR,1));

		free(eR); free(eN); free(m);
	}
	return 0;
}
