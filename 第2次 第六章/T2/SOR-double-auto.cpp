#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
using namespace std;
const double epsilon=0.000001;
double AError;
void matrix_make(double *h,int n){
	int n2=n*n; double *p=h;
	for(int i=0;i<n2*n2;i++){*p=0; p++;}
    p=h;
    *h=4;
	for(int i=1;i<n2;i++){
		p+=(n2+1); *p=4;
		if(i%n!=0) *(p-1)=*(p-n2)=-1;
	}
	for(int i=n,j=0;i<n2;i++,j++)
		*(h+i*n2+j)=*(h+j*n2+i)=-1;
	return ;
}
void vector_make(double *h,int n){
	*h=2;
	for(int i=1;i<n-1;i++){h++; *h=1;}
	h++; *h=2;
	for(int i=1;i<n-1;i++){
		h++; *h=1;
		for(int j=1;j<n-1;j++){h++; *h=0;}
		h++; *h=1;
	}
	h++; *h=2;
	for(int i=1;i<n-1;i++){h++; *h=1;}
	h++; *h=2;
	return ;
}
bool hault(int n,double *x){
    double *e=(double *)malloc(n*sizeof(double));
    for(int i=0;i<n;i++) *(e+i)=*(x+i)-1;
    AError=cblas_dnrm2(n,e,1);
    return (AError<epsilon);
}
int SOR(double *m,double *x,double *c,int n,double omega){
    double *y=(double *)malloc(n*sizeof(double));
    int cnt=0;
    do{
        cnt++;
        for(int i=0;i<n;i++){
            *(y+i)=*(c+i)-cblas_ddot(i,m+i*n,1,y,1);
            *(y+i)-=cblas_ddot(n-i-1,m+i*n+i+1,1,x+i+1,1);
            *(y+i)*=(omega/(*(m+i*n+i)));
            *(y+i)+=(1-omega)*(*(x+i));
        }
        cblas_dswap(n,x,1,y,1);
    }while(!hault(n,x));
    free(y);
    return cnt;
}
int _main(int n,double omega) {
	int n2=n*n,n4=n2*n2,cnt; //n should be at least 3!
	double *m=(double *)malloc(n4*sizeof(double));
	matrix_make(m,n);
	double *c=(double *)malloc(n2*sizeof(double));
	double *x=(double *)malloc(n2*sizeof(double));
	vector_make(c,n);
	for(int i=0;i<n2;i++) *(x+i)=0;

	clock_t s,e;
	s=clock();
	cnt=SOR(m,x,c,n2,omega);
	e=clock();

	printf("Size: %d\n",n);
	printf("Omega: %.2f\n",omega);
	printf("Count: %d\n",cnt);
	printf("Time(ms): %.3f\n",float(e-s)/CLOCKS_PER_SEC*1000);
	printf("\n"); //printf("A.Error: %.8f\nR.Error: %.8f\n\n",AError,AError/n);
	free(m); free(c); free(x);
	return cnt;
}
void _calc(int n){
    int min_cnt=1e9+7,cnt; double min_omega;
    for(double omega=0.11;omega<2;omega+=0.01){
        cnt=_main(n,omega);
        if(cnt<min_cnt){
            min_cnt=cnt; min_omega=omega;
        }
    }
    printf("--------\nSize: %d\n",n);
    printf("Min-Step-Count: %d\n",min_cnt);
    printf("Best-Omega: %.2f\n--------\n\n",min_omega);
    return ;
}
int main(){
    srand(time(0));
    //for(int n=5;n<30;n+=(rand()%4+1))
        //_calc(n);
    _calc(28); _calc(30);
    return 0;
}
