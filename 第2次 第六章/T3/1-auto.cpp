#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cblas.h>
using namespace std;
using namespace std::chrono;
const double epsilon=0.000001;
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
bool hault(int n,double *x,double *m,double *c,bool flag){
    double *e=(double *)malloc(n*sizeof(double)),AError,Residual;
    for(int i=0;i<n;i++) *(e+i)=*(x+i)-1;
    AError=cblas_dnrm2(n,e,1);
    if(flag==true){
        for(int i=0;i<n;i++) //复用一下变量e
            *(e+i)=cblas_ddot(n,m+i*n,1,x,1);
        cblas_daxpy(n,-1,c,1,e,1);
        Residual=cblas_dnrm2(n,e,1);
        printf("%.6f,%.6f\n",AError,Residual);
    }
    return (AError<epsilon);
}
int SOR(double *m,double *x,double *c,int n,double omega){
    double *y=(double *)malloc(n*sizeof(double));
    int cnt=0; bool flag;
    printf("0,"); hault(n,x,m,c,true);
    do{
        cnt++;
        for(int i=0;i<n;i++){
            *(y+i)=*(c+i)-cblas_ddot(i,m+i*n,1,y,1);
            *(y+i)-=cblas_ddot(n-i-1,m+i*n+i+1,1,x+i+1,1);
            *(y+i)*=(omega/(*(m+i*n+i)));
            *(y+i)+=(1-omega)*(*(x+i));
        }
        cblas_dswap(n,x,1,y,1);
        if(cnt%5==0){
            printf("%d,",cnt); flag=true;
        }else flag=false;
    }while(!hault(n,x,m,c,flag));
    printf("%d,",cnt); hault(n,x,m,c,true);
    free(y);
    return cnt;
}
int jordan(double *m,double *x,double *c,int n){
    double *y=(double *)malloc(n*sizeof(double)),*tmp;
    int cnt=0; bool flag;
    printf("0,"); hault(n,x,m,c,true);
    do{
        cnt++;
        for(int i=0;i<n;i++){
            *(y+i)=*(c+i)-cblas_ddot(n,x,1,m+i*n,1);
            *(y+i)+=*(x+i)*(*(m+i*n+i));
            *(y+i)/=*(m+i*n+i);
        }
        cblas_dswap(n,x,1,y,1);
        if(cnt%200==0){
            printf("%d,",cnt); flag=true;
        }else flag=false;
    }while(!hault(n,x,m,c,flag));
    printf("%d,",cnt); hault(n,x,m,c,true);
    free(y);
    return cnt;
}
int GS(double *m,double *x,double *c,int n){
    double *y=(double *)malloc(n*sizeof(double)),*tmp;
    int cnt=0; bool flag;
    printf("0,"); hault(n,x,m,c,true);
    do{
        cnt++;
        for(int i=0;i<n;i++){;
            *(y+i)=*(c+i)-cblas_ddot(i,m+i*n,1,y,1);
            *(y+i)-=cblas_ddot(n-i-1,m+i*n+i+1,1,x+i+1,1);
            *(y+i)/=*(m+i*n+i);
        }
        cblas_dswap(n,x,1,y,1);
        if(cnt%100==0){
            printf("%d,",cnt); flag=true;
        }else flag=false;
    }while(!hault(n,x,m,c,flag));
    printf("%d,",cnt); hault(n,x,m,c,true);
    free(y);
    return cnt;
}
int _main(int n,double omega) {
	int n2=n*n,n4=n2*n2,cnt;
	double *m=(double *)malloc(n4*sizeof(double));
	matrix_make(m,n);
	double *c=(double *)malloc(n2*sizeof(double));
	double *x=(double *)malloc(n2*sizeof(double));
	vector_make(c,n);
	time_point<steady_clock> s,e;

    printf("Size: 30\n");
	printf("Type: Jordan\nStep,AError,Residual\n");
	for(int i=0;i<n2;i++) *(x+i)=0;
	s=steady_clock::now();
	cnt=jordan(m,x,c,n2);
	e=steady_clock::now();

	printf("Type: GS\nStep,AError,Residual\n");
	for(int i=0;i<n2;i++) *(x+i)=0;
	s=steady_clock::now();
	cnt=GS(m,x,c,n2);
	e=steady_clock::now();

	printf("Type: SOR,Omega: %.2f\n",omega);
    printf("Step,AError,Residual\n");
	for(int i=0;i<n2;i++) *(x+i)=0;
	s=steady_clock::now();
	cnt=SOR(m,x,c,n2,omega);
	e=steady_clock::now();

	free(m); free(c); free(x);
	return 0;
}
int main(){
    _main(30,1.82);
    return 0;
}
