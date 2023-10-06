#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cblas.h>
#include <cmath>
using namespace std;
using namespace std::chrono;
void make_matrix_tb(double *a,double *b,double *c,int n){
    for(int i=0;i<n;i++){
       *(a+i)=*(c+i)=-1; *(b+i)=2;
    }
    double tmp=1;
    for(int i=0;i<n;i++){
        tmp/=2.0; *(a+i)*=tmp; *(b+i)*=tmp; *(c+i)*=tmp;
    }
    return ;
}
void make_vector(double *a,double *b,double *c,int n,double *y){
    for(int i=0;i<n;i++) *(y+i)=*(a+i)+*(b+i)+*(c+i);
    *y-=*a; *(y+n-1)-=*(c+n-1);
    return ;
}
void calc(double *a,double *b,double *c,int n,double *y,double *x){
    double *ap=a,*bp=b,*cp=c,*yp=y,*xp;
    *cp/=*bp; *yp/=*bp; *bp=1;
    for(int i=1;i<n;i++){
        ap++; bp++; *bp-=((*ap)*(*cp));
        yp++; *yp=(*yp-((*ap)*(*(yp-1))))/(*bp);
        cp++; *cp/=*bp; *bp=1;
    }
    cblas_dcopy(n,y,1,x,1);
    xp=x+n-1; cp=c+n-1;
    for(int i=1;i<n;i++){
        xp--; cp--; *xp-=(*cp)*(*(xp+1));
    }
    return ;
}
void _main(int n){
    time_point<steady_clock> s,e;
    double *a=(double *)malloc(n*sizeof(double));
    double *b=(double *)malloc(n*sizeof(double));
    double *c=(double *)malloc(n*sizeof(double));
    double *x=(double *)malloc(n*sizeof(double));
    double *y=(double *)malloc(n*sizeof(double)),error;
    long error_ratio=1e12;

    make_matrix_tb(a,b,c,n);
    make_vector(a,b,c,n,y);
    s=steady_clock::now();
    calc(a,b,c,n,y,x);
    e=steady_clock::now();
    printf("Type:2\nSize: %d\nTime(ms): %.3f\n",n,(e-s).count()/1000000.0);
    for(int i=0;i<n;i++) *(x+i)-=1;
    error=cblas_dnrm2(n,x,1);
    printf("A.Error: %.6f x10^(-12)\n",error*error_ratio);
    printf("R.Error: %.6f x10^(-12)\n\n",error/sqrt(n)*error_ratio);

    free(a); free(b); free(c); free(x); free(y);
    return ;
}
int main(){
    for(int n=5;n<1281;n*=2) _main(n);
    return 0;
}
