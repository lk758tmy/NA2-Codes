#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cblas.h>
#include <cmath>
using namespace std;
using namespace std::chrono;
void make_matrix_t(float *a,float *b,float *c,int n){
    for(int i=0;i<n;i++){
       *(a+i)=*(c+i)=-1; *(b+i)=2;
    }
    return ;
}
void make_vector(float *a,float *b,float *c,int n,float *y){
    for(int i=0;i<n;i++) *(y+i)=*(a+i)+*(b+i)+*(c+i);
    *y-=*a; *(y+n-1)-=*(c+n-1);
    return ;
}
void calc(float *a,float *b,float *c,int n,float *y,float *x){
    float *ap=a,*bp=b,*cp=c,*yp=y,*xp;
    *cp/=*bp; *yp/=*bp; *bp=1;
    for(int i=1;i<n;i++){
        ap++; bp++; *bp-=((*ap)*(*cp));
        yp++; *yp=(*yp-((*ap)*(*(yp-1))))/(*bp);
        cp++; *cp/=*bp; *bp=1;
    }
    cblas_scopy(n,y,1,x,1);
    xp=x+n-1; cp=c+n-1;
    for(int i=1;i<n;i++){
        xp--; cp--; *xp-=(*cp)*(*(xp+1));
    }
    return ;
}
void _main(int n){
    time_point<steady_clock> s,e;
    float *a=(float *)malloc(n*sizeof(float));
    float *b=(float *)malloc(n*sizeof(float));
    float *c=(float *)malloc(n*sizeof(float));
    float *x=(float *)malloc(n*sizeof(float));
    float *y=(float *)malloc(n*sizeof(float)),error;

    make_matrix_t(a,b,c,n);
    make_vector(a,b,c,n,y);
    s=steady_clock::now();
    calc(a,b,c,n,y,x);
    e=steady_clock::now();
    printf("Type:1\nSize: %d\nTime(ms): %.3f\n",n,(e-s).count()/1000000.0);
    for(int i=0;i<n;i++) *(x+i)-=1;
    error=cblas_snrm2(n,x,1);
    printf("A.Error: %.6f\nR.Error: %.6f\n\n",error,error/sqrt(n));

    free(a); free(b); free(c); free(x); free(y);
    return ;
}
int main(){
    for(int n=1000;n<1280001;n*=2) _main(n);
    return 0;
}
