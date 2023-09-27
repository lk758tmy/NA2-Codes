#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory.h>
#include <chrono>
#include <cblas.h>
using namespace std;
using namespace std::chrono;
void matrix_make(float *h,int n,float a,int type){
	memset(h,0,n*n*sizeof(float));
	for(int i=0;i<n;i++) *(h+i*n+i)=1;
	if(type==1)
        for(int i=1;i<n;i++) *(h+n*n-i-1)=*(h+i*n-1)=a;
    else
        for(int i=1;i<n;i++) *(h+i)=*(h+i*n)=a;
	return ;
}
void crout(float *m,int n){
    float s;
    for(int k=0;k<n;k++){
        for(int i=k;i<n;i++){
            s=cblas_sdot(k,m+i*n,1,m+k,n); //s=0;
            //for(int r=0;r<k;r++) s+=*(m+i*n+r)*(*(m+r*n+k));
            *(m+i*n+k)-=s;
        }
        for(int j=k+1;j<n;j++){
            s=cblas_sdot(k,m+k*n,1,m+j,n); //s=0;
            //for(int r=0;r<k;r++) s+=*(m+k*n+r)*(*(m+r*n+j));
            *(m+k*n+j)-=s; *(m+k*n+j)/=*(m+k*n+k);
        }
    }
    return ;
}
void _calc(float* m,int n,float a,int type){
    matrix_make(m,n,a,type);
    time_point<steady_clock> s,e;
	s=steady_clock::now();
    crout(m,n);
	e=steady_clock::now();
	printf("Type: %d\n",type);
	printf("Value of a: %f\n",a);
	printf("Size: %d\n",n);
	printf("Time(ms): %.3f\n\n",(e-s).count()/1000000.0);
	return ;
}
void _main(int n,int type) {
	int n2=n*n;
	float *m=(float *)malloc(n2*sizeof(float));
	_calc(m,n,0.001,type);
	_calc(m,n,1000,type);
	free(m);
	return ;
}
int main(){
    for(int n=250;n<2501;n+=250) _main(n,1);
    for(int n=250;n<2501;n+=250) _main(n,2);
    return 0;
}
