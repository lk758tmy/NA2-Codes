#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory.h>
#include <chrono>
#include <cblas.h>
using namespace std;
using namespace std::chrono;
void matrix_make(float *h,int n,float a){
	memset(h,0,n*n*sizeof(float));
	for(int i=0;i<n;i++) *(h+i*n+i)=1;
	for(int i=1;i<n;i++) *(h+i)=*(h+i*n)=a;
	return ;
}
void crout(float *m,int n){

    return ;
}
void _calc(float* m,int n,float a){
    matrix_make(m,n,a);
    time_point<steady_clock> s,e;
	s=steady_clock::now();
    crout(m,n);
	e=steady_clock::now();
	printf("Value of a: %f\n",a);
	printf("Size: %d\n",n);
	printf("Time(ms): %.3f\n\n",(e-s).count()/1000000.0);
	return ;
}
void _main(int n) {
	int n2=n*n;
	float *m=(float *)malloc(n2*sizeof(float));
	_calc(m,n,0.001);
	_calc(m,n,1000);
	free(m);
	return ;
}
int main(){
    //for(int n=500;n<10001;n+=500) _main(n);
    _main(10);
    return 0;
}
