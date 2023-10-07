#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cblas.h>
using namespace std;
using namespace std::chrono;
const float epsilon=0.000001;
float AError;
void matrix_make(float *h,int n){
	int n2=n*n; float *p=h;
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
void vector_make(float *h,int n){
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
bool hault(int n,float *x){
    float *e=(float *)malloc(n*sizeof(float));
    for(int i=0;i<n;i++) *(e+i)=*(x+i)-1;
    AError=cblas_snrm2(n,e,1);
    //printf("%.8f\n",AError);
    return (AError<epsilon);
}
int SOR(float *m,float *x,float *c,int n,float omega){
    float *y=(float *)malloc(n*sizeof(float));
    int cnt=0;
    do{
        cnt++;
        for(int i=0;i<n;i++){
            *(y+i)=*(c+i)-cblas_sdot(i,m+i*n,1,y,1);
            *(y+i)-=cblas_sdot(n-i-1,m+i*n+i+1,1,x+i+1,1);
            *(y+i)*=(omega/(*(m+i*n+i)));
            *(y+i)+=(1-omega)*(*(x+i));
        }
        cblas_sswap(n,x,1,y,1);
        //for(int i=0;i<n;i++) printf("%f ",*(x+i)); printf("\n");
    }while(!hault(n,x));
    free(y);
    return cnt;
}
int main() {
	int n,n2,n4,cnt; float omega;
	scanf("%d %f",&n,&omega); n2=n*n; n4=n2*n2; //n should be at least 3!
	float *m=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	vector_make(c,n);
	for(int i=0;i<n2;i++) *(x+i)=0;

	time_point<steady_clock> s,e;
	s=steady_clock::now();
	cnt=SOR(m,x,c,n2,omega);
    for(int i=0;i<n2;i++) printf("%f ",*(x+i)); printf("\n");
	e=steady_clock::now();

	printf("Size: %d\n",n);
	printf("Count: %d\n",cnt);
	printf("Time(ms): %.3f\n\n",(e-s).count()/1000000.0);
	printf("A.Error: %.8f\nR.Error: %.8f\n\n",AError,AError/n);
	free(m); free(c); free(x);
	return 0;
}
