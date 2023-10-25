#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cblas.h>
using namespace std;
using namespace std::chrono;
const float epsilon=0.000001;
float last_delta=0;
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
bool hault(int type,int n,float *x,float *y,float *m,float *c){
    if(type==1){
        float *r=(float *)malloc(n*sizeof(float));
        for(int i=0;i<n;i++)
            *(r+i)=cblas_sdot(n,m+i*n,1,x,1);
        cblas_saxpy(n,-1,c,1,r,1);
        if(cblas_snrm2(n,r,1)<epsilon){
            free(r); return true;
        }else{
            free(r); return false;
        }
    }else if(type==2){
        cblas_saxpy(n,-1,x,1,y,1);
        //printf("%.6f ",cblas_snrm2(n,y,1));
        return (cblas_snrm2(n,y,1)<epsilon);
    }else{
        cblas_saxpy(n,-1,x,1,y,1);
        float delta=cblas_snrm2(n,y,1);
        if(last_delta<0.0000001){
            last_delta=delta; return false;
        }
        //printf("%.6f ",abs(delta*delta/(last_delta-delta)));
        if(abs(delta*delta/(last_delta-delta))<epsilon){
            last_delta=delta; return true;
        }else{
            last_delta=delta; return false;
        }
    }
}
int jacobi(float *m,float *x,float *c,int n,int type){
    //type=1: 残量准则 2：相邻误差准则 3：后验误差准则
    float *y=(float *)malloc(n*sizeof(float)),*tmp;
    int cnt=0;
    do{
        cnt++;
        for(int i=0;i<n;i++){
            *(y+i)=*(c+i)-cblas_sdot(n,x,1,m+i*n,1);
            *(y+i)+=*(x+i)*(*(m+i*n+i));//用sdot多减了一项
            *(y+i)/=*(m+i*n+i);

        }
        cblas_sswap(n,x,1,y,1); //tmp=x; x=y; y=tmp;
    }while(!hault(type,n,x,y,m,c));
    free(y);
    return cnt;
}
int main() {
	int n,n2,n4,cnt,type;
	scanf("%d %d",&n,&type); n2=n*n; n4=n2*n2; //n should be at least 3!
	float *m=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	vector_make(c,n);
	for(int i=0;i<n2;i++) *(x+i)=0;

	time_point<steady_clock> s,e;
	s=steady_clock::now();
	last_delta=0;
	cnt=jacobi(m,x,c,n2,type);
    //for(int i=0;i<n2;i++) printf("%f ",*(x+i)); printf("\n");
	e=steady_clock::now();

	//for(int i=0;i<n2;i++) *(x+i)-=1;
	//float error=cblas_snrm2(n2,x,1);
	printf("Size: %d\n",n);
	printf("Count: %d\n",cnt);
	//printf("A.Error: %.6f\nR.Error: %.6f\n",error,error/n);
	printf("Time(ms): %.3f\n\n",(e-s).count()/1000000.0);
	free(m); free(c); free(x);
	return 0;
}
