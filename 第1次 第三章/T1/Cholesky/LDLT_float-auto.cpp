#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory.h>
#include <time.h>
#include <cblas.h>
using namespace std;
void matrix_make(float *h,int n){
	int n2=n*n; float *p=h;
	memset(h,0,n2*n2*sizeof(float));
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
void vector_print(float *h,int n){
	for(int i=0;i<n;i++){
		printf("%.5f ",*h); h++;
	} printf("\n\n");
	return ;
}
void matrix_print(float *h,int n){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			printf("%.4f ",*h);
			h++;
		}
		printf("\n");
	} printf("\n");
	return ;
}
void cholesky_llt(float *m,float *x,float *c,int n){
    float tmp,*y=(float *)malloc(n*sizeof(float));
	for(int i=0;i<n;i++){
        for(int j=0;j<i;j++)
            *(m+i*n+j)-=cblas_sdot(j,m+i*n,1,m+j*n,1);
        for(int j=0;j<i;j++){
            tmp=*(m+i*n+j); *(m+i*n+j)/=*(m+j*n+j);
            *(m+i*n+i)-=tmp*(*(m+i*n+j));
        }
        for(int j=i+1;j<n;j++) *(m+i*n+j)=0;
	} //Part I Ax=b => L : LL'x=b
	for(int i=0;i<n;i++)
        *(y+i)=*(c+i)-cblas_sdot(i,y,1,m+i*n,1);
    for(int i=0;i<n;i++) *(y+i)/=*(m+i*n+i);
	for(int i=n-1;i>-1;i--)
        *(x+i)=*(y+i)-cblas_sdot(n-1-i,x+i+1,1,m+(i+1)*n+i,n);
	//Part II Ly=b => y ; Dy*=y ; L'x=y* => x
	free(y); return ;
}
int auto_main(int n) {
	int n2=n*n,n4=n2*n2; //n should be at least 3!
	float *m=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	vector_make(c,n);
	clock_t s,e;
	s=clock();
	cholesky_llt(m,x,c,n2);
	e=clock();
	for(int i=0;i<n2;i++) *(x+i)-=1;
	printf("Error: %.6f\n",cblas_snrm2(n2,x,1));
	printf("Time(ms): %.3f\n",float(e-s)/CLOCKS_PER_SEC*1000);
	free(m); free(c); free(x);
	return 0;
}
int main(){
    for(int i=5;i<=120;i+=5){
        printf("Size: %d\n",i);
        auto_main(i); printf("\n");
    }
    return 0;
}
