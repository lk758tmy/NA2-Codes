#define _USE_MATH_DEFINES
#include <cstdio>
#include <utility>
#include <cmath>
#include <cblas.h>
using namespace std;
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
void precond_make(float *h,int n,float omega){
	int n2=n*n; float *p=h; omega/=2;
	for(int i=0;i<n2*n2;i++){*p=0; p++;}
	p=h;
	for(int i=0;i<n2;i++){
		*p=2;
		if(i%n!=0) *(p-1)=-omega;
		p+=(n2+1);
	}
	p=h+n*n2;
	for(int i=n;i<n2;i++){
		*p=-omega; p+=(n2+1);
	}
	return ;
}
void Cholskey(float *m,float *x,float *c,int n){
	float *y=(float *)malloc(n*sizeof(float));
	for(int i=0;i<n;i++){
		*(y+i)=*(c+i)-cblas_sdot(i,y,1,m+i*n,1);
		*(y+i)/=*(m+i*n+i);
	}
	for(int i=n-1;i>-1;i--){
		*(x+i)=*(y+i)-cblas_sdot(n-1-i,x+i+1,1,m+(i+1)*n+i,n);
		*(x+i)/=*(m+i*n+i);
	}
	free(y); return ;
}
int PCG(float *m,float *x,float *c,int n,float *r,float epsilon,float *q){
	float alpha,beta,*p=(float *)malloc(n*sizeof(float));
	float tmp_rdotz,*tmp_Amp=(float *)malloc(n*sizeof(float));
	float *z=(float *)malloc(n*sizeof(float));
	int cnt;
	
	for(int i=0;i<n;i++) *(r+i)=*(c+i);
	cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,x,1,-1,r,1);
	Cholskey(q,z,r,n);
	for(int i=0;i<n;i++) *(p+i)=-*(z+i);
	cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
	tmp_rdotz=cblas_sdot(n,r,1,z,1);
	alpha=tmp_rdotz/cblas_sdot(n,p,1,tmp_Amp,1);
	
	for(cnt=1;cnt<=n;cnt++){
		cblas_saxpy(n,alpha,p,1,x,1);
		cblas_saxpy(n,alpha,tmp_Amp,1,r,1);
		Cholskey(q,z,r,n);
		/*for(int j=0;j<n;j++) printf("%.2f ",*(x+j)); printf("\n");
		for(int j=0;j<n;j++) printf("%.2f ",*(z+j)); printf("\n");
		  for(int j=0;j<n;j++) printf("%.2f ",*(r+j)); printf("\n\n");*/
		beta=cblas_sdot(n,r,1,z,1)/tmp_rdotz;
		tmp_rdotz*=beta;
		if(cblas_snrm2(n,r,1)<epsilon) goto CG_END;
		cblas_saxpy(n,-1.0/beta,z,1,p,1);
		for(int j=0;j<n;j++) *(p+j)*=beta;
		cblas_sgemv(CblasRowMajor,CblasNoTrans,n,n,1,m,n,p,1,0,tmp_Amp,1);
		alpha=tmp_rdotz/cblas_sdot(n,p,1,tmp_Amp,1);
	}
	
	CG_END: free(p); free(tmp_Amp); free(z);
	return cnt;
}
int _main(int n,float omega){
	int n2=n*n,n4=n2*n2;
	float *m=(float *)malloc(n4*sizeof(float));
	float *q=(float *)malloc(n4*sizeof(float));
	matrix_make(m,n);
	precond_make(q,n,omega);
	
	float *c=(float *)malloc(n2*sizeof(float));
	float *x=(float *)malloc(n2*sizeof(float));
	float *r=(float *)malloc(n2*sizeof(float));
	for(int i=0;i<n;i++) *(x+i)=0;
	vector_make(c,n);
	
	printf("%d\t%.2f\t%d\n",n,omega,PCG(m,x,c,n2,r,1e-6,q));
	free(m); free(c); free(x); free(r); free(q);
    return 0;
}
int main(){
	int size;
	//scanf("%d",&size);
	for(size=3;size<=60;size+=3)
		for(float omega=1.01;omega<1.995;omega+=0.01) //0<omega<1不用尝试了
			_main(size,omega);
	return 0;
}
