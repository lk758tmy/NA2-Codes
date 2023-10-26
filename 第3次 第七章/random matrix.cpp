#include <cstdio>
#include <random>
#include <time.h>
#include <cblas.h>
#include <cmath>
#include <algorithm>
using namespace std;
void matrix_print(float *h,int n){
	for(int i=0;i<n;i++){
		//for(int j=0;j<n;j++){
			printf("%.2f\t",*h);
			h++;
		//}
		//printf("\n");
	} printf("\n");
	return ;
}
void rand_matrix(float *m,int n){
	//std::mt19937 gen(time(0));
	//std::uniform_int_distribution<int> dis_i(0, n-1);
	//std::uniform_real_distribution<float> dis_r(-pow(n,1-1.0/n),1);
	srand(time(0));
	for(int i=0;i<n*n;i++) *(m+i)=0;
	for(int i=0;i<n*n;i+=(n+1)) *(m+i)=10.0/sqrt(n);
	//int a,b,nmax=n*ceil(log(n));
	float alpha; int nmax=sqrt(n);
	float **cards=(float **)malloc(n*sizeof(float *)),**tmpp,*tmpi;
	for(int i=0;i<nmax;i++){
		for(int j=0;j<n;j++) *(cards+j)=m+n*j;
		for(int j=n-1;j>0;j--){
			tmpp=cards+rand()%j;
			tmpi=*(cards+j); *(cards+j)=*tmpp; *tmpp=tmpi;
		}
		for(int j=0;j<n-1;j+=2){
			alpha=((rand()%2==1)?1:-1);
			alpha*=(rand()%100+1)/100.0; //???????还是太大！
			cblas_saxpy(n,alpha,*(cards+j),1,*(cards+j+1),1);			
		}
	}
	/*for(int i=0;i<n*nmax;i++){
		b=a=dis_i(gen);
		while(b==a) b=dis_i(gen);
		alpha=pow(2,dis_r(gen));
		if(rand()%2==1) alpha=-alpha;
		if(rand()%2==1) cblas_saxpy(n,alpha,m+a*n,1,m+b*n,1);
		else cblas_saxpy(n,alpha,m+a,n,m+b,n);
	}*/
	/*for(int i=0;i<n;i++){
		b=a=dis_i(gen);
		while(b==a) b=dis_i(gen);
		if(rand()%2==1) cblas_sswap(n,m+a*n,1,m+b*n,1);
		else cblas_sswap(n,m+a,n,m+b,n);
	}*/
	return ;
}
int main(){
	//Suitable for n range from 5000 to 10000
	int n; scanf("%d",&n);
	float *m=(float *)malloc(n*n*sizeof(float));
	rand_matrix(m,n);
	//matrix_print(m,n);
	return 0;
}
