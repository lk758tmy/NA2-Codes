#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <chrono>
#include <algorithm>
using namespace std;
using namespace std::chrono;
const double epsilon=1e-8;
void matrix_make(double *m,int n){
	for(int i=0;i<n-1;i++){
		*(m+i*n+i)=2; *(m+i*n+i+1)=*(m+i*n+i+n)=-1;
		for(int j=2+i;j<n;j++)
			*(m+i*n+j)=*(m+j*n+i)=0;
	}
	*(m+n*n-1)=2;
	return ;
}
double sgn(double x){
	if(x>=0) return 1;
	else return -1;
}
int Cl_Jacobi(double *m,int n){
	int maxI,maxJ,cnt=0; double maxValue,c,s,t,kse;//,tmp;
	do{
		maxValue=0; //maxI=-1; maxJ=-1;
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){ //symmetric
				//if(i==j) continue;
				if(abs(*(m+i*n+j))>maxValue){
					maxValue=abs(*(m+i*n+j));
					maxI=i; maxJ=j;
				}
			}
		}
		if(maxValue<epsilon) break;
		kse=(*(m+n*maxI+maxI)-*(m+maxJ*n+maxJ))/(*(m+maxI*n+maxJ))/2;
		if(abs(kse)<11.4){ //tan(85)
			t=sgn(kse)/(abs(kse)+sqrt(1+kse*kse));
			c=1/sqrt(1+t*t); s=c*t;
		}else{
			t=0.25/kse;
			c=1/(1+t*t); s=c*t*2; c*=(1-t*t); //t=s/c;
		}
		//tmp=t*(*(m+n*maxI+maxJ));
		cblas_drot(n,m+n*maxI,1,m+n*maxJ,1,c,s);
		cblas_drot(n,m+maxI,n,m+maxJ,n,c,s);
		//cblas_dcopy(n,m+n*maxI,1,m+maxI,n);
		//cblas_dcopy(n,m+n*maxJ,1,m+maxJ,n);
		//*(m+n*maxI+maxI)+=tmp; *(m+n*maxJ+maxJ)-=tmp;
		//It doesn't work and I don't know why.
		cnt++;
	}while(1);
	printf("%.0f,",ceil(double(cnt)/(n*(n-1)/2)));
	return cnt;
}
int Cy_Jacobi(double *m,int n){
	int cnt=0,cntCycle=0; double c,s,t,kse; bool flag;
	do{
		flag=false;
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){
				if(abs(*(m+n*i+j))<epsilon) continue;
				flag=true;
				kse=(*(m+n*i+i)-*(m+j*n+j))/(*(m+i*n+j))/2;
				if(abs(kse)<11.4){
					t=sgn(kse)/(abs(kse)+sqrt(1+kse*kse));
					c=1/sqrt(1+t*t); s=c*t;
				}else{
					t=0.25/kse;
					c=1/(1+t*t); s=c*t*2; c*=(1-t*t);
				};
				cblas_drot(n,m+n*i,1,m+n*j,1,c,s);
				cblas_drot(n,m+i,n,m+j,n,c,s);
				cnt++;
			}
		}
		cntCycle++;
	}while(flag);
	printf("%d,",cntCycle);
	return cnt;
}
int T_Jacobi(double *m,int n,double delta){
	int cnt=0,cntCycle=0; double c,s,t,kse,th; bool flag;
	th=cblas_dnrm2(n*n,m,1); th-=cblas_dnrm2(n,m,n+1); th/=delta;
	do{
		flag=false;
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){
				if(abs(*(m+n*i+j))<th) continue;
				flag=true;
				kse=(*(m+n*i+i)-*(m+j*n+j))/(*(m+i*n+j))/2;
				if(abs(kse)<11.4){
					t=sgn(kse)/(abs(kse)+sqrt(1+kse*kse));
					c=1/sqrt(1+t*t); s=c*t;
				}else{
					t=0.25/kse;
					c=1/(1+t*t); s=c*t*2; c*=(1-t*t);
				};
				cblas_drot(n,m+n*i,1,m+n*j,1,c,s);
				cblas_drot(n,m+i,n,m+j,n,c,s);
				cnt++;
			}
		}
		cntCycle++;
		if(!flag){
			if(th<epsilon) break;
			else th/=delta;
		}
	}while(1);
	printf("%d,",cntCycle);
	return cnt;
}
int main(){
	double *m,*eigenNum,*eigenReal;
	time_point<steady_clock> s,e;
	//freopen("\data.csv","w",stdout);
	for(int n=100;n<102;n++){
	for(int q=0;q<3;q++){
		m=(double *)malloc(n*n*sizeof(double));
		eigenNum=(double *)malloc(n*sizeof(double));
		eigenReal=(double *)malloc(n*sizeof(double));
		for(int i=0;i<n;i++)
			*(eigenReal+i)=2-2*cos((i+1)*M_PI/(n+1));
		
		matrix_make(m,n);
		s=steady_clock::now();
		switch(q){
		case 0:
			printf("Classic,%d,%d",n,Cl_Jacobi(m,n)); break;
		case 1:
			printf("Cyclic,%d,%d",n,Cy_Jacobi(m,n)); break;
		case 2: 
			printf("Threshold,%d,%d",n,T_Jacobi(m,n,double(n+1)));
		}
		e=steady_clock::now();
		printf(",%.3f",(e-s).count()/1000000.0);
		cblas_dcopy(n,m,n+1,eigenNum,1);
		sort(eigenNum,eigenNum+n);
		cblas_daxpy(n,-1,eigenReal,1,eigenNum,1);
		printf(",%.3e\n",cblas_dnrm2(n,eigenNum,1));
		
		//matrix_make(m,n); power(m,n,v0);
		//matrix_make(m,n); power(m,n,v0);
		free(m); free(eigenNum); free(eigenReal);
    }}
	return 0;
}
