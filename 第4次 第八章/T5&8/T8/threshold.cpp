#include <cstdio>
#include <cmath>
#include <cblas.h>
using namespace std;
const double epsilon=1e-8;
double sgn(double x){
	if(x>=0) return 1;
	else return -1;
}
int T_Jacobi(double *m,int n,double delta){
	int cnt=0,cntCycle=0; double c,s,t,kse,th; bool flag;
	th=sqrt(2)*sqrt(m[7]*m[7]+m[6]*m[6]+m[5]*m[5]);
	th/=delta;
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
		//printf("%d %.3e\n",cntCycle,th);
		if(!flag){
			if(th<epsilon) break;
			else th/=delta;
		}
	}while(1);
	return cnt;
}
int main(){
	double m[9]={1e40,1e29,1e19,1e29,1e20,1e9,1e19,1e9,1.0};
	T_Jacobi(m,3,4.0);
	printf("%.5e %.5e %.5e\n",m[0],m[4],m[8]);
	free(m);
	return 0;
}
