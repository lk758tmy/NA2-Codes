#include <cstdio>
#include <cmath>
#include <cblas.h>
using namespace std;
const double epsilon=1e-8;
bool sturm_sameSign(double x,double y){
	if(x*y>0) return true;
	else if(x*y<0) return false;
	else if(x==0) return true;
	else return false;
}
int sturm_sameSignCnt(double *a,double *b,int n,double x){
	double p=1,q=*a-x,r; int cnt=0;
	if(sturm_sameSign(p,q)) cnt++;
	for(int i=1;i<n;i++){
		r=q*(*(a+i)-x)-p*(*(b+i))*(*(b+i));
		if(sturm_sameSign(q,r)) cnt++;
		p=q; q=r;
	}
	return cnt;
}
int sturm_dfs(double *a,double *b,int n,double l,double r){
	double m; int sl,sm,sr,dl,dr;
	do{
		m=(l+r)/2;
		sl=sturm_sameSignCnt(a,b,n,l);
		sm=sturm_sameSignCnt(a,b,n,m);
		sr=sturm_sameSignCnt(a,b,n,r);
		dl=sl-sm; dr=sm-sr;
		if(dl==0&&dr==0) return 0;
		else if(dl!=0&&dr!=0)
			return sturm_dfs(a,b,n,l,m)+sturm_dfs(a,b,n,m,r);
		else if(dl==0) l=m;
		else r=m;
		if(r-l<epsilon){
			printf("%.7f,",m); return dl+dr;
		}
	}while(1); return -1e9-7;
}
int sturm(double *a,double *b,int n,double l,double r){
	return sturm_dfs(a,b,n,l,r);
}
int main(){
	double *a,*b,*eR;
	for(int n=100;n<102;n++){
		a=(double *)malloc(n*sizeof(double *));
		b=(double *)malloc(n*sizeof(double *));
		eR=(double *)malloc(n*sizeof(double *));
		for(int i=0;i<n;i++){
			*(a+i)=2; *(b+i)=-1;
			*(eR+i)=2-2*cos((i+1)*M_PI/(n+1));
		}
		printf("\n%d\n\n",sturm(a,b,n,1,2));
		free(a); free(b); free(eR);
	}
	return 0;
}
