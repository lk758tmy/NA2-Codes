#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
#include <lapacke.h>
#include <stack>
using namespace std;
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
double maxbar(double *v,int n){
	double ans=0;
	for(int i=0;i<n;i++)
		ans=(abs(*(v+i))>abs(ans))?*(v+i):ans;
	return ans;
}
void gauss_colPivot(double *m,double *x,double *c,int n){
	stack<pair<int,int> > s;
	int maxRowId,rowId;
	double *mp=m,*cp=c,maxElement,tmp;
	for(int i=0;i<n;i++){
		rowId=i;
		mp=m+n*rowId+i; cp=c+rowId;
		maxElement=abs(*mp); maxRowId=rowId;
		for(int j=i+1;j<n;j++){
			mp+=n;
			if(abs(*mp)>maxElement){
				maxElement=abs(*mp); maxRowId=j;
			}
		}
		if(rowId!=maxRowId){
			cblas_dswap(n,m+n*rowId,1,m+n*maxRowId,1);
			s.push(make_pair(rowId,maxRowId));
			tmp=*cp; *cp=*(c+maxRowId); *(c+maxRowId)=tmp;	
		}
		mp=m+n*rowId+i;
		tmp=1.0/(*mp); (*cp)*=tmp;
		for(int j=i;j<n;j++){
			(*mp)*=tmp; *mp++;
		}
		mp=m+n*rowId+i; cp=c+rowId; tmp=*cp;
		for(int j=i+1;j<n;j++){
			cp++; mp+=n;
			*cp-=tmp*(*mp);
			cblas_daxpy(n,-*mp,m+n*rowId,1,m+n*j,1);
		}
	}
	double *xp=x+n-1,*xp2;
	*xp=*(c+n-1);
	for(int i=n-2;i>=0;i--){
		xp--; *xp=*(c+i);
		mp=m+n*(i+1)-1; xp2=x+n-1;
		for(int j=n-1;j>i;j--){
			*xp-=(*xp2)*(*mp);
			mp--; xp2--;
		}
	}
	while(!s.empty()){
		tmp=*(x+s.top().first);
		*(x+s.top().first)=*(x+s.top().second);
		*(x+s.top().second)=tmp;
		s.pop();
	}
	return ;
}
double inverse_power(double *a,int n,double &l,double *u,double *v0,double q){
	int cnt=0,*ipiv=(int *)malloc(n*sizeof(int));
	double *m=(double *)malloc(n*n*sizeof(double));
	double *v1=(double *)malloc(n*sizeof(double));
	double *v2=(double *)malloc(n*sizeof(double)),e,m0,m1=1e9;
	cblas_dcopy(n,v0,1,v1,1);
	for(int i=0;i<n;i++) *(a+i*n+i)-=q;
	
	e=cblas_ddot(n,v1,1,u,1)/cblas_dnrm2(n,v1,1);
	printf("%d\n%d,%.10e,%.10e,%.10e\n",n,cnt,abs(l-q),q,sqrt(1-e*e));
	do{
		cnt++; if(cnt==200) break;
		gauss_colPivot(a,v2,v1,n);
		//cblas_dcopy(n,v1,1,v2,1); cblas_dcopy(n*n,a,1,m,1);
		//LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,1,m,n,ipiv,v2,1);
		m0=m1; m1=maxbar(v2,n);
		//for(int i=0;i<n;i++) *(v1+ipiv[i]-1)=*(v2+i)/m1;
		for(int i=0;i<n;i++) *(v1+i)=*(v2+i)/m1;
		
		e=1/cblas_dnrm2(n,v1,1);
		e*=cblas_ddot(n,v1,1,u,1);
		printf("%d,%.10e,%.10e,%.10e\n",cnt,abs(l-q-1/m1),1/m1+q,sqrt(1-e*e));	
	}while(abs(1/m1-1/m0)>epsilon);
	free(v1); free(v2); free(ipiv); return q+1/m1;
}
int main(){
	double *a,*v0,*u,l;
	//freopen("\data.csv","w",stdout);
	for(int n=101;n<102;n++){
	for(double q=1.9999999925;q<2;q++){
		a=(double *)malloc(n*n*sizeof(double));
		v0=(double *)malloc(n*sizeof(double));
		u=(double *)malloc(n*sizeof(double));
		
		if(n==100){
				for(int i=0;i<n;i++) *(u+i)=sin((i+1)*50*M_PI/101);
				l=1/cblas_dnrm2(n,u,1);
				for(int i=0;i<n;i++) *(u+i)*=l;
				l=2-2*cos(50*M_PI/101);
		}if(n==101){
				for(int i=0;i<n;i++) *(u+i)=sin((i+1)*M_PI/2); // 51/102
				l=1/cblas_dnrm2(n,u,1);
				for(int i=0;i<n;i++) *(u+i)*=l;
				l=2;
		}
		srand(time(0));
		for(int i=0;i<n;i++) *(v0+i)=double(rand())/RAND_MAX;
		matrix_make(a,n); 
		if(n==100) inverse_power(a,n,l,u,v0,q);
		if(n==101) inverse_power(a,n,l,u,v0,q);
		free(u); free(a);		
	}}

	return 0;
}
