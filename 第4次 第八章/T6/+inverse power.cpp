#include <cstdio>
#include <cmath>
#include <cblas.h>
#include <time.h>
#include <random>
//#include <stack>
#include <lapacke.h>
using namespace std;
const double epsilon=1e-8;
bool sturm_sameSign(double x,double y){
	if(x*y>0) return true;
	else if(x*y<0) return false;
	else if(x==0) return true;
	else return false;
}
/*void gauss_colPivot(double *m,double *x,double *c,int n){
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
}*/
double maxbar(double *v,int n){
	double ans=0;
	for(int i=0;i<n;i++)
		ans=(abs(*(v+i))>abs(ans))?*(v+i):ans;
	return ans;
}
double sturm_inverse_1step(double *a,int n){
	double *v1=(double *)malloc(n*sizeof(double));
	double *v2=(double *)malloc(n*sizeof(double)),m0;
	int *ipiv=(int *)malloc(n*sizeof(int));
	srand(time(0));
	for(int i=0;i<n;i++) *(v1+i)=rand()/double(RAND_MAX);
	
	cblas_dcopy(n,v1,1,v2,1);
	LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,1,a,n,ipiv,v2,1);
	m0=maxbar(v2,n);
	for(int i=0;i<n;i++) *(v1+ipiv[i]-1)=*(v2+i)/m0;
	
	free(v1); free(v2); free(ipiv);
	return 1/m0;
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
			double *mm=(double *)malloc(n*n*sizeof(double));
			for(int i=0;i<n-1;i++){
				*(mm+i*n+i)=2-m; *(mm+i*n+i+1)=*(mm+i*n+i+n)=-1;
			}*(mm+n*n-1)=2-m;
			printf("%.10f,%.10f\n",m,m+sturm_inverse_1step(mm,n));
			free(mm); return dl+dr;
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
