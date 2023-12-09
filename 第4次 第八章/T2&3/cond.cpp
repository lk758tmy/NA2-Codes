#include <cstdio>
#include <cmath>
#include <lapacke.h>
using namespace std;
void matrix_make(double *m,int n,double q){
	for(int i=0;i<n-1;i++){
		*(m+i*n+i)=2-q; *(m+i*n+i+1)=*(m+i*n+i+n)=-1;
		for(int j=2+i;j<n;j++)
			*(m+i*n+j)=*(m+j*n+i)=0;
	}
	*(m+n*n-1)=2-q;
	
	double rcond,anorm;
	double *w1=(double *)malloc(4*n*sizeof(double));
	int *w2=(int *)malloc(n*n*sizeof(int)),lda=n,info;
	
	anorm=LAPACK_dlange("1",&n,&n,m,&lda,w1);
	LAPACK_dgetrf(&n,&n,m,&lda,w2,&info);
	LAPACK_dgecon("1",&n,m,&lda,&anorm,&rcond,w1,w2,&info);
	
	free(w1); free(w2);
	printf("Cond:%f\n",1.0/rcond);
	
	return ;
}
int main(){
	double *a=(double *)malloc(101*101*sizeof(double));
	matrix_make(a,100,2.001);
	matrix_make(a,100,3.001);
	matrix_make(a,101,2.001);
	matrix_make(a,101,3.001);
	return 0;
}
