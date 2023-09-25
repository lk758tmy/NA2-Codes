#include <cstdio>
#include <memory.h>
#include <lapack.h>
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
void _main(int n) {
	float rcond,anorm; int n2=n*n;
    //n should be at least 3!
	float *m=(float *)malloc(n2*n2*sizeof(float));
	float *w1=(float *)malloc(4*n2*sizeof(float));
	int *w2=(int *)malloc(n2*sizeof(int)),lda=n2,info;
	matrix_make(m,n);

    anorm=LAPACK_slange("1",&n2,&n2,m,&lda,w1);
	LAPACK_sgetrf(&n2,&n2,m,&lda,w2,&info);
    LAPACK_sgecon("1",&n2,m,&lda,&anorm,&rcond,w1,w2,&info);

    free(m); free(w1); free(w2);
    printf("Size:%d\tCond:%f\n",n,1.0/rcond);
	return;
}
int main(){
    for(int n=5;n<=120;n+=5) _main(n);
    return 0;
}
