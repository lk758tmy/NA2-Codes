#include <cstdio>
#include <cmath>
using namespace std;
float A[4],B[4],alpha,beta;
void multi(){
	A[1]=4*A[0]+A[1]*alpha;
	A[3]=4*A[2]+A[3]*alpha;
	A[0]*=alpha; A[2]*=alpha;
	B[0]*=alpha; B[3]*=beta;
	//B^m=((alpha^m,0),(0,beta^m));
	return ;
}
float norm2_A(){
	return sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+A[3]*A[3]);
}
float norm2_B(){
	return sqrt(B[0]*B[0]+B[3]*B[3]);
}
int main() {
	int N;
	scanf("%f%f%d",&alpha,&beta,&N);
	//¦Ñ(A)=alpha ¦Ñ(B)=max(alpha,beta)
	//Requires that alpha slightly smaller than beta
	A[0]=A[3]=alpha; A[1]=4; A[2]=0;
	B[0]=alpha; B[3]=beta; B[1]=B[2]=0;
	printf("%d\t%f\t%f\n",0,norm2_A(),norm2_B());
	for(int m=1;m<=N;m++){
		multi(); 
		if((m%25)!=0) continue;
		printf("%d\t%f\t%f\n",m,norm2_A(),norm2_B());
	}
	return 0;
}
