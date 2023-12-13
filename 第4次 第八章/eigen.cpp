#include <cstdio>
#include <cmath>
using namespace std;
int main(){
	freopen("eigen.csv","w",stdout);
	for(int n=100;n<102;n++){
		printf("%d\n",n);
		for(int i=0;i<n;i++)
			printf("%d,%.15f\n",i+1,2-2*cos((i+1)*M_PI/(n+1)));
		printf("\n");
	}
	printf("END");
	return 0;
}
