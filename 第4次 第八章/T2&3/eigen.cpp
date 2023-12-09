#include <cstdio>
#include <cmath>
using namespace std;
int main(){
	for(int n=100;n<102;n++){
		for(int i=0;i<n;i++)
			printf("%d %.6f\n",i+1,2-2*cos((i+1)*M_PI/(n+1)));
		printf("\n");
	}
	return 0;
}
