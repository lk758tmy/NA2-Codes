#include <cstdio>
#include <cmath>
using namespace std;
const double epsilon=1e-6;
double func(double x){
	return x*(-1+x*x/6)+sin(x); //Roots:5(5th)
}
double func_d(double x){
	return -1+x*x/2+cos(x);
}
int newton_SA(double &x,double(*f)(double),double(*d)(double),double real){
	int cnt=0; double y=x+1,t=f(x),s=d(x),z; int m;
	printf("%d,%.10f,%.10f\n",cnt,x,t);
	while(abs(t)>epsilon || abs(y-x)>epsilon){
		m=round(abs(1.0/(1-log(abs(s))/log(abs(t)))));
		cnt++; z=y; y=x; x-=m*t/s; t=f(x); s=d(x);
		printf("%d,%d,%.10f,%.10f",cnt,m,x,t);
		if(cnt==1) printf("\n");
		else printf(",%.10f\n",log(abs((x-real)/(y-real)))/log(abs((y-real)/(z-real))));
	}
	return cnt;
}
int main(){
	double x=10;
	newton_SA(x,func,func_d,0); printf("\n");
	return 0;
}

