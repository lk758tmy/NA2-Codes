#include <cstdio>
#include <cmath>
using namespace std;
const double epsilon=1e-6;
double func(double x){
	return 12+x*(-8+x*(-1+x)); //Roots:-3(1st),2(2nd)
}
double func_d(double x){
	return -8+x*(-2+3*x);
}
int newton(double &x,double(*f)(double),double(*d)(double),double real){
	int cnt=0; double y=x+1,t=f(x),z;
	printf("%d,%.10f,%.10f\n",cnt,x,t);
	while(abs(t)>epsilon || abs(y-x)>epsilon){
		cnt++; z=y; y=x; x-=t/d(x); t=f(x);
		printf("%d,%.10f,%.10f",cnt,x,t);
		if(cnt==1) printf("\n");
		else printf(",%.10f\n",log(abs((x-real)/(y-real)))/log(abs((y-real)/(z-real))));
	}
	return cnt; 
}
int main(){
	double x=5;
	newton(x,func,func_d,2); printf("\n");
	x=-5;
	newton(x,func,func_d,-3);
	return 0;
}

