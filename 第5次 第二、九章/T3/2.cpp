#include <cstdio>
#include <cmath>
using namespace std;
const double epsilon=1e-6;
double func(double x){
	return x*(-1+x*x/6)+sin(x); //Roots:5(5th)
}
int secant(double &x,double x2,double(*f)(double),double real){
	int cnt=0; double y=x2,s=f(x2),t=f(x),z;
	printf("%d,%.10f,%.10f\n",cnt,x,t);
	while(abs(t)>epsilon || abs(y-x)>epsilon){
		cnt++; z=y; y=x; x-=(y-z)*t/(t-s);
		s=t; t=f(x);
		printf("%d,%.10f,%.10f",cnt,x,t);
		if(cnt==1) printf("\n");
		else printf(",%.10f\n",log(abs((x-real)/(y-real)))/log(abs((y-real)/(z-real))));
	}
	return cnt; 
}
int main(){
	double x=5;
	secant(x,x+0.1,func,5); printf("\n");
	return 0;
}

