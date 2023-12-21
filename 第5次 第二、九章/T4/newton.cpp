#include <cstdio>
#include <cmath>
using namespace std;
const double epsilon=1e-6;
void func(double x1,double x2,double &y1,double &y2){
	y1=(x1+3)*(x2*x2-7)+18; y2=sin(x2*exp(x1)-1);
	return ;
}
void func_d(double x1,double x2,double* m){
	m[0]=x2*x2-7; m[1]=2*(x1+3)*x2; m[3]=exp(x1);
	m[3]*=cos(x2*m[3]-1); m[2]=m[3]*x2; //略微减少计算量
	return ;
}
int newton(double &x1,double &x2,void(*f)(double,double,double&,double&),
		void(*d)(double,double,double*),double r1,double r2){
	int cnt=0; double y1=x1+1,y2=x2+1,t1,t2,E=epsilon*epsilon;
	f(x1,x2,t1,t2); double z1,z2,m[4],tmpf,tmpx,tmpy,tmpz;
	printf("%d,%.10f,%.10f,%.10f,%.10f\n",cnt,x1,x2,t1,t2);
	while(t1*t1+t2*t2>E || (y1-x1)*(y1-x1)+(y2-x2)*(y2-x2)>E){
		cnt++; z1=y1; z2=y2; y1=x1; y2=x2;
		d(x1,x2,m); tmpf=1/(m[0]*m[3]-m[1]*m[2]);
		x1-=tmpf*(m[3]*t1-m[1]*t2); x2-=tmpf*(-m[2]*t1+m[0]*t2);
		f(x1,x2,t1,t2);
		printf("%d,%.10f,%.10f,%.10f,%.10f",cnt,x1,x2,t1,t2);
		if(cnt==1) printf("\n");
		else{
			tmpz=tmpy; tmpy=tmpx;
			tmpx=(x1-r1)*(x1-r1)+(x2-r2)*(x2-r2);
			printf(",%.10f\n",log(tmpx/tmpy)/log(tmpy/tmpz));
		}
	}
	return cnt; 
}
int main(){
	double x1=-0.15,x2=0.14;
	newton(x1,x2,func,func_d,0,1); printf("\n");
	return 0;
}

