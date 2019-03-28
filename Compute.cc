# include "iRRAM.h"
# include <functional>

using namespace iRRAM;

COMPLEX (*Tcoef)(int);
int Tk = 1;
int Td = 1;

COMPLEX seq(int k) {
	REAL a(1);
	if(k<=1)
		return COMPLEX(1);

	for(int i=2;i<=k;i++)
			a/=i;
	return COMPLEX(a);
}

COMPLEX evalHelper(int p, COMPLEX (*coef)(int), const COMPLEX& z, int k, int d) {

	COMPLEX result = coef(0);
	COMPLEX pow = z;
	COMPLEX cur;
	INTEGER a;
	int t = k-p + 32*d;
	for(int i = 1; i<=t;i++) {
		a=INTEGER(1);
		cur = coef(i)*pow;
		for(int j=d+i;j>i;j--) {
			a=a*INTEGER(j);
			cur = cur*COMPLEX(j);
		}
		cur=cur*COMPLEX(a);
		result = result + cur;
		pow =pow * z;
	}
	return result;
}





COMPLEX eval(int p,const COMPLEX& z) {

		
		return evalHelper(p,Tcoef,z,Tk,Td); 
}

void compute()
{
	Tcoef= seq;
	Tk = 1;
	Td = 0;
	COMPLEX z=pi();
	cout<<setRwidth(100);
	COMPLEX result = limit(eval,z);
	cout<<result._real<<" + "<<result._imag<<" i\n";
	return;
}
