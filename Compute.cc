# include "iRRAM.h"
# include "./Poly.h"
# include "./ANALYTIC.h"
# include <functional>

using namespace iRRAM;


COMPLEX seq(int k) {
	REAL a(1);
	if(k<=1)
		return COMPLEX(1);

	for(int i=2;i<=k;i++)
			a/=i;
	return COMPLEX(a);
}

COMPLEX evalHelper(int p, COMPLEX(*coef)(int), const COMPLEX& z, int k, int d) {
	COMPLEX result = coef(0);
	COMPLEX pow = z;
	COMPLEX cur;
	int t = k-p + 32*d;

	for(int i = 1; i<=t;i++) {
		cur = coef(i)*pow;
		for(int j=d+i;j>i;j--) {
			cur = cur*COMPLEX(j);
		}
		result = result + cur;
		pow =pow * z;
	}
	return result;
}

COMPLEX eval(COMPLEX(*coef)(int), const COMPLEX& z, int k, int d) {
	static COMPLEX (*Tcoef)(int) = coef;
	static int Tk = k;
	static int Td = d;
	
	COMPLEX(*lambda)(int, COMPLEX)  =  ([] (int p, COMPLEX z) -> COMPLEX {
		return evalHelper(p,Tcoef,z,Tk,Td); });
	
	return limit( (COMPLEX (*)(int, COMPLEX))lambda , z);
}

void compute()
{
		cout<<setRwidth(100);
		COMPLEX result = eval(seq,COMPLEX(REAL(0),pi()),1,0);
		cout<<result._real<<" + "<<result._imag<<" i\n";
		return;
}
