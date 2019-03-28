# include "iRRAM.h"
# include "./ANALYTIC.h"

using namespace iRRAM;

COMPLEX seq(int k) {
	static REAL prevSeq;
	static int prevK;
	REAL a(1);
	int i=2;

	if(k<=1)
		return COMPLEX(1);

	if(prevK<k) {
		a = prevSeq;
		i = prevK+1;
	}
	for(;i<=k;i++)
			a/=i;

	prevSeq = a;
	prevK = k;
	return COMPLEX(a);
}
/*
COMPLEX seq(int k) {
	REAL a(1);
	if(k<=1)
		return COMPLEX(1);
	for(int i=2;i<=k;i++)
		a/=i;
	return COMPLEX(a);
}
*/

COMPLEX evalHelper(int p, COMPLEX (*coef)(int), const COMPLEX& z, int k, int d) {
	COMPLEX result(0);
	COMPLEX pow(1);
	COMPLEX cur;
	INTEGER diffTerm(1);
	int t = k-p + 32*d;
	
	for(int i = 2;i<=d;i++)
		diffTerm*=i;

	for(int i = 0; i<=t;i++) {
		cur = (coef(i+d)*pow)*diffTerm;
		result = result + cur;

		diffTerm/=(i+1);
		diffTerm*=(i+d+1);
		pow =pow * z;
	}
	return result;
}

COMPLEX eval(COMPLEX(*coef)(int), const COMPLEX& z, int k, int d) {
	
	static COMPLEX (*Tcoef)(int) = coef;
	static int Tk = k;
	static int Td = d;
	
	COMPLEX(*lambda)(int, const COMPLEX&)  =  ([] (int p, const COMPLEX& z) -> COMPLEX {
			return evalHelper(p,Tcoef,z,Tk,Td); });

	return limit( lambda , z);
}

void compute()
{
	COMPLEX result;
	cout<<setRwidth(50);
	result = eval(seq,COMPLEX(REAL(0),pi()),1,1000);
	cout<<result._real<<" + "<<result._imag<<" i\n";
	return;
}
