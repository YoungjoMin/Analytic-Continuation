# include "ANALYTIC.h"

namespace iRRAM{

ANALYTIC::ANALYTIC(COMPLEX(*seq)(int), int k) : coef(seq), k(k) {}

COMPLEX ANALYTIC::evalHelper(int p, int d, const COMPLEX& z) {
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

		diffTerm*=(i+d+1);
		diffTerm/=(i+1);

		pow =pow * z;
	}
	return result;
}

/*
	 For |z| <= 1/2k,  where k given in the Constructor
*/
COMPLEX ANALYTIC::eval(int d, const COMPLEX& z) {
	static ANALYTIC * Tfp;
	static int Td;
	Tfp = this;
	Td = d;
	COMPLEX (*lambda)(int, const COMPLEX& z)  = ([] (int p, const COMPLEX& z) {
			return Tfp->evalHelper(p,Td, z);
	});

	return limit(lambda,z);
}


} //namespace iRRAM
