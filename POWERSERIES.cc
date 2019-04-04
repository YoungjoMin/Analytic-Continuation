# include "POWERSERIES.h"

namespace iRRAM{

POWERSERIES::POWERSERIES(COEF coef, int k) : coef(coef), k(k) {}
POWERSERIES::POWERSERIES(COMPLEX(*coef)(int), int k) : coef(coef), k(k) {}

COMPLEX POWERSERIES::evalHelper(int p, const COMPLEX& z, int d) {
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
COMPLEX POWERSERIES::eval(const COMPLEX& z, int d) {
	static POWERSERIES * Tfp;
	static int Td;
	Tfp = this;
	Td = d;
	COMPLEX (*lambda)(int, const COMPLEX& z)  = ([] (int p, const COMPLEX& z) {
			return Tfp->evalHelper(p,z, Td);
	});

	return limit(lambda,z);
}


} //namespace iRRAM
