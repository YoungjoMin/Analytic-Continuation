# include "ANALYTIC.h"

namespace iRRAM{

ANALYTIC::ANALYTIC(SEQUENCE seq, int k) {
	this->coef = seq;
	this->k = k;
}

REAL evalHelper(int p, const ANALYTIC& f, const REAL& z) {
	REAL result = f.coef(0);
	REAL pow = z;
	int t = f.k - p;
	for(int i = 1; i<=t ; i++) {
		result = result + (f.coef(i) * pow);
		pow = pow * z;
	}
	return result;
}

REAL ANALYTIC::operator() (const REAL& z) {
	return limit(evalHelper, *this, z);
}

} //namespace iRRAM
