# include "ANALYTIC.h"

namespace iRRAM{

ANALYTIC::ANALYTIC(REAL(*seq)(int), int k) {
	this->coef = seq;
	this->k = k;
	/*
	this->eval = 
		[=](int p , const REAL& r ) ->REAL {
	REAL result = coef(0);
	REAL pow = r;
	int t = k - p;
	for(int i = 1; i<=t ; i++) {
		result = result + (coef(i) * pow);
		pow = pow * r;
	}
	return result;
		};
	*/
}

REAL ANALYTIC::operator() (const REAL& z) {
	//return limit(this->eval, z);
	return REAL(1);
}

} //namespace iRRAM
