# include "POWERSERIES.h"

namespace iRRAM{

POWERSERIES::POWERSERIES() {
	coef = ([] (int j) -> COMPLEX { return COMPLEX(0);});
	k=1;
}
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


POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2) {
	COEF lambda = ([=] (int j) -> COMPLEX {
			return f1.coef(j)+f2.coef(j);
			});
	return POWERSERIES(lambda, std::max(f1.k,f2.k)+1);
}
POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2) {
	COEF lambda = ([=] (int j) -> COMPLEX {
			return f1.coef(j)-f2.coef(j);
			});
	return POWERSERIES(lambda, std::max(f1.k,f2.k)+1);
}
POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2) {
COEF lambda = ([=] (int j) -> COMPLEX {
					COMPLEX result(0);
					for(int i = 0;i<=j;i++)
						result= result + (f1.coef(i)*f2.coef(j-i));
					return result;
		}	);
	return POWERSERIES(lambda, f1.k+f2.k);
}

POWERSERIES& POWERSERIES::operator +=(const POWERSERIES& f) {
	*this = (*this)+f;
	return *this;
}
POWERSERIES& POWERSERIES::operator -=(const POWERSERIES& f) {
	*this = (*this)-f;
	return *this;
}
POWERSERIES& POWERSERIES::operator *=(const POWERSERIES& f) {
	*this = (*this)*f;
	return *this;
}


} //namespace iRRAM
