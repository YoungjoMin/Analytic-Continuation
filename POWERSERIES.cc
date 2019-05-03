# include "POWERSERIES.h"

namespace iRRAM{

POWERSERIES::POWERSERIES() {
	center = COMPLEX(0);
	coef = ([] (int j) -> COMPLEX { return COMPLEX(0);});
	k=1;
}
POWERSERIES::POWERSERIES(COEF coef, int k) : coef(coef), k(k) {}
POWERSERIES::POWERSERIES(COMPLEX(*coef)(int), int k) : coef(coef), k(k) {}
POWERSERIES::POWERSERIES(COMPLEX center, COEF coef, int k) : center(center), coef(coef), k(k) {}
POWERSERIES::POWERSERIES(COMPLEX center, COMPLEX(*coef)(int), int k) : center(center), coef(coef), k(k) {}

COMPLEX POWERSERIES::evalHelper(int p, const COMPLEX& z, int d) {
	COMPLEX result(0);
	COMPLEX pow(1);
	COMPLEX cur;
	COMPLEX dz = z - center;
	INTEGER diffTerm(1);
	int t = k-p + 32*d;
	
	for(int i = 2;i<=d;i++)
		diffTerm*=i;

	for(int i = 0; i<=t;i++) {
		cur = (coef(i+d)*pow)*diffTerm;
		result = result + cur;

		diffTerm*=(i+d+1);
		diffTerm/=(i+1);

		pow =pow * dz;
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


POWERSERIES POWERSERIES::differentiateHelper(int d) {
	COEF seq = ([=] (int j) -> COMPLEX {
			INTEGER diffTerm(1);
			for(int i = 1;i<=d;i++)
			    diffTerm*=(i+j);
			return (coef(j+d))*diffTerm;
			});
	int k1 = k + round((log(REAL(k)*REAL(d))/ln2())*d + REAL(0.5));
	int k2 =  round((exp(REAL(1))*k) + REAL(0.5));
	int newK = std::max(k1,k2);
	return POWERSERIES(seq,newK);
}

POWERSERIES POWERSERIES::integralHelper(int d) {
	COEF seq = ([=] (int j) -> COMPLEX {
			if(j<d)
				return COMPLEX(0);

			INTEGER intTerm(1);
			for(int i = 0;i<d;i++)
					intTerm*=(j-i);
			return (coef(j-d))/intTerm;
			});
	return POWERSERIES(seq,k);
}

POWERSERIES POWERSERIES::differentiate(int d) {
	if(d==0)
		return *this;
	if(d>0)
		return this->differentiateHelper(d);
	else
		return this->integralHelper(-d);
}
POWERSERIES POWERSERIES::continuation(const COMPLEX& z, int newK) {
	COEF seq = ([=] (int j) -> COMPLEX {
				INTEGER factorial(1);
				for(int i=2;i<=j;i++)
					factorial*=i;
				return this->eval(z,j)/factorial;
			});
	return POWERSERIES(z,seq,newK);
}
} //namespace iRRAM
