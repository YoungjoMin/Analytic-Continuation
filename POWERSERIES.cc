# include "POWERSERIES.h"

namespace iRRAM{
/**
 * @brief default constructor. assigns zero
 */
POWERSERIES::POWERSERIES() {
	memCount = -1;
	coefMemory = NULL;
	z0 = COMPLEX(0);
	coef = ([] (int j) -> COMPLEX { return COMPLEX(0);});
	k=1;
}
POWERSERIES::POWERSERIES(COEF coef, int k) : z0(0), coef(coef), k(k), memCount(-1), coefMemory(NULL){}
POWERSERIES::POWERSERIES(COMPLEX(*coef)(int), int k) : z0(0), coef(coef), k(k), memCount(-1), coefMemory(NULL) {}
POWERSERIES::POWERSERIES(COMPLEX z0, COEF coef, int k) : z0(z0), coef(coef), k(k), memCount(-1), coefMemory(NULL) {}
POWERSERIES::POWERSERIES(COMPLEX z0, COMPLEX(*coef)(int), int k) : z0(z0), coef(coef), k(k), memCount(-1), coefMemory(NULL) {}

POWERSERIES::POWERSERIES(const POWERSERIES& other) {
	k = other.k;
	coef = other.coef;
	z0 = other.z0;
	coefMemory = NULL;
	memCount = -1;
}
POWERSERIES::POWERSERIES(POWERSERIES&& other) {
	k = other.k;
	coef = std::move(other.coef);
	z0 = std::move(other.z0);
	coefMemory = other.coefMemory;
	memCount = other.memCount;
	other.memCount = -1;
	other.coefMemory = NULL;
}
POWERSERIES::~POWERSERIES() {
	if(coefMemory!=NULL)
		delete [] coefMemory;
}

int POWERSERIES::memorizeCoef(int count) {
	coefMemory = new COMPLEX[count+1];//memorize coef(0), ... coef(count)
	if(coefMemory==NULL)
		throw;

	for(int i = 0;i<=count;i++)
		coefMemory[i] = coef(i);
	return 1;
}

int POWERSERIES::findIterationCount(int p, int d) {
	int s = 0, e = k-p+64*d;
	int mid;
	double a, b;
	b = p-k + d*std::log2(k);
	while((e-s)>=2){
		mid = (e+s)/2;
		a = std::log2(mid+d+2)*d - mid;
		if(a<=b) e = mid;
		else s = mid;
	}

	return e;
}
/**
 * @brief this calculates value of the d-th derivative of the POWERSERIES at z with precison \f$2^p\f$.
 * @remarks To approximate Infinite sum, this use partial sum.
 * d-th derivative of the power series becomes  \f$ f^{d} = \displaystyle\sum_{j=0}^{\infty}{(j+1)(j+2) \cdots (j+d) a_{j+d} (z - z_0)^j} \f$ \n
 * the INTEGER diffTerm means the \f$ (j+1)(j+2) \cdots (j+d) \f$ terms.\n
 * the first for loop calculates first \f$ (j+1)(j+2) \cdots (j+d) \f$ term.\n
 * next loop calculates the partial sum. at each iteration it calculates i-th term (\f$ (j+1)(j+2) \cdots (i+d) a_{i+d} (z - z_0)^j\f$ )
 * and add to the result. \n\n\n
 * To explain why this algorithm iterates for k-p+32d times.\n
 * d-th derivative of the power series becomes  \f$ f^{d} = \displaystyle\sum_{j=0}^{\infty}{(j+1)(j+2) \cdots (j+d) a_{j+d} (z - z_0)^j} \f$ \n
 * Let partial sum of the power series \f$ f^{d}_{t} = \displaystyle\sum_{j=0}^{t}{(j+1)(j+2) \cdots (j+d) a_{j+d} (z - z_0)^j}  \f$\n
 * since coefficients satisfy \f$ |a_j| \leq 2^k k^j \f$ and domain is \f$ |z - z_0| \leq \frac{1}{2k}\f$ \n
 * By using Induction on d
 * \f{align*}{
 * |f^{d} - f^{d}_{t} | &\leq \displaystyle\sum_{j=t+1}^{\infty} {|(j+1)(j+2) \cdots (j+d)| |a_{j+d}| |z-z_0|^j }\\
 *  &\leq \displaystyle\sum_{j=t+1}^{\infty} {(j+1)(j+2) \cdots (j+d) 2^k k^{j+d} \frac{1}{2^j k^j} }\\
 *  &\leq 2^k k^d \displaystyle\sum_{j=t+1}^{\infty} {2^{-j} (j+1)(j+2) \cdots (j+d) }\\
 *  &\leq 2^k k^d 2^{-t} (t+d+2)^d \\
 *  &= 2^{k-t} k^d (t+d+2)^d
 * \f}
 * So, To make difference less than \f$2^p\f$ \n
 * \f$ |f^{d} - f^{d}_{t} | \leq 2^{k-t} k^d (t+d+2)^d \leq 2^p \f$ \n
 * So, \f$ 2^{k-t} k^d (t+d+2)^d \leq 2^p \f$ \n
 * By taking log, \f$ k-t + d\log_{2}{k} + d \log_{2}{(t+d+2)} \leq p \f$ \n
 * Then \f$ k - p + d\log_{2}{k(t+d+2)} \leq t \f$ \n
 * Supposing using int which is bounded by \f$ 2^{32} \f$ \n
 * By putting \f$2^{32}\f$ inside log, \f$ k - p + 32d \leq t\f$ \n
 * Therefore with \f$ t \geq k - p + 32d \f$, \f$ |f^{d} - f^{d}_{t} | \leq 2^p \f$ \n
 * With iterating more than k - p + 32d makes error less than \f$2^p\f$.
 * 
 * @warning given z should inside the domain \f$ |z - z_0| \leq  \frac{1}{2k} \f$ \n
 * @todo Instead of supposing that we are using int, use the inequality to get better iteration bound
 */
COMPLEX POWERSERIES::evalHelper(int p, const COMPLEX& z, int d) {
	COMPLEX result(0);
	COMPLEX pow(1);
	COMPLEX cur;
	COMPLEX dz = z - z0;
	INTEGER diffTerm(1);
	int t = findIterationCount(p,d);
	COEF CoefWithMem = ([=](int i)->COMPLEX{
		if(i<=memCount)
			return coefMemory[i];
		else
			return coef(i);
	});
	
	for(int i = 2;i<=d;i++)
		diffTerm*=i;// to calculate (j+1)(j+2)...(j+d)

	for(int i = 0; i<=t;i++) {
		cur = (CoefWithMem(i+d)*pow)*diffTerm;// add i-th term of the d-th derivative of the 
		result = result + cur;

		diffTerm*=(i+d+1);
		diffTerm/=(i+1);

		pow =pow * dz;
	}
	return result;
}

/**
 * @brief this calculates value of the d-th derivative of the POWERSERIES at z.
 * @see POWERSERIES::evalHelper
 * @warning given z should inside the domain \f$ |z - z_0| \leq  \frac{1}{2k} \f$
*/
COMPLEX POWERSERIES::eval(const COMPLEX& z, int d) {
	static POWERSERIES * Tfp;
	static int Td;//using capture is not able for iRRAM::limit. so, use static copies
	Tfp = this;
	Td = d;
	COMPLEX (*lambda)(int, const COMPLEX& z)  = ([] (int p, const COMPLEX& z) {
			return Tfp->evalHelper(p,z, Td);
	}); // to use iRRAM::limit

	return limit(lambda,z);// calculates lim p->inf evalHelper(p,z,d);
}

/**
 * let f1 is the power series with coefficient \f$ a_j \f$ and \f$ k_1 \f$ \n
 * and f2 is the power series with coefficient \f$ b_j \f$ and \f$ k_2 \f$ \n
 * then resulting power series has coefficient \f$ c_j = a_j + b_j \f$ \n
 * \f{align*}{
 * |c_j| &\leq |a_j| + |b_j|\\
 * &\leq 2^{k_1} {k_1}^j + 2^{k_2} {k_2}^j \\
 * &\leq 2^{\max(k_1,k_2)} {\max(k_1,k_2)}^j + 2^{\max(k_1,k_2)} {\max(k_1,k_2)}^j \\
 * &\leq 2^{\max(k_1,k_2) +1} {\max(k_1,k_2)}^j \\
 * &\leq 2^{\max(k_1,k_2)+1} (\max(k_1,k_2)+1)^j
 * \f}
 * therefore resulting POWERSERIES has \f$ max(k_1, k_2)+1 \f$ for k value.\n
 * @warning given two POWERSERIES should have same z0.
 */
POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2) {
	//suppose f1 and f2 has same z0
	COEF lambda = ([=] (int j) -> COMPLEX {
			return f1.coef(j)+f2.coef(j);
			}); //result coefficient sequence
	return POWERSERIES(f1.z0, lambda, std::max(f1.k,f2.k)+1);
}
/**
 * @see POWERSERIES::operator +
 * 
 * 
 * by the same reason of the operator +\n
 * this resulting POWERSERIES has \f$ max(k_1, k_2)+1 \f$ for k value.\n
 * @warning given two POWERSERIES should have same z0.
 */
POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2) {
	//suppose f1 and f2 has same z0
	COEF lambda = ([=] (int j) -> COMPLEX {
			return f1.coef(j)-f2.coef(j);
			});
	return POWERSERIES(f1.z0, lambda, std::max(f1.k,f2.k)+1);
}
/**
 * @warning given two POWERSERIES should have same z0.
 */
POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2) {
	//suppose f1 and f2 has same z0
COEF lambda = ([=] (int j) -> COMPLEX {
					COMPLEX result(0);
					for(int i = 0;i<=j;i++)
						result= result + (f1.coef(i)*f2.coef(j-i));
					return result;
		}	);
	return POWERSERIES(f1.z0, lambda, f1.k+f2.k);
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
POWERSERIES& POWERSERIES::operator=(const POWERSERIES& other) {
	if(this == (&other))
		return (*this);
	k = other.k;
	coef = other.coef;
	z0 = other.z0;
	coefMemory = NULL;
	memCount = -1;
	return (*this);
}
POWERSERIES& POWERSERIES::operator=(POWERSERIES&& other) {
	if(this == (&other))
		return (*this);
	k = other.k;
	coef = std::move(other.coef);
	z0 = std::move(other.z0);
	coefMemory = other.coefMemory;
	memCount = other.memCount;
	other.memCount = -1;
	other.coefMemory = NULL;
	return (*this);
}

POWERSERIES POWERSERIES::differentiateHelper(int d) {
	COEF seq = ([=] (int j) -> COMPLEX {
			INTEGER diffTerm(1);
			for(int i = 1;i<=d;i++)
			    diffTerm*=(i+j);
			return (coef(j+d))*diffTerm;
			});
	int k1 = k + round((log(REAL(k*d))/ln2())*d + REAL(0.5));
	int k2 =  round((exp(REAL(1))*k) + REAL(0.5));
	int newK = std::max(k1,k2);
	return POWERSERIES(z0, seq, newK);
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
	return POWERSERIES(z0, seq,k);
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
				static int prev=0;
				static INTEGER prevFactorial(0);
				INTEGER factorial(1);
				if( (prev+1)== j) {
					factorial = prevFactorial*j;
				}
				else {
					for(int i=2;i<=j;i++)
						factorial*=i;
				}
				prev = j;
				prevFactorial = std::move(factorial);
				return this->eval(z,j)/prevFactorial;
			});
	return POWERSERIES(z,seq, newK);
}
} //namespace iRRAM
