# pragma once

# include "iRRAM/lib.h"
# include <functional>
# include <algorithm>


namespace iRRAM {

using COEF = std::function<COMPLEX(int)>;
/**
 * @brief complex function of power series \f$ \displaystyle\sum_{j=0}^{\infty}{a_j (z - z_0)^j} \f$ with domain \f$ |z - z_0| \leq \frac{1}{2k}\f$
 * @remarks POWERSERIES is a class of power series function defined on complex plane.
 * \f[ 
 * \displaystyle\sum_{j=0}^{\infty}{a_j (z - z_0)^j}
 *  \f]
 * POWERSERIES has three member variable coef, centor, k.\n
 * coef is a function which gets nonnegative integer n and returns n-th coefficient \f$a_n\f$\n
 * COMPLEX z0 is centor of the power series \f$z_0 \f$.\n
 * k is positive integer about size of coefficient and the domain of the power series. it should be given when the POWERSERIES is initialized.\n
 * with given k, coefficient \f$ a_j \f$ should satisfy 
 * \f[ 
 * |a_j| \leq 2^k k^j
 * \f]
 * then POWERSERIES is function defined on a circle.
 * \f[ 
 * |z-z_0| \leq  \frac{1}{2k}
 *  \f]
 * 
 * @warning 
 * member variable function coef ((int) -> COMPLEX) should satisfy \f$ |a_j| \leq 2^k k^j \f$
 */
class POWERSERIES
{
public:
POWERSERIES();
POWERSERIES(COEF coef, int k);
POWERSERIES(COMPLEX(*coef)(int), int k);
POWERSERIES(COMPLEX z0, COEF coef, int k);
POWERSERIES(COMPLEX z0, COMPLEX(*coef)(int), int k);

friend POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2);

POWERSERIES& operator+=(const POWERSERIES& f);
POWERSERIES& operator-=(const POWERSERIES& f);
POWERSERIES& operator*=(const POWERSERIES& f);

COMPLEX eval(const COMPLEX& z, int d=0);//only for nonnegative d.
POWERSERIES differentiate(int d);//also for negative d. and return function has integral constant all 0.
POWERSERIES continuation(const COMPLEX& z); //POWERSEREIS at point z

private:
POWERSERIES differentiateHelper(int d);//for only differentiation (d>0)
POWERSERIES integralHelper(int d);// for only integral (d>0)
COMPLEX evalHelper(int p, const COMPLEX& z, int d);
/**
 * centor \f$ z_0 \f$ of the power series \f$ \displaystyle\sum_{j=0}^{\infty}{a_j (z - z_0)^j} \f$
 */
COMPLEX z0;
/**
 * coefficient function (int) -> COMPLEX\n
 * with input nonnegative integer n and returns n-th coefficient \f$a_n\f$ of the powerseries
 */
COEF coef;
/**
 * variable about size of coefficient and the domain of the powerseries
 */
int k;
};

} //namespace iRRAM
