# pragma once

# include "iRRAM/lib.h"
# include <functional>
# include <algorithm>


namespace iRRAM {

using COEF = std::function<COMPLEX(INTEGER)>;
/**
 * @brief complex function of power series \f$ \displaystyle\sum_{j=0}^{\infty}{a_j (z - w)^j} \f$ with domain \f$ |z - w| \leq \frac{1}{2k}\f$
 * @remarks POWERSERIES is a class of power series function defined on complex plane.
 * \f[ 
 * \displaystyle\sum_{j=0}^{\infty}{a_j (z - w)^j}
 *  \f]
 * POWERSERIES has three member variable coef, centor, k.\n
 * coef is a function which gets nonnegative integer n and returns n-th coefficient \f$a_n\f$\n
 * COMPLEX z0 is centor of the power series \f$w \f$.\n
 * k is positive integer about size of coefficient and the domain of the power series. it should be given when the POWERSERIES is initialized.\n
 * with given k, coefficient \f$ a_j \f$ should satisfy 
 * \f[ 
 * |a_j| \leq 2^k k^j
 * \f]
 * then POWERSERIES is function defined on a circle.
 * \f[ 
 * |z-w| \leq  \frac{1}{2k}
 *  \f]
 * 
 * @warning 
 * member variable function coef ((int) -> COMPLEX) should satisfy \f$ |a_j| \leq 2^k k^j \f$
 */
class POWERSERIES
{
public:
POWERSERIES();
POWERSERIES(COEF coef, INTEGER K, INTEGER k, COMPLEX w);
POWERSERIES(COMPLEX(*coef)(INTEGER), INTEGER K, INTEGER k, COMPLEX w);
POWERSERIES(const POWERSERIES& other);

friend POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2);

POWERSERIES& operator+=(const POWERSERIES& f);
POWERSERIES& operator-=(const POWERSERIES& f);
POWERSERIES& operator*=(const POWERSERIES& f);

POWERSERIES& operator=(const POWERSERIES& other);

COMPLEX eval(const COMPLEX& z, INTEGER d);//only for nonnegative d.
POWERSERIES differentiate(INTEGER d);//also for negative d. and return function has integral constant all 0.
private:

COMPLEX evalHelper(int p, const COMPLEX& z);
POWERSERIES differentiateHelper(INTEGER d);//for only differentiation (d>0)
POWERSERIES integralHelper(INTEGER d);// for only integral (d>0)
/**
 * centor \f$ w \f$ of the power series \f$ \displaystyle\sum_{j=0}^{\infty}{a_j (z - w)^j} \f$
 */
COMPLEX w;
/**
 * coefficient function (INTEGER) -> COMPLEX\n
 * with input nonnegative integer n and returns n-th coefficient \f$a_n\f$ of the powerseries
 */
COEF coef;
/**
 * variable about size of coefficient and the domain of the powerseries
 */
INTEGER K, k;
};

}//namespace iRRAM

