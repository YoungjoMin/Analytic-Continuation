# pragma once

# include "iRRAM/lib.h"
# include <functional>
# include <algorithm>


namespace iRRAM {

using COEF = std::function<COMPLEX(INTEGER)>;
class ANALYTIC
{
public:
ANALYTIC(COEF coef, INTEGER L, INTEGER l, REAL x);
ANALYTIC(COMPLEX(*coef)(INTEGER), INTEGER L, INTEGER l, REAL x);
ANALYTIC(const ANALYTIC& other);

friend ANALYTIC operator +(const ANALYTIC& f1, const ANALYTIC& f2);
friend ANALYTIC operator -(const ANALYTIC& f1, const ANALYTIC& f2);
friend ANALYTIC operator *(const ANALYTIC& f1, const ANALYTIC& f2);

ANALYTIC& operator+=(const ANALYTIC& f);
ANALYTIC& operator-=(const ANALYTIC& f);
ANALYTIC& operator*=(const ANALYTIC& f);

ANALYTIC& operator=(const ANALYTIC& other);

COMPLEX eval(const COMPLEX& z, INTEGER d);//only for nonnegative d.
ANALYTIC differentiate(INTEGER d);//also for negative d. and return function has integral constant all 0.
private:

COMPLEX evalHelper(int p, const COMPLEX& z);
ANALYTIC differentiateHelper(INTEGER d);//for only differentiation (d>0)
ANALYTIC integralHelper(INTEGER d);// for only integral (d>0)

REAL x;
COEF coef;
INTEGER K, k;
};

}//namespace iRRAM

