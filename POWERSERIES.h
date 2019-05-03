# pragma once

# include "iRRAM/lib.h"
# include <functional>
# include <algorithm>


namespace iRRAM {

using COEF = std::function<COMPLEX(int)>;

class POWERSERIES
{
public:
POWERSERIES();
POWERSERIES(COEF coef, int k);
POWERSERIES(COMPLEX(*coef)(int), int k);
POWERSERIES(COMPLEX center, COEF coef, int k);
POWERSERIES(COMPLEX center, COMPLEX(*coef)(int), int k);

friend POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2);

POWERSERIES& operator+=(const POWERSERIES& f);
POWERSERIES& operator-=(const POWERSERIES& f);
POWERSERIES& operator*=(const POWERSERIES& f);


//	 for |z| <= 1/2k
COMPLEX eval(const COMPLEX& z, int d=0);
POWERSERIES differentiate(int d);//also for negative d constant term = 0
POWERSERIES continuation(const COMPLEX& z, int k); //POWERSEREIS at point z

private:
POWERSERIES differentiateHelper(int d);//for only differentiation (d>0)
POWERSERIES integralHelper(int d);// for only integral (d>0)
COMPLEX evalHelper(int p, const COMPLEX& z, int d);

COMPLEX center;
COEF coef;
int k;
};

} //namespace iRRAM
