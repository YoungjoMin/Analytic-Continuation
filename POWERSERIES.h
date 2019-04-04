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

friend POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2);//TODO
friend POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2);//TODO
friend POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2);//TODO

POWERSERIES& operator+=(const POWERSERIES& f);
POWERSERIES& operator-=(const POWERSERIES& f);
POWERSERIES& operator*=(const POWERSERIES& f);


COMPLEX eval(const COMPLEX& z, int d=0);
POWERSERIES differentiate(int d);//TODO also for negative d constant term = 0
POWERSERIES continuation(const COMPLEX& z, int k); //POWERSEREIS at point z

private:
COMPLEX evalHelper(int p, const COMPLEX& z, int d);
COEF coef;
int k;
};

} //namespace iRRAM
