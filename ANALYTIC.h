# pragma once

# include "iRRAM/lib.h"

namespace iRRAM {
class ANALYTIC
{
public:
ANALYTIC(COMPLEX (*seq)(int), int k);
COMPLEX eval(int d, const COMPLEX& z);

private:
COMPLEX evalHelper(int p, int d, const COMPLEX& z);

COMPLEX (*coef)(int);
int k;
};

} //namespace iRRAM
