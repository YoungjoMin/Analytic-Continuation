# pragma once

# include "iRRAM/lib.h"

namespace iRRAM {
class ANALYTIC
{
public:
ANALYTIC(REAL (*seq)(int), int k);
REAL operator () (const REAL& z);

private:
REAL (*eval)(int, REAL);
REAL (*coef)(int);
int k;
};

} //namespace iRRAM
