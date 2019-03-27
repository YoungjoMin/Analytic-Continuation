# pragma once

# include "iRRAM/lib.h"
# include <functional>

namespace iRRAM {

using SEQUENCE = std::function<REAL(int)>;
//how about using recurrence relation

class ANALYTIC
{
public:
ANALYTIC(SEQUENCE seq, int k);
REAL operator () (const REAL& z);

friend REAL evalHelper(int p, const ANALYTIC& f, const REAL& z);

private:

SEQUENCE coef;
int k;
};

} //namespace iRRAM
