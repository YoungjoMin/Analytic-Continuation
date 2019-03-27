# include "iRRAM.h"
# include "./Poly.h"
# include "./ANALYTIC.h"

using namespace iRRAM;


COMPLEX seq(int k) {
	REAL a = 1;
	if(k<=1)
		return COMPLEX(1);

	for(int i=2;i<=k;i++)
			a/=i;
	return COMPLEX(a);
}

void compute()
{
		ANALYTIC f(seq, 1);
		cout<<setRwidth(100);
		COMPLEX result = f(COMPLEX(0.5));
		cout<<result._real<<" + "<<result._imag<<" i\n";
		return;
}
