# include "iRRAM.h"
# include "./POWERSERIES.h"

using namespace iRRAM;

COMPLEX invXSeq(int k) {
	return COMPLEX(1);
}

COMPLEX seq(int k) {
	static REAL prevSeq;
	static int prevK;
	REAL a(1);
	int i=2;

	if(k<=1)
		return COMPLEX(1);

	if(prevK<k) {
		a = prevSeq;
		i = prevK+1;
	}
	for(;i<=k;i++)
			a/=i;

	prevSeq = a;
	prevK = k;
	return COMPLEX(a);
}
/*
COMPLEX seq(int k) {
	REAL a(1);
	if(k<=1)
		return COMPLEX(1);
	for(int i=2;i<=k;i++)
		a/=i;
	return COMPLEX(a);
}
*/

void compute2()
{
	COMPLEX pos[5]= {COMPLEX(0.5), COMPLEX(-0.5), COMPLEX(REAL(0),REAL(0.5)),
									COMPLEX(0, REAL(-0.5)), COMPLEX(pi()/8)};
	POWERSERIES f(invXSeq,1);
	COMPLEX result;
	for(int i=0;i<5;i++) {
		result = f.eval(pos[i]);
		cout<<"using taylor series\n";
	  cout<<result._real<<" + "<<result._imag<<" i\n";
		result = COMPLEX(1)/(COMPLEX(1)-pos[i]);
		cout<<"using 1/1-x\n";
	  cout<<result._real<<" + "<<result._imag<<" i\n\n\n";
	}
	for(int i=0;i<5;i++) {
		result= f.eval(pos[i],1);
		cout<<"using taylor series\n";
	  cout<<result._real<<" + "<<result._imag<<" i\n";

		result = COMPLEX(1)-pos[i];
		result = result*result;
		result = COMPLEX(1)/result;
		cout<<"using 1/(1-x)^2\n";
	  cout<<result._real<<" + "<<result._imag<<" i\n\n\n";

	}


}

void compute1()
{

	POWERSERIES f(seq,1);
	COMPLEX result;
	cout<<"f = e^x\n";
	cout<<setRwidth(50);
	result = f.eval(COMPLEX(REAL(0),pi()),100);
	cout<<"D^100 f(0+pi i) = \n";
	cout<<result._real<<" + "<<result._imag<<" i\n\n";

	result = f.eval(COMPLEX(REAL(0),REAL(0)));
	cout<<"f(0) = \n";
	cout<<result._real<<" + "<<result._imag<<" i\n\n";

	result = f.eval(COMPLEX(REAL(0),pi()));
	cout<<"f(0+pi i) = \n";
	cout<<result._real<<" + "<<result._imag<<" i\n\n";

	result = f.eval(COMPLEX(REAL(0.25),REAL(0.25)), 1);
	cout<<"D^1 f(0.25+0.25i) = \n";
	cout<<result._real<<" + "<<result._imag<<" i\n\n";

	result = f.eval(COMPLEX(REAL(0),REAL(0.5)), 1000);
	cout<<"D^1000 f(0.5 i) = \n";
	cout<<result._real<<" + "<<result._imag<<" i\n\n";
	return;
}
void compute()
{
	compute1();
	compute2();
}
