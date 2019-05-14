# include "iRRAM.h"
# include "./POWERSERIES.h"
# include <chrono>

using namespace iRRAM;

COMPLEX fiboSeq(int k) {
	if(k<=1)
		return COMPLEX(1);
	INTEGER a, b(1), c(1);
	int i;
	for(i=2;i<=k;i++) {
		a = b;
		b = c;
		c = a+b;
	}
	return COMPLEX(c);
}

COMPLEX IdentitySeq(int k) {
	return (k==1 ? COMPLEX(1) : COMPLEX(0));
}
COMPLEX InvXSeq(int k) {
	return COMPLEX(1);
}

COMPLEX ExpSeq(int k) {
	INTEGER a(1);
	if(k<=1)
		return COMPLEX(1);
	for(int i=2;i<=k;i++)
		a*=i;
	return COMPLEX(REAL(1)/REAL(a));
}
void myPrint(COMPLEX z) {
	cout<<z._real<<" + "<<z._imag<<" i\n";
}

void IntTest()
{
	COMPLEX ret1, ret2, z;
	POWERSERIES f(InvXSeq,1);
	POWERSERIES g(IdentitySeq,1);
	POWERSERIES h(ExpSeq,1);


	POWERSERIES k = f.differentiate(-1);
	cout<<"integral 1/1-x = -log(1-x)\n";
	cout<<" integral 1/1-x at x = 1/2\n";
	z = COMPLEX(REAL(1)/REAL(2));
	ret1 = k.eval(z);
	myPrint(ret1);
	cout<<" -log(1-x) at x = 1/2"<<"\n";
	ret2 = -log(REAL(1)-z);
	myPrint(ret2);

	cout<<"\n integral 1/1-x at x = 1/2i\n";
	z = COMPLEX(REAL(0), REAL(1)/REAL(2));
	ret1 = k.eval(z);
	myPrint(ret1);
	cout<<" -log(1-x) at x = 1/2i"<<"\n";
	ret2 = -log(COMPLEX(1) - z);
	myPrint(ret2);
	cout<<"\n\n";


	k = g.differentiate(-3);
	cout<< "3 times integral x  = 1/24 x^4\n";
	cout<< "integral x at x = 0.1\n";
	z = COMPLEX(REAL(1)/REAL(10));
	ret1 = k.eval(z);
	myPrint(ret1);
	cout<< "using 1/24 x^4 at x = 0.1\n";
	ret2 = (z*z*z*z)/REAL(24);
	myPrint(ret2);
	cout<<"\n\n";


	k = h.differentiate(-2);
	cout<<" 2 times integral e^x = e^x - 1 -x\n";
	cout<< "integral x at x = 0.1 + 0.1i\n";
	z = COMPLEX(REAL(1)/REAL(10), REAL(1)/REAL(10));
	ret1 = k.eval(z);
	myPrint(ret1);
	cout<<"using e^x -1 - x\n";
	ret2 = h.eval(z) - COMPLEX(1) - z;
	myPrint(ret2);
	cout<<"\n\n";
	return;
}

void DiffTest()
{
	COMPLEX ret1, ret2, z;
	POWERSERIES f(InvXSeq,1);
	POWERSERIES g(IdentitySeq,1);
	POWERSERIES h(ExpSeq,1);


	POWERSERIES k = f.differentiate(3);
	cout<<"3 times differentiate 1/1-x = 3! / (1-x)^4\n";
	cout<<"3 diff 1/1-x at x = 1/2\n";
	z = COMPLEX(REAL(1)/REAL(2));
	ret1 = k.eval(z);
	myPrint(ret1);
	cout<<" 3! / (1-x)^4 at x = 1/2"<<"\n";
	ret2 = COMPLEX(6) / ((1-z)*(1-z)*(1-z)*(1-z));
	myPrint(ret2);

	cout<<"\n 3 diff 1/1-x at x = 1/2i\n";
	z = COMPLEX(REAL(0), REAL(1)/REAL(2));
	ret1 = k.eval(z);
	myPrint(ret1);
	cout<<" 3! / (1-x)^4 at x = 1/2"<<"\n";
	ret2 = COMPLEX(6) / ((1-z)*(1-z)*(1-z)*(1-z));
	myPrint(ret2);
	cout<<"\n\n";


	k = h.differentiate(10);
	cout<<" 10 times differentiate e^x = e^x\n";
	cout<< "10 diff x at x = 0.1 + 0.1i\n";
	z = COMPLEX(REAL(1)/REAL(10), REAL(1)/REAL(10));
	ret1 = k.eval(z);
	myPrint(ret1);
	cout<<"using e^x\n";
	ret2 = h.eval(z);
	myPrint(ret2);
	cout<<"\n\n";
	return;
}

void ACTest()
{
	
	COMPLEX ret1, ret, z;
	POWERSERIES f(ExpSeq,1);

	cout<<"continuation of f= exp(x) at x= 0.1 then calc value of f at 0.1+0.1i\n";
	POWERSERIES f1 = f.continuation(COMPLEX(REAL(1)/REAL(10)));
	z = COMPLEX(REAL(1)/REAL(10), REAL(1)/REAL(10));
	auto s1 = std::chrono::system_clock::now();
	ret1 = f1.eval(z);
	auto e1 = std::chrono::system_clock::now();
	auto s2 = std::chrono::system_clock::now();
	ret  = f.eval(z);
	auto e2 = std::chrono::system_clock::now();
	cout<<"values\n";
	myPrint(ret1);
	myPrint(ret);
	std::cout<<"continuation"<<(e1-s1).count()<<'\n';
	std::cout<<"evaluation"<<(e2-s2).count()<<'\n';



}

void compute()
{
	//IntTest();
	//DiffTest();
	ACTest();
	
}
