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

void EvalCompTime()
{
	
	COMPLEX  z(REAL(0.1),REAL(0.1));
	COMPLEX approx, real = exp(z);
	int k[6] = {1, 5, 10, 50, 100, 500};
	int d[6] = {0, 5, 10, 50, 100, 500};
	for(int i = 0;i<6;i++) {
		for(int j = 0;j<6;j++) {
			POWERSERIES f(ExpSeq, k[i]);
			auto s = std::chrono::system_clock::now();
			real = f.evalHelper(0, z, d[j]);
			auto e = std::chrono::system_clock::now();
			std::chrono::duration<double> diff = e-s;
			std::cout<<"with i = "<<i<<" , j = "<<j<<" and time : "<< diff.count()<<"s\n";
		}
	}

}
void ContiCompTime()
{
	COMPLEX z;
	int n[2] = {0, 1};
	int p[4] = {0, -10, -100, -1000};
		POWERSERIES g[5];
		POWERSERIES f(InvXSeq, 1);
		g[0] = f;
		g[1] = g[0].continuation(z,1);
		g[2] = g[1].continuation(z,1);
		g[3] = g[2].continuation(z,1);
		g[4] = g[3].continuation(z,1);
	for(int i = 0;i<4;i++) {
		int j=2;
		for(int j = 0;j<5;j++) {
			f = f.continuation(z,1);
			auto s = std::chrono::system_clock::now();
			COMPLEX real = g[j].evalHelper(p[i], z, 0);
			auto e = std::chrono::system_clock::now();
			std::chrono::duration<double> diff = e-s;
			std::cout<<"with i = "<<i<<" , j = "<<j<<" and time : "<< diff.count()<<"s\n";
		}
	}
}
void compute()
{
	//EvalCompTime();
	ContiCompTime();
}
