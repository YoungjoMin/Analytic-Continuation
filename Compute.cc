# include "iRRAM.h"
# include "./POWERSERIES.h"
# include <chrono>

using namespace iRRAM;
COMPLEX logSeq(int k) {
	if(k==0)
		return COMPLEX(0);
	else
		return COMPLEX(REAL(1)/REAL(k));
}

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
	POWERSERIES g[6];
	POWERSERIES f(InvXSeq, 1);
	g[0] = f;
	g[1] = g[0].continuation(z,1);
	g[1].memorizeCoef(5000);
	std::cout<<"memorize g1 OK\n";
	g[2] = g[1].continuation(z,1);
	g[2].memorizeCoef(1500);
	std::cout<<"memorize g2 OK\n";
	g[3] = g[2].continuation(z,1);
	g[3].memorizeCoef(200);
	std::cout<<"memorize g2 OK\n";
	g[4] = g[3].continuation(z,1);
	//g[4].memorizeCoef(500);
	g[5] = g[4].continuation(z,1);
	for(int j = 0;j<6;j++) {
		f = f.continuation(z,1);
		auto s = std::chrono::system_clock::now();
		COMPLEX real = g[j].evalHelper(0, z, 0);
		auto e = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = e-s;
		std::cout<<"with j = "<<j<<" and time : "<< diff.count()<<"s\n";
	}
}

void iterCnt()
{
	int cnt[6];
	POWERSERIES f;
	cnt[0] = f.findIterationCount(0,0);
	cnt[1] = f.findIterationCount(0,cnt[0]);
	cnt[2] = f.findIterationCount(0,cnt[1]);
	cnt[3] = f.findIterationCount(0,cnt[2]);
	cnt[4] = f.findIterationCount(0,cnt[3]);
	cnt[5] = f.findIterationCount(0,cnt[4]);
	std::cout<<cnt[0]<<'\n'
		<<cnt[1]<<'\n'
		<<cnt[2]<<'\n'
		<<cnt[3]<<'\n'
		<<cnt[4]<<'\n'
		<<cnt[5]<<'\n';
}

void compute()
{
	REAL ang = pi()/6;

	COMPLEX next(cos(ang),sin(ang));
	POWERSERIES f(COMPLEX(1),logSeq,1);
	POWERSERIES g = f.continuation(next,1);
	COMPLEX val0 = f.eval(COMPLEX(1));
	cout<<val0._real<<"\n";
	cout<<val0._imag<<"\n";
	COMPLEX val1 = f.eval(next);
	cout<<val1._real<<"\n";
	cout<<val1._imag<<"\n";
	COMPLEX val2 = g.eval(COMPLEX(1));
	cout<<val2._real<<"\n";
	cout<<val2._imag<<"\n";
	COMPLEX val3 = g.eval(next);
	cout<<val3._real<<"\n";
	cout<<val3._imag<<"\n";
}
