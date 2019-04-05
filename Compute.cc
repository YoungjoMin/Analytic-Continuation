# include "iRRAM.h"
# include "./POWERSERIES.h"

using namespace iRRAM;


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

void valueTest() {
	int numOfFunc = 3;
	COEF coefList[] = {IdentitySeq, InvXSeq, ExpSeq};
	int  kList[]    = {1, 1, 1};

	int numOfComplex = 2;//5;
	COMPLEX complexList[] = {COMPLEX(0.1), COMPLEX(0.1,0.1), COMPLEX(0,-0.2),COMPLEX(0), 
		COMPLEX(-0.1,-0.1)};
	int numOfD = 2;//4;
	int DList[] = {1, 5, 10, 100};

	cout<<setRwidth(40);

	POWERSERIES * f = new POWERSERIES[numOfFunc];
	for(int i = 0;i<numOfFunc;i++) {
		f[i] = POWERSERIES(coefList[i], kList[i]);
	}

	for(int i = 0;i<numOfFunc;i++) {
		for(int j = 0;j<numOfComplex;j++) {
			for(int k = 0;k<numOfD;k++) {
				POWERSERIES &g = f[i];
				COMPLEX &x = complexList[j];
				int curD = DList[k];
				
				COMPLEX intgx = g.eval(x,-curD);
				COMPLEX diffgx = g.eval(x,curD);
				COMPLEX diffgx2 = (g.differentiate(curD)).eval(x);
				COMPLEX intgx2 = (g.differentiate(-curD)).eval(x);
				
				if(!bound(abs(intgx - intgx2),-50)) {
					std::cout<<"error at integrate (i,j,k) = ("<<i<<", "<<j<<", "<<k<<")\n";
					cout<<abs(intgx - intgx2)<<"\n";
					cout<<(intgx)._real <<"\n"<<(intgx)._imag<<" i\n";
					cout<<(intgx2)._real <<"\n"<<(intgx2)._imag<<" i\n";
				}
				if(!bound(abs(diffgx - diffgx2),-50)) {
					std::cout<<"error at differentiate(i,j,k) = ("<<i<<", "<<j<<", "<<k<<")\n";
					cout<<abs(diffgx - diffgx2)<<"\n";
					cout<<(diffgx)._real <<"\n"<<(diffgx)._imag<<" i\n";
					cout<<(diffgx2)._real <<"\n"<<(diffgx2)._imag<<" i\n";
				}
				std::cout<<"("<<i<<", "<<j<<", "<<k<<") endd\n";
			}
		}
	}
	
	std::cout<<"testEnd"<<std::endl;
}

void removeTest() {

	COMPLEX z(REAL(0.5));
	POWERSERIES * f = new POWERSERIES(ExpSeq,1);
	POWERSERIES * g = new POWERSERIES(InvXSeq,1);

	POWERSERIES sum = (*f)+(*g);
	POWERSERIES sub = (*f)-(*g);
	POWERSERIES mul = (*f)*(*g);

	COMPLEX beforeSum = sum.eval(z);
	COMPLEX beforeSub = sub.eval(z);
	COMPLEX beforeMul = mul.eval(z);

	delete f;
	delete g;

	COMPLEX afterSum = sum.eval(z);
	COMPLEX afterSub = sub.eval(z);
	COMPLEX afterMul = mul.eval(z);

	cout<<beforeSum._real<<"\n";
	cout<<afterSum._real<<"\n";
	cout<<beforeSub._real<<"\n";
	cout<<afterSub._real<<"\n";
	cout<<beforeMul._real<<"\n";
	cout<<afterMul._real<<"\n";

}

void modifyTest() {

	COMPLEX z(REAL(0.5));
	POWERSERIES f(ExpSeq,1);
	POWERSERIES g(InvXSeq,1);

	POWERSERIES sum = f+g;
	POWERSERIES sub = f-g;
	POWERSERIES mul = f*g;

	COMPLEX beforeSum = sum.eval(z);
	COMPLEX beforeSub = sub.eval(z);
	COMPLEX beforeMul = mul.eval(z);

	f = POWERSERIES(IdentitySeq,1);
	g = POWERSERIES(InvXSeq,1);

	COMPLEX afterSum = sum.eval(z);
	COMPLEX afterSub = sub.eval(z);
	COMPLEX afterMul = mul.eval(z);

	cout<<beforeSum._real<<"\n";
	cout<<afterSum._real<<"\n";
	cout<<beforeSub._real<<"\n";
	cout<<afterSub._real<<"\n";
	cout<<beforeMul._real<<"\n";
	cout<<afterMul._real<<"\n";
}


void compute()
{
	//modifyTest();
	//removeTest();
	valueTest();
}
