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

	int numOfComplex = 5;
	COMPLEX complexList[] = {COMPLEX(0.1), COMPLEX(0.1,0.1), COMPLEX(0,-0.2),COMPLEX(0), 
		COMPLEX(-0.1,-0.1)};

	cout<<setRwidth(40);

	POWERSERIES * f = new POWERSERIES[numOfFunc];
	for(int i = 0;i<numOfFunc;i++) {
		f[i] = POWERSERIES(coefList[i], kList[i]);
	}

	for(int i = 0;i<numOfFunc;i++) {
		for(int j = 0;j<numOfFunc;j++) {
			for(int k = 0;k<numOfComplex;k++) {
				POWERSERIES &g = f[i];
				POWERSERIES &h = f[j];
				COMPLEX &x = complexList[k];

				COMPLEX gx = g.eval(x);
				COMPLEX hx = h.eval(x);
				COMPLEX addghx = (g+h).eval(x);
				COMPLEX subghx = (g-h).eval(x);
				COMPLEX mulghx = (g*h).eval(x);

				if(!bound(abs(gx+hx - addghx),-50)) {
					std::cout<<"error at (i,j,k) = ("<<i<<", "<<j<<", "<<k<<")\n";
					cout<<abs(gx+hx-addghx)<<"\n";
					cout<<(gx+hx)._real <<"\n"<<(gx+hx)._imag<<"\n";
					cout<<(addghx)._real <<"\n"<<(addghx)._imag<<"\n";
				}
				if(!bound(abs(gx-hx - subghx),-50)) {
					std::cout<<"error at (i,j,k) = ("<<i<<", "<<j<<", "<<k<<")\n";
					cout<<abs(gx-hx-subghx)<<"\n";
					cout<<(gx-hx)._real <<"\n"<<(gx-hx)._imag<<"\n";
					cout<<(subghx)._real <<"\n"<<(subghx)._imag<<"\n";
				}
				if(!bound(abs(gx*hx - mulghx),-50)) {
					std::cout<<"error at (i,j,k) = ("<<i<<", "<<j<<", "<<k<<")\n";
					cout<<abs(gx*hx-mulghx)<<"\n";
					cout<<(gx*hx)._real <<"\n"<<(gx*hx)._imag<<"\n";
					cout<<(mulghx)._real <<"\n"<<(mulghx)._imag<<"\n";
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

void arithmeticAssignmentTest() {
	COMPLEX z(0.5);
	POWERSERIES f(ExpSeq,1);
	POWERSERIES g(InvXSeq,1);

	cout<< "f(z) = "<<f.eval(z)._real<<"\n";
	cout<< "g(z) = "<<g.eval(z)._real<<"\n";

	f +=g;
	cout<< "f += g\n";
	cout<< "f(z) = "<<f.eval(z)._real<<"\n";	

	f -=g;
	cout<< "f -= g\n";
	cout<< "f(z) = "<<f.eval(z)._real<<"\n";
	
	f *=g;
	cout<< "f *= g\n";
	cout<< "f(z) = "<<f.eval(z)._real<<"\n";
	cout<< "f(z)/2 = "<<(f.eval(z)/2)._real<<"\n";
}

void compute()
{
	modifyTest();
	removeTest();
	valueTest();
	arithmeticAssignmentTest();
}
