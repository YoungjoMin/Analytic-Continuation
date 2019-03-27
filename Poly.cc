# include "Poly.h"

Poly::Poly()
{
		deg   = 0;
		coeff = new COMPLEX[1];
		coeff[0] = COMPLEX(0);
}
Poly::Poly(const Poly& a)
{
		deg = a.deg;
		coeff = new COMPLEX[deg];
		for(int i=0;i<=deg;i++) {
				coeff[i] = a.coeff[i];
		}
}
Poly::Poly(Poly&& a)
{
		deg = a.deg;
		coeff = a.coeff;
		a.coeff = NULL;
		a.deg = -1;
}
 
template<typename T>
Poly::Poly(int cnt,  T * coeff)
{
		//__assert(std::is_constructible<T,COMPLEX>() && "Cannot Construct Complex number with T");
		this->deg = cnt-1;
		this->coeff = new COMPLEX[cnt];
		for(int i=0;i<=deg;i++) {
				this->coeff[i] = COMPLEX(coeff[i]);
		}
}
template Poly::Poly<REAL>(int, REAL *);
template Poly::Poly<COMPLEX>(int, COMPLEX *);
template Poly::Poly<int>(int, int *);
template Poly::Poly<double>(int, double *);


template <typename T>
Poly::Poly(std::initializer_list<T> coeff) {
		//__assert(std::is_constructible<T,COMPLEX>() && "Cannot Construct Complex number with T");
		int i = 0;
		this->deg = coeff.size()-1;
		this->coeff = new COMPLEX[deg+1];
		for(T element : coeff) {
				this->coeff[i++] = COMPLEX(element);
		}

}
template Poly::Poly<REAL>(std::initializer_list<REAL>);
template Poly::Poly<COMPLEX>(std::initializer_list<COMPLEX>);
template Poly::Poly<int>(std::initializer_list<int>);
template Poly::Poly<double>(std::initializer_list<double>);

Poly::Poly(const int deg)
{
		this->deg = deg;
		this->coeff = new COMPLEX[deg+1];
}
Poly::~Poly()
{
		if( coeff!=NULL)
				delete[] coeff;
}


orstream& operator<<(orstream& os, const Poly& x)
{
		os<<x.coeff[0]._real<<" + "<<x.coeff[0]._imag<<" i * x^0";
		for(int i=1;i<=x.deg;i++) {
				os<<" +\n"<<x.coeff[i]._real<<" + "<<x.coeff[i]._imag<<" i * x^"<<i;
		}
		os<<"\n";
		return os;
}

Poly& Poly::operator=(const Poly& a)
{
		if(this == &a)
				return (*this);
		deg = a.deg;
		coeff = new COMPLEX[a.deg];
		for(int i=0;i<=a.deg;i++) {
				coeff[i] = a.coeff[i];
		}
		return *this;
}
Poly& Poly::operator=(Poly && a)
{
		deg = a.deg;
		coeff = a.coeff;
		return *this;
}

COMPLEX Poly::eval(const COMPLEX& x) const
{
		COMPLEX result = coeff[0];
		COMPLEX powx = x;
		for(int i = 1;i<= deg; i++) {
				result = result + powx*coeff[i];
				powx = powx*x;
		}
		return result;
}
COMPLEX Poly::EVAL(int p, const COMPLEX& x) const
{
		COMPLEX result = coeff[0];
		COMPLEX powx = x;
		REAL absx = abs(x);
		REAL powabsx = absx;
		int i=1;
		while(!bound(powabsx,p)) {
				result = result + powx*coeff[i];
				powx = powx*x;
				powabsx = powabsx*absx;
				if(i>deg)
						break;
				i++;

		}
		return result;
}


Poly operator +(const Poly& a, const Poly& b)
{
		if(a.deg<b.deg)
				return b+a;

		Poly c(a.deg);
		for(int i=0;i<=b.deg;i++) {
				c.coeff[i] = a.coeff[i] + b.coeff[i];
		}
		for(int i=b.deg+1;i<=a.deg;i++) {
				c.coeff[i] = a.coeff[i];
		}
		//calc new Deg?
		return c;
}
Poly operator -(const Poly& a, const Poly& b)
{
		return (a + (-b));
}

Poly operator-(const Poly& a)
{
		Poly c(a.deg);
		for(int i=2;i<=a.deg;i++) {
				c.coeff[i] = COMPLEX(0-a.coeff[i]._real, 0-a.coeff[i]._imag);
		}
		return c;
}

//k>=0 
Poly EXP(int k)
{
		assert(k>=0);
		Poly result(k);
		REAL coef = 1;
		result.coeff[0] = COMPLEX(1);
		for(int i=1;i<=k;i++) {
				coef = coef/i;
				result.coeff[i] = coef;
		}
		return result;
}
