
#ifndef iRAM_POLI_H
#define iRAM_POLI_H

# include "iRRAM/lib.h"
# include <initializer_list>
# include <type_traits>
# include <algorithm>
# include <iostream>
# include <cassert>

using namespace iRRAM;
class Poly
{
public:
		Poly();
		Poly(const Poly& a);
		Poly(Poly && a);
		template<typename T> 	Poly(int deg, T * coeff);
		template<typename T>  Poly(std::initializer_list<T> coeff);
		~Poly();

		COMPLEX eval(const COMPLEX& x) const;
		COMPLEX EVAL(int p, const COMPLEX& x) const;

		friend orstream& operator<<(orstream& os, const Poly& x);

		Poly& operator=(const Poly& a);
		Poly& operator=(Poly&& a);

		friend Poly operator+(const Poly & a, const Poly & b);
		friend Poly operator-(const Poly & a, const Poly & b);

		friend Poly operator-(const Poly & a);

		//assume k>=0 calculate exp(x) Taylor poly at 0 with degree k
		friend Poly EXP(int k);

private:
	 Poly(const int deg);

	 int deg;
	 COMPLEX * coeff;

};
Poly EXP(int k);
#endif
