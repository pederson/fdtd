#include <fdtd.h>

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>


using namespace fdtd;

// compile with:
// 			clang++-3.8 -std=c++14 -stdlib=libc++ -I../ testTEM.cpp -o test


int main(int argc, char * argv[]){

	// time-stepping parameters
	double cfl = 0.9;
	double dx = 1e-6;
	double dt=cfl*dx/c0;

	// std::cout << typeid(y).name() << std::endl;

	// YeeCell typedefs
	typedef YeeFieldsTEM<double, std::array> yeefieldT;
	static bool pmlx = true;
	static bool pmly = false;
	static bool pmlz = false;
	typedef PML<TEM, true, false, false> pmlT;
	typedef YeeCell<yeefieldT, VacuumPolarization, VacuumMagnetization, Neighbor1, pmlT> YeeCellVacTEM;
	
	// define the domain data structure and connect neighbors
	std::vector<YeeCellVacTEM> cells1(50);
	for (auto i=1; i<49; i++){
		cells1[i].setNeighborMin(0, cells1[i-1]);
		cells1[i].setNeighborMax(0, cells1[i+1]);
	}

	
	// setup the PML parameters
	PMLParameterModel p(dx);
	int nPML = 10;
	for (auto i=1; i<nPML+2; i++){
		double x = static_cast<double>(nPML-i+1)/static_cast<double>(nPML);
		double xm = x - 0.5/static_cast<double>(nPML);
		if (pmlx){
			cells1[i].setPMLParametersEx(p.K(x), p.S(x), p.A(x), dt);
			cells1[i].setPMLParametersHx(p.K(xm), p.S(xm), p.A(xm), dt);

			cells1[50-i].setPMLParametersEx(p.K(x), p.S(x), p.A(x), dt);
			cells1[49-i].setPMLParametersHx(p.K(xm), p.S(xm), p.A(xm), dt);
		}
	}

	for(auto j=0; j<50; j++){
			std::cout << " EKx: " << cells1[j].pmlEKx() ;
			std::cout << " ESx: " << cells1[j].pmlESx() ;
			std::cout << " EAx: " << cells1[j].pmlEAx() ;
			std::cout << " EBx: " << cells1[j].pmlEBx() ;
			std::cout << " ECx: " << cells1[j].pmlECx() ;
			std::cout << " HKx: " << cells1[j].pmlHKx() ;
			std::cout << " HSx: " << cells1[j].pmlHSx() ;
			std::cout << " HAx: " << cells1[j].pmlHAx() ;
			std::cout << " HBx: " << cells1[j].pmlHBx() ;
			std::cout << " HCx: " << cells1[j].pmlHCx() ;
			std::cout << std::endl;
	}
	throw -1;


	// start time-stepping
	for (auto t=0; t<200; t++){
		std::for_each(++cells1.begin(), --cells1.end(), YeeUpdateD<TEM>(dt,dx));
		// if (t<19) cells1[25].Dz() = sin(2*pi*0.05*t);
		cells1[25].Dz() += exp(-(t-10)*(t-10)*0.05);
		std::for_each(++cells1.begin(), --cells1.end(), UpdatePMLD<TEM>(dt,dx));
		std::for_each(++cells1.begin(), --cells1.end(), ConstantUpdateE<TEM>(1));
	
		std::for_each(++cells1.begin(), --cells1.end(), YeeUpdateB<TEM>(dt,dx));
		std::for_each(++cells1.begin(), --cells1.end(), UpdatePMLB<TEM>(dt,dx));
		std::for_each(++cells1.begin(), --cells1.end(), ConstantUpdateH<TEM>(1));
	
		// std::cout << "Ez:" ;
		for(auto i=0; i<50; i++) std::cout << ", " << cells1[i].Ez() ;
		std::cout << std::endl;
	}

	return 0;
}