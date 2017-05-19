#include <fdtd.h>

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>


using namespace fdtd;

// compile with:
// 			clang++-3.8 -std=c++14 -stdlib=libc++ -I../ test.cpp -o test


int main(int argc, char * argv[]){

	typedef YeeCell<YeeFields3D<double, std::array>, VacuumPolarization, VacuumMagnetization, Neighbor3, PML<ThreeD, true, true, true>> YeeCellVac3D;
	YeeCellVac3D y;

	// std::cout << typeid(y).name() << std::endl;

	typedef YeeFieldsTEM<double, std::array> yeefieldT;
	typedef PML<TEM, true, false, false> pmlT;
	typedef YeeCell<yeefieldT, VacuumPolarization, VacuumMagnetization, Neighbor1, pmlT> YeeCellVacTEM;
	std::vector<YeeCellVacTEM> cells1(50);

	for (auto i=1; i<49; i++){
		cells1[i].setNeighborMin(0, cells1[i-1]);
		cells1[i].setNeighborMax(0, cells1[i+1]);
	}

	double dt=1;
	double dx=c0;
	PMLParameterModel p(c0);
	int nPML = 10;
	for (auto i=0; i<nPML; i++){
		double x = double(i)/double(nPML);
		double xm = x + 0.5/double(nPML);
		cells1[i].setPMLParametersE(p.K(x), p.S(x), p.A(x), dt);
		cells1[49-i].setPMLParametersE(p.K(x), p.S(x), p.A(x), dt);
		cells1[i].setPMLParametersH(p.K(xm), p.S(xm), p.A(xm), dt);
		cells1[48-i].setPMLParametersH(p.K(xm), p.S(xm), p.A(xm), dt);

	}


	for (auto t=0; t<200; t++){
		std::for_each(++cells1.begin(), --cells1.end(), YeeUpdateD<TEM>(dt,dx));
		// if (t<19) cells1[25].Dz() = sin(2*pi*0.05*t);
		cells1[25].Dz() = exp(-(t-10)*(t-10)*0.05);
		std::for_each(++cells1.begin(), --cells1.end(), UpdatePMLD<TEM>(dt,dx));
		std::for_each(++cells1.begin(), --cells1.end(), ConstantUpdateE<TEM>(1));
	
		std::for_each(++cells1.begin(), --cells1.end(), YeeUpdateB<TEM>(dt,dx));
		std::for_each(++cells1.begin(), --cells1.end(), ConstantUpdateH<TEM>(1));
	
		// std::cout << "Ez:" ;
		for(auto i=0; i<50; i++) std::cout << ", " << cells1[i].Ez() ;
		std::cout << std::endl;
	}

	return 0;
}