#include <fdtd.h>

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>


using namespace fdtd;

// compile with:
// 			clang++-3.8 -std=c++14 -stdlib=libc++ -I../ testTM.cpp -o test


int main(int argc, char * argv[]){

	// time-stepping parameters
	double cfl = 0.9;
	double dx = 1e-6;
	double dt=cfl*dx/(c0*fdtd::sqrt(2.0));

	// std::cout << typeid(y).name() << std::endl;

	// YeeCell typedefs
	typedef YeeFieldsTM<double, std::array> yeefieldT;
	static bool pmlx = true;
	static bool pmly = true;
	static bool pmlz = false;
	typedef PML<TM, true, true, false> pmlT;
	typedef YeeCell<yeefieldT, VacuumPolarization, VacuumMagnetization, Neighbor2, pmlT> YeeCellVacTM;
	
	// define the domain data structure and connect neighbors
	std::vector<std::vector<YeeCellVacTM>> cells1(50);
	for (auto i=0; i<50; i++) cells1[i].resize(50);
	for (auto i=1; i<49; i++){
		for (auto j=1; j<49; j++){
			cells1[i][j].setNeighborMin(0, cells1[i-1][j]);
			cells1[i][j].setNeighborMax(0, cells1[i+1][j]);
			cells1[i][j].setNeighborMin(1, cells1[i][j-1]);
			cells1[i][j].setNeighborMax(1, cells1[i][j+1]);
		}
	}

	
	// setup the PML parameters
	PMLParameterModel p(dx);
	int nPML = 10;
	if (pmlx){
		for (auto i=1; i<nPML+2; i++){
			double x = static_cast<double>(nPML-i+1)/static_cast<double>(nPML);
			double xm = x - 0.5/static_cast<double>(nPML);
			for (auto j=0; j<50; j++){
				cells1[i][j].setPMLParametersEx(p.K(x), p.S(x), p.A(x), dt);
				cells1[i][j].setPMLParametersHx(p.K(xm), p.S(xm), p.A(xm), dt);

				cells1[50-i][j].setPMLParametersEx(p.K(x), p.S(x), p.A(x), dt);
				cells1[49-i][j].setPMLParametersHx(p.K(xm), p.S(xm), p.A(xm), dt);
			}
		}
	}
	if (pmly){
		for (auto i=1; i<nPML+2; i++){
			double x = static_cast<double>(nPML-i+1)/static_cast<double>(nPML);
			double xm = x - 0.5/static_cast<double>(nPML);
			for (auto j=0; j<50; j++){
				cells1[j][i].setPMLParametersEy(p.K(x), p.S(x), p.A(x), dt);
				cells1[j][i].setPMLParametersHy(p.K(xm), p.S(xm), p.A(xm), dt);

				cells1[j][50-i].setPMLParametersEy(p.K(x), p.S(x), p.A(x), dt);
				cells1[j][49-i].setPMLParametersHy(p.K(xm), p.S(xm), p.A(xm), dt);
			}
		}
	}


	for(auto i=0; i<50; i++){
		for(auto j=0; j<50; j++){
			std::cout << ", " << cells1[i][j].pmlHCy() ;
		}
	}
	std::cout << std::endl;
	throw -1;


	// std::cout << "am here" << std::endl;
	// start time-stepping
	for (auto t=0; t<200; t++){
		for (auto i=1; i<49; i++) std::for_each(++cells1[i].begin(), --cells1[i].end(), YeeUpdateD<TM>(dt,dx));
		// if (t<19) cells1[25].Dz() = sin(2*pi*0.05*t);
		cells1[25][25].Dz() += exp(-(t-10)*(t-10)*0.05);
		for (auto i=1; i<49; i++) std::for_each(++cells1[i].begin(), --cells1[i].end(), UpdatePMLD<TM>(dt,dx));
		for (auto i=1; i<49; i++) std::for_each(++cells1[i].begin(), --cells1[i].end(), ConstantUpdateE<TM>(1));
	
		// std::cout << "here" << std::endl;

		for (auto i=1; i<49; i++) std::for_each(++cells1[i].begin(), --cells1[i].end(), YeeUpdateB<TM>(dt,dx));
		for (auto i=1; i<49; i++) std::for_each(++cells1[i].begin(), --cells1[i].end(), UpdatePMLB<TM>(dt,dx));
		for (auto i=1; i<49; i++) std::for_each(++cells1[i].begin(), --cells1[i].end(), ConstantUpdateH<TM>(1));
	
		// std::cout << "here1" << std::endl;
		// std::cout << "Ez:" ;
		for(auto i=0; i<50; i++){
			for(auto j=0; j<50; j++){
				std::cout << ", " << cells1[i][j].Ez() ;
			}
		}
		std::cout << std::endl;
	}

	return 0;
}