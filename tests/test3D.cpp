#include <fdtd.h>

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>

#include <cmath>


using namespace fdtd;

// compile with:
// 			clang++-3.8 -std=c++14 -stdlib=libc++ -I../ test3D.cpp -o test3D


int main(int argc, char * argv[]){

	// time-stepping parameters
	const double cfl = 0.99;
	const double dx = 1e-6;
	const double dt=cfl*dx/(c0*fdtd::sqrt(3.0));

	// YeeCell typedefs
	typedef YeeFields3D<double, std::array> yeefieldT;
	static bool pmlx = true;
	static bool pmly = true;
	static bool pmlz = true;
	typedef PML<ThreeD, true, true, true> pmlT;
	typedef YeeCell<yeefieldT, VacuumPolarization, VacuumMagnetization, Neighbor3, pmlT> YeeCellVac3D;
	
	// define the domain data structure and connect neighbors
	std::vector<std::vector<std::vector<YeeCellVac3D>>> cells1(50);
	for (auto i=0; i<50; i++) {
		cells1[i].resize(50);
		for (auto j=0; j<50; j++) cells1[i][j].resize(50);
	}
	for (auto i=1; i<49; i++){
		for (auto j=1; j<49; j++){
			for (auto k=1; k<49; k++){
				cells1[i][j][k].setNeighborMin(0, cells1[i-1][j][k]);
				cells1[i][j][k].setNeighborMax(0, cells1[i+1][j][k]);
				cells1[i][j][k].setNeighborMin(1, cells1[i][j-1][k]);
				cells1[i][j][k].setNeighborMax(1, cells1[i][j+1][k]);
				cells1[i][j][k].setNeighborMin(2, cells1[i][j][k-1]);
				cells1[i][j][k].setNeighborMax(2, cells1[i][j][k+1]);
			}
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
				for (auto k=0; k<50; k++){
					if (pmlx){
						cells1[i][j][k].setPMLParametersEx(p.K(x), p.S(x), p.A(x), dt);
						cells1[i][j][k].setPMLParametersHx(p.K(xm), p.S(xm), p.A(xm), dt);

						cells1[50-i][j][k].setPMLParametersEx(p.K(x), p.S(x), p.A(x), dt);
						cells1[49-i][j][k].setPMLParametersHx(p.K(xm), p.S(xm), p.A(xm), dt);
					}
					if (pmly){
						cells1[j][i][k].setPMLParametersEy(p.K(x), p.S(x), p.A(x), dt);
						cells1[j][i][k].setPMLParametersHy(p.K(xm), p.S(xm), p.A(xm), dt);

						cells1[j][50-i][k].setPMLParametersEy(p.K(x), p.S(x), p.A(x), dt);
						cells1[j][49-i][k].setPMLParametersHy(p.K(xm), p.S(xm), p.A(xm), dt);
					}
					if (pmlz){
						cells1[j][k][i].setPMLParametersEz(p.K(x), p.S(x), p.A(x), dt);
						cells1[j][k][i].setPMLParametersHz(p.K(xm), p.S(xm), p.A(xm), dt);

						cells1[j][k][50-i].setPMLParametersEz(p.K(x), p.S(x), p.A(x), dt);
						cells1[j][k][49-i].setPMLParametersHz(p.K(xm), p.S(xm), p.A(xm), dt);
					}
				}
			}
		}
	}
	// if (pmly){
	// 	for (auto i=1; i<nPML+2; i++){
	// 		double x = static_cast<double>(nPML-i+1)/static_cast<double>(nPML);
	// 		double xm = x - 0.5/static_cast<double>(nPML);
	// 		for (auto j=0; j<50; j++){
	// 			for (auto k=0; k<50; k++){
	// 				cells1[j][i][k].setPMLParametersEy(p.K(x), p.S(x), p.A(x), dt);
	// 				cells1[j][i][k].setPMLParametersHy(p.K(xm), p.S(xm), p.A(xm), dt);

	// 				cells1[j][50-i][k].setPMLParametersEy(p.K(x), p.S(x), p.A(x), dt);
	// 				cells1[j][49-i][k].setPMLParametersHy(p.K(xm), p.S(xm), p.A(xm), dt);
	// 			}
	// 		}
	// 	}
	// }
	// if (pmlz){
	// 	for (auto i=1; i<nPML+2; i++){
	// 		double x = static_cast<double>(nPML-i+1)/static_cast<double>(nPML);
	// 		double xm = x - 0.5/static_cast<double>(nPML);
	// 		for (auto j=0; j<50; j++){
	// 			for (auto k=0; k<50; k++){
	// 				cells1[j][k][i].setPMLParametersEz(p.K(x), p.S(x), p.A(x), dt);
	// 				cells1[j][k][i].setPMLParametersHz(p.K(xm), p.S(xm), p.A(xm), dt);

	// 				cells1[j][k][50-i].setPMLParametersEz(p.K(x), p.S(x), p.A(x), dt);
	// 				cells1[j][k][49-i].setPMLParametersHz(p.K(xm), p.S(xm), p.A(xm), dt);
	// 			}
	// 		}
	// 	}
	// }


	// // for(auto i=0; i<50; i++){
	// for(auto j=0; j<50; j++){
	// 		std::cout << " EKx: " << cells1[j][25].pmlEKx() ;
	// 		std::cout << " ESx: " << cells1[j][25].pmlESx() ;
	// 		std::cout << " EAx: " << cells1[j][25].pmlEAx() ;
	// 		std::cout << " EBx: " << cells1[j][25].pmlEBx() ;
	// 		std::cout << " ECx: " << cells1[j][25].pmlECx() ;
	// 		std::cout << " HKx: " << cells1[j][25].pmlHKx() ;
	// 		std::cout << " HSx: " << cells1[j][25].pmlHSx() ;
	// 		std::cout << " HAx: " << cells1[j][25].pmlHAx() ;
	// 		std::cout << " HBx: " << cells1[j][25].pmlHBx() ;
	// 		std::cout << " HCx: " << cells1[j][25].pmlHCx() ;
	// 		std::cout << std::endl;
	// }
	// 	for(auto j=0; j<50; j++){
	// 		std::cout << " EBy: " << cells1[25][j].pmlEBy() ;
	// 		std::cout << " ECy: " << cells1[25][j].pmlECy() ;
	// 		std::cout << " HBy: " << cells1[25][j].pmlHBy() ;
	// 		std::cout << " HCy: " << cells1[25][j].pmlHCy() ;
	// 		std::cout << std::endl;
	// 	}
	// // }
	// std::cout << std::endl;
	// throw -1;


	// std::cout << "am here" << std::endl;
	// start time-stepping
	for (auto t=0; t<200; t++){
		std::cout << "t: " << t << "\r" << std::flush;
		for (auto i=1; i<49; i++){
			for (auto j=1; j<49; j++) std::for_each(++cells1[i][j].begin(), --cells1[i][j].end(), YeeUpdateD<ThreeD>(dt,dx));
		}
		cells1[25][25][25].Dz() += sin(2*pi*c0/(20*dx)*static_cast<double>(t)*dt);
		// cells1[25][25][25].Dz() = exp(-(t-10)*(t-10)*dt*dt*c0/(20*dx));
		// cells1[25][25][26].Dz() = 0;
		// cells1[25][25][24].Dz() = 0;
		for (auto i=1; i<49; i++){
			for (auto j=1; j<49; j++) std::for_each(++cells1[i][j].begin(), --cells1[i][j].end(), UpdatePMLD<ThreeD>(dt,dx));
		}
		for (auto i=1; i<49; i++){
			for (auto j=1; j<49; j++) std::for_each(++cells1[i][j].begin(), --cells1[i][j].end(), ConstantUpdateE<ThreeD>(1));
		}
		// std::cout << "here" << std::endl;

		for (auto i=1; i<49; i++){
			for (auto j=1; j<49; j++) std::for_each(++cells1[i][j].begin(), --cells1[i][j].end(), YeeUpdateB<ThreeD>(dt,dx));
		}
		for (auto i=1; i<49; i++){
			for (auto j=1; j<49; j++) std::for_each(++cells1[i][j].begin(), --cells1[i][j].end(), UpdatePMLB<ThreeD>(dt,dx));
		}
		for (auto i=1; i<49; i++){
			for (auto j=1; j<49; j++) std::for_each(++cells1[i][j].begin(), --cells1[i][j].end(), ConstantUpdateH<ThreeD>(1));
		}

		// for(auto i=0; i<50; i++){
		// 	for(auto j=0; j<50; j++){
		// 		std::cout << ", " << cells1[i][j][25].Ez() ;
		// 	}
		// }
		// std::cout << std::endl;
	}

	return 0;
}