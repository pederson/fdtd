// compile with:
// 			clang++-3.9 -std=c++14 -I../ testBoundaries.cpp -o testBoundaries


#include <fdtd.h>
#include "../include/BoundaryUpdates.hpp"

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>


using namespace fdtd;

typedef double 				scalar_type;
typedef ThreeD 				Mode;
typedef YeeCell<YeeFields<Mode, scalar_type, std::array>,
				VacuumPolarization<scalar_type>,
				VacuumMagnetization<scalar_type>
				>			yee_type;

struct cell_with_neighb : public yee_type{
	decltype(auto) getNeighborMin(Dir d){
		return *this;
	};
};

typedef cell_with_neighb	cell_type;

int main(int argc, char * argv[]){

	// vector of cells that represents a boundary
	std::vector<cell_type> v(50);
	// std::vector<cell_type> ghost(50);

	BoundaryData bd = make_pec_boundary<Mode, Dir::X, Orientation::MIN>(v.begin(), v.end());


	return 0;
}