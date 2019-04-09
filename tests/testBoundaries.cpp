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

	BoundaryData bd_pec = make_pec_boundary<Mode, Dir::X, Orientation::MIN>(v.begin(), v.end());
	bd_pec.print_summary();

	BoundaryData bd_pmc = make_pmc_boundary<Mode, Dir::Y, Orientation::MIN>(v.begin(), v.end());
	bd_pmc.print_summary();

	BoundaryData bd_symm = make_symmetric_boundary<Mode, Dir::X, Orientation::MAX>(v.begin(), v.end());
	bd_symm.print_summary();

	BoundaryData bd_asymm = make_antisymmetric_boundary<Mode, Dir::Y, Orientation::MAX>(v.begin(), v.end());
	bd_asymm.print_summary();


	// periodic boundary
	std::vector<cell_type> v_comm(50);

	// bloch-periodic boundary


	// parallel boundary

	return 0;
}