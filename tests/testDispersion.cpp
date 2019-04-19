// compile with:
// 			clang++-3.9 -std=c++14 -I../ testDispersion.cpp -o testDispersion

// This test file is really just to test the mechanics of the
// dispersion updaters, and is NOT a test for correctness


#include <fdtd.h>
#include "../include/DispersiveMaterials.hpp"

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

	decltype(auto) getNeighborMax(Dir d){
		return *this;
	};
};



typedef cell_with_neighb	cell_type;


int main(int argc, char * argv[]){

	int ncells = 5;

	// vector of cells that represents a boundary
	std::vector<cell_type> v(ncells);
	for (auto it=v.begin(); it!=v.end(); it++){
		(*it).Ex() = it - v.begin();
		(*it).Ey() = it - v.begin() + 1;
		(*it).Ez() = it - v.begin() + 2;
	}
	// std::vector<cell_type> ghost(10);



	std::cout << "******************************** CONSTANT DISPERSION " << std::endl;
	std::cout << "******************************** ----------------- " << std::endl;
	matmap["constant"] = make_material_pair<cell_type>(ConstantUpdate<Mode, FieldType::Electric>(4.0), ConstantUpdate<Mode, FieldType::Magnetic>(2.0));

	// BoundaryData bd_pec = make_pec_boundary<Mode, Dir::X, Orientation::MIN>(v.begin(), v.end());
	// bd_pec.print_summary();

	// std::cout << "******************************** CONSTANT DISPERSION " << std::endl;
	// std::cout << "******************************** ----------------- " << std::endl;
	// BoundaryData bd_pmc = make_pmc_boundary<Mode, Dir::Y, Orientation::MIN>(v.begin(), v.end());
	// bd_pmc.print_summary();
	// bd_pmc.updateE();

	// std::cout << "******************************** CONSTANT DISPERSION " << std::endl;
	// std::cout << "******************************** ----------------- " << std::endl;
	// BoundaryData bd_symm = make_symmetric_boundary<Mode, Dir::X, Orientation::MAX>(v.begin(), v.end());
	// bd_symm.print_summary();

	// std::cout << "******************************** CONSTANT DISPERSION " << std::endl;
	// std::cout << "******************************** ----------------- " << std::endl;
	// BoundaryData bd_asymm = make_antisymmetric_boundary<Mode, Dir::Y, Orientation::MAX>(v.begin(), v.end());
	// bd_asymm.print_summary();


	// // periodic boundary
	// std::cout << "******************************** CONSTANT DISPERSION " << std::endl;
	// std::cout << "******************************** ----------------- " << std::endl;
	// std::vector<cell_type> v_comm(ncells);
	// for (auto it=v_comm.begin(); it!=v_comm.end(); it++) {
	// 	(*it).Ex() = 1.0;
	// 	(*it).Ey() = 1.0;
	// 	(*it).Ez() = 1.0;
	// }
	// BoundaryData bd_per = make_periodic_boundary<Mode, Dir::X, Orientation::MIN>(v.begin(), v.end(), vector_periodic_fctor(v.begin(), v_comm));
	// bd_per.print_summary();



	return 0;
}