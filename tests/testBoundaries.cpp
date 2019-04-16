// compile with:
// 			clang++-3.9 -std=c++14 -I../ testBoundaries.cpp -o testBoundaries


#include <fdtd.h>
#include "../include/BoundaryUpdates.hpp"

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>


using namespace fdtd;

typedef std::complex<double> 				scalar_type;
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



struct vector_periodic_fctor{
	typedef std::vector<cell_type> 				VectorType;
	typedef typename VectorType::iterator 		ItType;

	ItType 					mBegin;
	VectorType * 			mPeriodic;

	vector_periodic_fctor() {};
	vector_periodic_fctor(ItType beg, VectorType & v) : mPeriodic(&v), mBegin(beg) {};

	// const cell_type & operator()(ItType it) const {
	// 	return (*mPeriodic)[it - mBegin];
	// }

	cell_type & operator()(ItType it) {
		return (*mPeriodic)[it - mBegin];
	}
};


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
	BoundaryData bd_per = make_periodic_boundary<Mode, Dir::X, Orientation::MIN>(v.begin(), v.end(), vector_periodic_fctor(v.begin(), v_comm));
	bd_per.print_summary();

	// bloch-periodic boundary
	BoundaryData bd_bper = make_bloch_periodic_boundary<Mode, Dir::Z, Orientation::MAX>(v.begin(), v.end(), vector_periodic_fctor(v.begin(), v_comm), std::complex<double>(1.0, 1.0));
	bd_bper.print_summary();


	// parallel boundary

	return 0;
}