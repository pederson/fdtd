// compile with:
// 			clang++-3.9 -std=c++14 -I../ testBoundaries.cpp -o testBoundaries

#include <mpi.h>

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

	int ncells = 5;

	// vector of cells that represents a boundary
	std::vector<cell_type> v(ncells);
	for (auto it=v.begin(); it!=v.end(); it++){
		(*it).Ex() = it - v.begin();
		(*it).Ey() = it - v.begin() + 1;
		(*it).Ez() = it - v.begin() + 2;
	}
	// std::vector<cell_type> ghost(10);

	std::cout << "******************************** PEC BOUNDARY " << std::endl;
	std::cout << "******************************** ----------------- " << std::endl;
	BoundaryData bd_pec = make_pec_boundary<Mode, Dir::X, Orientation::MIN>(v.begin(), v.end());
	bd_pec.print_summary();

	std::cout << "******************************** PMC BOUNDARY " << std::endl;
	std::cout << "******************************** ----------------- " << std::endl;
	BoundaryData bd_pmc = make_pmc_boundary<Mode, Dir::Y, Orientation::MIN>(v.begin(), v.end());
	bd_pmc.print_summary();
	bd_pmc.updateE();

	std::cout << "******************************** SYMMETRIC BOUNDARY " << std::endl;
	std::cout << "******************************** ----------------- " << std::endl;
	BoundaryData bd_symm = make_symmetric_boundary<Mode, Dir::X, Orientation::MAX>(v.begin(), v.end());
	bd_symm.print_summary();

	std::cout << "******************************** ANTISYMMETRIC BOUNDARY " << std::endl;
	std::cout << "******************************** ----------------- " << std::endl;
	BoundaryData bd_asymm = make_antisymmetric_boundary<Mode, Dir::Y, Orientation::MAX>(v.begin(), v.end());
	bd_asymm.print_summary();


	// periodic boundary
	std::cout << "******************************** PERIODIC BOUNDARY " << std::endl;
	std::cout << "******************************** ----------------- " << std::endl;
	std::vector<cell_type> v_comm(ncells);
	for (auto it=v_comm.begin(); it!=v_comm.end(); it++) {
		(*it).Ex() = 1.0;
		(*it).Ey() = 1.0;
		(*it).Ez() = 1.0;
	}
	BoundaryData bd_per = make_periodic_boundary<Mode, Dir::X, Orientation::MIN>(v.begin(), v.end(), vector_periodic_fctor(v.begin(), v_comm));
	bd_per.print_summary();

	// // bloch-periodic boundary
	// std::cout << "******************************** BLOCH-PERIODIC BOUNDARY " << std::endl;
	// std::cout << "******************************** ----------------- " << std::endl;
	// BoundaryData bd_bper = make_bloch_periodic_boundary<Mode, Dir::Z, Orientation::MAX>(v.begin(), v.end(), vector_periodic_fctor(v.begin(), v_comm), std::complex<double>(1.0, 1.0));
	// bd_bper.print_summary();
	// // bd_bper.updateE();


	// parallel boundary
	std::cout << "******************************** PARALLEL BOUNDARY " << std::endl;
	std::cout << "******************************** ----------------- " << std::endl;
	MPI_Init(&argc, &argv);
	int rank, neighb;

	std::cout << "initialized" << std::endl;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::cout << "got rank" << std::endl;

	if (rank == 0) 	neighb = 1;
	else 			neighb = 0;


	for (auto it=v.begin(); it!=v.end(); it++){
		(*it).Ex() = rank;
		(*it).Ey() = rank;
	}

	std::cout << "I'm rank " << rank << " and neighbor is " << neighb << std::endl;
	BoundaryData bd_para = make_parallel_boundary<Mode, Dir::X, Orientation::MAX>(v.begin(), v.end(), neighb);
	bd_para.print_summary();
	bd_para.updateE();

	std::cout << "I'm rank " << rank << " and neighbor is " << neighb << std::endl;
	for (auto it=v.begin(); it!=v.end(); it++){
		std::cout << (*it).Ex() << std::endl;
	}

	MPI_Finalize();

	return 0;
}