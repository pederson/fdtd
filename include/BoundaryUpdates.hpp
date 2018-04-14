#ifndef _BOUNDARYUPDATES_H
#define _BOUNDARYUPDATES_H

#include "FDTDConstants.hpp"
#include "YeeUpdates.hpp"

namespace fdtd{


enum class Boundary : char {
	NONE,
	Periodic,
	BlochPeriodic,
	PEC,
	PMC,
	Parallel
};


template <Boundary b, class Mode, Dir d, Orientation o, template <typename> class FieldPolicy = GetField>
struct UpdateBoundaryD{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};

template <Boundary b, class Mode, Dir d, Orientation o, template <typename> class FieldPolicy = GetField>
struct UpdateBoundaryB{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// // specialization for 3D, periodic
// template <>
// struct UpdateBoundary<Periodic, ThreeD, Dir::X, Orientation::MIN>{
// 	double dt, dx;

// 	UpdateBoundary<Periodic, ThreeD, Dir::X>(double deltat, double deltax): dt(deltat), dx(deltax) {};

// 	template <class YeeCell>
// 	void operator()(YeeCell & f,
// 					double kx, double ky, double kz,
// 					double Lx, double Ly, double Lz){
// 		update(f, dt, dx, kx, ky, kz, Lx, Ly, Lz);
// 	};

// 	template <class YeeCell>
// 	static void update(YeeCell & f, double delt, double delx,
// 					   double kx, double ky, double kz,
// 					   double Lx, double Ly, double Lz){

// 		std::complex<double> ph = exp(1i*(kx*Lx + ky*Ly + kz*Lz));

// 		// x direction
// 		f.pmlHIxz() = f.pmlEBx()*f.pmlHIxz() + f.pmlECx()/delx*(f.Hz() - f.getNeighborMin(0).Hz());
// 		f.Dy() -= f.pmlHIxz()*delt;

// 		f.pmlHIxy() = f.pmlEBx()*f.pmlHIxy() + f.pmlECx()/delx*(f.Hy() - f.getNeighborMin(0).Hy());
// 		f.Dz() += f.pmlHIxy()*delt;
// 	};
// };





//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// Define difference operators for boundaries
template <typename EMField, Dir d, 
		  typename FieldGetter = GetField<EMField>, 
		  typename NeighborGetter = GetNeighbor<d, EMField::neighb_side>,
		  FieldType F = EMField::field_type,
		  Boundary b = Boundary::NONE,
		  typename ...Args>
struct BoundaryDifferenceOperator{
	static_assert(std::is_same<Field, EMField>::value, "BoundaryDifferenceOperator needs a valid EMField");
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// Define PEC difference operator


////////////// E //////////////
template <typename EMField, Dir d, typename FieldGetter, typename NeighborGetter>
struct BoundaryDifferenceOperator<EMField, d, FieldGetter, NeighborGetter, FieldType::Electric, Boundary::PEC>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {
		return 0.0;
	}
};


////////////// H //////////////
template <typename EMField, Dir d, typename FieldGetter, typename NeighborGetter>
struct BoundaryDifferenceOperator<EMField, d, FieldGetter, NeighborGetter, FieldType::Magnetic, Boundary::PEC>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {
		return 0.0;
	}
};

template <typename EMField, Dir d>
struct PECDifferenceOperatorTypedef{
	typedef BoundaryDifferenceOperator<EMField, d, 
							   GetField<EMField>, 
							   GetNeighbor<d, EMField::neighb_side>, 
							   EMField::field_type, 
							   Boundary::PEC> type;
};
template <typename EMField, Dir d>
using PECDifferenceOperator = typename PECDifferenceOperatorTypedef<EMField, d>::type;



// define the temporal policy for PEC
struct TemporalPEC : public TemporalScheme {
	static constexpr double curl_coeff = 0.0;

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get(YeeCell & f){
		return 0.0;
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void increment(YeeCell & f, ValueT hold){
		return;
	}
};


// define the D update for PEC boundary
template <class Mode, Dir d, Orientation o, template <typename> class FieldPolicy>
struct UpdateBoundaryD<Boundary::PEC, Mode, d, o, FieldPolicy>{
	double dt, dx;

	template <typename YeeCell>
	void operator()(YeeCell & f){
		YeeUpdateD<Mode, TemporalPEC, PMLCoeff<false>, FieldPolicy, PECDifferenceOperator>::update(f, dt, dx);
	}
};

// define the update for PEC boundary
template <class Mode, Dir d, Orientation o, template <typename> class FieldPolicy>
struct UpdateBoundaryB<Boundary::PEC, Mode, d, o, FieldPolicy>{
	double dt, dx;

	template <typename YeeCell>
	void operator()(YeeCell & f){
		YeeUpdateB<Mode, TemporalPEC, PMLCoeff<false>, FieldPolicy, PECDifferenceOperator>::update(f, dt, dx);
	}
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// Define PMC difference operator


// ////////////// E //////////////
// template <typename EMField, Dir d, typename FieldGetter, typename NeighborGetter>
// struct BoundaryDifferenceOperator<EMField, d, FieldGetter, NeighborGetter, FieldType::Electric, Boundary::PMC>{
// 	template <class YeeCell>
// 	static decltype(auto) get(YeeCell & f) {
// 		return 0.0;
// 	}
// };


// ////////////// H //////////////
// template <typename EMField, Dir d, typename FieldGetter, typename NeighborGetter>
// struct BoundaryDifferenceOperator<EMField, d, FieldGetter, NeighborGetter, FieldType::Magnetic, Boundary::PMC>{
// 	template <class YeeCell>
// 	static decltype(auto) get(YeeCell & f) {
// 		return 0.0;
// 	}
// };

// template <typename EMField, Dir d>
// struct PMCDifferenceOperatorTypedef{
// 	typedef BoundaryDifferenceOperator<EMField, d, 
// 							   GetField<EMField>, 
// 							   GetNeighbor<d, EMField::neighb_side>, 
// 							   EMField::field_type, 
// 							   Boundary::PMC> type;
// };
// template <typename EMField, Dir d>
// using PMCDifferenceOperator = typename PMCDifferenceOperatorTypedef<EMField, d>::type;



// // define the temporal policy for PMC
// struct TemporalPMC : public TemporalScheme {
// 	static constexpr double curl_coeff = 0.0;

// 	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
// 	static decltype(auto) get(YeeCell & f){
// 		return 0.0;
// 	}

// 	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
// 	static void increment(YeeCell & f, ValueT hold){
// 		return;
// 	}
// };


// // define the D update for PMC boundary
// template <class Mode, Dir d, template <typename> class FieldPolicy>
// struct UpdateBoundaryD<Boundary::PMC, Mode, d, Orientation::MIN, FieldPolicy>{
// 	double dt, dx;

// 	template <typename YeeCell>
// 	void operator()(YeeCell & f){
// 		// do nothing
// 	}
// };

// template <class Mode, Dir d, template <typename> class FieldPolicy>
// struct UpdateBoundaryD<Boundary::PMC, Mode, d, Orientation::MAX, FieldPolicy>{
// 	double dt, dx;

// 	template <typename YeeCell>
// 	void operator()(YeeCell & f){
// 		// do nothing
// 	}
// };

// // define the update for PMC boundary
// template <class Mode, Dir d, Orientation o, template <typename> class FieldPolicy>
// struct UpdateBoundaryB<Boundary::PMC, Mode, d, o, FieldPolicy>{
// 	double dt, dx;

// 	template <typename YeeCell>
// 	void operator()(YeeCell & f){
// 		YeeUpdateB<Mode, TemporalPMC, PMLCoeff<false>, FieldPolicy, PMCDifferenceOperator>::update(f, dt, dx);
// 	}
// };




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// Define periodic difference operator




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// Define Bloch-periodic difference operator


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// Define parallel difference operator


}// end namespace fdtd

#endif