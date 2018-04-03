#ifndef _YEEUPDATES_H
#define _YEEUPDATES_H

#include "FDTDConstants.hpp"
#include "DefaultInterfaces.hpp"

namespace fdtd{




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// Define difference operators on YeeCell objects
// which will be used to construct the YeeAlgorithm for
// any given mode
template <typename EMField, Dir d, 
		  typename FieldGetter = GetField<EMField>, 
		  typename  NeighborGetter = GetNeighbor<d, EMField::neighb_side>,
		  FieldType F = EMField::field_type>
struct DifferenceOperator{
	static_assert(std::is_same<Field, EMField>::value, "DifferenceOperator needs a valid EMField");


};


////////////// E //////////////
template <typename EMField, Dir d, typename FieldGetter, typename NeighborGetter>
struct DifferenceOperator<EMField, d, FieldGetter, NeighborGetter, FieldType::Electric>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {
		return FieldGetter::get(NeighborGetter::get(f)) - FieldGetter::get(f);
	}
};


////////////// H //////////////
template <typename EMField, Dir d, typename FieldGetter, typename NeighborGetter>
struct DifferenceOperator<EMField, d, FieldGetter, NeighborGetter, FieldType::Magnetic>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {
		return FieldGetter::get(f) - FieldGetter::get(NeighborGetter::get(f));
	}
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



struct TemporalScheme {};

struct BD1 : public TemporalScheme {
	static constexpr double curl_coeff = 1.0;

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get(YeeCell & f){
		return FieldGetter::get(f);
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void increment(YeeCell & f, ValueT hold){
		return;
	}
};


struct BD3 : public TemporalScheme {
	static constexpr double curl_coeff = 24.0/23.0;

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get(YeeCell & f){
		return 21.0/23.0*FieldGetter::get(f) 
			  +3.0/23.0*FieldGetter::get(f.BD(0))
			  -1.0/23.0*FieldGetter::get(f.BD(1));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void increment(YeeCell & f, ValueT hold){
		FieldGetter::get(f.BD(1)) = FieldGetter::get(f.BD(0));
		FieldGetter::get(f.BD(0)) = hold;
		return;
	}
};


struct BD4 : public TemporalScheme {
	static constexpr double curl_coeff = 24.0/22.0;

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get(YeeCell & f){
		return 17.0/22.0*FieldGetter::get(f) 
			  +9.0/22.0*FieldGetter::get(f.BD(0))
			  -5.0/22.0*FieldGetter::get(f.BD(1))
			  +1.0/22.0*FieldGetter::get(f.BD(2));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void increment(YeeCell & f, ValueT hold){
		FieldGetter::get(f.BD(2)) = FieldGetter::get(f.BD(1));
		FieldGetter::get(f.BD(1)) = FieldGetter::get(f.BD(0));
		FieldGetter::get(f.BD(0)) = hold;
		return;
	}
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************





template <class Mode>
struct NullFieldUpdate{
	static_assert(std::is_same<EMMode, Mode>::value, "NullUpdate needs a valid Mode");
	
	double dt, dx;

	NullFieldUpdate(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************





// these updates require PML values available
// i.e.  pmlEKx(), pmlEKy(), and pmlEKz()
template <class Mode, class TimePolicy = BD1, class PMLCoeffPolicy = PMLCoeff<true>>
struct YeeUpdateD{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	static_assert(std::is_same<TemporalScheme, TimePolicy>::value, "YeeUpdate needs a valid TemporalScheme");
};


// specialization for 3D
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateD<ThreeD, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateD<ThreeD, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hDx = f.Dx();
		f.Dx() = TimePolicy::template get<fdtd::Dx>
				+TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlEKy(f)*fdtd::DifferenceOperator<Hz, Dir::Y>::get(f)
											 - 1.0/PMLCoeffPolicy::pmlEKz(f)*fdtd::DifferenceOperator<Hy, Dir::Z>::get(f));
		auto hDy = f.Dy();
		f.Dy() = TimePolicy::template get<fdtd::Dy>
				+TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlEKz(f)*fdtd::DifferenceOperator<Hx, Dir::Z>::get(f)
											 - 1.0/PMLCoeffPolicy::pmlEKx(f)*fdtd::DifferenceOperator<Hz, Dir::X>::get(f));
		auto hDz = f.Dz();
		f.Dz() = TimePolicy::template get<fdtd::Dz>
				+TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlEKx(f)*fdtd::DifferenceOperator<Hy, Dir::X>::get(f)
											 - 1.0/PMLCoeffPolicy::pmlEKy(f)*fdtd::DifferenceOperator<Hx, Dir::Y>::get(f));
	
		TimePolicy::template increment<fdtd::Dx>(f, hDx);
		TimePolicy::template increment<fdtd::Dy>(f, hDy);
		TimePolicy::template increment<fdtd::Dz>(f, hDz);
	};
};

// specialization for TE
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateD<TE, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateD<TE, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hDx = f.Dx();
		f.Dx() = TimePolicy::template get<fdtd::Dx>
				+TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlEKy(f)*fdtd::DifferenceOperator<Hz, Dir::Y>::get(f));
		auto hDy = f.Dy();
		f.Dy() = TimePolicy::template get<fdtd::Dy>
				+TimePolicy::curl_coeff*dt/dx*(-1.0/PMLCoeffPolicy::pmlEKx(f)*fdtd::DifferenceOperator<Hz, Dir::X>::get(f));
	
		TimePolicy::template increment<fdtd::Dx>(f, hDx);
		TimePolicy::template increment<fdtd::Dy>(f, hDy);
	};
};


// specialization for TM
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateD<TM, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateD<TM, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hDz = f.Dz();
		f.Dz() = TimePolicy::template get<fdtd::Dz>(f) 
			   + TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlEKx(f)*fdtd::DifferenceOperator<Hy, Dir::X>::get(f)
					  						 - 1.0/PMLCoeffPolicy::pmlEKy(f)*fdtd::DifferenceOperator<Hx, Dir::Y>::get(f));
		
		TimePolicy::template increment<fdtd::Dz>(f, hDz);
	};
};



// specialization for TEM
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateD<TEM, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateD<TEM, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hDz = f.Dz();
		f.Dz() = TimePolicy::template get<fdtd::Dz>
				+TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlEKx(f)*fdtd::DifferenceOperator<Hy, Dir::X>::get(f));

		TimePolicy::template increment<fdtd::Dz>(f, hDz);
	};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <class Mode, class TimePolicy = BD1, class PMLCoeffPolicy = PMLCoeff<true>>
struct YeeUpdateB{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	static_assert(std::is_same<TemporalScheme, TimePolicy>::value, "YeeUpdate needs a valid TemporalScheme");
};



// specialization for 3D
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateB<ThreeD, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateB<ThreeD, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hBx = f.Bx();
		f.Bx() = TimePolicy::template get<fdtd::Bx>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKy(f)*DifferenceOperator<fdtd::Ez, Dir::Y>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKz(f)*DifferenceOperator<fdtd::Ey, Dir::Z>::get(f));
		auto hBy = f.By();
		f.By() = TimePolicy::template get<fdtd::By>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKz(f)*DifferenceOperator<fdtd::Ex, Dir::Z>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKx(f)*DifferenceOperator<fdtd::Ez, Dir::X>::get(f));
		auto hBz = f.Bz();
		f.Bz() = TimePolicy::template get<fdtd::Bz>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKx(f)*DifferenceOperator<fdtd::Ey, Dir::X>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKy(f)*DifferenceOperator<fdtd::Ex, Dir::Y>::get(f));

		TimePolicy::template increment<fdtd::Bx>(f, hBx);
		TimePolicy::template increment<fdtd::By>(f, hBy);
		TimePolicy::template increment<fdtd::Bz>(f, hBz);

	};
};


// specialization for TE
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateB<TE, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateB<TE, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hBz = f.Bz();
		f.Bz() = TimePolicy::template get<fdtd::Bz>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKx(f)*DifferenceOperator<fdtd::Ey, Dir::X>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKy(f)*DifferenceOperator<fdtd::Ex, Dir::Y>::get(f));

		TimePolicy::template increment<fdtd::Bz>(f, hBz);
	};
};


// specialization for TM
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateB<TM, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateB<TM, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hBx = f.Bx();
		f.Bx() = TimePolicy::template get<fdtd::Bx>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKy(f)*DifferenceOperator<fdtd::Ez, Dir::Y>::get(f));
		auto hBy = f.By();
		f.By() = TimePolicy::template get<fdtd::By>(f)
				-TimePolicy::curl_coeff*dt/dx*(-1.0/PMLCoeffPolicy::pmlHKx(f)*DifferenceOperator<fdtd::Ez, Dir::X>::get(f));

		TimePolicy::template increment<fdtd::Bx>(f, hBx);
		TimePolicy::template increment<fdtd::By>(f, hBy);
	};
};


// specialization for TEM
template <typename TimePolicy, typename PMLCoeffPolicy>
struct YeeUpdateB<TEM, TimePolicy, PMLCoeffPolicy>{
	double dt, dx;

	YeeUpdateB<TEM, TimePolicy, PMLCoeffPolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hBy = f.By();
		f.By() = TimePolicy::template get<fdtd::By>(f)
				-TimePolicy::curl_coeff*dt/dx*(-1.0/PMLCoeffPolicy::pmlHKx(f)*DifferenceOperator<fdtd::Ez, Dir::X>::get(f));
		TimePolicy::template increment<fdtd::By>(f, hBy);
	};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <typename FieldType>
inline void updateDx3D(FieldType & Dx, const FieldType & Hy, const FieldType & Hz,
					 const FieldType & Hy_ZMin, const FieldType & Hz_YMin,
					 double dt, double dx){
	Dx += dt/dx*(Hz - Hz_YMin - Hy + Hy_ZMin);
}

template <typename FieldType>
inline void updateDy3D(FieldType & Dy, const FieldType & Hx, const FieldType & Hz,
					 const FieldType & Hx_ZMin, const FieldType & Hz_XMin,
					 double dt, double dx){
	Dy += dt/dx*(Hx - Hx_ZMin - Hz + Hz_XMin);
}

template <typename FieldType>
inline void updateDz3D(FieldType & Dz, const FieldType & Hx, const FieldType & Hy,
					 const FieldType & Hx_YMin, const FieldType & Hy_XMin,
					 double dt, double dx){
	Dz += dt/dx*(Hy - Hy_XMin - Hx + Hx_YMin);
}





}// end namespace fdtd

#endif