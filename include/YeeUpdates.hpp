#ifndef _YEEUPDATES_H
#define _YEEUPDATES_H

#include "FDTDConstants.hpp"
#include "DefaultInterfaces.hpp"

namespace fdtd{


template <typename FieldGetter,
		  typename RightGetter,
		  typename LeftGetter>
struct DifferenceOp{
	template <class Object>
	static decltype(auto) get(Object && o){
		return (FieldGetter::get(RightGetter::get(o)) - FieldGetter::get(LeftGetter::get(o)));
	}
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************





template <Dir d, typename FieldGetter, 
		  template<Dir, Orientation> class NeighborGetter = GetNeighbor>
struct CentralDifferenceTypedef{
	typedef DifferenceOp<FieldGetter, 
					     NeighborGetter<d, Orientation::MAX>, 
					     NeighborGetter<d, Orientation::MIN>> type;
};

template <Dir d, typename FieldGetter, template<Dir, Orientation> class NeighborGetter = GetNeighbor>
using CentralDifferenceOperator = typename CentralDifferenceTypedef<d, FieldGetter, NeighborGetter>::type;


template <Dir d, typename FieldGetter, 
		  template<Dir, Orientation> class NeighborGetter = GetNeighbor>
struct ForwardDifferenceTypedef{
	typedef DifferenceOp<FieldGetter, 
					     NeighborGetter<d, Orientation::MAX>, 
					     GetSelf> type;
};

template <Dir d, typename FieldGetter, template<Dir, Orientation> class NeighborGetter = GetNeighbor>
using ForwardDifferenceOperator = typename ForwardDifferenceTypedef<d, FieldGetter, NeighborGetter>::type;



template <Dir d, typename FieldGetter, 
		  template<Dir, Orientation> class NeighborGetter = GetNeighbor>
struct BackwardDifferenceTypedef{
	typedef DifferenceOp<FieldGetter, 
					     GetSelf, 
					     NeighborGetter<d, Orientation::MIN>> type;
};

template <Dir d, typename FieldGetter, template<Dir, Orientation> class NeighborGetter = GetNeighbor>
using BackwardDifferenceOperator = typename BackwardDifferenceTypedef<d, FieldGetter, NeighborGetter>::type;




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
		  FieldType F = EMField::field_type,
		  typename ...Args>
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


template <typename EMField, Dir d>
struct DefaultDifferenceOperatorTypedef{
	typedef DifferenceOperator<EMField, d, 
							   GetField<EMField>, 
							   GetNeighbor<d, EMField::neighb_side>,
							   EMField::field_type> 	type;
};
template <typename EMField, Dir d>
using DefaultDifferenceOperator = typename DefaultDifferenceOperatorTypedef<EMField, d>::type;





//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// Define components of the curl operator
template <Dir I, Dir J, 
		  typename FieldGetter,
		  typename DeltaXPolicy, 
		  template<Dir, typename> class DifferenceOperator>
struct CurlOperator{
	// Diagonal components are zero 0
	template <typename Object>
	static constexpr std::enable_if_t<I==J, decltype(FieldGetter::get(std::declval<Object>()))>
	get(Object && o) {return 0.0;};

	template <typename Object>
	static constexpr std::enable_if_t<I!=J, decltype(FieldGetter::get(std::declval<Object>()))>
	get(Object && o) {
		return LeviCivita<MutuallyOrthogonal<I,J>::value,I,J>::value*DifferenceOperator<MutuallyOrthogonal<I,J>::value,FieldGetter>::get(o)/DeltaXPolicy::get();
	};
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// template <typename EMField, Dir I, Dir J>
// struct FDTDCurlOperatorTypedef
// {
// 	typedef CurlOperator<I,J,GetField<EMField>,PMLDeltaX,DefaultDifferenceOperator<EMField, I>> type;
// };
// template <typename EMField, Dir I, Dir J>
// using FDTDCurlOperator = typename FDTDCurlOperatorTypedef<EMField, I, J>::type;


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



struct TemporalScheme {};

struct BD1 : public TemporalScheme {
	static constexpr double curl_coeff = 1.0;
	static constexpr double last_coeff = 1.0;

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get(YeeCell & f){
		return FieldGetter::get(f);
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void increment(YeeCell && f, ValueT hold){
		return;
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void decrement(YeeCell && f, ValueT hold){
		return;
	}
};


struct BD3 : public TemporalScheme {
	static constexpr double curl_coeff = 24.0/23.0;
	static constexpr double last_coeff = -1.0/23.0;

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get(const YeeCell & f){
		return 21.0/23.0*FieldGetter::get(f) 
			  +3.0/23.0*FieldGetter::get(f.BD(0))
			  -1.0/23.0*FieldGetter::get(f.BD(1));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get_reverse(const YeeCell & f){
		return +21.0*FieldGetter::get(f.BD(0)) 
			  +3.0*FieldGetter::get(f.BD(1))
			  -23.0*FieldGetter::get(f);
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get_last(const YeeCell & f){
		return FieldGetter::get(f.BD(1));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get_first(const YeeCell & f){
		return FieldGetter::get(f.BD(0));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void increment(YeeCell && f, ValueT hold){
		FieldGetter::get(f.BD(1)) = FieldGetter::get(f.BD(0));
		FieldGetter::get(f.BD(0)) = hold;
		return;
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void decrement(YeeCell && f, ValueT hold){
		FieldGetter::get(f.BD(0)) = FieldGetter::get(f.BD(1));
		FieldGetter::get(f.BD(1)) = hold;
		return;
	}
};


struct BD4 : public TemporalScheme {
	static constexpr double curl_coeff = 24.0/22.0;
	static constexpr double last_coeff = 1.0/22.0;

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get(YeeCell & f){
		return 17.0/22.0*FieldGetter::get(f) 
			  +9.0/22.0*FieldGetter::get(f.BD(0))
			  -5.0/22.0*FieldGetter::get(f.BD(1))
			  +1.0/22.0*FieldGetter::get(f.BD(2));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get_reverse(YeeCell & f){
		return -17.0*FieldGetter::get(f.BD(0)) 
			  -9.0*FieldGetter::get(f.BD(1))
			  +5.0*FieldGetter::get(f.BD(2))
			  -22.0*FieldGetter::get(f);
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get_last(const YeeCell & f){
		return FieldGetter::get(f.BD(2));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell>
	static decltype(auto) get_first(const YeeCell & f){
		return FieldGetter::get(f.BD(0));
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void increment(YeeCell && f, ValueT hold){
		FieldGetter::get(f.BD(2)) = FieldGetter::get(f.BD(1));
		FieldGetter::get(f.BD(1)) = FieldGetter::get(f.BD(0));
		FieldGetter::get(f.BD(0)) = hold;
		return;
	}

	template <typename EMField, typename FieldGetter = GetField<EMField>, typename YeeCell, typename ValueT>
	static void decrement(YeeCell && f, ValueT hold){
		FieldGetter::get(f.BD(0)) = FieldGetter::get(f.BD(1));
		FieldGetter::get(f.BD(1)) = FieldGetter::get(f.BD(2));
		FieldGetter::get(f.BD(2)) = hold;
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


namespace yu_details{
	template <FieldType F>
	struct Coeff{static constexpr double value = 1.0;};
	template<> struct Coeff<FieldType::Magnetic>{static constexpr double value = -1.0;};

	template <FieldType F>
	struct CurlType{static constexpr FieldType value = FieldType::Magnetic;};
	template<> struct CurlType<FieldType::Magnetic>{static constexpr FieldType value = FieldType::Electric;};
}

template <class Mode, 
		  class EMField,
		  class TimePolicy 			= BD1, 
		  class PMLCoeffPolicy 		= GetPML,
		  template <typename> class FieldPolicy = GetField,
		  template <typename,Dir> class DifferencePolicy = DefaultDifferenceOperator>
struct YeeUpdate : public TimePolicy{
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	static_assert(std::is_base_of<TemporalScheme, TimePolicy>::value, "YeeUpdate needs a valid TemporalScheme");
private:
	double dt, dx;

	typedef FieldPolicy<EMField> 	FP;
	typedef TimePolicy 				TP;
	typedef PMLCoeffPolicy 			PML;
	// typedef DifferencePolicy 		DP;
	static constexpr Dir 			I = FieldDir<EMField>::value;
	static constexpr FieldType 		FT = EMField::field_type;

	// this is the curl update for one direction
	template <Dir d>
	struct atomic_update{
		static constexpr Dir 			J = d;
		static constexpr Dir 			K = MutuallyOrthogonal<I,J>::value;
		static constexpr FieldType 		CurlType =  yu_details::CurlType<FT>::value;
		typedef typename FieldComponent<CurlType, K>::type 	CurlField;
		typedef std::conditional_t<CurlType == FieldType::Electric,
								   typename FieldComponents<Mode>::electric,
								   typename FieldComponents<Mode>::magnetic> 	tuple_type;

		template <typename YeeCell>
		static constexpr std::enable_if_t<Detail::tuple_contains_type<CurlField, tuple_type>::value, std::remove_reference_t<decltype(FP::get(std::declval<YeeCell>()))>>
		get(YeeCell && f, double delx){
			return LeviCivita<I, J, K>::value*GetPML::F<FT, J>(f)/(delx)*DifferencePolicy<CurlField, J>::get(f);
		}

		template <typename YeeCell>
		static constexpr std::enable_if_t<!Detail::tuple_contains_type<CurlField, tuple_type>::value, std::remove_reference_t<decltype(FP::get(std::declval<YeeCell>()))>>
		get(YeeCell && f, double delx){
			return 0.0;
		}
	};

public:
	YeeUpdate(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell && f){
		update_all(f, dt, dx);
	};

	// update with the curl term given, with increment of the time policy
	template <class YeeCell>
	static void	update_all(YeeCell && f, double delt, double delx){
		static constexpr Dir d1 = (I == Dir::X ? Dir::Y : (I == Dir::Y ? Dir::Z : Dir::X));
		static constexpr Dir d2 = MutuallyOrthogonal<I, d1>::value;

		auto hold = FP::get(f);
		FP::get(f) = TimePolicy::template get<EMField, FP>(f)
					+yu_details::Coeff<FT>::value*TimePolicy::curl_coeff*delt*(atomic_update<d1>::get(f, delx) + atomic_update<d2>::get(f, delx));

		TimePolicy::template increment<EMField, FP>(f, hold);
	}	

	// update with the curl term given, with increment of the time policy
	template <class YeeCell, typename ValueType>
	static void	update(YeeCell && f, double delt, ValueType curl){
		auto hold = FP::get(f);
		FP::get(f) = TimePolicy::template get<EMField, FP>(f)
					+yu_details::Coeff<FT>::value*TimePolicy::curl_coeff*delt*curl;

		TimePolicy::template increment<EMField, FP>(f, hold);
	}	

	// update in a single direction, with increment of the time policy
	template <Dir d, class YeeCell>
	static void	update(YeeCell && f, double delt, double delx){
		static constexpr Dir 			J = d;
		static constexpr Dir 			K = MutuallyOrthogonal<I,J>::value;

		auto hold = FP::get(f);
		FP::get(f) = TimePolicy::template get<EMField, FP>(f)
					+PML::F<FT, K>(f)*yu_details::Coeff<FT>::value*LeviCivita<I, J, K>::value*TimePolicy::curl_coeff*delt/delx*DifferencePolicy<typename FieldComponent<yu_details::CurlType<FT>::value, K>::type, J>::get(f);

		TimePolicy::template increment<EMField, FP>(f, hold);
	}

	// update with left and right values provided as input
	template <Dir d, class YeeCell, typename ValueType>
	static void	update(YeeCell && f, double delt, double delx, ValueType right, ValueType left){
		static constexpr Dir 			J = d;
		static constexpr Dir 			K = MutuallyOrthogonal<I,J>::value;
		
		// auto hold = FP::get(f);
		FP::get(f) = TimePolicy::template get<EMField, FP>(f)
					+yu_details::Coeff<FT>::value*LeviCivita<I, J, K>::value*TimePolicy::curl_coeff*delt/delx*(right - left);

		// TimePolicy::template increment<EMField, FP>(f, hold);
	}

	template <typename YeeCell, typename ValueType>
	static void increment(YeeCell && f, ValueType && val){
		TimePolicy::template increment<EMField, FP>(std::forward<YeeCell>(f), std::forward<ValueType>(val));
	}




	///////////// REVERSE
	// reverse_update in a single direction, with decrement of the time policy
	template <class YeeCell, typename ValueType, typename T = EMField>
	static void	reverse_update(YeeCell && f, double delt, ValueType curl){

		auto hold = TimePolicy::template get_first<EMField, FP>(f);
		auto val = TimePolicy::template get_reverse<EMField, FP>(f)
					-yu_details::Coeff<FT>::value*TimePolicy::curl_coeff/TimePolicy::last_coeff*delt*curl;

		TimePolicy::template decrement<EMField, FP>(f, val);
		FP::get(f) = hold;
	}

	// // reverse_update with left and right values provided as input
	// template <Dir d, class YeeCell, typename ValueType, typename T = EMField>
	// static void	reverse_update(YeeCell && f, double delt, double delx, ValueType right, ValueType left){
	// 	static constexpr Dir 			J = d;
	// 	static constexpr Dir 			K = MutuallyOrthogonal<I,J>::value;
		
	// 	// auto hold = FP::get(f);
	// 	FP::get(f) = TimePolicy::template get<EMField, FP>(f)
	// 				+yu_details::Coeff<FT>::value*LeviCivita<I, J, K>::value*TimePolicy::curl_coeff*delt/delx*(right - left);

	// 	// TimePolicy::template decrement<EMField, FP>(f, hold);
	// }
	
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// Most Generalized version of the FDTD update


// template <class Mode, 
// 		  class TimePolicy 			= BD1, 
// 		  class PMLCoeffPolicy 		= PMLCoeff<true>,
// 		  template <typename> class FieldPolicy = GetField,
// 		  template <typename,Dir> class DifferencePolicy = DefaultDifferenceOperator>
// struct YeeUpdate{
// 	// template <typename, typename, typename, typename, typename> class DifferencePolicy = DifferenceOperator
// 	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
// 	static_assert(std::is_same<TemporalScheme, TimePolicy>::value, "YeeUpdate needs a valid TemporalScheme");
// private:
// 	double dt, dx;

// 	template <typename EMField>
// 	struct atomic_update{
// 		static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
// 		// make sure this is only B or D
// 		// static_assert(EMField::)

// 		typedef FieldPolicy<EMField> 	FP;
// 		typedef TimePolicy 				TP;
// 		typedef DifferencePolicy 		DP;

// 		// ***** NEED TO INCORPORATE PML POLICY!!
// 		template <class YeeCell>
// 		static void get(YeeCell && f, double dt, double dx){
// 			auto hold = FP::get(f);
// 			FP::get(f) = TP::template get<EMField, FP>(f) + TP::curl_coeff*dt/dx*DP<EMField, Dir::Y>::get(f);
// 		}
// 	};

// public:
// 	YeeUpdate(double deltat, double deltax): dt(deltat), dx(deltax) {};

// 	template <class YeeCell>
// 	void operator()(YeeCell && f){
// 		update(f, dt, dx);
// 	};


// 	template <class YeeCell>
// 	static void update(YeeCell & f, double delt, double delx){
// 		auto hDx = FieldPolicy<fdtd::Dx>::get(f);
// 		FieldPolicy<fdtd::Dx>::get(f) = TimePolicy::template get<fdtd::Dx, FieldPolicy<fdtd::Dx>>(f)
// 				+TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKy(f)*DifferencePolicy<Hz, Dir::Y>::get(f)
// 											     - 1.0/PMLCoeffPolicy::pmlEKz(f)*DifferencePolicy<Hy, Dir::Z>::get(f));
// 		auto hDy = FieldPolicy<fdtd::Dy>::get(f);
// 		FieldPolicy<fdtd::Dy>::get(f) = TimePolicy::template get<fdtd::Dy, FieldPolicy<fdtd::Dy>>(f)
// 				+TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKz(f)*DifferencePolicy<Hx, Dir::Z>::get(f)
// 											     - 1.0/PMLCoeffPolicy::pmlEKx(f)*DifferencePolicy<Hz, Dir::X>::get(f));
// 		auto hDz = FieldPolicy<fdtd::Dz>::get(f);
// 		FieldPolicy<fdtd::Dz>::get(f) = TimePolicy::template get<fdtd::Dz, FieldPolicy<fdtd::Dz>>(f)
// 				+TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKx(f)*DifferencePolicy<Hy, Dir::X>::get(f)
// 											     - 1.0/PMLCoeffPolicy::pmlEKy(f)*DifferencePolicy<Hx, Dir::Y>::get(f));
	
// 		TimePolicy::template increment<fdtd::Dx, FieldPolicy<fdtd::Dx>>(f, hDx);
// 		TimePolicy::template increment<fdtd::Dy, FieldPolicy<fdtd::Dy>>(f, hDy);
// 		TimePolicy::template increment<fdtd::Dz, FieldPolicy<fdtd::Dz>>(f, hDz);

// 	}

	
// };



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <class Mode, 
		  class TimePolicy 			= BD1, 
		  class PMLCoeffPolicy 		= PMLCoeff<true>,
		  template <typename> class FieldPolicy = GetField,
		  template <typename,Dir> class DifferencePolicy = DefaultDifferenceOperator>
struct YeeUpdateD{
	// template <typename, typename, typename, typename, typename> class DifferencePolicy = DifferenceOperator
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	static_assert(std::is_same<TemporalScheme, TimePolicy>::value, "YeeUpdate needs a valid TemporalScheme");
};


// specialization for 3D
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateD<ThreeD, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateD<ThreeD, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell && f){
		update(f, dt, dx);
	};


	template <class YeeCell>
	static void update(YeeCell & f, double delt, double delx){
		auto hDx = FieldPolicy<fdtd::Dx>::get(f);
		FieldPolicy<fdtd::Dx>::get(f) = TimePolicy::template get<fdtd::Dx, FieldPolicy<fdtd::Dx>>(f)
				+TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKy(f)*DifferencePolicy<Hz, Dir::Y>::get(f)
											     - 1.0/PMLCoeffPolicy::pmlEKz(f)*DifferencePolicy<Hy, Dir::Z>::get(f));
		auto hDy = FieldPolicy<fdtd::Dy>::get(f);
		FieldPolicy<fdtd::Dy>::get(f) = TimePolicy::template get<fdtd::Dy, FieldPolicy<fdtd::Dy>>(f)
				+TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKz(f)*DifferencePolicy<Hx, Dir::Z>::get(f)
											     - 1.0/PMLCoeffPolicy::pmlEKx(f)*DifferencePolicy<Hz, Dir::X>::get(f));
		auto hDz = FieldPolicy<fdtd::Dz>::get(f);
		FieldPolicy<fdtd::Dz>::get(f) = TimePolicy::template get<fdtd::Dz, FieldPolicy<fdtd::Dz>>(f)
				+TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKx(f)*DifferencePolicy<Hy, Dir::X>::get(f)
											     - 1.0/PMLCoeffPolicy::pmlEKy(f)*DifferencePolicy<Hx, Dir::Y>::get(f));
	
		TimePolicy::template increment<fdtd::Dx, FieldPolicy<fdtd::Dx>>(f, hDx);
		TimePolicy::template increment<fdtd::Dy, FieldPolicy<fdtd::Dy>>(f, hDy);
		TimePolicy::template increment<fdtd::Dz, FieldPolicy<fdtd::Dz>>(f, hDz);

	}
};

// specialization for TE
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateD<TE, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateD<TE, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell && f){
		update(f, dt, dx);
	};

	template <class YeeCell>
	static void update(YeeCell & f, double delt, double delx){
		auto hDx = FieldPolicy<fdtd::Dx>::get(f);
		FieldPolicy<fdtd::Dx>::get(f) = TimePolicy::template get<fdtd::Dx>(f)
				+TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKy(f)*DifferencePolicy<Hz, Dir::Y>::get(f));
		auto hDy = FieldPolicy<fdtd::Dy>::get(f);
		FieldPolicy<fdtd::Dy>::get(f) = TimePolicy::template get<fdtd::Dy>(f)
				+TimePolicy::curl_coeff*delt/delx*(-1.0/PMLCoeffPolicy::pmlEKx(f)*DifferencePolicy<Hz, Dir::X>::get(f));
	
		TimePolicy::template increment<fdtd::Dx>(f, hDx);
		TimePolicy::template increment<fdtd::Dy>(f, hDy);
	};
};


// specialization for TM
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateD<TM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateD<TM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell && f){
		update(f, dt, dx);
	};

	template <class YeeCell>
	static void update(YeeCell & f, double delt, double delx){
		auto hDz = FieldPolicy<fdtd::Dz>::get(f);
		FieldPolicy<fdtd::Dz>::get(f) = TimePolicy::template get<fdtd::Dz>(f) 
			   						  + TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKx(f)*DifferencePolicy<Hy, Dir::X>::get(f)
					  												    - 1.0/PMLCoeffPolicy::pmlEKy(f)*DifferencePolicy<Hx, Dir::Y>::get(f));
		
		TimePolicy::template increment<fdtd::Dz>(f, hDz);
	};
};



// specialization for TEM
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateD<TEM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateD<TEM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell && f){
		update(f, dt, dx);
	};

	template <class YeeCell>
	static void update(YeeCell & f, double delt, double delx){
		// std::cout << "here 1" << std::endl;
		auto hDz = FieldPolicy<fdtd::Dz>::get(f);
		// std::cout << "here 2" << std::endl;
		FieldPolicy<fdtd::Dz>::get(f) = TimePolicy::template get<fdtd::Dz>(f)
									  + TimePolicy::curl_coeff*delt/delx*(1.0/PMLCoeffPolicy::pmlEKx(f)*DifferencePolicy<Hy, Dir::X>::get(f));

		// std::cout << "here 3" << std::endl;
		TimePolicy::template increment<fdtd::Dz>(f, hDz);
		// std::cout << "here 4" << std::endl;
	};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <class Mode, 
		  class TimePolicy 			= BD1, 
		  class PMLCoeffPolicy 		= PMLCoeff<true>,
		  template <typename> class FieldPolicy = GetField,
		  template <typename,Dir> class DifferencePolicy = DefaultDifferenceOperator>
struct YeeUpdateB{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	static_assert(std::is_same<TemporalScheme, TimePolicy>::value, "YeeUpdate needs a valid TemporalScheme");
};



// specialization for 3D
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateB<ThreeD, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateB<ThreeD, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell && f){
		auto hBx = FieldPolicy<fdtd::Bx>::get(f);
		FieldPolicy<fdtd::Bx>::get(f) = TimePolicy::template get<fdtd::Bx>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKy(f)*DifferencePolicy<fdtd::Ez, Dir::Y>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKz(f)*DifferencePolicy<fdtd::Ey, Dir::Z>::get(f));
		auto hBy = FieldPolicy<fdtd::By>::get(f);
		FieldPolicy<fdtd::By>::get(f) = TimePolicy::template get<fdtd::By>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKz(f)*DifferencePolicy<fdtd::Ex, Dir::Z>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKx(f)*DifferencePolicy<fdtd::Ez, Dir::X>::get(f));
		auto hBz = FieldPolicy<fdtd::Bz>::get(f);
		FieldPolicy<fdtd::Bz>::get(f) = TimePolicy::template get<fdtd::Bz>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKx(f)*DifferencePolicy<fdtd::Ey, Dir::X>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKy(f)*DifferencePolicy<fdtd::Ex, Dir::Y>::get(f));

		TimePolicy::template increment<fdtd::Bx>(f, hBx);
		TimePolicy::template increment<fdtd::By>(f, hBy);
		TimePolicy::template increment<fdtd::Bz>(f, hBz);

	};
};


// specialization for TE
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateB<TE, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateB<TE, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell && f){
		auto hBz = FieldPolicy<fdtd::Bz>::get(f);
		FieldPolicy<fdtd::Bz>::get(f) = TimePolicy::template get<fdtd::Bz>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKx(f)*DifferencePolicy<fdtd::Ey, Dir::X>::get(f)
					 						  -1.0/PMLCoeffPolicy::pmlHKy(f)*DifferencePolicy<fdtd::Ex, Dir::Y>::get(f));

		TimePolicy::template increment<fdtd::Bz>(f, hBz);
	};
};


// specialization for TM
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateB<TM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateB<TM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell && f){
		auto hBx = FieldPolicy<fdtd::Bx>::get(f);
		FieldPolicy<fdtd::Bx>::get(f) = TimePolicy::template get<fdtd::Bx>(f)
				-TimePolicy::curl_coeff*dt/dx*(1.0/PMLCoeffPolicy::pmlHKy(f)*DifferencePolicy<fdtd::Ez, Dir::Y>::get(f));
		auto hBy = FieldPolicy<fdtd::By>::get(f);
		FieldPolicy<fdtd::By>::get(f) = TimePolicy::template get<fdtd::By>(f)
				-TimePolicy::curl_coeff*dt/dx*(-1.0/PMLCoeffPolicy::pmlHKx(f)*DifferencePolicy<fdtd::Ez, Dir::X>::get(f));

		TimePolicy::template increment<fdtd::Bx>(f, hBx);
		TimePolicy::template increment<fdtd::By>(f, hBy);
	};
};


// specialization for TEM
template <typename TimePolicy, 
		  typename PMLCoeffPolicy, 
		  template <typename> class FieldPolicy,
		  template <typename,Dir> class DifferencePolicy>
struct YeeUpdateB<TEM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>{
	double dt, dx;

	YeeUpdateB<TEM, TimePolicy, PMLCoeffPolicy, FieldPolicy, DifferencePolicy>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell && f){
		auto hBy = FieldPolicy<fdtd::By>::get(f);
		FieldPolicy<fdtd::By>::get(f) = TimePolicy::template get<fdtd::By>(f)
				-TimePolicy::curl_coeff*dt/dx*(-1.0/PMLCoeffPolicy::pmlHKx(f)*DifferencePolicy<fdtd::Ez, Dir::X>::get(f));
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



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



}// end namespace fdtd

#endif