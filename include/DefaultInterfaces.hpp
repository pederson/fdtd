#ifndef _DEFAULTINTERFACES_H
#define _DEFAULTINTERFACES_H

#include "FDTDConstants.hpp"

namespace fdtd{



struct GetSelf{
	template <typename Object>
	static decltype(auto) get(Object && o){return o;};
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// define the default field-getting behavior
template <typename EMField>
struct GetField{
	static_assert(std::is_same<Field, EMField>::value, "GetField needs a valid EMField");
};

template <>
struct GetField<fdtd::Ex>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Ex();}
};

template <>
struct GetField<fdtd::Ey>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Ey();}
};

template <>
struct GetField<fdtd::Ez>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Ez();}
};



template <>
struct GetField<fdtd::Dx>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Dx();}
};

template <>
struct GetField<fdtd::Dy>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Dy();}
};

template <>
struct GetField<fdtd::Dz>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Dz();}
};


template <>
struct GetField<fdtd::Bx>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Bx();}
};

template <>
struct GetField<fdtd::By>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.By();}
};

template <>
struct GetField<fdtd::Bz>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Bz();}
};


template <>
struct GetField<fdtd::Hx>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Hx();}
};

template <>
struct GetField<fdtd::Hy>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Hy();}
};

template <>
struct GetField<fdtd::Hz>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Hz();}
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



// Define the default neighbor-getting behavior
template <Dir d, Orientation o>
struct GetNeighbor{};

template <Dir d>
struct GetNeighbor <d, Orientation::MIN>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.getNeighborMin(static_cast<int>(d));};
};

template <Dir d>
struct GetNeighbor <d, Orientation::MAX>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.getNeighborMax(static_cast<int>(d));};
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// Define the PML coefficient depending on whether the user 
// wants to include it or not
template <bool pml_option>
struct PMLCoeff{
};

template<>
struct PMLCoeff<true>{
	template <class YeeCell>
	static double pmlEKx(YeeCell & f) {return f.pmlEKx();};

	template <class YeeCell>
	static double pmlEKy(YeeCell & f) {return f.pmlEKy();};

	template <class YeeCell>
	static double pmlEKz(YeeCell & f) {return f.pmlEKz();};


	template <class YeeCell>
	static double pmlHKx(YeeCell & f) {return f.pmlHKx();};

	template <class YeeCell>
	static double pmlHKy(YeeCell & f) {return f.pmlHKy();};

	template <class YeeCell>
	static double pmlHKz(YeeCell & f) {return f.pmlHKz();};
};



template<>
struct PMLCoeff<false>{
	template <class YeeCell>
	static double pmlEKx(YeeCell & f) {return 1.0;};

	template <class YeeCell>
	static double pmlEKy(YeeCell & f) {return 1.0;};

	template <class YeeCell>
	static double pmlEKz(YeeCell & f) {return 1.0;};


	template <class YeeCell>
	static double pmlHKx(YeeCell & f) {return 1.0;};

	template <class YeeCell>
	static double pmlHKy(YeeCell & f) {return 1.0;};

	template <class YeeCell>
	static double pmlHKz(YeeCell & f) {return 1.0;};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



}// end namespace fdtd

#endif