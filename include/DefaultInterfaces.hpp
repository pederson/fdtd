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


#define FDTD_GET_FIELD_DEF(name) 	\
	template <>						\
	struct GetField<fdtd::name>{		\
		template <class YeeCell>	\
		static decltype(auto) get(YeeCell && f) {return f.name();}		\
	};

FDTD_GET_FIELD_DEF(Ex);
FDTD_GET_FIELD_DEF(Ey);
FDTD_GET_FIELD_DEF(Ez);

FDTD_GET_FIELD_DEF(Dx);
FDTD_GET_FIELD_DEF(Dy);
FDTD_GET_FIELD_DEF(Dz);

FDTD_GET_FIELD_DEF(Hx);
FDTD_GET_FIELD_DEF(Hy);
FDTD_GET_FIELD_DEF(Hz);

FDTD_GET_FIELD_DEF(Bx);
FDTD_GET_FIELD_DEF(By);
FDTD_GET_FIELD_DEF(Bz);



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

// default PML interface
struct GetPML{
	template <FieldType ft, Dir I, Dir J, typename YeeCell>
	static decltype(auto) integrator(YeeCell && f){return f.template pmlI<ft, I, J>();};

	template <FieldType ft, Dir I, typename YeeCell>
	static decltype(auto) B(YeeCell && f){return f.template pmlB<ft, I>();};

	template <FieldType ft, Dir I, typename YeeCell>
	static decltype(auto) C(YeeCell && f){return f.template pmlC<ft, I>();};

	template <FieldType ft, Dir I, typename YeeCell>
	static decltype(auto) F(YeeCell && f){return f.template pmlF<ft, I>();};

	template <FieldType ft, Dir I, typename YeeCell>
	static decltype(auto) G(YeeCell && f){return f.template pmlG<ft, I>();};
};

//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



}// end namespace fdtd

#endif