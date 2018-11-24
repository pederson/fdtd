#ifndef _STATICFIELDS_H
#define _STATICFIELDS_H

#include "FDTDConstants.hpp"
// #include "YeeFields.hpp"

#include <type_traits>

namespace fdtd{


template<typename scalar_type, template <class,std::size_t> class ContainerT>
struct StaticFields{
private:
	struct Efields{
	private:
		ContainerT<scalar_type, 3>		E;

	public:
		constexpr Efields(){
			E.fill(0.0);
		};

		scalar_type & Ex() {return std::get<0>(E);};
		scalar_type & Ey() {return std::get<1>(E);};
		scalar_type & Ez() {return std::get<2>(E);};

		// const scalar_type & Ex() {return std::get<0>(E);};
		// const scalar_type & Ey() {return std::get<1>(E);};
		// const scalar_type & Ez() {return std::get<2>(E);};
	};

	struct Bfields{
	private:
		ContainerT<scalar_type, 3>		B;
	public:
		constexpr Bfields(){
			B.fill(0.0);
		};

		scalar_type & Bx() {return std::get<0>(B);};
		scalar_type & By() {return std::get<1>(B);};
		scalar_type & Bz() {return std::get<2>(B);};

		// const scalar_type & Bx() {return std::get<0>(B);};
		// const scalar_type & By() {return std::get<1>(B);};
		// const scalar_type & Bz() {return std::get<2>(B);};
	};
	
	Bfields mB;
	Efields mE;

public:
	Bfields & static_magnetic_field() {return mB;};
	// const Bfields & static_magnetic_field() {return mB;};

	Efields & static_electric_field() {return mE;};
	// const Efields & static_electric_field() {return mE;};
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

}// end namespace fdtd

#endif