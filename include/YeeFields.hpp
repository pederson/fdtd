#ifndef _YEEFIELDS_H
#define _YEEFIELDS_H

#include "FDTDConstants.hpp"

#include <type_traits>

namespace fdtd{


// the ContainerT must be able to overload std::get<I>()
template<typename Mode,
		 typename FieldType,
		 std::size_t NumE,
		 std::size_t NumH,
		 template<class, std::size_t> class ContainerT>
struct YeeFields{
	typedef FieldType type;
	typedef Mode mode;

	static_assert(std::is_base_of<EMMode, Mode>::value, "Mode must be of a predefined type");

	constexpr YeeFields(){
		D.fill(0.0); E.fill(0.0);
		B.fill(0.0); H.fill(0.0);
	}


	ContainerT<FieldType, NumE>		D;
	ContainerT<FieldType, NumH>		B;
	ContainerT<FieldType, NumE>		E;
	ContainerT<FieldType, NumH>		H;

};



template<typename FieldType, template <class,std::size_t> class ContainerT>
struct YeeFields3D: public YeeFields<ThreeD, FieldType, 3, 3, ContainerT>{
	typedef YeeFields<ThreeD, FieldType, 3, 3, ContainerT> BaseT;
	
	FieldType & Ex() {return std::get<0>(BaseT::E);};
	FieldType & Ey() {return std::get<1>(BaseT::E);};
	FieldType & Ez() {return std::get<2>(BaseT::E);};

	FieldType & Dx() {return std::get<0>(BaseT::D);};
	FieldType & Dy() {return std::get<1>(BaseT::D);};
	FieldType & Dz() {return std::get<2>(BaseT::D);};

	FieldType & Hx() {return std::get<0>(BaseT::H);};
	FieldType & Hy() {return std::get<1>(BaseT::H);};
	FieldType & Hz() {return std::get<2>(BaseT::H);};

	FieldType & Bx() {return std::get<0>(BaseT::B);};
	FieldType & By() {return std::get<1>(BaseT::B);};
	FieldType & Bz() {return std::get<2>(BaseT::B);};



	const FieldType & Ex() const {return std::get<0>(BaseT::E);};
	const FieldType & Ey() const {return std::get<1>(BaseT::E);};
	const FieldType & Ez() const {return std::get<2>(BaseT::E);};

	const FieldType & Dx() const {return std::get<0>(BaseT::D);};
	const FieldType & Dy() const {return std::get<1>(BaseT::D);};
	const FieldType & Dz() const {return std::get<2>(BaseT::D);};

	const FieldType & Hx() const {return std::get<0>(BaseT::H);};
	const FieldType & Hy() const {return std::get<1>(BaseT::H);};
	const FieldType & Hz() const {return std::get<2>(BaseT::H);};

	const FieldType & Bx() const {return std::get<0>(BaseT::B);};
	const FieldType & By() const {return std::get<1>(BaseT::B);};
	const FieldType & Bz() const {return std::get<2>(BaseT::B);};
};



template<typename FieldType, template <class,std::size_t> class ContainerT>
struct YeeFieldsTE: public YeeFields<TE, FieldType, 2, 1, ContainerT>{
	typedef YeeFields<TE, FieldType, 2, 1, ContainerT> BaseT;
	
	FieldType & Ex() {return std::get<0>(BaseT::E);};
	FieldType & Ey() {return std::get<1>(BaseT::E);};

	FieldType & Dx() {return std::get<0>(BaseT::D);};
	FieldType & Dy() {return std::get<1>(BaseT::D);};

	FieldType & Hz() {return std::get<0>(BaseT::H);};

	FieldType & Bz() {return std::get<0>(BaseT::B);};


	const FieldType & Ex() const {return std::get<0>(BaseT::E);};
	const FieldType & Ey() const {return std::get<1>(BaseT::E);};

	const FieldType & Dx() const {return std::get<0>(BaseT::D);};
	const FieldType & Dy() const {return std::get<1>(BaseT::D);};

	const FieldType & Hz() const {return std::get<0>(BaseT::H);};

	const FieldType & Bz() const {return std::get<0>(BaseT::B);};
};



template<typename FieldType, template <class,std::size_t> class ContainerT>
struct YeeFieldsTM: public YeeFields<TM, FieldType, 1, 2, ContainerT>{
	typedef YeeFields<TM, FieldType, 1, 2, ContainerT> BaseT;
	
	FieldType & Ez() {return std::get<0>(BaseT::E);};

	FieldType & Dz() {return std::get<0>(BaseT::D);};

	FieldType & Hx() {return std::get<0>(BaseT::H);};
	FieldType & Hy() {return std::get<1>(BaseT::H);};

	FieldType & Bx() {return std::get<0>(BaseT::B);};
	FieldType & By() {return std::get<1>(BaseT::B);};


	const FieldType & Ez() const {return std::get<0>(BaseT::E);};

	const FieldType & Dz() const {return std::get<0>(BaseT::D);};

	const FieldType & Hx() const {return std::get<0>(BaseT::H);};
	const FieldType & Hy() const {return std::get<1>(BaseT::H);};

	const FieldType & Bx() const {return std::get<0>(BaseT::B);};
	const FieldType & By() const {return std::get<1>(BaseT::B);};
};



template<typename FieldType, template <class,std::size_t> class ContainerT>
struct YeeFieldsTEM: public YeeFields<TEM, FieldType, 1, 1, ContainerT>{
	typedef YeeFields<TEM, FieldType, 1, 1, ContainerT> BaseT;
	
	FieldType & Ez() {return std::get<0>(BaseT::E);};

	FieldType & Dz() {return std::get<0>(BaseT::D);};

	FieldType & Hy() {return std::get<0>(BaseT::H);};

	FieldType & By() {return std::get<0>(BaseT::B);};



	const FieldType & Ez() const {return std::get<0>(BaseT::E);};

	const FieldType & Dz() const {return std::get<0>(BaseT::D);};

	const FieldType & Hy() const {return std::get<0>(BaseT::H);};

	const FieldType & By() const {return std::get<0>(BaseT::B);};
};

}// end namespace fdtd

#endif