#ifndef _YEEFIELDS_H
#define _YEEFIELDS_H

#include "FDTDConstants.hpp"

#include <type_traits>

namespace fdtd{


// the ContainerT must be able to overload std::get<I>()
template<typename Mode,
		 typename scalar_type,
		 std::size_t NumE,
		 std::size_t NumH,
		 template<class, std::size_t> class ContainerT>
struct YeeFieldData{
	typedef scalar_type type;
	typedef Mode mode;

	static_assert(std::is_base_of<EMMode, Mode>::value, "Mode must be of a predefined type");

	constexpr YeeFieldData(){
		D.fill(0.0); E.fill(0.0);
		B.fill(0.0); H.fill(0.0);
	}


	ContainerT<scalar_type, NumE>		D;
	ContainerT<scalar_type, NumH>		B;
	ContainerT<scalar_type, NumE>		E;
	ContainerT<scalar_type, NumH>		H;

};


template<typename Mode, typename scalar_type, template <class,std::size_t> class ContainerT>
struct YeeFields: public YeeFieldData<Mode, scalar_type, Mode::numE, Mode::numH, ContainerT>{
	typedef YeeFieldData<Mode, scalar_type, Mode::numE, Mode::numH, ContainerT> BaseT;
};

template<typename scalar_type, template <class,std::size_t> class ContainerT>
struct YeeFields<ThreeD, scalar_type, ContainerT>: public YeeFieldData<ThreeD, scalar_type, ThreeD::numE, ThreeD::numH, ContainerT>{
	typedef YeeFieldData<ThreeD, scalar_type, ThreeD::numE, ThreeD::numH, ContainerT> BaseT;

	scalar_type & Ex() {return std::get<0>(BaseT::E);};
	scalar_type & Ey() {return std::get<1>(BaseT::E);};
	scalar_type & Ez() {return std::get<2>(BaseT::E);};

	scalar_type & Dx() {return std::get<0>(BaseT::D);};
	scalar_type & Dy() {return std::get<1>(BaseT::D);};
	scalar_type & Dz() {return std::get<2>(BaseT::D);};

	scalar_type & Hx() {return std::get<0>(BaseT::H);};
	scalar_type & Hy() {return std::get<1>(BaseT::H);};
	scalar_type & Hz() {return std::get<2>(BaseT::H);};

	scalar_type & Bx() {return std::get<0>(BaseT::B);};
	scalar_type & By() {return std::get<1>(BaseT::B);};
	scalar_type & Bz() {return std::get<2>(BaseT::B);};



	const scalar_type & Ex() const {return std::get<0>(BaseT::E);};
	const scalar_type & Ey() const {return std::get<1>(BaseT::E);};
	const scalar_type & Ez() const {return std::get<2>(BaseT::E);};

	const scalar_type & Dx() const {return std::get<0>(BaseT::D);};
	const scalar_type & Dy() const {return std::get<1>(BaseT::D);};
	const scalar_type & Dz() const {return std::get<2>(BaseT::D);};

	const scalar_type & Hx() const {return std::get<0>(BaseT::H);};
	const scalar_type & Hy() const {return std::get<1>(BaseT::H);};
	const scalar_type & Hz() const {return std::get<2>(BaseT::H);};

	const scalar_type & Bx() const {return std::get<0>(BaseT::B);};
	const scalar_type & By() const {return std::get<1>(BaseT::B);};
	const scalar_type & Bz() const {return std::get<2>(BaseT::B);};


	// templated call
	template <Dir d>
	scalar_type & E() const {return std::get<static_cast<char>(d)>(BaseT::E);};
	template <Dir d>
	scalar_type & D() const {return std::get<static_cast<char>(d)>(BaseT::D);};
	template <Dir d>
	scalar_type & H() const {return std::get<static_cast<char>(d)>(BaseT::H);};
	template <Dir d>
	scalar_type & B() const {return std::get<static_cast<char>(d)>(BaseT::B);};
};



template<typename scalar_type, template <class,std::size_t> class ContainerT>
struct YeeFields<TE, scalar_type, ContainerT>: public YeeFieldData<TE, scalar_type, TE::numE, TE::numH, ContainerT>{
	typedef YeeFieldData<TE, scalar_type, TE::numE, TE::numH, ContainerT> BaseT;

	scalar_type & Ex() {return std::get<0>(BaseT::E);};
	scalar_type & Ey() {return std::get<1>(BaseT::E);};

	scalar_type & Dx() {return std::get<0>(BaseT::D);};
	scalar_type & Dy() {return std::get<1>(BaseT::D);};

	scalar_type & Hz() {return std::get<0>(BaseT::H);};

	scalar_type & Bz() {return std::get<0>(BaseT::B);};


	const scalar_type & Ex() const {return std::get<0>(BaseT::E);};
	const scalar_type & Ey() const {return std::get<1>(BaseT::E);};

	const scalar_type & Dx() const {return std::get<0>(BaseT::D);};
	const scalar_type & Dy() const {return std::get<1>(BaseT::D);};

	const scalar_type & Hz() const {return std::get<0>(BaseT::H);};

	const scalar_type & Bz() const {return std::get<0>(BaseT::B);};

};



template<typename scalar_type, template <class,std::size_t> class ContainerT>
struct YeeFields<TM, scalar_type, ContainerT>: public YeeFieldData<TM, scalar_type, TM::numE, TM::numH, ContainerT>{
	typedef YeeFieldData<TM, scalar_type, TM::numE, TM::numH, ContainerT> BaseT;

	scalar_type & Ez() {return std::get<0>(BaseT::E);};

	scalar_type & Dz() {return std::get<0>(BaseT::D);};

	scalar_type & Hx() {return std::get<0>(BaseT::H);};
	scalar_type & Hy() {return std::get<1>(BaseT::H);};

	scalar_type & Bx() {return std::get<0>(BaseT::B);};
	scalar_type & By() {return std::get<1>(BaseT::B);};


	const scalar_type & Ez() const {return std::get<0>(BaseT::E);};

	const scalar_type & Dz() const {return std::get<0>(BaseT::D);};

	const scalar_type & Hx() const {return std::get<0>(BaseT::H);};
	const scalar_type & Hy() const {return std::get<1>(BaseT::H);};

	const scalar_type & Bx() const {return std::get<0>(BaseT::B);};
	const scalar_type & By() const {return std::get<1>(BaseT::B);};

};



template<typename scalar_type, template <class,std::size_t> class ContainerT>
struct YeeFields<TEM, scalar_type, ContainerT>: public YeeFieldData<TEM, scalar_type, TEM::numE, TEM::numH, ContainerT>{
	typedef YeeFieldData<TEM, scalar_type, TEM::numE, TEM::numH, ContainerT> BaseT;

	scalar_type & Ez() {return std::get<0>(BaseT::E);};

	scalar_type & Dz() {return std::get<0>(BaseT::D);};

	scalar_type & Hy() {return std::get<0>(BaseT::H);};

	scalar_type & By() {return std::get<0>(BaseT::B);};



	const scalar_type & Ez() const {return std::get<0>(BaseT::E);};

	const scalar_type & Dz() const {return std::get<0>(BaseT::D);};

	const scalar_type & Hy() const {return std::get<0>(BaseT::H);};

	const scalar_type & By() const {return std::get<0>(BaseT::B);};

};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <int n, typename Mode, typename scalar_type>
struct BDFields {
	typedef YeeFields<Mode, scalar_type, std::array> 	Fields;
private: 
	Fields mFields[n];
	unsigned int mCurr;
public:
	BDFields<n, Mode, scalar_type>() : mCurr(0) {};

	BDFields & BD() {return *this;};
	Fields & BD(int i) {return mFields[(mCurr+i)%n];};
	void increment() {mCurr = (mCurr==0 ? n-1 : mCurr-1);};
};



}// end namespace fdtd

#endif