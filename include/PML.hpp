#ifndef _PML_H
#define _PML_H

#include <array>
#include <algorithm>

#include "FDTDConstants.hpp"

namespace fdtd{



template<typename Mode> class StoredPMLx;
template<typename Mode> class StoredPMLy;
template<typename Mode> class StoredPMLz;	
class NoPMLx;
class NoPMLy;
class NoPMLz;

template <typename Mode,
		  bool X, bool Y, bool Z,
		  class PMLTypeX = StoredPMLx<Mode>,
		  class PMLTypeY = StoredPMLy<Mode>,
		  class PMLTypeZ = StoredPMLz<Mode>
		  >
class PML : public std::conditional<X, PMLTypeX, NoPMLx>::type
		  , public std::conditional<Y, PMLTypeY, NoPMLy>::type
		  , public std::conditional<Z, PMLTypeZ, NoPMLz>::type
{
public:

};








//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <class Mode>
struct UpdatePMLD{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};



// specialization for 3D
template<>
struct UpdatePMLD<ThreeD>{
	double dt, dx;

	UpdatePMLD<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		// x direction
		f.pmlHIxz() = f.pmlEBx()*f.pmlHIxz() + f.pmlECx()/dx*(f.Hz() - f.getNeighborMin(2).Hz());
		f.Dy() -= f.pmlHIxz()*dt;

		f.pmlHIxy() = f.pmlEBx()*f.pmlHIxy() + f.pmlECx()/dx*(f.Hy() - f.getNeighborMin(0).Hy());
		f.Dz() += f.pmlHIxy()*dt;

		// y direction
		f.pmlHIyz() = f.pmlEBy()*f.pmlHIyz() + f.pmlECy()/dx*(f.Hz() - f.getNeighborMin(1).Hz());	
		f.Dx() += f.pmlHIyz()*dt;

		f.pmlHIyx() = f.pmlEBy()*f.pmlHIyx() + f.pmlECy()/dx*(f.Hx() - f.getNeighborMin(1).Hx());	
		f.Dz() -= f.pmlHIyx()*dt;

		// z direction
		f.pmlHIzy() = f.pmlEBz()*f.pmlHIzy() + f.pmlECz()/dx*(f.Hy() - f.getNeighborMin(2).Hy());	
		f.Dx() -= f.pmlHIzy()*dt;

		f.pmlHIzx() = f.pmlEBz()*f.pmlHIzx() + f.pmlECz()/dx*(f.Hx() - f.getNeighborMin(2).Hx());	
		f.Dy() += f.pmlHIzx()*dt;
	};
};



// specialization for TE
template<>
struct UpdatePMLD<TE>{
	double dt, dx;

	UpdatePMLD<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		// x direction
		f.pmlHIxz() = f.pmlEBx()*f.pmlHIxz() + f.pmlECx()/dx*(f.Hz() - f.getNeighborMin(2).Hz());
		f.Dy() -= f.pmlHIxz()*dt;

		// y direction
		f.pmlHIyz() = f.pmlEBy()*f.pmlHIyz() + f.pmlECy()/dx*(f.Hz() - f.getNeighborMin(1).Hz());	
		f.Dx() += f.pmlHIyz()*dt;
	};
};


// specialization for TM
template<>
struct UpdatePMLD<TM>{
	double dt, dx;

	UpdatePMLD<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		// x direction
		f.pmlHIxy() = f.pmlEBx()*f.pmlHIxy() + f.pmlECx()/dx*(f.Hy() - f.getNeighborMin(0).Hy());
		f.Dz() += f.pmlHIxy()*dt;

		// y direction
		f.pmlHIyx() = f.pmlEBy()*f.pmlHIyx() + f.pmlECy()/dx*(f.Hx() - f.getNeighborMin(1).Hx());	
		f.Dz() -= f.pmlHIyx()*dt;

	};
};


// specialization for TEM
template<>
struct UpdatePMLD<TEM>{
	double dt, dx;

	UpdatePMLD<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.pmlHIxy() = f.pmlEBx()*f.pmlHIxy() + f.pmlECx()/dx*(f.Hy() - f.getNeighborMin(0).Hy());	
		f.Dz() += f.pmlHIxy()*dt;
	};
};






//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <class Mode>
struct UpdatePMLB{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};


// specialization for 3D
template<>
struct UpdatePMLB<ThreeD>{
	double dt, dx;

	UpdatePMLB<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		// x direction
		f.pmlEIxz() = f.pmlHBx()*f.pmlEIxz() + f.pmlHCx()/dx*(f.Ez() - f.getNeighborMin(2).Ez());
		f.By() += f.pmlEIxz()*dt;

		f.pmlEIxy() = f.pmlHBx()*f.pmlEIxy() + f.pmlHCx()/dx*(f.Ey() - f.getNeighborMin(0).Ey());
		f.Bz() -= f.pmlEIxy()*dt;

		// y direction
		f.pmlEIyz() = f.pmlHBy()*f.pmlEIyz() + f.pmlHCy()/dx*(f.Ez() - f.getNeighborMin(1).Ez());	
		f.Bx() -= f.pmlEIyz()*dt;

		f.pmlEIyx() = f.pmlHBy()*f.pmlEIyx() + f.pmlHCy()/dx*(f.Ex() - f.getNeighborMin(1).Ex());	
		f.Bz() += f.pmlEIyx()*dt;

		// z direction
		f.pmlEIzy() = f.pmlHBz()*f.pmlEIzy() + f.pmlHCz()/dx*(f.Ey() - f.getNeighborMin(2).Ey());	
		f.Bx() += f.pmlEIzy()*dt;

		f.pmlEIzx() = f.pmlHBz()*f.pmlEIzx() + f.pmlHCz()/dx*(f.Ex() - f.getNeighborMin(2).Ex());	
		f.By() -= f.pmlHIzx()*dt;
	};
};


// specialization for TE
template<>
struct UpdatePMLB<TE>{
	double dt, dx;

	UpdatePMLB<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		// x direction
		f.pmlEIxy() = f.pmlHBx()*f.pmlEIxy() + f.pmlHCx()/dx*(f.Ey() - f.getNeighborMin(0).Ey());
		f.Bz() -= f.pmlEIxy()*dt;

		// y direction
		f.pmlEIyx() = f.pmlHBy()*f.pmlEIyx() + f.pmlHCy()/dx*(f.Ex() - f.getNeighborMin(1).Ex());	
		f.Bz() += f.pmlEIyx()*dt;
	};
};


// specialization for TM
template<>
struct UpdatePMLB<TM>{
	double dt, dx;

	UpdatePMLB<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		// x direction
		f.pmlEIxz() = f.pmlHBx()*f.pmlEIxz() + f.pmlHCx()/dx*(f.Ez() - f.getNeighborMin(2).Ez());
		f.By() += f.pmlEIxz()*dt;

		// y direction
		f.pmlEIyz() = f.pmlHBy()*f.pmlEIyz() + f.pmlHCy()/dx*(f.Ez() - f.getNeighborMin(1).Ez());	
		f.Bx() -= f.pmlEIyz()*dt;
	};
};


// specialization for TEM
template<>
struct UpdatePMLB<TEM>{
	double dt, dx;

	UpdatePMLB<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.pmlEIxz() = f.pmlHBx()*f.pmlEIxz() + f.pmlHCx()/dx*(f.getNeighborMax(0).Ez() - f.Ez());	
		f.By() += f.pmlEIxz()*dt;
	};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************





struct NoPMLx{



	// Electric PML parameters
	constexpr double pmlEKx() const {return 1.0;};

	constexpr double pmlESx() const {return 0.0;};

	constexpr double pmlEAx() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlEBx() const {return 1.0;};

	constexpr double pmlECx() const {return 0.0;};


	// convolution terms
	constexpr double pmlEIxx() const {return 0.0;};
	constexpr double pmlEIxy() const {return 0.0;};
	constexpr double pmlEIxz() const {return 0.0;};





	// Magnetic PML parameters
	constexpr double pmlHKx() const {return 1.0;};

	constexpr double pmlHSx() const {return 0.0;};

	constexpr double pmlHAx() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlHBx() const {return 1.0;};

	constexpr double pmlHCx() const {return 0.0;};


	// convolution terms
	constexpr double pmlHIxx() const {return 0.0;};
	constexpr double pmlHIxy() const {return 0.0;};
	constexpr double pmlHIxz() const {return 0.0;};
};




struct NoPMLy{

	// PML parameters
	constexpr double pmlEKy() const {return 1.0;};

	constexpr double pmlESy() const {return 0.0;};

	constexpr double pmlEAy() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlEBy() const {return 1.0;};

	constexpr double pmlECy() const {return 0.0;};


	// convolution terms
	constexpr double pmlEIyx() const {return 0.0;};
	constexpr double pmlEIyy() const {return 0.0;};
	constexpr double pmlEIyz() const {return 0.0;};








	// Magnetic PML parameters
	constexpr double pmlHKy() const {return 1.0;};

	constexpr double pmlHSy() const {return 0.0;};

	constexpr double pmlHAy() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlHBy() const {return 1.0;};

	constexpr double pmlHCy() const {return 0.0;};


	// convolution terms
	constexpr double pmlHIyx() const {return 0.0;};
	constexpr double pmlHIyy() const {return 0.0;};
	constexpr double pmlHIyz() const {return 0.0;};
};




struct NoPMLz{

	// PML parameters
	constexpr double pmlEKz() const {return 1.0;};

	constexpr double pmlESz() const {return 0.0;};

	constexpr double pmlEAz() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlEBz() const {return 1.0;};

	constexpr double pmlECz() const {return 0.0;};


	// convolution terms
	constexpr double pmlEIzx() const {return 0.0;};
	constexpr double pmlEIzy() const {return 0.0;};
	constexpr double pmlEIzz() const {return 0.0;};






	// PML parameters
	constexpr double pmlHKz() const {return 1.0;};

	constexpr double pmlHSz() const {return 0.0;};

	constexpr double pmlHAz() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlHBz() const {return 1.0;};

	constexpr double pmlHCz() const {return 0.0;};


	// convolution terms
	constexpr double pmlHIzx() const {return 0.0;};
	constexpr double pmlHIzy() const {return 0.0;};
	constexpr double pmlHIzz() const {return 0.0;};
};












//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <typename Mode>
struct PMLIx{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};

template <>
struct PMLIx<ThreeD>{
	// E convolution terms
	double EIxx;
	double EIxy;
	double EIxz;

	// H convolution terms
	double HIxx;
	double HIxy;
	double HIxz;

	// accessors
	double & pmlEIxx() {return EIxx;};
	double & pmlEIxy() {return EIxy;};
	double & pmlEIxz() {return EIxz;};
	double & pmlHIxx() {return HIxx;};
	double & pmlHIxy() {return HIxy;};
	double & pmlHIxz() {return HIxz;};
};



template <>
struct PMLIx<TE>{
	// E convolution terms
	double EIxx;
	double EIxy;

	// H convolution terms
	double HIxz;

	// accessors
	double & pmlEIxx() {return EIxx;};
	double & pmlEIxy() {return EIxy;};
	double & pmlHIxz() {return HIxz;};
};



template <>
struct PMLIx<TM>{
	// E convolution terms
	double EIxz;

	// H convolution terms
	double HIxx;
	double HIxy;

	// accessors
	double & pmlEIxz() {return EIxz;};
	double & pmlHIxx() {return HIxx;};
	double & pmlHIxy() {return HIxy;};
};


template <>
struct PMLIx<TEM>{
	// E convolution terms
	double EIxz;

	// H convolution terms
	double HIxy;

	// accessors
	double & pmlEIxz() {return EIxz;};
	double & pmlHIxy() {return HIxy;};
};




// x PML Stored class
template <typename Mode>
struct StoredPMLx : public PMLIx<Mode>{

	// PML parameters
	double EKx;
	double ESx;
	double EAx;

	// derived PML parameters
	double EBx;
	double ECx;


	// PML parameters
	double HKx;
	double HSx;
	double HAx;

	// derived PML parameters
	double HBx;
	double HCx;


	StoredPMLx()
	: EKx(1.0), ESx(0.0), EAx(0.0)
	, EBx(1.0), ECx(0.0)
	, HKx(1.0), HSx(0.0), HAx(0.0)
	, HBx(1.0), HCx(0.0) {};

	// StoredPMLx<ThreeD>(double K, double S, double A, double dt)
	// : EKx(K), ESx(S), EAx(A)
	// , EBx(1.0), ECx(0.0)
	// , EIxx(0.0), EIxy(0.0), EIxz(0.0) {setPMLParametersE(EKx, ESx, EAx, dt);
	// 								   setPMLParametersH(HKx, HSx, HAx, dt);};


	void setPMLParametersEx(double K, double S, double A, double dt){
		EKx = K; ESx = S; EAx = A;
		EBx = exp(-dt/eps0*(ESx/EKx + EAx));
		ECx = (S==0 && A == 0) ? 0.0 : ESx/EKx*1.0/(ESx+EKx*EAx)*(EBx-1.0);
	}

	void setPMLParametersHx(double K, double S, double A, double dt){
		HKx = K; HSx = S; HAx = A;
		HBx = exp(-dt/eps0*(HSx/HKx + HAx));
		HCx = (S==0 && A == 0) ? 0.0 : HSx/HKx*1.0/(HSx+HKx*HAx)*(HBx-1.0);
	}


	// PML parameters
	constexpr double & pmlEKx() {return EKx;};
	constexpr double & pmlESx() {return ESx;};
	constexpr double & pmlEAx() {return EAx;};

	constexpr double & pmlHKx() {return HKx;};
	constexpr double & pmlHSx() {return HSx;};
	constexpr double & pmlHAx() {return HAx;};

	// derived PML parameters
	constexpr double & pmlEBx() {return EBx;};
	constexpr double & pmlECx() {return ECx;};

	constexpr double & pmlHBx() {return HBx;};
	constexpr double & pmlHCx() {return HCx;};
};






//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <typename Mode>
struct PMLIy{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};

template <>
struct PMLIy<ThreeD>{
	// E convolution terms
	double EIyx;
	double EIyy;
	double EIyz;

	// H convolution terms
	double HIyx;
	double HIyy;
	double HIyz;

	// accessors
	double & pmlEIyx() {return EIyx;};
	double & pmlEIyy() {return EIyy;};
	double & pmlEIyz() {return EIyz;};
	double & pmlHIyx() {return HIyx;};
	double & pmlHIyy() {return HIyy;};
	double & pmlHIyz() {return HIyz;};
};



template <>
struct PMLIy<TE>{
	// E convolution terms
	double EIyx;
	double EIyy;

	// H convolution terms
	double HIyz;

	// accessors
	double & pmlEIyx() {return EIyx;};
	double & pmlEIyy() {return EIyy;};
	double & pmlHIyz() {return HIyz;};
};




template <>
struct PMLIy<TM>{
	// E convolution terms
	double EIyz;

	// H convolution terms
	double HIyx;
	double HIyy;

	// accessors
	double & pmlEIyz() {return EIyz;};
	double & pmlHIyx() {return HIyx;};
	double & pmlHIyy() {return HIyy;};
};




template <>
struct PMLIy<TEM>{
	// E convolution terms
	double EIyz;

	// H convolution terms
	double HIyy;

	// accessors
	double & pmlEIyz() {return EIyz;};
	double & pmlHIyy() {return HIyy;};
};




// y PML Stored class
template <typename Mode>
struct StoredPMLy : public PMLIy<Mode>{

	// PML parameters
	double EKy;
	double ESy;
	double EAy;

	// derived PML parameters
	double EBy;
	double ECy;


	// PML parameters
	double HKy;
	double HSy;
	double HAy;

	// derived PML parameters
	double HBy;
	double HCy;


	StoredPMLy()
	: EKy(1.0), ESy(0.0), EAy(0.0)
	, EBy(1.0), ECy(0.0)
	, HKy(1.0), HSy(0.0), HAy(0.0)
	, HBy(1.0), HCy(0.0) {};

	// StoredPMLy<ThreeD>(double K, double S, double A, double dt)
	// : EKy(K), ESy(S), EAy(A)
	// , EBy(1.0), ECy(0.0)
	// , EIyy(0.0), EIyy(0.0), EIyz(0.0) {setPMLParametersE(EKy, ESy, EAy, dt);
	// 								   setPMLParametersH(HKy, HSy, HAy, dt);};


	void setPMLParametersEy(double K, double S, double A, double dt){
		EKy = K; ESy = S; EAy = A;
		EBy = exp(-dt/eps0*(ESy/EKy + EAy));
		ECy = ESy/EKy*1.0/(ESy+EKy*EAy)*(EBy-1);
	}

	void setPMLParametersHy(double K, double S, double A, double dt){
		HKy = K; HSy = S; HAy = A;
		HBy = exp(-dt/eps0*(HSy/HKy + HAy));
		HCy = HSy/HKy*1.0/(HSy+HKy*HAy)*(HBy-1);
	}


	// PML parameters
	constexpr double & pmlEKy() {return EKy;};
	constexpr double & pmlESy() {return ESy;};
	constexpr double & pmlEAy() {return EAy;};

	constexpr double & pmlHKy() {return HKy;};
	constexpr double & pmlHSy() {return HSy;};
	constexpr double & pmlHAy() {return HAy;};

	// derived PML parameters
	constexpr double & pmlEBy() {return EBy;};
	constexpr double & pmlECy() {return ECy;};

	constexpr double & pmlHBy() {return HBy;};
	constexpr double & pmlHCy() {return HCy;};
};







//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************










template <typename Mode>
struct PMLIz{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};

template <>
struct PMLIz<ThreeD>{
	// E convolution terms
	double EIzx;
	double EIzy;
	double EIzz;

	// H convolution terms
	double HIzx;
	double HIzy;
	double HIzz;

	// accessors
	double & pmlEIzx() {return EIzx;};
	double & pmlEIzy() {return EIzy;};
	double & pmlEIzz() {return EIzz;};
	double & pmlHIzx() {return HIzx;};
	double & pmlHIzy() {return HIzy;};
	double & pmlHIzz() {return HIzz;};
};


template <>
struct PMLIz<TE>{
	// E convolution terms
	double EIzx;
	double EIzy;

	// H convolution terms
	double HIzz;

	// accessors
	double & pmlEIzx() {return EIzx;};
	double & pmlEIzy() {return EIzy;};
	double & pmlHIzz() {return HIzz;};
};


template <>
struct PMLIz<TM>{
	// E convolution terms
	double EIzz;

	// H convolution terms
	double HIzx;
	double HIzy;

	// accessors
	double & pmlEIzz() {return EIzz;};
	double & pmlHIzx() {return HIzx;};
	double & pmlHIzy() {return HIzy;};
};




template <>
struct PMLIz<TEM>{
	// E convolution terms
	double EIzz;

	// H convolution terms
	double HIzy;

	// accessors
	double & pmlEIzz() {return EIzz;};
	double & pmlHIzy() {return HIzy;};
};




// y PML Stored class
template <typename Mode>
struct StoredPMLz : public PMLIz<Mode>{

	// PML parameters
	double EKz;
	double ESz;
	double EAz;

	// derived PML parameters
	double EBz;
	double ECz;


	// PML parameters
	double HKz;
	double HSz;
	double HAz;

	// derived PML parameters
	double HBz;
	double HCz;


	StoredPMLz()
	: EKz(1.0), ESz(0.0), EAz(0.0)
	, EBz(1.0), ECz(0.0)
	, HKz(1.0), HSz(0.0), HAz(0.0)
	, HBz(1.0), HCz(0.0) {};

	// StoredPMLz<ThreeD>(double K, double S, double A, double dt)
	// : EKz(K), ESz(S), EAz(A)
	// , EBz(1.0), ECz(0.0)
	// , EIz(0.0), EIz(0.0), EIzz(0.0) {setPMLParametersE(EKz, ESz, EAz, dt);
	// 								   setPMLParametersH(HKz, HSz, HAz, dt);};


	void setPMLParametersEz(double K, double S, double A, double dt){
		EKz = K; ESz = S; EAz = A;
		EBz = exp(-dt/eps0*(ESz/EKz + EAz));
		// EBz = exp(-dt*(ESz/EKz + EAz));
		ECz = ESz/EKz*1.0/(ESz+EKz*EAz)*(EBz-1);
	}

	void setPMLParametersHz(double K, double S, double A, double dt){
		HKz = K; HSz = S; HAz = A;
		HBz = exp(-dt/eps0*(HSz/HKz + HAz));
		// HBz = exp(-dt*(HSz/HKz + HAz));
		HCz = HSz/HKz*1.0/(HSz+HKz*HAz)*(HBz-1);
	}


	// PML parameters
	constexpr double & pmlEKz() {return EKz;};
	constexpr double & pmlESz() {return ESz;};
	constexpr double & pmlEAz() {return EAz;};

	constexpr double & pmlHKz() {return HKz;};
	constexpr double & pmlHSz() {return HSz;};
	constexpr double & pmlHAz() {return HAz;};

	// derived PML parameters
	constexpr double & pmlEBz() {return EBz;};
	constexpr double & pmlECz() {return ECz;};

	constexpr double & pmlHBz() {return HBz;};
	constexpr double & pmlHCz() {return HCz;};
};





//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



struct PMLParameterModel{
	double mM, mMa;
	double msMax, mkMax, maMax;

	PMLParameterModel(double dx)
	: mM(3), mMa(1), msMax(0.8*(mM+1)/(imp0*dx)), mkMax(5), maMax(0.05) {};

	PMLParameterModel(double m, double ma, double sMax, double kMax, double aMax)
	: mM(m), mMa(ma), msMax(sMax), mkMax(kMax), maMax(aMax) {};
	
	double K(double x) const {
		if (x<0.0) return 1.0;
		return (x > 1.0 ? 1.0 : 1.0+pow(x,mM)*(mkMax-1.0));
	};

	double S(double x) const {
		if (x<0.0) return 0.0;
		return (x > 1.0 ? 0.0 : pow(x,mM)*msMax);
	};

	double A(double x) const {
		if (x<0.0) return 0.0;
		return (x > 1.0 ? 0.0 : pow(1.0-x, mMa)*maMax);
	};

};








}// end namespace fdtd

#endif