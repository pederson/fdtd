#ifndef _PML_H
#define _PML_H

#include <array>

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
struct UpdatePMLD<TEM>{
	double dt, dx;

	UpdatePMLD<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.pmlEIxy() = f.pmlEBx()*f.pmlEIxy() + f.pmlECx()/dx*(f.Hy() - f.getNeighborMin(0).Hy());	
		f.Dz() += f.pmlEIxy();
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
struct PMLIx<TM>{
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
struct PMLIx<TEM>{
	// E convolution terms
	double EIxy;

	// H convolution terms
	double HIxz;

	// accessors
	double & pmlEIxy() {return EIxy;};
	double & pmlHIxz() {return HIxz;};
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


	void setPMLParametersE(double K, double S, double A, double dt){
		EKx = K; ESx = S; EAx = A;
		EBx = exp(-dt/eps0*(ESx/EKx + EAx));
		ECx = ESx/EKx*1.0/(ESx+EKx*EAx)*(EBx-1);
	}

	void setPMLParametersH(double K, double S, double A, double dt){
		HKx = K; HSx = S; HAx = A;
		HBx = exp(-dt/eps0*(HSx/HKx + HAx));
		HCx = HSx/HKx*1.0/(HSx+HKx*HAx)*(HBx-1);
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
struct PMLIy<TM>{
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
struct PMLIy<TEM>{
	// E convolution terms
	double EIyy;

	// H convolution terms
	double HIyz;

	// accessors
	double & pmlEIyy() {return EIyy;};
	double & pmlHIyz() {return HIyz;};
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


	void setPMLParametersE(double K, double S, double A, double dt){
		EKy = K; ESy = S; EAy = A;
		EBy = eyp(-dt/eps0*(ESy/EKy + EAy));
		ECy = ESy/EKy*1.0/(ESy+EKy*EAy)*(EBy-1);
	}

	void setPMLParametersH(double K, double S, double A, double dt){
		HKy = K; HSy = S; HAy = A;
		HBy = eyp(-dt/eps0*(HSy/HKy + HAy));
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
struct PMLIz<TM>{
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
struct PMLIz<TEM>{
	// E convolution terms
	double EIzy;

	// H convolution terms
	double HIzz;

	// accessors
	double & pmlEIzy() {return EIzy;};
	double & pmlHIzz() {return HIzz;};
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


	void setPMLParametersE(double K, double S, double A, double dt){
		EKz = K; ESz = S; EAz = A;
		EBz = ezp(-dt/eps0*(ESz/EKz + EAz));
		ECz = ESz/EKz*1.0/(ESz+EKz*EAz)*(EBz-1);
	}

	void setPMLParametersH(double K, double S, double A, double dt){
		HKz = K; HSz = S; HAz = A;
		HBz = ezp(-dt/eps0*(HSz/HKz + HAz));
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
	
	double K(double x) const {return 1.0+pow(x,mM)*(mkMax-1.0);};
	double S(double x) const {return pow(x,mM)*msMax;};
	double A(double x) const {return pow(1.0-x, mMa)*maMax;};

};








}// end namespace fdtd

#endif