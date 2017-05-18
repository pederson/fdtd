#ifndef _PML_H
#define _PML_H

#include <array>

namespace fdtd{




template <class PMLxPolicy,
		  class PMLyPolicy,
		  class PMLzPolicy>
class PML : public PMLxPolicy, public PMLyPolicy, public PMLzPolicy
{
public:
	typedef PMLxPolicy PMLxT;
	typedef PMLyPolicy PMLyT;
	typedef PMLzPolicy PMLzT;
};










//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <class Mode>
struct UpdatePML{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};


// specialization for 3D
template<>
struct UpdatePML<TEM>{
	double dt, dx;

	UpdatePML<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.pmlIxy() = f.pmlBx()*f.pmlIxy() + f.pmlCx()/dx*(f.Hy() - f.getNeighborMin(0).Hy());	
	};
};







//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************





struct NoPMLx{

	// PML parameters
	constexpr double pmlKx() const {return 1.0;};

	constexpr double pmlSx() const {return 0.0;};

	constexpr double pmlAx() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlBx() const {return 1.0;};

	constexpr double pmlCx() const {return 0.0;};


	// convolution terms
	constexpr double pmlIxx() const {return 0.0;};
	constexpr double pmlIxy() const {return 0.0;};
	constexpr double pmlIxz() const {return 0.0;};
};




struct NoPMLy{

	// PML parameters
	constexpr double pmlKy() const {return 1.0;};

	constexpr double pmlSy() const {return 0.0;};

	constexpr double pmlAy() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlBy() const {return 1.0;};

	constexpr double pmlCy() const {return 0.0;};


	// convolution terms
	constexpr double pmlIyx() const {return 0.0;};
	constexpr double pmlIyy() const {return 0.0;};
	constexpr double pmlIyz() const {return 0.0;};
};




struct NoPMLz{

	// PML parameters
	constexpr double pmlKz() const {return 1.0;};

	constexpr double pmlSz() const {return 0.0;};

	constexpr double pmlAz() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlBz() const {return 1.0;};

	constexpr double pmlCz() const {return 0.0;};


	// convolution terms
	constexpr double pmlIzx() const {return 0.0;};
	constexpr double pmlIzy() const {return 0.0;};
	constexpr double pmlIzz() const {return 0.0;};
};








//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************





template <typename Mode>
struct StoredPMLx{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};



// specialization for 3D
template <>
struct StoredPMLx<ThreeD>{

	// PML parameters
	double Kx;
	double Sx;
	double Ax;

	// derived PML parameters
	double Bx;
	double Cx;

	// convolution terms
	double Ixx;
	double Ixy;
	double Ixz;


	StoredPMLx<ThreeD>()
	: Kx(1.0), Sx(0.0), Ax(0.0)
	, Bx(1.0), Cx(0.0)
	, Ixx(0.0), Ixy(0.0), Ixz(0.0) {};

	StoredPMLx<ThreeD>(double K, double S, double A, double dt)
	: Kx(K), Sx(S), Ax(A)
	, Bx(1.0), Cx(0.0)
	, Ixx(0.0), Ixy(0.0), Ixz(0.0) {setPMLParameters(Kx, Sx, Ax, dt);};


	void setPMLParameters(double K, double S, double A, double dt){
		Kx = K; Sx = S; Ax = A;
		Bx = exp(-dt/eps0*(Sx/Kx + Ax));
		Cx = Sx/Kx*1.0/(Sx+Kx*Ax)*(Bx-1);
	}


	// PML parameters
	constexpr double & pmlKx() {return Kx;};

	constexpr double & pmlSx() {return Sx;};

	constexpr double & pmlAx() {return Ax;};

	// derived PML parameters
	constexpr double & pmlBx() {return Bx;};

	constexpr double & pmlCx() {return Cx;};


	// convolution terms
	double & pmlIxx() {return Ixx;};
	double & pmlIxy() {return Ixy;};
	double & pmlIxz() {return Ixz;};
};





// specialization for TE
template <>
struct StoredPMLx<TE>{

	// PML parameters
	double Kx;
	double Sx;
	double Ax;

	// derived PML parameters
	double Bx;
	double Cx;

	// convolution terms
	double Ixz;


	StoredPMLx<TE>()
	: Kx(1.0), Sx(0.0), Ax(0.0)
	, Bx(1.0), Cx(0.0)
	, Ixz(0.0) {};

	StoredPMLx<TE>(double K, double S, double A, double dt)
	: Kx(K), Sx(S), Ax(A)
	, Bx(1.0), Cx(0.0)
	, Ixz(0.0) {setPMLParameters(Kx, Sx, Ax, dt);};


	void setPMLParameters(double K, double S, double A, double dt){
		Kx = K; Sx = S; Ax = A;
		Bx = exp(-dt/eps0*(Sx/Kx + Ax));
		Cx = Sx/Kx*1.0/(Sx+Kx*Ax)*(Bx-1);
	}



	// PML parameters
	constexpr double & pmlKx() {return Kx;};

	constexpr double & pmlSx() {return Sx;};

	constexpr double & pmlAx() {return Ax;};

	// derived PML parameters
	constexpr double & pmlBx() {return Bx;};

	constexpr double & pmlCx() {return Cx;};


	// convolution terms
	double & pmlIxz() {return Ixz;};
};




// specialization for TM
template <>
struct StoredPMLx<TM>{

	// PML parameters
	double Kx;
	double Sx;
	double Ax;

	// derived PML parameters
	double Bx;
	double Cx;

	// convolution terms
	double Ixx;
	double Ixy;

	StoredPMLx<TM>()
	: Kx(1.0), Sx(0.0), Ax(0.0)
	, Bx(1.0), Cx(0.0)
	, Ixx(0.0), Ixy(0.0) {};

	StoredPMLx<TM>(double K, double S, double A, double dt)
	: Kx(K), Sx(S), Ax(A)
	, Bx(1.0), Cx(0.0)
	, Ixx(0.0), Ixy(0.0) {setPMLParameters(Kx, Sx, Ax, dt);};


	void setPMLParameters(double K, double S, double A, double dt){
		Kx = K; Sx = S; Ax = A;
		Bx = exp(-dt/eps0*(Sx/Kx + Ax));
		Cx = Sx/Kx*1.0/(Sx+Kx*Ax)*(Bx-1);
	}


	// PML parameters
	constexpr double & pmlKx() {return Kx;};

	constexpr double & pmlSx() {return Sx;};

	constexpr double & pmlAx() {return Ax;};

	// derived PML parameters
	constexpr double & pmlBx() {return Bx;};

	constexpr double & pmlCx() {return Cx;};


	// convolution terms
	double & pmlIxx() {return Ixx;};
	double & pmlIxy() {return Ixy;};
};



// specialization for TEM
template <>
struct StoredPMLx<TEM>{

	// PML parameters
	double Kx;
	double Sx;
	double Ax;

	// derived PML parameters
	double Bx;
	double Cx;

	// convolution terms
	double Ixy;


	StoredPMLx<TEM>()
	: Kx(1.0), Sx(0.0), Ax(0.0)
	, Bx(1.0), Cx(0.0)
	, Ixy(0.0) {};

	StoredPMLx<TEM>(double K, double S, double A, double dt)
	: Kx(K), Sx(S), Ax(A)
	, Bx(1.0), Cx(0.0)
	, Ixy(0.0) {setPMLParameters(Kx, Sx, Ax, dt);};


	void setPMLParameters(double K, double S, double A, double dt){
		Kx = K; Sx = S; Ax = A;
		Bx = exp(-dt/eps0*(Sx/Kx + Ax));
		Cx = Sx/Kx*1.0/(Sx+Kx*Ax)*(Bx-1);
	}


	// PML parameters
	constexpr double & pmlKx() {return Kx;};

	constexpr double & pmlSx() {return Sx;};

	constexpr double & pmlAx() {return Ax;};

	// derived PML parameters
	constexpr double & pmlBx() {return Bx;};

	constexpr double & pmlCx() {return Cx;};


	// convolution terms
	double & pmlIxy() {return Ixy;};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <typename Mode>
struct StoredPMLy{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template <>
struct StoredPMLy<ThreeD>{

	// PML parameters
	double Ky;
	double Sy;
	double Ay;

	// derived PML parameters
	double By;
	double Cy;

	// convolution terms
	double Iyx;
	double Iyy;
	double Iyz;



	// PML parameters
	constexpr double & pmlKy() {return Ky;};

	constexpr double & pmlSy() {return Sy;};

	constexpr double & pmlAy() {return Ay;};

	// derived PML parameters
	constexpr double & pmlBy() {return By;};

	constexpr double & pmlCy() {return Cy;};


	// convolution terms
	double & pmlIyx() {return Iyx;};
	double & pmlIyy() {return Iyy;};
	double & pmlIyz() {return Iyz;};
};





// specialization for TE
template <>
struct StoredPMLy<TE>{

	// PML parameters
	double Ky;
	double Sy;
	double Ay;

	// derived PML parameters
	double By;
	double Cy;

	// convolution terms
	double Iyz;



	// PML parameters
	constexpr double & pmlKy() {return Ky;};

	constexpr double & pmlSy() {return Sy;};

	constexpr double & pmlAy() {return Ay;};

	// derived PML parameters
	constexpr double & pmlBy() {return By;};

	constexpr double & pmlCy() {return Cy;};


	// convolution terms
	double & pmlIyz() {return Iyz;};
};




// specialization for TM
template <>
struct StoredPMLy<TM>{

	// PML parameters
	double Ky;
	double Sy;
	double Ay;

	// derived PML parameters
	double By;
	double Cy;

	// convolution terms
	double Iyx;
	double Iyy;



	// PML parameters
	constexpr double & pmlKy() {return Ky;};

	constexpr double & pmlSy() {return Sy;};

	constexpr double & pmlAy() {return Ay;};

	// derived PML parameters
	constexpr double & pmlBy() {return By;};

	constexpr double & pmlCy() {return Cy;};


	// convolution terms
	double & pmlIyx() {return Iyx;};
	double & pmlIyy() {return Iyy;};
};



// specialization for TEM
template <>
struct StoredPMLy<TEM>{

	// PML parameters
	double Ky;
	double Sy;
	double Ay;

	// derived PML parameters
	double By;
	double Cy;

	// convolution terms
	double Iyy;



	// PML parameters
	constexpr double & pmlKy() {return Ky;};

	constexpr double & pmlSy() {return Sy;};

	constexpr double & pmlAy() {return Ay;};

	// derived PML parameters
	constexpr double & pmlBy() {return By;};

	constexpr double & pmlCy() {return Cy;};


	// convolution terms
	double & pmlIyy() {return Iyy;};
};






//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <typename Mode>
struct StoredPMLz{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template <>
struct StoredPMLz<ThreeD>{

	// PML parameters
	double Kz;
	double Sz;
	double Az;

	// derived PML parameters
	double Bz;
	double Cz;

	// convolution terms
	double Izx;
	double Izy;
	double Izz;



	// PML parameters
	constexpr double & pmlKz() {return Kz;};

	constexpr double & pmlSz() {return Sz;};

	constexpr double & pmlAz() {return Az;};

	// derived PML parameters
	constexpr double & pmlBz() {return Bz;};

	constexpr double & pmlCz() {return Cz;};


	// convolution terms
	double & pmlIzx() {return Izx;};
	double & pmlIzy() {return Izy;};
	double & pmlIzz() {return Izz;};
};





// specialization for TE
template <>
struct StoredPMLz<TE>{

	// PML parameters
	double Kz;
	double Sz;
	double Az;

	// derived PML parameters
	double Bz;
	double Cz;

	// convolution terms
	double Izz;



	// PML parameters
	constexpr double & pmlKz() {return Kz;};

	constexpr double & pmlSz() {return Sz;};

	constexpr double & pmlAz() {return Az;};

	// derived PML parameters
	constexpr double & pmlBz() {return Bz;};

	constexpr double & pmlCz() {return Cz;};


	// convolution terms
	double & pmlIzz() {return Izz;};
};




// specialization for TM
template <>
struct StoredPMLz<TM>{

	// PML parameters
	double Kz;
	double Sz;
	double Az;

	// derived PML parameters
	double Bz;
	double Cz;

	// convolution terms
	double Izx;
	double Izy;



	// PML parameters
	constexpr double & pmlKz() {return Kz;};

	constexpr double & pmlSz() {return Sz;};

	constexpr double & pmlAz() {return Az;};

	// derived PML parameters
	constexpr double & pmlBz() {return Bz;};

	constexpr double & pmlCz() {return Cz;};


	// convolution terms
	double & pmlIzz() {return Izx;};
	double & pmlIzy() {return Izy;};
};



// specialization for TEM
template <>
struct StoredPMLz<TEM>{

	// PML parameters
	double Kz;
	double Sz;
	double Az;

	// derived PML parameters
	double Bz;
	double Cz;

	// convolution terms
	double Izy;



	// PML parameters
	constexpr double & pmlKz() {return Kz;};

	constexpr double & pmlSz() {return Sz;};

	constexpr double & pmlAz() {return Az;};

	// derived PML parameters
	constexpr double & pmlBz() {return Bz;};

	constexpr double & pmlCz() {return Cz;};


	// convolution terms
	double & pmlIzy() {return Izy;};
};


}// end namespace fdtd

#endif