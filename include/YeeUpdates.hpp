#ifndef _YEEUPDATES_H
#define _YEEUPDATES_H

#include "FDTDConstants.hpp"

namespace fdtd{


struct TemporalScheme {};
struct BD2 : public TemporalScheme {};
struct BD3 : public TemporalScheme {};
struct BD4 : public TemporalScheme {};



template <int n, typename Mode, typename FieldType>
struct BDFields {
	typedef YeeFields<Mode, FieldType, std::array> 	Fields;
private: 
	Fields mFields[n];
	unsigned int mCurr;
public:
	BDFields<n, Mode, FieldType>() : mCurr(2) {};

	BDFields & BD() {return *this;};
	Fields & BD(int i) {return mFields[(n-1-mCurr+i)%n];};
	void increment() {mCurr = std::min(n-1, (mCurr-1));};
};

// Define difference operators on YeeCell objects
// which will be used to construct the YeeAlgorithm for
// any given mode

template <typename EMField, Dir d>
struct DifferenceOperator{
	static_assert(std::is_same<Field, EMField>::value, "DifferenceOperator needs a valid EMField");
};


////////////// E //////////////

template <>
struct DifferenceOperator<fdtd::Ex, Dir::Y>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.getNeighborMax(1).Ex() - f.Ex();
	}
};

template <>
struct DifferenceOperator<fdtd::Ex, Dir::Z>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.getNeighborMax(2).Ex() - f.Ex();
	}
};


template <>
struct DifferenceOperator<fdtd::Ey, Dir::X>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.getNeighborMax(0).Ey() - f.Ey();
	}
};

template <>
struct DifferenceOperator<fdtd::Ey, Dir::Z>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.getNeighborMax(2).Ey() - f.Ey();
	}
};

template <>
struct DifferenceOperator<fdtd::Ez, Dir::X>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.getNeighborMax(0).Ez() - f.Ez();
	}
};

template <>
struct DifferenceOperator<fdtd::Ez, Dir::Y>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.getNeighborMax(1).Ez() - f.Ez();
	}
};





//////////////// H ////////////////


template <>
struct DifferenceOperator<fdtd::Hx, Dir::Y>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.Hx() - f.getNeighborMin(1).Hx();
	}
};

template <>
struct DifferenceOperator<fdtd::Hx, Dir::Z>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.Hx() - f.getNeighborMin(2).Hx();
	}
};


template <>
struct DifferenceOperator<fdtd::Hy, Dir::X>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.Hy() - f.getNeighborMin(0).Hy();
	}
};

template <>
struct DifferenceOperator<fdtd::Hy, Dir::Z>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.Hy() - f.getNeighborMin(2).Hy();
	}
};

template <>
struct DifferenceOperator<fdtd::Hz, Dir::X>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.Hz() - f.getNeighborMin(0).Hz();
	}
};

template <>
struct DifferenceOperator<fdtd::Hz, Dir::Y>{
	template <class YeeCell>
	static double get(YeeCell & f) {
		return f.Hz() - f.getNeighborMin(1).Hz();
	}
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
template <class Mode, class TimeScheme = BD2>
struct YeeUpdateD{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	static_assert(std::is_same<TemporalScheme, TimeScheme>::value, "YeeUpdate needs a valid TemporalScheme");
};


// specialization for 3D
template<>
struct YeeUpdateD<ThreeD>{
	double dt, dx;

	YeeUpdateD<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dx() += dt/dx*(1.0/f.pmlEKy()*(f.Hz() - f.getNeighborMin(1).Hz()) - 1.0/f.pmlEKz()*(f.Hy() - f.getNeighborMin(2).Hy()));
		f.Dy() += dt/dx*(1.0/f.pmlEKz()*(f.Hx() - f.getNeighborMin(2).Hx()) - 1.0/f.pmlEKx()*(f.Hz() - f.getNeighborMin(0).Hz()));
		f.Dz() += dt/dx*(1.0/f.pmlEKx()*(f.Hy() - f.getNeighborMin(0).Hy()) - 1.0/f.pmlEKy()*(f.Hx() - f.getNeighborMin(1).Hx()));
	};
};

// specialization for TE
template<>
struct YeeUpdateD<TE>{
	double dt, dx;

	YeeUpdateD<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dx() += dt/dx*( 1.0/f.pmlEKy()*(f.Hz() - f.getNeighborMin(1).Hz()));
		f.Dy() += dt/dx*(-1.0/f.pmlEKx()*(f.Hz() - f.getNeighborMin(0).Hz()));
	};
};


// specialization for TM
template<>
struct YeeUpdateD<TM>{
	double dt, dx;

	YeeUpdateD<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dz() += dt/dx*(1.0/f.pmlEKx()*fdtd::DifferenceOperator<Hy, Dir::X>::get(f)
					   - 1.0/f.pmlEKy()*fdtd::DifferenceOperator<Hx, Dir::Y>::get(f));
		// f.Dz() += dt/dx*(1.0/f.pmlEKx()*(f.Hy() - f.getNeighborMin(0).Hy()) - 1.0/f.pmlEKy()*(f.Hx() - f.getNeighborMin(1).Hx()));
	};
};



// specialization for TEM
template<>
struct YeeUpdateD<TEM>{
	double dt, dx;

	YeeUpdateD<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dz() += dt/dx*(1.0/f.pmlEKx()*(f.Hy() - f.getNeighborMin(0).Hy()));
	};
};


// ************ THESE ARE THE BD4 UPDATES
// specialization for 3D
template<>
struct YeeUpdateD<ThreeD, BD4>{
	double dt, dx;

	YeeUpdateD<ThreeD, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double hDx = f.Dx();
		f.Dx() = 17.0/22.0*f.Dx() 
				+9.0/22.0*f.BD(0).Dx()
				-5.0/22.0*f.BD(1).Dx()
				+1.0/22.0*f.BD(2).Dx()
				+24.0/22.0*dt/dx*(1.0/f.pmlEKy()*fdtd::DifferenceOperator<fdtd::Hz, Dir::Y>::get(f) - 1.0/f.pmlEKz()*fdtd::DifferenceOperator<fdtd::Hy, Dir::Z>::get(f));

		double hDy = f.Dy();
		f.Dy() = 17.0/22.0*f.Dy() 
				+9.0/22.0*f.BD(0).Dy()
				-5.0/22.0*f.BD(1).Dy()
				+1.0/22.0*f.BD(2).Dy()
				+24.0/22.0*dt/dx*(1.0/f.pmlEKz()*fdtd::DifferenceOperator<fdtd::Hx, Dir::Z>::get(f) - 1.0/f.pmlEKx()*fdtd::DifferenceOperator<fdtd::Hz, Dir::X>::get(f));
		
		double hDz = f.Dz();
		f.Dz() = 17.0/22.0*f.Dz() 
				+9.0/22.0*f.BD(0).Dz()
				-5.0/22.0*f.BD(1).Dz()
				+1.0/22.0*f.BD(2).Dz()
				+24.0/22.0*dt/dx*(1.0/f.pmlEKx()*fdtd::DifferenceOperator<fdtd::Hy, Dir::X>::get(f) - 1.0/f.pmlEKy()*fdtd::DifferenceOperator<fdtd::Hx, Dir::Y>::get(f));
	
		f.BD().increment();
		f.BD(0).Dx() = hDx;
		f.BD(0).Dy() = hDy;
		f.BD(0).Dz() = hDz;
	};
};

// specialization for TE
template<>
struct YeeUpdateD<TE, BD4>{
	double dt, dx;

	YeeUpdateD<TE, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dx() += dt/dx*( 1.0/f.pmlEKy()*(f.Hz() - f.getNeighborMin(1).Hz()));
		f.Dy() += dt/dx*(-1.0/f.pmlEKx()*(f.Hz() - f.getNeighborMin(0).Hz()));
	
		double hDx = f.Dx();
		f.Dx() = 17.0/22.0*f.Dx() 
				+9.0/22.0*f.BD(0).Dx()
				-5.0/22.0*f.BD(1).Dx()
				+1.0/22.0*f.BD(2).Dx()
				+24.0/22.0*dt/dx*(1.0/f.pmlEKy()*fdtd::DifferenceOperator<fdtd::Hz, Dir::Y>::get(f));

		double hDy = f.Dy();
		f.Dy() = 17.0/22.0*f.Dy() 
				+9.0/22.0*f.BD(0).Dy()
				-5.0/22.0*f.BD(1).Dy()
				+1.0/22.0*f.BD(2).Dy()
				+24.0/22.0*dt/dx*(-1.0/f.pmlEKx()*fdtd::DifferenceOperator<fdtd::Hz, Dir::X>::get(f));
		
		f.BD().increment();
		f.BD(0).Dx() = hDx;
		f.BD(0).Dy() = hDy;
	};
};


// specialization for TM
template<>
struct YeeUpdateD<TM, BD4>{
	double dt, dx;

	YeeUpdateD<TM, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){

		double hDz = f.Dz();
		f.Dz() = 17.0/22.0*f.Dz() 
				+9.0/22.0*f.BD(0).Dz()
				-5.0/22.0*f.BD(1).Dz()
				+1.0/22.0*f.BD(2).Dz()
				+24.0/22.0*dt/dx*(1.0/f.pmlEKx()*fdtd::DifferenceOperator<fdtd::Hy, Dir::X>::get(f) - 1.0/f.pmlEKy()*fdtd::DifferenceOperator<fdtd::Hx, Dir::Y>::get(f));
	
		f.BD().increment();
		f.BD(0).Dz() = hDz;
	};
};



// specialization for TEM
template<>
struct YeeUpdateD<TEM, BD4>{
	double dt, dx;

	YeeUpdateD<TEM, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){

		double hDz = f.Dz();
		f.Dz() = 17.0/22.0*f.Dz() 
				+9.0/22.0*f.BD(0).Dz()
				-5.0/22.0*f.BD(1).Dz()
				+1.0/22.0*f.BD(2).Dz()
				+24.0/22.0*dt/dx*(1.0/f.pmlEKx()*fdtd::DifferenceOperator<fdtd::Hy, Dir::X>::get(f));
	
		f.BD().increment();
		f.BD(0).Dz() = hDz;
	};
};








//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************





// ordinary Yee updates without PML
template <class Mode>
struct PlainYeeUpdateD{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};


// specialization for 3D
template<>
struct PlainYeeUpdateD<ThreeD>{
	double dt, dx;

	PlainYeeUpdateD<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dx() += dt/dx*((f.Hz() - f.getNeighborMin(1).Hz()) - (f.Hy() - f.getNeighborMin(2).Hy()));
		f.Dy() += dt/dx*((f.Hx() - f.getNeighborMin(2).Hx()) - (f.Hz() - f.getNeighborMin(0).Hz()));
		f.Dz() += dt/dx*((f.Hy() - f.getNeighborMin(0).Hy()) - (f.Hx() - f.getNeighborMin(1).Hx()));
	};
};

// specialization for TE
template<>
struct PlainYeeUpdateD<TE>{
	double dt, dx;

	PlainYeeUpdateD<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dx() += dt/dx*( (f.Hz() - f.getNeighborMin(1).Hz()));
		f.Dy() += dt/dx*(-(f.Hz() - f.getNeighborMin(0).Hz()));
	};
};


// specialization for TM
template<>
struct PlainYeeUpdateD<TM>{
	double dt, dx;

	PlainYeeUpdateD<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dz() += dt/dx*((f.Hy() - f.getNeighborMin(0).Hy()) - (f.Hx() - f.getNeighborMin(1).Hx()));
	};
};



// specialization for TEM
template<>
struct PlainYeeUpdateD<TEM>{
	double dt, dx;

	PlainYeeUpdateD<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dz() += dt/dx*((f.Hy() - f.getNeighborMin(0).Hy()));
	};
};








//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <class Mode, class TimeScheme = BD2>
struct YeeUpdateB{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	static_assert(std::is_same<TemporalScheme, TimeScheme>::value, "YeeUpdate needs a valid TemporalScheme");
};



// specialization for 3D
template<>
struct YeeUpdateB<ThreeD>{
	double dt, dx;

	YeeUpdateB<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bx() -= dt/dx*(1.0/f.pmlHKy()*(f.getNeighborMax(1).Ez() - f.Ez()) - 1.0/f.pmlHKz()*(f.getNeighborMax(2).Ey() - f.Ey()));
		f.By() -= dt/dx*(1.0/f.pmlHKz()*(f.getNeighborMax(2).Ex() - f.Ex()) - 1.0/f.pmlHKx()*(f.getNeighborMax(0).Ez() - f.Ez()));
		f.Bz() -= dt/dx*(1.0/f.pmlHKx()*(f.getNeighborMax(0).Ey() - f.Ey()) - 1.0/f.pmlHKy()*(f.getNeighborMax(1).Ex() - f.Ex()));
	};
};


// specialization for TE
template<>
struct YeeUpdateB<TE>{
	double dt, dx;

	YeeUpdateB<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bz() -= dt/dx*(1.0/f.pmlHKx()*(f.getNeighborMax(0).Ey() - f.Ey()) - 1.0/f.pmlHKy()*(f.getNeighborMax(1).Ex() - f.Ex()));
	};
};


// specialization for TM
template<>
struct YeeUpdateB<TM>{
	double dt, dx;

	YeeUpdateB<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bx() -= dt/dx*( 1.0/f.pmlHKy()*(f.getNeighborMax(1).Ez() - f.Ez()));
		f.By() -= dt/dx*(-1.0/f.pmlHKx()*(f.getNeighborMax(0).Ez() - f.Ez()));
	};
};


// specialization for TEM
template<>
struct YeeUpdateB<TEM>{
	double dt, dx;

	YeeUpdateB<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.By() -= dt/dx*(-1.0/f.pmlHKx()*(f.getNeighborMax(0).Ez() - f.Ez()));
	};
};





// *********** THIS IS BD4 SCHEME
// specialization for 3D
template<>
struct YeeUpdateB<ThreeD, BD4>{
	double dt, dx;

	YeeUpdateB<ThreeD, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		auto hBx = f.Bx();
		f.Bx() = 17.0/22.0*f.Bx() 
				+9.0/22.0*f.BD(0).Bx()
				-5.0/22.0*f.BD(1).Bx()
				+1.0/22.0*f.BD(2).Bx()
				-24.0/22.0*dt/dx*(1.0/f.pmlHKy()*fdtd::DifferenceOperator<fdtd::Ez, Dir::Y>::get(f) - 1.0/f.pmlHKz()*fdtd::DifferenceOperator<fdtd::Ey, Dir::Z>::get(f));
		
		auto hBy = f.By();
		f.By() = 17.0/22.0*f.By() 
				+9.0/22.0*f.BD(0).By()
				-5.0/22.0*f.BD(1).By()
				+1.0/22.0*f.BD(2).By()
				-24.0/22.0*dt/dx*(1.0/f.pmlHKz()*fdtd::DifferenceOperator<fdtd::Ex, Dir::Z>::get(f) - 1.0/f.pmlHKx()*fdtd::DifferenceOperator<fdtd::Ez, Dir::X>::get(f));
		
		auto hBz = f.Bz();
		f.Bz() = 17.0/22.0*f.Bz() 
				+9.0/22.0*f.BD(0).Bz()
				-5.0/22.0*f.BD(1).Bz()
				+1.0/22.0*f.BD(2).Bz()
				-24.0/22.0*dt/dx*(1.0/f.pmlHKx()*fdtd::DifferenceOperator<fdtd::Ey, Dir::X>::get(f) - 1.0/f.pmlHKy()*fdtd::DifferenceOperator<fdtd::Ex, Dir::Y>::get(f));
	
		f.BD().increment();
		f.BD(0).Bx() = hBx;
		f.BD(0).By() = hBy;
		f.BD(0).Bz() = hBz;
	};
};


// specialization for TE
template<>
struct YeeUpdateB<TE, BD4>{
	double dt, dx;

	YeeUpdateB<TE, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bz() -= dt/dx*(1.0/f.pmlHKx()*(f.getNeighborMax(0).Ey() - f.Ey()) - 1.0/f.pmlHKy()*(f.getNeighborMax(1).Ex() - f.Ex()));
		auto hBz = f.Bz();
		f.Bz() = 17.0/22.0*f.Bz() 
				+9.0/22.0*f.BD(0).Bz()
				-5.0/22.0*f.BD(1).Bz()
				+1.0/22.0*f.BD(2).Bz()
				-24.0/22.0*dt/dx*(1.0/f.pmlHKx()*fdtd::DifferenceOperator<fdtd::Ey, Dir::X>::get(f) - 1.0/f.pmlHKy()*fdtd::DifferenceOperator<fdtd::Ex, Dir::Y>::get(f));
	
		f.BD().increment();
		f.BD(0).Bz() = hBz;
	};
};


// specialization for TM
template<>
struct YeeUpdateB<TM, BD4>{
	double dt, dx;

	YeeUpdateB<TM, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bx() -= dt/dx*( 1.0/f.pmlHKy()*(f.getNeighborMax(1).Ez() - f.Ez()));
		f.By() -= dt/dx*(-1.0/f.pmlHKx()*(f.getNeighborMax(0).Ez() - f.Ez()));
	
		auto hBx = f.Bx();
		f.Bx() = 17.0/22.0*f.Bx() 
				+9.0/22.0*f.BD(0).Bx()
				-5.0/22.0*f.BD(1).Bx()
				+1.0/22.0*f.BD(2).Bx()
				-24.0/22.0*dt/dx*(1.0/f.pmlHKy()*fdtd::DifferenceOperator<fdtd::Ez, Dir::Y>::get(f));
		
		auto hBy = f.By();
		f.By() = 17.0/22.0*f.By() 
				+9.0/22.0*f.BD(0).By()
				-5.0/22.0*f.BD(1).By()
				+1.0/22.0*f.BD(2).By()
				-24.0/22.0*dt/dx*(-1.0/f.pmlHKx()*fdtd::DifferenceOperator<fdtd::Ez, Dir::X>::get(f));
		
		f.BD().increment();
		f.BD(0).Bx() = hBx;
		f.BD(0).By() = hBy;
	};
};


// specialization for TEM
template<>
struct YeeUpdateB<TEM, BD4>{
	double dt, dx;

	YeeUpdateB<TEM, BD4>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.By() -= dt/dx*(-1.0/f.pmlHKx()*(f.getNeighborMax(0).Ez() - f.Ez()));

		auto hBy = f.By();
		f.By() = 17.0/22.0*f.By() 
				+9.0/22.0*f.BD(0).By()
				-5.0/22.0*f.BD(1).By()
				+1.0/22.0*f.BD(2).By()
				-24.0/22.0*dt/dx*(-1.0/f.pmlHKx()*fdtd::DifferenceOperator<fdtd::Ez, Dir::X>::get(f));
		

		f.BD().increment();
		f.BD(0).By() = hBy;
	};
};







//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************







template <class Mode>
struct PlainYeeUpdateB{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};



// specialization for 3D
template<>
struct PlainYeeUpdateB<ThreeD>{
	double dt, dx;

	PlainYeeUpdateB<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bx() -= dt/dx*(f.getNeighborMax(1).Ez() - f.Ez() - f.getNeighborMax(2).Ey() + f.Ey());
		f.By() -= dt/dx*(f.getNeighborMax(2).Ex() - f.Ex() - f.getNeighborMax(0).Ez() + f.Ez());
		f.Bz() -= dt/dx*(f.getNeighborMax(0).Ey() - f.Ey() - f.getNeighborMax(1).Ex() + f.Ex());
	};
};


// specialization for TE
template<>
struct PlainYeeUpdateB<TE>{
	double dt, dx;

	PlainYeeUpdateB<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bz() -= dt/dx*(f.getNeighborMax(0).Ey() - f.Ey() - f.getNeighborMax(1).Ex() + f.Ex());
	};
};


// specialization for TM
template<>
struct PlainYeeUpdateB<TM>{
	double dt, dx;

	PlainYeeUpdateB<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bx() -= dt/dx*( f.getNeighborMax(1).Ez() - f.Ez());
		f.By() -= dt/dx*(- f.getNeighborMax(0).Ez() + f.Ez());
	};
};


// specialization for TEM
template<>
struct PlainYeeUpdateB<TEM>{
	double dt, dx;

	PlainYeeUpdateB<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.By() -= dt/dx*(- f.getNeighborMax(0).Ez() + f.Ez());
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