#ifndef _YEEUPDATES_H
#define _YEEUPDATES_H

#include "FDTDConstants.hpp"

namespace fdtd{




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











template <class Mode>
struct YeeUpdateD{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};


// specialization for 3D
template<>
struct YeeUpdateD<ThreeD>{
	double dt, dx;

	YeeUpdateD<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dx() += dt/dx*(f.Hz() - f.getNeighborMin(1).Hz() - f.Hy() + f.getNeighborMin(2).Hy());
		f.Dy() += dt/dx*(f.Hx() - f.getNeighborMin(2).Hx() - f.Hz() + f.getNeighborMin(0).Hz());
		f.Dz() += dt/dx*(f.Hy() - f.getNeighborMin(0).Hy() - f.Hx() + f.getNeighborMin(1).Hx());
	};
};

// specialization for TE
template<>
struct YeeUpdateD<TE>{
	double dt, dx;

	YeeUpdateD<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dx() += dt/dx*( f.Hz() - f.getNeighborMin(1).Hz());
		f.Dy() += dt/dx*(-f.Hz() + f.getNeighborMin(0).Hz());
	};
};


// specialization for TM
template<>
struct YeeUpdateD<TM>{
	double dt, dx;

	YeeUpdateD<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dz() += dt/dx*(f.Hy() - f.getNeighborMin(0).Hy() - f.Hx() + f.getNeighborMin(1).Hx());
	};
};



// specialization for TEM
template<>
struct YeeUpdateD<TEM>{
	double dt, dx;

	YeeUpdateD<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Dz() += dt/dx*(f.Hy() - f.getNeighborMin(0).Hy());
	};
};








//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************








template <class Mode>
struct YeeUpdateB{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};



// specialization for 3D
template<>
struct YeeUpdateB<ThreeD>{
	double dt, dx;

	YeeUpdateB<ThreeD>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bx() -= dt/dx*(f.getNeighborMax(1).Ez() - f.Ez() - f.getNeighborMax(2).Ey() + f.Ey());
		f.By() -= dt/dx*(f.getNeighborMax(2).Ex() - f.Ex() - f.getNeighborMax(0).Ez() + f.Ez());
		f.Bz() -= dt/dx*(f.getNeighborMax(0).Ey() - f.Ey() - f.getNeighborMax(1).Ex() + f.Ex());
	};
};


// specialization for TE
template<>
struct YeeUpdateB<TE>{
	double dt, dx;

	YeeUpdateB<TE>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bz() -= dt/dx*(f.getNeighborMax(0).Ey() - f.Ey() - f.getNeighborMax(1).Ex() + f.Ex());
	};
};


// specialization for TM
template<>
struct YeeUpdateB<TM>{
	double dt, dx;

	YeeUpdateB<TM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Bx() -= dt/dx*( f.getNeighborMax(1).Ez() - f.Ez());
		f.By() -= dt/dx*(- f.getNeighborMax(0).Ez() + f.Ez());
	};
};


// specialization for TEM
template<>
struct YeeUpdateB<TEM>{
	double dt, dx;

	YeeUpdateB<TEM>(double deltat, double deltax): dt(deltat), dx(deltax) {};

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