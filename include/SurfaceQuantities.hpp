#ifndef _SURFACEQUANTITIES_H
#define _SURFACEQUANTITIES_H

#include <array>
#include <algorithm>

#include "FDTDConstants.hpp"

namespace fdtd{




template <typename Mode>
struct SurfaceQuantities{
	static_assert(std::is_same<EMMode, Mode>::value, "SurfaceQuantities needs a valid Mode");
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <> struct SurfaceQuantities<ThreeD>{
	double mJex, mJey, mJez;
	double mJhx, mJhy, mJhz;
	double mRho;

	constexpr SurfaceQuantities<ThreeD>()
	: mJex(0.0), mJey(0.0), mJez(0.0)
	, mJhx(0.0), mJhy(0.0), mJhz(0.0)
	, mRho(0.0) {};

	double & rho() {return mRho;};
	double & Jex() {return mJex;};
	double & Jey() {return mJey;};
	double & Jez() {return mJez;};
	double & Jhx() {return mJhx;};
	double & Jhy() {return mJhy;};
	double & Jhz() {return mJhz;};

	constexpr const double & rho() const {return mRho;};
	constexpr const double & Jex() const {return mJex;};
	constexpr const double & Jey() const {return mJey;};
	constexpr const double & Jez() const {return mJez;};
	constexpr const double & Jhx() const {return mJhx;};
	constexpr const double & Jhy() const {return mJhy;};
	constexpr const double & Jhz() const {return mJhz;};
};






//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <> struct SurfaceQuantities<TE>{
	double mJex, mJey;
	double mJhz;
	double mRho;

	constexpr SurfaceQuantities<TE>()
	: mJex(0.0), mJey(0.0)
	, mJhz(0.0)
	, mRho(0.0) {};

	double & rho() {return mRho;};
	double & Jex() {return mJex;};
	double & Jey() {return mJey;};
	double & Jhz() {return mJhz;};

	constexpr const double & rho() const {return mRho;};
	constexpr const double & Jex() const {return mJex;};
	constexpr const double & Jey() const {return mJey;};
	constexpr const double & Jhz() const {return mJhz;};
};





//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <> struct SurfaceQuantities<TM>{
	double mJez;
	double mJhx, mJhy;
	double mRho;

	constexpr SurfaceQuantities<TM>()
	: mJez(0.0)
	, mJhx(0.0), mJhy(0.0)
	, mRho(0.0) {};

	double & rho() {return mRho;};
	double & Jez() {return mJez;};
	double & Jhx() {return mJhx;};
	double & Jhy() {return mJhy;};

	constexpr const double & rho() const {return mRho;};
	constexpr const double & Jez() const {return mJez;};
	constexpr const double & Jhx() const {return mJhx;};
	constexpr const double & Jhy() const {return mJhy;};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <> struct SurfaceQuantities<TEM>{
	double mJez;
	double mJhy;
	double mRho;

	constexpr SurfaceQuantities<TEM>()
	: mJez(0.0)
	, mJhy(0.0)
	, mRho(0.0) {};

	double & rho() {return mRho;};
	double & Jez() {return mJez;};
	double & Jhy() {return mJhy;};

	constexpr const double & rho() const {return mRho;};
	constexpr const double & Jez() const {return mJez;};
	constexpr const double & Jhy() const {return mJhy;};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <class Mode>
struct UpdateSurfaceChargeDensity{
	static_assert(std::is_same<EMMode, Mode>::value, "UpdateSurfaceChargeDensity needs a valid Mode");
	
};


template <>
struct UpdateSurfaceChargeDensity<ThreeD>{
	double dx;
	UpdateSurfaceChargeDensity<ThreeD>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.rho() = fdtd::eps0/dx*(f.Ex() + f.Ey() + f.Ez() - f.getNeighborMin(0).Ex() - f.getNeighborMin(1).Ey() - f.getNeighborMin(2).Ez());
	}
};


template <>
struct UpdateSurfaceChargeDensity<TE>{
	double dx;
	UpdateSurfaceChargeDensity<TE>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.rho() = fdtd::eps0/dx*(f.Ex() + f.Ey() - f.getNeighborMin(0).Ex() - f.getNeighborMin(1).Ey());
	}
};


template <>
struct UpdateSurfaceChargeDensity<TM>{
	double dx;
	UpdateSurfaceChargeDensity<TM>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){}
};


template <>
struct UpdateSurfaceChargeDensity<TEM>{
	double dx;
	UpdateSurfaceChargeDensity<TEM>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){}
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <class Mode>
struct UpdateSurfaceCurrentE{
	static_assert(std::is_same<EMMode, Mode>::value, "UpdateSurfaceCurrentE needs a valid Mode");
	
};


template <>
struct UpdateSurfaceCurrentE<ThreeD>{
	double dx;
	UpdateSurfaceCurrentE<ThreeD>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jex() = 1.0/dx*(f.Hz() - f.getNeighborMin(1).Hz() - f.Hy() + f.getNeighborMin(2).Hy());
		f.Jey() = 1.0/dx*(f.Hx() - f.getNeighborMin(2).Hx() - f.Hz() + f.getNeighborMin(0).Hz());
		f.Jez() = 1.0/dx*(f.Hy() - f.getNeighborMin(0).Hy() - f.Hx() + f.getNeighborMin(1).Hx());
	}
};



template <>
struct UpdateSurfaceCurrentE<TE>{
	double dx;
	UpdateSurfaceCurrentE<TE>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jex() = 1.0/dx*(f.Hz() - f.getNeighborMin(1).Hz());
		f.Jey() = 1.0/dx*(- f.Hz() + f.getNeighborMin(0).Hz());
	}
};


template <>
struct UpdateSurfaceCurrentE<TM>{
	double dx;
	UpdateSurfaceCurrentE<TM>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jez() = 1.0/dx*(f.Hy() - f.getNeighborMin(0).Hy() - f.Hx() + f.getNeighborMin(1).Hx());
	}
};


template <>
struct UpdateSurfaceCurrentE<TEM>{
	double dx;
	UpdateSurfaceCurrentE<TEM>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jez() = 1.0/dx*(f.Hy() - f.getNeighborMin(0).Hy());
	}
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <class Mode>
struct UpdateSurfaceCurrentH{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	
};



template <>
struct UpdateSurfaceCurrentH<ThreeD>{
	double dx;
	UpdateSurfaceCurrentH<ThreeD>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jhx() = -1.0/dx*(f.getNeighborMax(1).Ez() - f.Ez() - f.getNeighborMax(2).Ey() + f.Ey());
		f.Jhy() = -1.0/dx*(f.getNeighborMax(2).Ex() - f.Ex() - f.getNeighborMax(0).Ez() + f.Ez());
		f.Jhz() = -1.0/dx*(f.getNeighborMax(0).Ey() - f.Ey() - f.getNeighborMax(1).Ex() + f.Ex());
	}
};



template <>
struct UpdateSurfaceCurrentH<TE>{
	double dx;
	UpdateSurfaceCurrentH<TE>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jhz() = -1.0/dx*(f.getNeighborMax(0).Ey() - f.Ey() - f.getNeighborMax(1).Ex() + f.Ex());
	}
};


template <>
struct UpdateSurfaceCurrentH<TM>{
	double dx;
	UpdateSurfaceCurrentH<TM>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jhx() = -1.0/dx*( f.getNeighborMax(1).Ez() - f.Ez());
		f.Jhy() = -1.0/dx*(-f.getNeighborMax(0).Ez() + f.Ez());
	}
};


template <>
struct UpdateSurfaceCurrentH<TEM>{
	double dx;
	UpdateSurfaceCurrentH<TEM>(double delta_x) : dx(delta_x){};

	template <class YeeCell>
	void operator()(YeeCell & f){
		f.Jhy() = -1.0/dx*(-f.getNeighborMax(0).Ez() + f.Ez());
	}
};






}// end namespace fdtd

#endif