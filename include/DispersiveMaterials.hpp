#ifndef _DISPERSIVEMATERIALS_H
#define _DISPERSIVEMATERIALS_H

#include "FDTDConstants.hpp"

namespace fdtd{


template <typename scalar_type = double>
struct VacuumPolarization{
	constexpr double permittivity_r() const {return 1.0;};
	constexpr scalar_type Px() const {return 0;};
	constexpr scalar_type Py() const {return 0;};
	constexpr scalar_type Pz() const {return 0;};
};

template <typename scalar_type = double>
struct VacuumMagnetization{
	constexpr double permeability_r() const {return 1.0;};
	constexpr scalar_type Mx() const {return 0;};
	constexpr scalar_type My() const {return 0;};
	constexpr scalar_type Mz() const {return 0;};
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <typename scalar_type = double>
struct SinglePolarization{
	scalar_type mPx, mPy, mPz;
	double mPermittivityR;

	SinglePolarization()
	: mPermittivityR(1.0)
	, mPx(0.0), mPy(0.0), mPz(0.0) {};

	const double & permittivity_r() const {return mPermittivityR;};
	const scalar_type & Px() const {return mPx;};
	const scalar_type & Py() const {return mPy;};
	const scalar_type & Pz() const {return mPz;};


	double & permittivity_r() {return mPermittivityR;}
	scalar_type & Px() {return mPx;};
	scalar_type & Py() {return mPy;};
	scalar_type & Pz() {return mPz;};
};

template <typename scalar_type = double>
struct DoublePolarization : public SinglePolarization<scalar_type>{
	scalar_type mJx, mJy, mJz;

	DoublePolarization()
	: SinglePolarization<scalar_type>(), mJx(0.0), mJy(0.0), mJz(0.0) {};


	const scalar_type & Jx() const {return mJx;};
	const scalar_type & Jy() const {return mJy;};
	const scalar_type & Jz() const {return mJz;};


	// double & permittivity_r() {return mPermittivityR;}
	scalar_type & Jx() {return mJx;};
	scalar_type & Jy() {return mJy;};
	scalar_type & Jz() {return mJz;};
};

template <typename scalar_type = double>
struct SingleMagnetization{
	scalar_type mMx, mMy, mMz;
	double mPermeabilityR;

	SingleMagnetization()
	: mPermeabilityR(1.0)
	, mMx(0.0), mMy(0.0), mMz(0.0) {};

	const double & permeability_r() const {return mPermeabilityR;};
	const scalar_type & Mx() const {return mMx;};
	const scalar_type & My() const {return mMy;};
	const scalar_type & Mz() const {return mMz;};

	double & permeability_r() {return mPermeabilityR;}
	scalar_type & Mx() {return mMx;};
	scalar_type & My() {return mMy;};
	scalar_type & Mz() {return mMz;};
};

template <typename scalar_type = double>
struct DoubleMagnetization : public SingleMagnetization<scalar_type>{
	scalar_type mKx, mKy, mKz;

	DoubleMagnetization()
	: SingleMagnetization<scalar_type>(), mKx(0.0), mKy(0.0), mKz(0.0) {};

	const scalar_type & Kx() const {return mKx;};
	const scalar_type & Ky() const {return mKy;};
	const scalar_type & Kz() const {return mKz;};

	scalar_type & Kx() {return mKx;};
	scalar_type & Ky() {return mKy;};
	scalar_type & Kz() {return mKz;};
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

    struct Jx : public EType{
      static constexpr offset_type    off = {0.0, -0.5, -0.5};
      static constexpr const char *   name = "Jx";
      static constexpr Dir            direction = Dir::X;
  };
    struct Jy : public EType{
      static constexpr offset_type    off = {-0.5, 0.0, -0.5};
      static constexpr const char *   name = "Jy";
      static constexpr Dir            direction = Dir::Y;
  };
    struct Jz : public EType{
      static constexpr offset_type    off = {-0.5, -0.5, 0.0};
      static constexpr const char *   name = "Jz";
      static constexpr Dir            direction = Dir::Z;
  };
    struct Px : public EType{
      static constexpr offset_type    off = Jx::off;
      static constexpr const char *   name = "Px";
      static constexpr Dir            direction = Dir::X;
    };
    struct Py : public EType{
      static constexpr offset_type    off = Jy::off;
      static constexpr const char *   name = "Py";
      static constexpr Dir            direction = Dir::Y;
    };
    struct Pz : public EType{
      static constexpr offset_type    off = Jz::off;
      static constexpr const char *   name = "Pz";
      static constexpr Dir            direction = Dir::Z;
    };


    struct Kx : public HType{
      static constexpr offset_type    off = {-0.5, 0.0, 0.0};
      static constexpr const char *   name = "Kx";
      static constexpr Dir            direction = Dir::X;
  };
    struct Ky : public HType{
      static constexpr offset_type    off = {0.0, -0.5, 0.0};
      static constexpr const char *   name = "Ky";
      static constexpr Dir            direction = Dir::Y;
  };
    struct Kz : public HType{
      static constexpr offset_type    off = {0.0, 0.0, -0.5};
      static constexpr const char *   name = "Kz";
      static constexpr Dir            direction = Dir::Z;
  };
    struct Mx : public HType{
      static constexpr offset_type    off = Kx::off;
      static constexpr const char *   name = "Mx";
      static constexpr Dir            direction = Dir::X;
    };
    struct My : public HType{
      static constexpr offset_type    off = Ky::off;
      static constexpr const char *   name = "My";
      static constexpr Dir            direction = Dir::Y;
    };
    struct Mz : public HType{
      static constexpr offset_type    off = Kz::off;
      static constexpr const char *   name = "Mz";
      static constexpr Dir            direction = Dir::Z;
    };


template <FieldType ft, Dir d>
struct CurrentComponent{};

template <> struct CurrentComponent<FieldType::Electric, Dir::X>{typedef Jx type;};
template <> struct CurrentComponent<FieldType::Electric, Dir::Y>{typedef Jy type;};
template <> struct CurrentComponent<FieldType::Electric, Dir::Z>{typedef Jz type;};

template <> struct CurrentComponent<FieldType::Magnetic, Dir::X>{typedef Kx type;};
template <> struct CurrentComponent<FieldType::Magnetic, Dir::Y>{typedef Ky type;};
template <> struct CurrentComponent<FieldType::Magnetic, Dir::Z>{typedef Kz type;};


template <FieldType ft, Dir d>
struct PolarizationComponent{};

template <> struct PolarizationComponent<FieldType::Electric, Dir::X>{typedef Px type;};
template <> struct PolarizationComponent<FieldType::Electric, Dir::Y>{typedef Py type;};
template <> struct PolarizationComponent<FieldType::Electric, Dir::Z>{typedef Pz type;};

template <> struct PolarizationComponent<FieldType::Magnetic, Dir::X>{typedef Mx type;};
template <> struct PolarizationComponent<FieldType::Magnetic, Dir::Y>{typedef My type;};
template <> struct PolarizationComponent<FieldType::Magnetic, Dir::Z>{typedef Mz type;};



template <>
struct GetField<fdtd::Jx>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Jx();}
};

template <>
struct GetField<fdtd::Jy>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Jy();}
};

template <>
struct GetField<fdtd::Jz>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Jz();}
};



template <>
struct GetField<fdtd::Px>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Px();}
};

template <>
struct GetField<fdtd::Py>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Py();}
};

template <>
struct GetField<fdtd::Pz>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Pz();}
};


template <>
struct GetField<fdtd::Kx>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Kx();}
};

template <>
struct GetField<fdtd::Ky>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Ky();}
};

template <>
struct GetField<fdtd::Kz>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Kz();}
};


template <>
struct GetField<fdtd::Mx>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Mx();}
};

template <>
struct GetField<fdtd::My>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.My();}
};

template <>
struct GetField<fdtd::Mz>{
	template <class YeeCell>
	static decltype(auto) get(YeeCell & f) {return f.Mz();}
};

//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

struct StoredValue{
	double mVal;

	StoredValue(double K) : mVal(K) {};

	template <typename T>
	static double get(T && f, StoredValue & s){return s.mVal;};
};


template <typename Tipo>
struct Stored{
	double mVal;

	Stored(double K) : mVal(K) {};

	template <typename T>
	static double get(T && f, Stored & s){return s.mVal;};
};

//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


namespace constant{
	struct ConstantUpdate{


		template <typename T>
		static void explicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double eps_0){
			D = eps_rel*eps_0*E;
		}


	//////////////// Verlet (2nd order) ////////////////////////

		template <typename T>
		static void implicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double eps_0){

			E = D/(eps_rel*eps_0);

		}

	};

	template <typename EMField>
	struct ConstantImplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double eps_0){
			ConstantUpdate::implicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
											  eps_r, eps_0); 
		}
	};

	template <typename EMField>
	struct ConstantExplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double eps_0){
			ConstantUpdate::explicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
											  eps_r, eps_0); 
		}
	};
}// end namespace constant






// This is very generalized and pretty good, but it would be convenient to have a type
// that could store the material parameters within the struct itself
template <class Mode, 
		  FieldType ftype = FieldType::Electric,
		  bool forward = false,
		  class StaticValue = StoredValue>
struct ConstantUpdate : public StaticValue
{
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");

	static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);

	ConstantUpdate() {};
	ConstantUpdate(double K) : StaticValue(K) {};


	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							constant::ConstantImplicit>(f, StaticValue::get(f, *this), val);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							constant::ConstantExplicit>(f, StaticValue::get(f, *this), val);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f){
		get(f);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f, double delta_t){
		get(f);
	};

};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

namespace conductive{
	struct ConductiveUpdate{


		template <typename T>
		static void explicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double eps_0, 
								   double cond_freq, 
								   double dt){
			P += dt*cond_freq*E;
			D = eps_rel*eps_0*E + P;
		}


	//////////////// Verlet (2nd order) ////////////////////////

		template <typename T>
		static void implicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double eps_0, 
								   double cond_freq, 
								   double dt){

			double factor = dt*cond_freq/(eps_rel*eps_0);

			// first update P to the level of D
			P += factor*D/(1.0+factor);
			// then update E to the level of D
			E = (D - P)/(eps_rel*eps_0);

		}

	};

	template <typename EMField>
	struct ConductiveImplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double eps_0, double wp, double dt){
			ConductiveUpdate::implicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
											  eps_r, eps_0, wp, dt); 
		}
	};

	template <typename EMField>
	struct ConductiveExplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double eps_0, double wp, double dt){
			ConductiveUpdate::explicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
											  eps_r, eps_0, wp, dt); 
		}
	};
}// end namespace conductive






// This is very generalized and pretty good, but it would be convenient to have a type
// that could store the material parameters within the struct itself
template <class Mode,
		  FieldType ftype = FieldType::Electric,
		  bool forward = false, 
		  class StaticValue = StoredValue,
		  class ConductiveFreq = Stored<void>>
struct ConductiveUpdate : public StaticValue, public ConductiveFreq
{
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");

	static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
	double dt;

	ConductiveUpdate(double deltat): dt(deltat) {};
	ConductiveUpdate(double K, double w0, double deltat): StaticValue(K), ConductiveFreq(w0), dt(deltat) {};


	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							conductive::ConductiveImplicit>(f, StaticValue::get(f, *this), val, ConductiveFreq::get(f, *this), dt);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							conductive::ConductiveExplicit>(f, StaticValue::get(f, *this), val, ConductiveFreq::get(f, *this), dt);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f){
		get(f);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f, double delta_t){
		dt = delta_t;
		get(f);
	};

};






//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

namespace lorentz{
	struct LorentzUpdate{


		template <typename T>
		static void explicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double delta,
								   double eps_0, 
								   double lorentz_freq, 
								   double gamma,
								   double dt){
			static_assert(std::is_same<T,void>::value, "EXPLICIT UPDATE FOR LORENTZ MATERIAL IS INCOMPLETE!");
		}


	//////////////// Verlet (2nd order) ////////////////////////

		template <typename T>
		static void implicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double delta,
								   double eps_0, 
								   double lorentz_freq, 
								   double gamma,
								   double dt){

			double b = 2.0*gamma*dt;
			double c = dt*lorentz_freq*lorentz_freq*(1.0+delta/eps_rel);

			// first update P to the level of D
			P += 1.0/(b+dt*c)*(b*P + dt*J);

			// finally update J a half step past D
			J += 1.0/(b+dt*c)*(-c*P + J);

			// then update E to the level of D
			E = (D - P)/(eps_rel*eps_0);
		}

	};

	template <typename EMField>
	struct LorentzImplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double delta, double eps_0, double wp, double gamma, double dt){
			LorentzUpdate::implicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
										 eps_r, delta, eps_0, wp, gamma, dt); 
		}
	};

	template <typename EMField>
	struct LorentzExplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double delta, double eps_0, double wp, double gamma, double dt){
			LorentzUpdate::explicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
										 eps_r, delta, eps_0, wp, gamma, dt); 
		}
	};
}// end namespace lorentz






// This is very generalized and pretty good, but it would be convenient to have a type
// that could store the material parameters within the struct itself
template <class Mode,
		  FieldType ftype = FieldType::Electric,
		  bool forward = false, 
		  class StaticValue = StoredValue,
		  class Delta 		= Stored<void>,
		  class LorentzFreq = Stored<int>,
		  class Gamma 		= Stored<double>
		  >
struct LorentzUpdate : public StaticValue, public Delta, public LorentzFreq, public Gamma
{
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");

	static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
	double dt;

	LorentzUpdate(double deltat): dt(deltat) {};
	LorentzUpdate(double K, double del, double w0, double gamma, double deltat): StaticValue(K), Delta(del), LorentzFreq(w0), Gamma(gamma), dt(deltat) {};


	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							lorentz::LorentzImplicit>(f, StaticValue::get(f, *this), Delta::get(f, *this), val, LorentzFreq::get(f, *this), Gamma::get(f, *this), dt);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							lorentz::LorentzExplicit>(f, StaticValue::get(f, *this), Delta::get(f, *this), val, LorentzFreq::get(f, *this), Gamma::get(f, *this), dt);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f){
		get(f);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f, double delta_t){
		dt = delta_t;
		get(f);
	};

};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


namespace drude{
	struct DrudeUpdate{

	//////////////// implicit Euler (1st order) ////////////////////////
		// // using implicit Euler (1st-order in time)
		// template <typename T>
		// static void implicit_update(T& D, T& E, T& P, T& J,
		// 						   double eps_rel,
		// 						   double eps_0, 
		// 						   double drude_freq, 
		// 						   double gamma,
		// 						   double dt){
		// 	double b = gamma*dt;
		// 	double c = dt*drude_freq*drude_freq/eps_rel;

		// 	J = (J + c*(D - P))/(1.0+b+dt*c);
		// 	P += dt*J;
		// 	E = (D - P)/(eps_rel*eps_0);
		// }

		template <typename T>
		static void explicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double eps_0, 
								   double drude_freq, 
								   double gamma,
								   double dt){
			double b = gamma*dt;
			double c = dt*drude_freq*drude_freq*eps_0;

			P = (-dt*J + (1.0+b)*P)/(1.0+b);
			J = (J - c*E)/(1.0+b);
			D = eps_rel*eps_0*E + P;
		}


	//////////////// Crank-Nicolson (2nd order) ////////////////////////

		// // using Crank-Nicolson (2nd-order in time)
		// template <typename T>
		// static void implicit_update(T& D, T& E, T& P, T& J,
		// 						   double eps_rel,
		// 						   double eps_0, 
		// 						   double drude_freq, 
		// 						   double gamma,
		// 						   double dt){
		// 	double b = gamma*dt*0.5;
		// 	double c = 0.5*dt*drude_freq*drude_freq/eps_rel;
		// 	double d = 0.5*dt*drude_freq*drude_freq*eps_0;

		// 	auto K1 = (1.0+b)*J + d*E + c*D;
		// 	auto K2 = 0.5*dt*J + P;

		// 	double a = 1+b+0.5*dt*c;

		// 	J = 1.0/a * (K1 - c*K2);
		// 	P = 1.0/a * (0.5*dt*K1 + (1.0+b)*K2);
		// 	E = (D - P)/(eps_rel*eps_0);
		// }

		// template <typename T>
		// static void explicit_update(T& D, T& E, T& P, T& J,
		// 						   double eps_rel,
		// 						   double eps_0, 
		// 						   double drude_freq, 
		// 						   double gamma,
		// 						   double dt){
		// 	double b = gamma*dt;
		// 	double c = dt*drude_freq*drude_freq*eps_0;

		// 	P = (-dt*J + (1.0+b)*P)/(1.0+b);
		// 	J = (J - c*E)/(1.0+b);
		// 	D = eps_rel*eps_0*E + P;
		// }


	//////////////// Verlet (2nd order) ////////////////////////

		template <typename T>
		static void implicit_update(T& D, T& E, T& P, T& J,
								   double eps_rel,
								   double eps_0, 
								   double drude_freq, 
								   double gamma,
								   double dt){
			double b = gamma*dt*0.5;
			double w = dt*drude_freq; // normalized plasma freq

			// first update P to the level of D
			P += dt*J;

			// then update E to the level of D
			E = (D - P)/(eps_rel*eps_0);

			// finally update J a half step past D
			J = J*(1.0-b)/(1.0+b) + w*w/dt*eps_0*E;
			
		}


	//////////////// Recursive Convolution (1st order) i.e. exponential time-stepping ////////////////////////

		// template <typename T>
		// static void implicit_update(T& D, T& E, T& P, T& J,
		// 						   double eps_rel,
		// 						   double eps_0, 
		// 						   double drude_freq, 
		// 						   double gamma,
		// 						   double dt){
			
		// 	// P = I_1
		// 	// J = I_2

		// 	// normalized quantities
		// 	double G = gamma*dt;
		// 	double W = drude_freq*dt;

		// 	double X_1 = eps_0*W*W/G;
		// 	double X_2 = eps_0*(W/G)*(W/G)*(exp(-G)-1.0);

		// 	double R_1 = 1.0;
		// 	double R_2 = exp(-G);

		// 	// update new J and P
		// 	P = E*X_1 + R_1*P;
		// 	J = E*X_2 + R_2*J;

		// 	// update E with new values of J, P
		// 	E = (D - R_1*P - R_2*J)/(eps_rel*eps_0 + X_1 + X_2);		

		// }

	};

	template <typename EMField>
	struct DrudeImplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double eps_0, double wp, double gamma, double dt){
			DrudeUpdate::implicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
										 eps_r, eps_0, wp, gamma, dt); 
		}
	};

	template <typename EMField>
	struct DrudeExplicit{
		static constexpr Dir d = FieldDir<EMField>::value; 
		static constexpr FieldType ft = EMField::field_type;
		typedef typename FluxComponent<ft, d>::type 		DType;
		typedef typename FieldComponent<ft, d>::type 	EType;
		typedef typename CurrentComponent<ft, d>::type 	JType;
		typedef typename PolarizationComponent<ft, d>::type 	PType;

		template <typename T>
		static void get(T && f, double eps_r, double eps_0, double wp, double gamma, double dt){
			DrudeUpdate::explicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
										 eps_r, eps_0, wp, gamma, dt); 
		}
	};
}// end namespace drude



// This is very generalized and pretty good, but it would be convenient to have a type
// that could store the material parameters within the struct itself
template <class Mode,
		  FieldType ftype = FieldType::Electric,
		  bool forward = false, 
		  class StaticValue = StoredValue,
		  class DrudeFreq 	= Stored<void>,
		  class Gamma 		= Stored<int>
		  >
struct DrudeUpdate : public StaticValue, public DrudeFreq, public Gamma
{
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");

	static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
	double dt;

	DrudeUpdate(double deltat): dt(deltat) {};
	DrudeUpdate(double K, double w0, double gamma, double deltat): StaticValue(K), DrudeFreq(w0), Gamma(gamma), dt(deltat) {};


	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							drude::DrudeImplicit>(f, StaticValue::get(f, *this), val, DrudeFreq::get(f, *this), Gamma::get(f, *this), dt);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	get(YeeCell && f){
		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							drude::DrudeExplicit>(f, StaticValue::get(f, *this), val, DrudeFreq::get(f, *this), Gamma::get(f, *this), dt);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f){
		get(f);
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f, double delta_t){
		dt = delta_t;
		get(f);
	};

};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

namespace drude{
	static constexpr double m_e = 9.11e-31;
	static constexpr double q_e = 1.6e-19;
} // end namespace drude




// drude update with a magnetic field term included
// The magnetic field could be self-consistent with the wave, or it could
// be imposed externally. This behavior is controlled through the template
// class 'MagneticField' which is required to have a static function ::get
// which returns a struct with accessors .Bx(), .By(), and .Bz() 
template <class Mode, 
			class StaticValue,
			class DrudeFreq,
			class Gamma,
			class MagneticField>
struct MagnetizedDrudeUpdateParametrized{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template<class StaticValue,
			class DrudeFreq,
			class Gamma,
			class MagneticField>
struct MagnetizedDrudeUpdateParametrized<ThreeD, StaticValue, DrudeFreq, Gamma, MagneticField>{
	double dt;



	MagnetizedDrudeUpdateParametrized<ThreeD, StaticValue, DrudeFreq, Gamma, MagneticField>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell && f){


		// all updates together (impicitly)
		double wp = dt*DrudeFreq::get(f);
		double vc = dt*Gamma::get(f);
		double vx = dt*MagneticField::get(f).Bx()*drude::q_e/drude::m_e;
		double vy = dt*MagneticField::get(f).By()*drude::q_e/drude::m_e;
		double vz = dt*MagneticField::get(f).Bz()*drude::q_e/drude::m_e;
		double oneplus = (1.0+vc+wp*wp);
		double denom = oneplus*(oneplus*oneplus + (vx*vx + vy*vy + vz*vz));

		// A matrix
		// first col
		double a11 = oneplus*oneplus + vx*vx;
		double a21 = vz*oneplus+vx*vy;
		double a31 = -vy*oneplus+vx*vz;
		
		// second col
		double a12 = -vz*oneplus+vx*vy;
		double a22 = oneplus*oneplus + vy*vy;
		double a32 = vx*oneplus+vy*vz;

		// third col
		double a13 = vy*oneplus+vx*vz;
		double a23 = -vx*oneplus+vy*vz;
		double a33 = oneplus*oneplus + vz*vz;


		auto Jxh = 1.0/denom*(a11*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a12*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) + a13*(f.Jz()-wp*wp/dt*f.Pz()+wp*wp/dt*f.Dz()));
		auto Jyh = 1.0/denom*(a21*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a22*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) + a23*(f.Jz()-wp*wp/dt*f.Pz()+wp*wp/dt*f.Dz()));
		auto Jzh = 1.0/denom*(a31*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a32*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) + a33*(f.Jz()-wp*wp/dt*f.Pz()+wp*wp/dt*f.Dz()));

		auto Pxh = f.Px() + 1.0/denom*(a11*(dt*f.Jx()+wp*wp*f.Dx()-wp*wp*f.Px()) + a12*(dt*f.Jy()+wp*wp*f.Dy()-wp*wp*f.Py()) + a13*(dt*f.Jz()+wp*wp*f.Dz()-wp*wp*f.Pz()));
		auto Pyh = f.Py() + 1.0/denom*(a21*(dt*f.Jx()+wp*wp*f.Dx()-wp*wp*f.Px()) + a22*(dt*f.Jy()+wp*wp*f.Dy()-wp*wp*f.Py()) + a23*(dt*f.Jz()+wp*wp*f.Dz()-wp*wp*f.Pz()));
		auto Pzh = f.Pz() + 1.0/denom*(a31*(dt*f.Jx()+wp*wp*f.Dx()-wp*wp*f.Px()) + a32*(dt*f.Jy()+wp*wp*f.Dy()-wp*wp*f.Py()) + a33*(dt*f.Jz()+wp*wp*f.Dz()-wp*wp*f.Pz()));

		f.Jx() = Jxh;
		f.Jy() = Jyh;
		f.Jz() = Jzh;

		f.Px() = Pxh;
		f.Py() = Pyh;
		f.Pz() = Pzh;

		f.Ex() = (f.Dx() - f.Px())/(eps0*StaticValue::get(f));
		f.Ey() = (f.Dy() - f.Py())/(eps0*StaticValue::get(f));
		f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));

	};

	

};


// specialization for TM
template<class StaticValue,
			class DrudeFreq,
			class Gamma,
			class MagneticField>
struct MagnetizedDrudeUpdateParametrized<TM, StaticValue, DrudeFreq, Gamma, MagneticField>{
	double dt;



	MagnetizedDrudeUpdateParametrized<TM, StaticValue, DrudeFreq, Gamma, MagneticField>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell && f){
		double b = Gamma::get(f)*dt;
		double c = dt*DrudeFreq::get(f)*DrudeFreq::get(f)*(1.0/StaticValue::get(f));


		// all updates together (impicitly)
		double wp = dt*DrudeFreq::get(f);
		double vc = dt*Gamma::get(f);
		double vx = dt*real(MagneticField::get(f).Bx())*drude::q_e/drude::m_e;
		double vy = dt*real(MagneticField::get(f).By())*drude::q_e/drude::m_e;
		double vz = dt*real(MagneticField::get(f).Bz())*drude::q_e/drude::m_e;
		double oneplus = (1.0+vc+wp*wp);
		double denom = oneplus*(oneplus*oneplus + (vx*vx + vy*vy + vz*vz));

		// A matrix
		// first col
		double a11 = oneplus*oneplus + vx*vx;
		double a21 = vz*oneplus+vx*vy;
		double a31 = -vy*oneplus+vx*vz;
		
		// second col
		double a12 = -vz*oneplus+vx*vy;
		double a22 = oneplus*oneplus + vy*vy;
		double a32 = vx*oneplus+vy*vz;

		// third col
		double a13 = vy*oneplus+vx*vz;
		double a23 = -vx*oneplus+vy*vz;
		double a33 = oneplus*oneplus + vz*vz;


		// auto Jxh = 1.0/denom*(a11*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a12*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) + a13*(f.Jz()-wp*wp/dt*f.Pz()+wp*wp/dt*f.Dz()));
		// auto Jyh = 1.0/denom*(a21*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a22*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) + a23*(f.Jz()-wp*wp/dt*f.Pz()+wp*wp/dt*f.Dz()));
		auto Jzh = 1.0/denom*( a33*(f.Jz()-wp*wp/dt*f.Pz()+wp*wp/dt*f.Dz()));

		// auto Pxh = f.Px() + 1.0/denom*(a11*(dt*f.Jx()+wp*wp*f.Dx()-wp*wp*f.Px()) + a12*(dt*f.Jy()+wp*wp*f.Dy()-wp*wp*f.Py()) + a13*(dt*f.Jz()+wp*wp*f.Dz()-wp*wp*f.Pz()));
		// auto Pyh = f.Py() + 1.0/denom*(a21*(dt*f.Jx()+wp*wp*f.Dx()-wp*wp*f.Px()) + a22*(dt*f.Jy()+wp*wp*f.Dy()-wp*wp*f.Py()) + a23*(dt*f.Jz()+wp*wp*f.Dz()-wp*wp*f.Pz()));
		auto Pzh = f.Pz() + 1.0/denom*( a33*(dt*f.Jz()+wp*wp*f.Dz()-wp*wp*f.Pz()));

		// f.Jx() = Jxh;
		// f.Jy() = Jyh;
		f.Jz() = Jzh;

		// f.Px() = Pxh;
		// f.Py() = Pyh;
		f.Pz() = Pzh;

		// f.Ex() = (f.Dx() - f.Px())/(eps0*StaticValue::get(f));
		// f.Ey() = (f.Dy() - f.Py())/(eps0*StaticValue::get(f));
		f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));

	};

};


// specialization for TE
template<class StaticValue,
			class DrudeFreq,
			class Gamma,
			class MagneticField>
struct MagnetizedDrudeUpdateParametrized<TE, StaticValue, DrudeFreq, Gamma, MagneticField>{
	double dt;



	MagnetizedDrudeUpdateParametrized<TE, StaticValue, DrudeFreq, Gamma, MagneticField>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell && f){


		// FIXME: should incorporate static value into these constants somehow...
		// I did the math without it, but really should have incorporated it in the first place
		// I think it just changes the oneplus value from 1+... to eps+...


		// all updates together (impicitly)
		double wp = dt*DrudeFreq::get(f);
		double vc = dt*Gamma::get(f);
		double vx = dt*real(MagneticField::get(f).Bx())*drude::q_e/drude::m_e;
		double vy = dt*real(MagneticField::get(f).By())*drude::q_e/drude::m_e;
		double vz = dt*real(MagneticField::get(f).Bz())*drude::q_e/drude::m_e;
		double oneplus = (1.0+vc+wp*wp);
		double denom = oneplus*(oneplus*oneplus + (vx*vx + vy*vy + vz*vz));

		// A matrix
		// first col
		double a11 = oneplus*oneplus + vx*vx;
		double a21 = vz*oneplus+vx*vy;
		double a31 = -vy*oneplus+vx*vz;
		
		// second col
		double a12 = -vz*oneplus+vx*vy;
		double a22 = oneplus*oneplus + vy*vy;
		double a32 = vx*oneplus+vy*vz;

		// third col
		double a13 = vy*oneplus+vx*vz;
		double a23 = -vx*oneplus+vy*vz;
		double a33 = oneplus*oneplus + vz*vz;




		auto Jxh = 1.0/denom*(a11*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a12*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) );
		auto Jyh = 1.0/denom*(a21*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a22*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) );
		// auto Jzh = 1.0/denom*(a31*(f.Jx()-wp*wp/dt*f.Px()+wp*wp/dt*f.Dx()) + a32*(f.Jy()-wp*wp/dt*f.Py()+wp*wp/dt*f.Dy()) + a33*(f.Jz()-wp*wp/dt*f.Pz()+wp*wp/dt*f.Dz()));

		auto Pxh = f.Px() + 1.0/denom*(a11*(dt*f.Jx()+wp*wp*f.Dx()-wp*wp*f.Px()) + a12*(dt*f.Jy()+wp*wp*f.Dy()-wp*wp*f.Py()) );
		auto Pyh = f.Py() + 1.0/denom*(a21*(dt*f.Jx()+wp*wp*f.Dx()-wp*wp*f.Px()) + a22*(dt*f.Jy()+wp*wp*f.Dy()-wp*wp*f.Py()) );
		// auto Pzh = f.Pz() + 1.0/denom*(a31*(dt*f.Jx()+wp*wp/dt*f.Dx()-wp*wp*f.Px()) + a32*(dt*f.Jy()+wp*wp/dt*f.Dy()-wp*wp*f.Py()) + a33*(dt*f.Jz()+wp*wp/dt*f.Dz()-wp*wp*f.Pz()));


		f.Jx() = Jxh;
		f.Jy() = Jyh;
		// f.Jz() = Jzh;

		f.Px() = Pxh;
		f.Py() = Pyh;
		// f.Pz() = Pzh;

		f.Ex() = (f.Dx() - f.Px())/(eps0*StaticValue::get(f));
		f.Ey() = (f.Dy() - f.Py())/(eps0*StaticValue::get(f));
		// f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));
	};

};


// specialization for 1D TEM
template<class StaticValue,
			class DrudeFreq,
			class Gamma,
			class MagneticField>
struct MagnetizedDrudeUpdateParametrized<TEM, StaticValue, DrudeFreq, Gamma, MagneticField>{
	double dt;



	MagnetizedDrudeUpdateParametrized<TEM, StaticValue, DrudeFreq, Gamma, MagneticField>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell && f){
		double b = Gamma::get(f)*dt;
		double c = dt*DrudeFreq::get(f)*DrudeFreq::get(f)*(1.0/StaticValue::get(f));

		auto Pzhold = f.Pz();
		f.Pz() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Pz()) + dt*(f.Jz()+c*f.Dz()));
		f.Jz() = 1.0/(1.0+b+dt*c)*(-c*(Pzhold) + (f.Jz()+c*f.Dz()));

		f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));
	};

};


}// end namespace fdtd

#endif