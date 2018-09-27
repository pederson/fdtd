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



template <class FieldType>
inline void ConstantUpdate(FieldType & Out, const FieldType & In, double c){Out=In/c;}


template <class FieldType>
inline void ConstantCalculate(FieldType & Out, const FieldType & In, double c){Out=In*c;}


template <class Mode, bool forward = false>
struct ConstantUpdateE{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template<>
struct ConstantUpdateE<ThreeD>{
	double eps;

	ConstantUpdateE<ThreeD>(double c): eps(eps0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Ex(), f.Dx(), eps);
		ConstantUpdate(f.Ey(), f.Dy(), eps);
		ConstantUpdate(f.Ez(), f.Dz(), eps);
	};


	template<class YeeCell>
	void calculate(YeeCell & f){
		ConstantCalculate(f.Dx(), f.Ex(), eps);
		ConstantCalculate(f.Dy(), f.Ey(), eps);
		ConstantCalculate(f.Dz(), f.Ez(), eps);
	};
};

// specialization for TE
template<>
struct ConstantUpdateE<TE>{
	double eps;

	ConstantUpdateE<TE>(double c): eps(eps0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Ex(), f.Dx(), eps);
		ConstantUpdate(f.Ey(), f.Dy(), eps);
	};


	template<class YeeCell>
	void calculate(YeeCell & f){
		ConstantCalculate(f.Dx(), f.Ex(), eps);
		ConstantCalculate(f.Dy(), f.Ey(), eps);
	};
};


// specialization for TM
template<>
struct ConstantUpdateE<TM>{
	double eps;

	ConstantUpdateE<TM>(double c): eps(eps0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Ez(), f.Dz(), eps);
	};

	template<class YeeCell>
	void calculate(YeeCell & f){
		ConstantCalculate(f.Dz(), f.Ez(), eps);
	};
};

// specialization for TEM
template <bool forward>
struct ConstantUpdateE<TEM, forward>{
	double eps;

	ConstantUpdateE<TEM, forward>(double c): eps(eps0*c) {};

	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	operator()(YeeCell && f){
		ConstantUpdate(f.Ez(), f.Dz(), eps);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	operator() (YeeCell && f){ConstantCalculate(f.Dz(), f.Ez(), eps);};


	template<class YeeCell>
	void calculate(YeeCell && f){
		ConstantCalculate(f.Dz(), f.Ez(), eps);
	};
};









//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************












template <class Mode, bool forward = false>
struct ConstantUpdateH{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template<>
struct ConstantUpdateH<ThreeD>{
	double mu;

	ConstantUpdateH<ThreeD>(double c): mu(mu0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Hx(), f.Bx(), mu);
		ConstantUpdate(f.Hy(), f.By(), mu);
		ConstantUpdate(f.Hz(), f.Bz(), mu);
	};
};


// specialization for TE
template<>
struct ConstantUpdateH<TE>{
	double mu;

	ConstantUpdateH<TE>(double c): mu(mu0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Hz(), f.Bz(), mu);
	};
};


// specialization for TM
template<>
struct ConstantUpdateH<TM>{
	double mu;

	ConstantUpdateH<TM>(double c): mu(mu0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Hx(), f.Bx(), mu);
		ConstantUpdate(f.Hy(), f.By(), mu);
	};
};


// specialization for TEM
// template<>
// struct ConstantUpdateH<TEM>{
// 	double mu;

// 	ConstantUpdateH<TEM>(double c): mu(mu0*c) {};

// 	template<class YeeCell>
// 	void operator()(YeeCell && f){
// 		ConstantUpdate(f.Hy(), f.By(), mu);
// 	};
// };

// specialization for TEM
template <bool forward>
struct ConstantUpdateH<TEM, forward>{
	double mu;

	ConstantUpdateH<TEM, forward>(double c): mu(mu0*c) {};

	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	operator()(YeeCell && f){
		ConstantUpdate(f.Hy(), f.By(), mu);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	operator() (YeeCell && f){ConstantCalculate(f.By(), f.Hy(), mu);};


	template<class YeeCell>
	void calculate(YeeCell && f){
		ConstantCalculate(f.By(), f.Hy(), mu);
	};
};







//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <class Mode>
struct ConductiveUpdateE{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template<>
struct ConductiveUpdateE<ThreeD>{
	double eps_r;
	double w0;
	double dt;

	double factor;

	ConductiveUpdateE<ThreeD>(double c, double omega0, double deltat): eps_r(c), w0(omega0), dt(deltat), factor(deltat*omega0/(c*eps0)) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Px() = (f.Px() + factor*f.Dx())/(1.0+factor);
		f.Py() = (f.Py() + factor*f.Dy())/(1.0+factor);
		f.Pz() = (f.Pz() + factor*f.Dz())/(1.0+factor);

		f.Ex() = (f.Dx() - f.Px())/(eps0*eps_r);
		f.Ey() = (f.Dy() - f.Py())/(eps0*eps_r);
		f.Ez() = (f.Dz() - f.Pz())/(eps0*eps_r);
	};

};

// specialization for TE
template<>
struct ConductiveUpdateE<TE>{
	double eps_r;
	double w0;
	double dt;

	double factor;

	ConductiveUpdateE<TE>(double c, double omega0, double deltat): eps_r(c), w0(omega0), dt(deltat), factor(deltat*omega0/(c*eps0)) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Px() = (f.Px() + factor*f.Dx())/(1.0+factor);
		f.Py() = (f.Py() + factor*f.Dy())/(1.0+factor);

		f.Ex() = (f.Dx() - f.Px())/(eps0*eps_r);
		f.Ey() = (f.Dy() - f.Py())/(eps0*eps_r);
	};

};


// specialization for TM
template<>
struct ConductiveUpdateE<TM>{
	double eps_r;
	double w0;
	double dt;

	double factor;

	ConductiveUpdateE<TM>(double c, double omega0, double deltat): eps_r(c), w0(omega0), dt(deltat), factor(deltat*omega0/(c*eps0)) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Pz() = (f.Pz() + factor*f.Dz())/(1.0+factor);

		f.Ez() = (f.Dz() - f.Pz())/(eps0*eps_r);
	};

};

// specialization for TEM
template<>
struct ConductiveUpdateE<TEM>{
	double eps_r;
	double w0;
	double dt;

	double factor;

	ConductiveUpdateE<TEM>(double c, double omega0, double deltat): eps_r(c), w0(omega0), dt(deltat), factor(deltat*omega0/(c*eps0)) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		f.Pz() = (f.Pz() + factor*f.Dz())/(1.0+factor);

		f.Ez() = (f.Dz() - f.Pz())/(eps0*eps_r);
	};

};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <class Mode, 
			class StaticValue,
			class Delta, 
			class LorentzFreq,
			class Gamma>
struct LorentzUpdateParametrized{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template<class StaticValue,
			class Delta, 
			class LorentzFreq,
			class Gamma>
struct LorentzUpdateParametrized<ThreeD, StaticValue, Delta, LorentzFreq, Gamma>{
	double eps_r;
	double dt;



	LorentzUpdateParametrized<ThreeD, StaticValue, Delta, LorentzFreq, Gamma>(double c, double deltat): eps_r(c), dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double b = 2.0*Gamma::get(f)*dt;
		double c = dt*LorentzFreq::get(f)*LorentzFreq::get(f)*(1.0+Delta::get(f)/eps_r);

		f.Px() = (f.Px() + 1.0/(b+dt*c)*(b*f.Px() + dt*f.Jx()));
		f.Jx() = (f.Jx() + 1.0/(b+dt*c)*(-c*f.Px() + f.Jx()));

		f.Py() = (f.Py() + 1.0/(b+dt*c)*(b*f.Py() + dt*f.Jy()));
		f.Jy() = (f.Jy() + 1.0/(b+dt*c)*(-c*f.Py() + f.Jy()));

		f.Pz() = (f.Pz() + 1.0/(b+dt*c)*(b*f.Pz() + dt*f.Jz()));
		f.Jz() = (f.Jz() + 1.0/(b+dt*c)*(-c*f.Pz() + f.Jz()));

		f.Ex() = (f.Dx() - f.Px())/(eps0*eps_r);
		f.Ey() = (f.Dy() - f.Py())/(eps0*eps_r);
		f.Ez() = (f.Dz() - f.Pz())/(eps0*eps_r);
	};

};


// specialization for TM
template<class StaticValue,
			class Delta, 
			class LorentzFreq,
			class Gamma>
struct LorentzUpdateParametrized<TM, StaticValue, Delta, LorentzFreq, Gamma>{
	double eps_r;
	double dt;



	LorentzUpdateParametrized<TM, StaticValue, Delta, LorentzFreq, Gamma>(double c, double deltat): eps_r(c), dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double b = 2.0*Gamma::get(f)*dt;
		double c = dt*LorentzFreq::get(f)*LorentzFreq::get(f)*(1.0+Delta::get(f)/eps_r);


		f.Pz() = (f.Pz() + 1.0/(b+dt*c)*(b*f.Pz() + dt*f.Jz()));
		f.Jz() = (f.Jz() + 1.0/(b+dt*c)*(-c*f.Pz() + f.Jz()));

		f.Ez() = (f.Dz() - f.Pz())/(eps0*eps_r);
	};

};


// specialization for TE
template<class StaticValue,
			class Delta, 
			class LorentzFreq,
			class Gamma>
struct LorentzUpdateParametrized<TE, StaticValue, Delta, LorentzFreq, Gamma>{
	double eps_r;
	double dt;



	LorentzUpdateParametrized<TE, StaticValue, Delta, LorentzFreq, Gamma>(double c, double deltat): eps_r(c), dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double b = 2.0*Gamma::get(f)*dt;
		double c = dt*LorentzFreq::get(f)*LorentzFreq::get(f)*(1.0+Delta::get(f)/eps_r);

		f.Px() = (f.Px() + 1.0/(b+dt*c)*(b*f.Px() + dt*f.Jx()));
		f.Jx() = (f.Jx() + 1.0/(b+dt*c)*(-c*f.Px() + f.Jx()));

		f.Py() = (f.Py() + 1.0/(b+dt*c)*(b*f.Py() + dt*f.Jy()));
		f.Jy() = (f.Jy() + 1.0/(b+dt*c)*(-c*f.Py() + f.Jy()));

		f.Ex() = (f.Dx() - f.Px())/(eps0*eps_r);
		f.Ey() = (f.Dy() - f.Py())/(eps0*eps_r);
	};

};


// specialization for 1D TEM
template<class StaticValue,
			class Delta, 
			class LorentzFreq,
			class Gamma>
struct LorentzUpdateParametrized<TEM, StaticValue, Delta, LorentzFreq, Gamma>{
	double eps_r;
	double dt;



	LorentzUpdateParametrized<TEM, StaticValue, Delta, LorentzFreq, Gamma>(double c, double deltat): eps_r(c), dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double b = 2.0*Gamma::get(f)*dt;
		double c = dt*LorentzFreq::get(f)*LorentzFreq::get(f)*(1.0+Delta::get(f)/eps_r);

		f.Pz() = (f.Pz() + 1.0/(b+dt*c)*(b*f.Pz() + dt*f.Jz()));
		f.Jz() = (f.Jz() + 1.0/(b+dt*c)*(-c*f.Pz() + f.Jz()));

		f.Ez() = (f.Dz() - f.Pz())/(eps0*eps_r);
	};

};





//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************




template <class Mode, 
			class StaticValue,
			class DrudeFreq,
			class Gamma>
struct DrudeUpdateParametrized{
	static_assert(std::is_same<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
};


// specialization for 3D
template<class StaticValue,
			class DrudeFreq,
			class Gamma>
struct DrudeUpdateParametrized<ThreeD, StaticValue, DrudeFreq, Gamma>{
	double dt;



	DrudeUpdateParametrized<ThreeD, StaticValue, DrudeFreq, Gamma>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double b = Gamma::get(f)*dt;
		double c = dt*DrudeFreq::get(f)*DrudeFreq::get(f)*(1.0/StaticValue::get(f));

		auto Pxhold = f.Px();
		f.Px() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Px()) + dt*(f.Jx()+c*f.Dx()));
		f.Jx() = 1.0/(1.0+b+dt*c)*(-c*(Pxhold) + (f.Jx()+c*f.Dx()));
		auto Pyhold = f.Py();
		f.Py() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Py()) + dt*(f.Jy()+c*f.Dy()));
		f.Jy() = 1.0/(1.0+b+dt*c)*(-c*(Pyhold) + (f.Jy()+c*f.Dy()));
		auto Pzhold = f.Pz();
		f.Pz() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Pz()) + dt*(f.Jz()+c*f.Dz()));
		f.Jz() = 1.0/(1.0+b+dt*c)*(-c*(Pzhold) + (f.Jz()+c*f.Dz()));

		f.Ex() = (f.Dx() - f.Px())/(eps0*StaticValue::get(f));
		f.Ey() = (f.Dy() - f.Py())/(eps0*StaticValue::get(f));
		f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));

	};

};


// specialization for TM
template<class StaticValue,
			class DrudeFreq,
			class Gamma>
struct DrudeUpdateParametrized<TM, StaticValue, DrudeFreq, Gamma>{
	double dt;



	DrudeUpdateParametrized<TM, StaticValue, DrudeFreq, Gamma>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double b = Gamma::get(f)*dt;
		double c = dt*DrudeFreq::get(f)*DrudeFreq::get(f)*(1.0/StaticValue::get(f));


		double Pzhold = f.Pz();
		f.Pz() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Pz()) + dt*(f.Jz()+c*f.Dz()));
		f.Jz() = 1.0/(1.0+b+dt*c)*(-c*(Pzhold) + (f.Jz()+c*f.Dz()));

		f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));
	};

};


// specialization for TE
template<class StaticValue,
			class DrudeFreq,
			class Gamma>
struct DrudeUpdateParametrized<TE, StaticValue, DrudeFreq, Gamma>{
	double dt;



	DrudeUpdateParametrized<TE, StaticValue, DrudeFreq, Gamma>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		double b = Gamma::get(f)*dt;
		double c = dt*DrudeFreq::get(f)*DrudeFreq::get(f)*(1.0/StaticValue::get(f));

		double Pxhold = f.Px();
		f.Px() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Px()) + dt*(f.Jx()+c*f.Dx()));
		f.Jx() = 1.0/(1.0+b+dt*c)*(-c*(Pxhold) + (f.Jx()+c*f.Dx()));
		double Pyhold = f.Py();
		f.Py() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Py()) + dt*(f.Jy()+c*f.Dy()));
		f.Jy() = 1.0/(1.0+b+dt*c)*(-c*(Pyhold) + (f.Jy()+c*f.Dy()));

		f.Ex() = (f.Dx() - f.Px())/(eps0*StaticValue::get(f));
		f.Ey() = (f.Dy() - f.Py())/(eps0*StaticValue::get(f));
	};

};


// specialization for 1D TEM
template<class StaticValue,
			class DrudeFreq,
			class Gamma>
struct DrudeUpdateParametrized<TEM, StaticValue, DrudeFreq, Gamma>{
	double dt;



	DrudeUpdateParametrized<TEM, StaticValue, DrudeFreq, Gamma>(double deltat): dt(deltat) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
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