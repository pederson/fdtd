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
template <bool forward>
struct ConstantUpdateE<TE, forward>{
	double eps;

	ConstantUpdateE<TE, forward>(double c): eps(eps0*c) {};

	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	operator()(YeeCell && f){
		ConstantUpdate(f.Ex(), f.Dx(), eps);
		ConstantUpdate(f.Ey(), f.Dy(), eps);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	operator()(YeeCell && f){
		ConstantCalculate(f.Dx(), f.Ex(), eps);
		ConstantCalculate(f.Dy(), f.Ey(), eps);
	};

	template <class YeeCell>
	void calculate(YeeCell && f){
		ConstantCalculate(f.Dx(), f.Ex(), eps);
		ConstantCalculate(f.Dy(), f.Ey(), eps);
	};
};


// specialization for TM
template <bool forward>
struct ConstantUpdateE<TM, forward>{
	double eps;

	ConstantUpdateE<TM, forward>(double c): eps(eps0*c) {};

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


	template <class YeeCell>
	void calculate(YeeCell && f){
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
template <bool forward>
struct ConstantUpdateH<TE, forward>{
	double mu;

	ConstantUpdateH<TE, forward>(double c): mu(mu0*c) {};

	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	operator()(YeeCell && f){
		ConstantUpdate(f.Hz(), f.Bz(), mu);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	operator() (YeeCell && f){
		ConstantCalculate(f.Hz(), f.Bz(), mu);
	};


	template<class YeeCell>
	void calculate(YeeCell && f){
		ConstantCalculate(f.Hz(), f.Bz(), mu);
	};
};


// specialization for TM
template <bool forward>
struct ConstantUpdateH<TM, forward>{
	double mu;

	ConstantUpdateH<TM, forward>(double c): mu(mu0*c) {};

	// use SFINAE to enable this only when forward = false;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<!T, void>::type
	operator()(YeeCell && f){
		ConstantUpdate(f.Hx(), f.Bx(), mu);
		ConstantUpdate(f.Hy(), f.By(), mu);
	};

	// use SFINAE to enable this only when forward = true;
	template <class YeeCell, bool T = forward>
	typename std::enable_if<T, void>::type
	operator() (YeeCell && f){
		ConstantCalculate(f.Bx(), f.Hx(), mu);
		ConstantCalculate(f.By(), f.Hy(), mu);
	};


	template<class YeeCell>
	void calculate(YeeCell && f){
		ConstantCalculate(f.Bx(), f.Hx(), mu);
		ConstantCalculate(f.By(), f.Hy(), mu);
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
	void operator()(YeeCell & f){

		// proceeds in two steps
		// 1) magnetic updated (implicit)
		// 2) regular drude update (implicit)

		// // 1) MAGNETIC UPDATE
		// auto Jxh = drude::q_e/drude::m_e*(-MagneticField::get(f).Bz()*f.Jy() + MagneticField::get(f).By()*f.Jz());
		// auto Jyh = drude::q_e/drude::m_e*( MagneticField::get(f).Bz()*f.Jx() - MagneticField::get(f).Bx()*f.Jz());
		// auto Jzh = drude::q_e/drude::m_e*(-MagneticField::get(f).By()*f.Jx() + MagneticField::get(f).Bx()*f.Jy());

		// f.Jx() = Jxh;
		// f.Jy() = Jyh;
		// f.Jz() = Jzh;

		// // 2) DRUDE UPDATE
		// double b = Gamma::get(f)*dt;
		// double c = dt*DrudeFreq::get(f)*DrudeFreq::get(f)*(1.0/StaticValue::get(f));

		// auto Pxhold = f.Px();
		// f.Px() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Px()) + dt*(f.Jx()+c*f.Dx()));
		// f.Jx() = 1.0/(1.0+b+dt*c)*(-c*(Pxhold) + (f.Jx()+c*f.Dx()));
		// auto Pyhold = f.Py();
		// f.Py() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Py()) + dt*(f.Jy()+c*f.Dy()));
		// f.Jy() = 1.0/(1.0+b+dt*c)*(-c*(Pyhold) + (f.Jy()+c*f.Dy()));
		// auto Pzhold = f.Pz();
		// f.Pz() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Pz()) + dt*(f.Jz()+c*f.Dz()));
		// f.Jz() = 1.0/(1.0+b+dt*c)*(-c*(Pzhold) + (f.Jz()+c*f.Dz()));

		// f.Ex() = (f.Dx() - f.Px())/(eps0*StaticValue::get(f));
		// f.Ey() = (f.Dy() - f.Py())/(eps0*StaticValue::get(f));
		// f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));


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
	void operator()(YeeCell & f){
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
	void operator()(YeeCell & f){


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