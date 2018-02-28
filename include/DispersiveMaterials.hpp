#ifndef _DISPERSIVEMATERIALS_H
#define _DISPERSIVEMATERIALS_H

#include "FDTDConstants.hpp"

namespace fdtd{



struct VacuumPolarization{
	constexpr double permittivity_r() const {return 1.0;};
	constexpr double Px() const {return 0;};
	constexpr double Py() const {return 0;};
	constexpr double Pz() const {return 0;};
};


struct VacuumMagnetization{
	constexpr double permeability_r() const {return 1.0;};
	constexpr double Mx() const {return 0;};
	constexpr double My() const {return 0;};
	constexpr double Mz() const {return 0;};
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



struct SinglePolarization{
	double mPx, mPy, mPz, mPermittivityR;

	SinglePolarization()
	: mPermittivityR(1.0)
	, mPx(0.0), mPy(0.0), mPz(0.0) {};

	const double & permittivity_r() const {return mPermittivityR;};
	const double & Px() const {return mPx;};
	const double & Py() const {return mPy;};
	const double & Pz() const {return mPz;};


	double & permittivity_r() {return mPermittivityR;}
	double & Px() {return mPx;};
	double & Py() {return mPy;};
	double & Pz() {return mPz;};
};

struct DoublePolarization : public SinglePolarization{
	double mJx, mJy, mJz;

	DoublePolarization()
	: SinglePolarization(), mJx(0.0), mJy(0.0), mJz(0.0) {};


	const double & Jx() const {return mJx;};
	const double & Jy() const {return mJy;};
	const double & Jz() const {return mJz;};


	double & permittivity_r() {return mPermittivityR;}
	double & Jx() {return mJx;};
	double & Jy() {return mJy;};
	double & Jz() {return mJz;};
};


struct SingleMagnetization{
	double mMx, mMy, mMz, mPermeabilityR;

	SingleMagnetization()
	: mPermeabilityR(1.0)
	, mMx(0.0), mMy(0.0), mMz(0.0) {};

	const double & permeability_r() const {return mPermeabilityR;};
	const double & Mx() const {return mMx;};
	const double & My() const {return mMy;};
	const double & Mz() const {return mMz;};

	double & permeability_r() {return mPermeabilityR;}
	double & Mx() {return mMx;};
	double & My() {return mMy;};
	double & Mz() {return mMz;};
};


struct DoubleMagnetization : public SingleMagnetization{
	double mKx, mKy, mKz;

	DoubleMagnetization()
	: mKx(0.0), mKy(0.0), mKz(0.0) {};

	const double & Kx() const {return mKx;};
	const double & Ky() const {return mKy;};
	const double & Kz() const {return mKz;};

	double & Kx() {return mKx;};
	double & Ky() {return mKy;};
	double & Kz() {return mKz;};
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


template <class Mode>
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
template<>
struct ConstantUpdateE<TEM>{
	double eps;

	ConstantUpdateE<TEM>(double c): eps(eps0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Ez(), f.Dz(), eps);
	};

	template<class YeeCell>
	void calculate(YeeCell & f){
		ConstantCalculate(f.Dz(), f.Ez(), eps);
	};
};









//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************












template <class Mode>
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
template<>
struct ConstantUpdateH<TEM>{
	double mu;

	ConstantUpdateH<TEM>(double c): mu(mu0*c) {};

	template<class YeeCell>
	void operator()(YeeCell & f){
		ConstantUpdate(f.Hy(), f.By(), mu);
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

		double Pxhold = f.Px();
		f.Px() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Px()) + dt*(f.Jx()+c*f.Dx()));
		f.Jx() = 1.0/(1.0+b+dt*c)*(-c*(Pxhold) + (f.Jx()+c*f.Dx()));
		double Pyhold = f.Py();
		f.Py() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Py()) + dt*(f.Jy()+c*f.Dy()));
		f.Jy() = 1.0/(1.0+b+dt*c)*(-c*(Pyhold) + (f.Jy()+c*f.Dy()));
		double Pzhold = f.Pz();
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

		double Pzhold = f.Pz();
		f.Pz() = 1.0/(1.0+b+dt*c)*((1.0+b)*(f.Pz()) + dt*(f.Jz()+c*f.Dz()));
		f.Jz() = 1.0/(1.0+b+dt*c)*(-c*(Pzhold) + (f.Jz()+c*f.Dz()));

		f.Ez() = (f.Dz() - f.Pz())/(eps0*StaticValue::get(f));
	};

};


}// end namespace fdtd

#endif