#ifndef _DISPERSIVEMATERIALS_H
#define _DISPERSIVEMATERIALS_H

#include "FDTDConstants.hpp"
#include "DefaultInterfaces.hpp"

namespace fdtd{


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************		


// restrict the set of possible
// constant coefficient dispersions
enum class Dispersion : char {
	Vacuum,
	Constant,
	Conductive,
	Drude,
	Lorentz,
	MagnetizedDrude,
	MultiCoefficient,
	Anisotropic
};
template <> struct NameArray<Dispersion>{
static constexpr std::array<const char *, 8> value = {"Vacuum", 		// translational symmetry
												   "Constant", 	// translational symmetry
												   "Conductive",				// reflection symmetry
												   "Drude",				// reflection symmetry
												   "Lorentz",			// reflection symmetry
												   "MagnetizedDrude",		// reflection symmetry
												   "MultiCoefficient", 			
												   "Anisotropic"};};
constexpr std::array<const char *, 8> NameArray<Dispersion>::value;



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************		


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
	double mPermittivityR;
	const double & permittivity_r() const {return mPermittivityR;};
	double & permittivity_r() {return mPermittivityR;}

	FDTD_DECLARE_MEMBER(scalar_type, Px);
	FDTD_DECLARE_MEMBER(scalar_type, Py);
	FDTD_DECLARE_MEMBER(scalar_type, Pz);


	SinglePolarization()
	: mPermittivityR(1.0)
	, mPx(0.0), mPy(0.0), mPz(0.0) {};
};

template <typename scalar_type = double>
struct DoublePolarization : public SinglePolarization<scalar_type>{

	FDTD_DECLARE_MEMBER(scalar_type, Jx);
	FDTD_DECLARE_MEMBER(scalar_type, Jy);
	FDTD_DECLARE_MEMBER(scalar_type, Jz);

	DoublePolarization()
	: SinglePolarization<scalar_type>(), mJx(0.0), mJy(0.0), mJz(0.0) {};


};

template <typename scalar_type = double>
struct SingleMagnetization{
	double mPermeabilityR;
	const double & permeability_r() const {return mPermeabilityR;};
	double & permeability_r() {return mPermeabilityR;}

	FDTD_DECLARE_MEMBER(scalar_type, Mx);
	FDTD_DECLARE_MEMBER(scalar_type, My);
	FDTD_DECLARE_MEMBER(scalar_type, Mz);

	SingleMagnetization()
	: mPermeabilityR(1.0)
	, mMx(0.0), mMy(0.0), mMz(0.0) {};
};

template <typename scalar_type = double>
struct DoubleMagnetization : public SingleMagnetization<scalar_type>{
	
	FDTD_DECLARE_MEMBER(scalar_type, Kx);
	FDTD_DECLARE_MEMBER(scalar_type, Ky);
	FDTD_DECLARE_MEMBER(scalar_type, Kz);

	DoubleMagnetization()
	: SingleMagnetization<scalar_type>(), mKx(0.0), mKy(0.0), mKz(0.0) {};

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

					

FDTD_GET_FIELD_DEF(Jx);
FDTD_GET_FIELD_DEF(Jy);
FDTD_GET_FIELD_DEF(Jz);
FDTD_GET_FIELD_DEF(Px);
FDTD_GET_FIELD_DEF(Py);
FDTD_GET_FIELD_DEF(Pz);
FDTD_GET_FIELD_DEF(Kx);
FDTD_GET_FIELD_DEF(Ky);
FDTD_GET_FIELD_DEF(Kz);
FDTD_GET_FIELD_DEF(Mx);
FDTD_GET_FIELD_DEF(My);
FDTD_GET_FIELD_DEF(Mz);


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


// this is a generalized version of the single-component material update
template <typename UpdateStruct, bool implicit, typename EMField>
struct DispersiveAtomicUpdate{
	static constexpr Dir d 				= FieldDir<EMField>::value; 
	static constexpr FieldType ft 		= EMField::field_type;
	typedef typename FluxComponent<ft, d>::type 		DType;
	typedef typename FieldComponent<ft, d>::type 		EType;
	typedef typename CurrentComponent<ft, d>::type 		JType;
	typedef typename PolarizationComponent<ft, d>::type PType;

	template <typename T, typename ... Args, bool I = implicit>
	static std::enable_if_t<implicit, void> 
	get(T && f, Args... args){
		UpdateStruct::implicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
										  args...); 
	}

	template <typename T, typename ... Args, bool I = implicit>
	static std::enable_if_t<!implicit, void> 
	get(T && f, Args... args){
		UpdateStruct::explicit_update(GetField<DType>::get(f), GetField<EType>::get(f), GetField<PType>::get(f), GetField<JType>::get(f),
										  args...); 
	}
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// This is very generalized and is designed to delegate 
// the data storage and get(...) functionality to the template params
// this way, the data storage policy can handle any read/write of that
// data if necessary
template <class 		Mode, 
		  FieldType 	ftype,
		  bool 			forward,
		  class 		StoragePolicy,
		  class 		GetterPolicy>
struct DispersiveUpdate : public StoragePolicy, public GetterPolicy
{
private:
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
public:
	DispersiveUpdate() {};

	// any extra constructor arguments are passed to the StoragePolicy constructor
	template <typename... Args>
	DispersiveUpdate(Args... args) : StoragePolicy(args...) {};

	DispersiveUpdate(const DispersiveUpdate & d) : StoragePolicy(d) {};

	// the get(...) function is defined by the GetterPolicy
	// using GetterPolicy

	// pass just a cell and use the stored time-step
	template <class YeeCell>
	void operator()(YeeCell && f){
		GetterPolicy::get(f, *static_cast<StoragePolicy*>(this));
	};

	// allows the user to pass in a variable time-step
	template <class YeeCell>
	void operator()(YeeCell && f, double delta_t){
		GetterPolicy::get(f, *static_cast<StoragePolicy*>(this), delta_t);
	};

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



	// call the constant get(...) function for a generic StoragePolicy
	template <typename Mode, FieldType ftype, bool forward, typename StoragePolicy>
	struct ConstantCall {
	private:
		template <typename EMField>
		using UpdateType = DispersiveAtomicUpdate<constant::ConstantUpdate, !forward, EMField>;
		static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
	public:
		template <typename YeeCell, typename... Args>
		static void get(YeeCell && f, StoragePolicy & sp, Args... args){
			Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							UpdateType>(f, StoragePolicy::get(sp, f, args...), val);
		}
	};



	// class that defines a constant-coefficient storage policy
	// and the interface to get(...) those values in order to 
	// interface with the ConstantCall function
	struct ConstantStorage{
	private:
		double mK;
	public:
		ConstantStorage(double K) : mK(K) {};

		double & K() {return mK;};
		const double & K() const {return mK;};

		template <typename ... Args>
		static decltype(auto) get(ConstantStorage & c, Args... args){return c.K();}
	
		void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "<Constant>" << std::endl;
				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Constant>" <<  mK << "</Constant>" << std::endl;
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "</Constant>" << std::endl;
		}
	};

	

}// end namespace constant




template <typename Mode, FieldType ftype, bool forward = false>
using ConstantUpdate = DispersiveUpdate<Mode, ftype, forward, 
										   constant::ConstantStorage, 
										   constant::ConstantCall<Mode, ftype, forward, constant::ConstantStorage>>;




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


	// call the conductive get(...) function for a generic StoragePolicy
	template <typename Mode, FieldType ftype, bool forward, typename StoragePolicy>
	struct ConductiveCall {
	private:
		template <typename EMField>
		using UpdateType = DispersiveAtomicUpdate<conductive::ConductiveUpdate, !forward, EMField>;
		static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
	public:
		template <typename YeeCell, typename... Args>
		static void get(YeeCell && f, StoragePolicy & sp, Args... args){
			Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							UpdateType>(f, StoragePolicy::getK(sp, f, args...), val, 
										   StoragePolicy::getFreq(sp, f, args...), 
										   StoragePolicy::getDt(sp, f, args...));
		}

		// call with a specified time-step
		template <typename YeeCell, typename... Args>
		static void get(YeeCell && f, StoragePolicy & sp, double dt, Args... args){
			Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							UpdateType>(f, StoragePolicy::getK(sp, f, args...), val, 
										   StoragePolicy::getFreq(sp, f, args...), 
										   dt);
		}
	};



	// class that defines a conductive-coefficient storage policy
	// and the interface to get(...) those values in order to 
	// interface with the ConductiveCall function
	struct ConductiveStorage{
		FDTD_DECLARE_MEMBER(double, K);
		FDTD_DECLARE_MEMBER(double, Freq);
		FDTD_DECLARE_MEMBER(double, dt);
	public:
		ConductiveStorage(double K, double Freq, double delt) : mK(K), mFreq(Freq), mdt(delt) {};
		ConductiveStorage(double K, double Freq) : mK(K), mFreq(Freq) {};
		template <typename ... Args>
		static decltype(auto) getK(ConductiveStorage & c, Args... args){return c.K();}
		template <typename ... Args>
		static decltype(auto) getFreq(ConductiveStorage & c, Args... args){return c.Freq();}
		template <typename ... Args>
		static decltype(auto) getDt(ConductiveStorage & c, Args... args){return c.dt();}

		void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "<Conductive>" << std::endl;
				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Constant>" <<  mK << "</Constant>" << std::endl;

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Freq>" <<  mFreq << "</Freq>" << std::endl;

			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "</Conductive>" << std::endl;
		}
	};
}// end namespace conductive


template <typename Mode, FieldType ftype, bool forward = false>
using ConductiveUpdate = DispersiveUpdate<Mode, ftype, forward, 
										   conductive::ConductiveStorage, 
										   conductive::ConductiveCall<Mode, ftype, forward, conductive::ConductiveStorage>>;




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

	// call the lorentz get(...) function for a generic StoragePolicy
	template <typename Mode, FieldType ftype, bool forward, typename StoragePolicy>
	struct LorentzCall {
	private:
		template <typename EMField>
		using UpdateType = DispersiveAtomicUpdate<lorentz::LorentzUpdate, !forward, EMField>;
		static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
	public:
		template <typename YeeCell, typename... Args>
		static void get(YeeCell && f, StoragePolicy & sp, Args... args){
			Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							UpdateType>(f, StoragePolicy::getK(sp, f, args...),
										   StoragePolicy::getDelta(sp, f, args...), val, 
										   StoragePolicy::getFreq(sp, f, args...),
										   StoragePolicy::getGamma(sp, f, args...), 
										   StoragePolicy::getDt(sp, f, args...));
		}

		// call with a specified time-step
		template <typename YeeCell, typename... Args>
		static void get(YeeCell && f, StoragePolicy & sp, double dt, Args... args){
			Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							UpdateType>(f, StoragePolicy::getK(sp, f, args...),
										   StoragePolicy::getDelta(sp, f, args...), val, 
										   StoragePolicy::getFreq(sp, f, args...),
										   StoragePolicy::getGamma(sp, f, args...), 
										   dt);
		}
	};



	// class that defines a conductive-coefficient storage policy
	// and the interface to get(...) those values in order to 
	// interface with the ConductiveCall function
	struct LorentzStorage{
		FDTD_DECLARE_MEMBER(double, K);
		FDTD_DECLARE_MEMBER(double, Delta);
		FDTD_DECLARE_MEMBER(double, Freq);
		FDTD_DECLARE_MEMBER(double, Gamma);
		FDTD_DECLARE_MEMBER(double, dt);
	public:
		LorentzStorage(double K, double Delta, double Freq, double Gamma, double dt) : mK(K), mDelta(Delta), mFreq(Freq), mGamma(Gamma), mdt(dt) {};
		LorentzStorage(double K, double Delta, double Freq, double Gamma) : mK(K), mDelta(Delta), mFreq(Freq), mGamma(Gamma) {};
		template <typename ... Args>
		static decltype(auto) getK(LorentzStorage & c, Args... args){return c.K();}
		template <typename ... Args>
		static decltype(auto) getDelta(LorentzStorage & c, Args... args){return c.Delta();}
		template <typename ... Args>
		static decltype(auto) getFreq(LorentzStorage & c, Args... args){return c.Freq();}
		template <typename ... Args>
		static decltype(auto) getGamma(LorentzStorage & c, Args... args){return c.Gamma();}
		template <typename ... Args>
		static decltype(auto) getDt(LorentzStorage & c, Args... args){return c.dt();}
	
		void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "<Lorentz>" << std::endl;
				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Constant>" <<  mK << "</Constant>" << std::endl;

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Delta>" <<  mDelta << "</Delta>" << std::endl;

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Freq>" <<  mFreq << "</Freq>" << std::endl;

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Gamma>" <<  mGamma << "</Gamma>" << std::endl;
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "</Lorentz>" << std::endl;
		}
	};

}// end namespace lorentz


template <typename Mode, FieldType ftype, bool forward = false>
using LorentzUpdate = DispersiveUpdate<Mode, ftype, forward, 
										   lorentz::LorentzStorage, 
										   lorentz::LorentzCall<Mode, ftype, forward, lorentz::LorentzStorage>>;





//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


namespace drude{
	static constexpr double m_e = 9.11e-31;
	static constexpr double q_e = 1.6e-19;

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


		//////////////// Matrix Exponential ////////////////////////

		// implementation here

	};


	// call the drude get(...) function for a generic StoragePolicy
	template <typename Mode, FieldType ftype, bool forward, typename StoragePolicy>
	struct DrudeCall {
	private:
		template <typename EMField>
		using UpdateType = DispersiveAtomicUpdate<drude::DrudeUpdate, !forward, EMField>;
		static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
	public:
		template <typename YeeCell, typename... Args>
		static void get(YeeCell && f, StoragePolicy & sp, Args... args){
			Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							UpdateType>(f, StoragePolicy::getK(sp, f, args...), val, 
										   StoragePolicy::getFreq(sp, f, args...),
										   StoragePolicy::getGamma(sp, f, args...), 
										   StoragePolicy::getDt(sp, f, args...));
		}

		// call with a specified time-step
		template <typename YeeCell, typename... Args>
		static void get(YeeCell && f, StoragePolicy & sp, double dt, Args... args){
			Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
													   typename FieldComponents<Mode>::electric, 
													   typename FieldComponents<Mode>::magnetic>, 
							UpdateType>(f, StoragePolicy::getK(sp, f, args...), val, 
										   StoragePolicy::getFreq(sp, f, args...),
										   StoragePolicy::getGamma(sp, f, args...), 
										   dt);
		}
	};



	// class that defines a conductive-coefficient storage policy
	// and the interface to get(...) those values in order to 
	// interface with the ConductiveCall function
	struct DrudeStorage{
		FDTD_DECLARE_MEMBER(double, K);
		FDTD_DECLARE_MEMBER(double, Freq);
		FDTD_DECLARE_MEMBER(double, Gamma);
		FDTD_DECLARE_MEMBER(double, dt);
	public:
		DrudeStorage(double K, double Freq, double Gamma, double dt) : mK(K), mFreq(Freq), mGamma(Gamma), mdt(dt) {};
		DrudeStorage(double K, double Freq, double Gamma) : mK(K), mFreq(Freq), mGamma(Gamma) {};
		template <typename ... Args>
		static decltype(auto) getK(DrudeStorage & c, Args... args){return c.K();}
		template <typename ... Args>
		static decltype(auto) getFreq(DrudeStorage & c, Args... args){return c.Freq();}
		template <typename ... Args>
		static decltype(auto) getGamma(DrudeStorage & c, Args... args){return c.Gamma();}
		template <typename ... Args>
		static decltype(auto) getDt(DrudeStorage & c, Args... args){return c.dt();}
	

		void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "<Drude>" << std::endl;
				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Constant>" <<  mK << "</Constant>" << std::endl;

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Freq>" <<  mFreq << "</Freq>" << std::endl;

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				os << "<Gamma>" <<  mGamma << "</Gamma>" << std::endl;
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "</Drude>" << std::endl;
		}
	};


	// interface with the ConductiveCall function
	struct FluidDrudeStorage{
	private:
		static constexpr double Kval = fdtd::sqrt(drude::q_e*drude::q_e/(drude::m_e*fdtd::eps0));
		FDTD_DECLARE_MEMBER(double, dt);
	public:
		FluidDrudeStorage(double dt) : mdt(dt) {};
		FluidDrudeStorage() {};
		template <typename ... Args>
		static constexpr decltype(auto) getK(FluidDrudeStorage & c, Args... args){return 1.0;}
		template <typename CellType, typename ... Args>
		static decltype(auto) getFreq(FluidDrudeStorage & c, CellType && f, Args... args){return Kval*std::sqrt(f.ne());}
		template <typename CellType, typename ... Args>
		static decltype(auto) getGamma(FluidDrudeStorage & c, CellType && f, Args... args){return f.nu_m();}
		template <typename ... Args>
		static decltype(auto) getDt(FluidDrudeStorage & c, Args... args){return c.dt();}
	

		void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "<FluidDrude>" << std::endl;
			for (auto i=0; i<ntabs; i++) os << "\t" ;
			os << "</FluidDrude>" << std::endl;
		}
	};

}// end namespace drude


template <typename Mode, FieldType ftype, bool forward = false>
using DrudeUpdate = DispersiveUpdate<Mode, ftype, forward, 
										   drude::DrudeStorage, 
										   drude::DrudeCall<Mode, ftype, forward, drude::DrudeStorage>>;


template <typename Mode, FieldType ftype, bool forward = false>
using FluidDrudeUpdate = DispersiveUpdate<Mode, ftype, forward, 
										   drude::FluidDrudeStorage, 
										   drude::DrudeCall<Mode, ftype, forward, drude::FluidDrudeStorage>>;



// // This is very generalized and pretty good, but it would be convenient to have a type
// // that could store the material parameters within the struct itself
// template <class Mode,
// 		  FieldType ftype = FieldType::Electric,
// 		  bool forward = false, 
// 		  class StaticValue = StoredValue,
// 		  class DrudeFreq 	= Stored<void>,
// 		  class Gamma 		= Stored<int>
// 		  >
// struct DrudeUpdate : public StaticValue, public DrudeFreq, public Gamma
// {
// private:
// 	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
// 	static constexpr double val = (ftype == FieldType::Electric ? fdtd::eps0 : fdtd::mu0);
// 	double dt;

// 	template <typename EMField>
// 	using UpdateType = DispersiveAtomicUpdate<drude::DrudeUpdate, !forward, EMField>;

// public:
// 	DrudeUpdate(double deltat): dt(deltat) {};
// 	DrudeUpdate(double K, double w0, double gamma, double deltat): StaticValue(K), DrudeFreq(w0), Gamma(gamma), dt(deltat) {};

// 	template <class YeeCell>
// 	void get(YeeCell && f){
// 		Detail::for_each_tuple_type<std::conditional_t<ftype == FieldType::Electric, 
// 													   typename FieldComponents<Mode>::electric, 
// 													   typename FieldComponents<Mode>::magnetic>, 
// 							UpdateType>(f, StaticValue::get(f, *this), val, DrudeFreq::get(f, *this), Gamma::get(f, *this), dt);
// 	};

// 	// allows the user to pass in a variable time-step
// 	template <class YeeCell>
// 	void operator()(YeeCell && f){
// 		get(f);
// 	};

// 	// allows the user to pass in a variable time-step
// 	template <class YeeCell>
// 	void operator()(YeeCell && f, double delta_t){
// 		dt = delta_t;
// 		get(f);
// 	};

// };



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



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