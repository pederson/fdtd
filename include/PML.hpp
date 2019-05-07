#ifndef _PML_H
#define _PML_H

#include <array>
#include <algorithm>
#include <sstream>

#include "FDTDConstants.hpp"
#include "YeeUpdates.hpp"

namespace fdtd{




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

template <typename Mode, typename scalar_type = double>
struct NewPML{
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
private:
	std::array<std::array<std::array<scalar_type, 3>, 3>,2>	mI;

	std::array<std::array<double, 3>,2>	mK;
	std::array<std::array<double, 3>,2>	mS;
	std::array<std::array<double, 3>,2>	mA;

	std::array<std::array<double, 3>,2>	mB;
	std::array<std::array<double, 3>,2>	mC;
	std::array<std::array<double, 3>,2>	mF;
	std::array<std::array<double, 3>,2>	mG;

public:

	constexpr NewPML() {
		std::array<double, 3> z = {0.0, 0.0, 0.0};
		std::array<double, 3> o = {1.0, 1.0, 1.0};

		std::array<scalar_type, 3> zI = {0.0, 0.0, 0.0};
		std::array<std::array<scalar_type, 3>, 3> zzI = {zI, zI, zI};
		mI.fill(zzI);

		mK.fill(o);
		mS.fill(z);
		mA.fill(z);
		mB.fill(z);
		mC.fill(z);
		mF.fill(o);
		mG.fill(z);
	}

	// integrators
	template <FieldType ft, Dir I, Dir J>
	constexpr scalar_type & pmlI(){return mI[static_cast<char>(ft)][static_cast<char>(I)][static_cast<char>(J)];};

	// constants
	template <FieldType ft, Dir I>
	constexpr double & pmlK(){return mK[static_cast<char>(ft)][static_cast<char>(I)];};

	template <FieldType ft, Dir I>
	constexpr double & pmlS(){return mS[static_cast<char>(ft)][static_cast<char>(I)];};

	template <FieldType ft, Dir I>
	constexpr double & pmlA(){return mA[static_cast<char>(ft)][static_cast<char>(I)];};

	// derived constants
	template <FieldType ft, Dir I>
	constexpr double & pmlB(){return mB[static_cast<char>(ft)][static_cast<char>(I)];};

	template <FieldType ft, Dir I>
	constexpr double & pmlC(){return mC[static_cast<char>(ft)][static_cast<char>(I)];};

	template <FieldType ft, Dir I>
	constexpr double & pmlF(){return mF[static_cast<char>(ft)][static_cast<char>(I)];};

	template <FieldType ft, Dir I>
	constexpr double & pmlG(){return mG[static_cast<char>(ft)][static_cast<char>(I)];};

	template <FieldType ft, Dir I>
	void setPMLParameters(double K, double S, double A, double dt){
		pmlK<ft, I>() = K;
		pmlS<ft, I>() = S;
		pmlA<ft, I>() = A;

		double nu = 1.0/eps0*(S/K + A);
		pmlB<ft, I>() = exp(-dt*nu);
		pmlC<ft, I>() = (nu==0) ? 0.0 : S/K*1.0/(S+K*A)*(pmlB<ft, I>()-1.0);
		
		// convolutional PML
		// pmlF<ft, I>() = 1.0/K ;
		// pmlG<ft, I>() = 1.0;

		// matrix exponential time-stepping PML
		double u = -S/(eps0*K*K);
		pmlF<ft, I>() = (nu==0) ? 1.0/K : 1.0/K - 1.0/dt*(1.0-pmlB<ft, I>() - nu*dt)*u/(nu*nu);
		pmlG<ft, I>() = (nu==0) ? 0.0	: 1.0/dt*(1.0-pmlB<ft, I>())/nu;
	}
};


template <typename Mode, typename scalar_type = double>
struct EmptyPML{
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");

	// integrators
	template <FieldType ft, Dir I, Dir J>
	constexpr scalar_type pmlI(){return 0.0;};

	// constants
	template <FieldType ft, Dir I>
	constexpr scalar_type pmlK(){return 1.0;};

	template <FieldType ft, Dir I>
	constexpr scalar_type pmlS(){return 0.0;};

	template <FieldType ft, Dir I>
	constexpr scalar_type pmlA(){return 0.0;};

	// derived constants
	template <FieldType ft, Dir I>
	constexpr scalar_type pmlB(){return 1.0;};

	template <FieldType ft, Dir I>
	constexpr scalar_type pmlC(){return 0.0;};

	template <FieldType ft, Dir I>
	constexpr scalar_type pmlF(){return 1.0;};

	template <FieldType ft, Dir I>
	constexpr scalar_type pmlG(){return 0.0;};
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



struct PMLOptions{
private:
	bool mConds[3][2];
	std::size_t mNPML;

public:

	PMLOptions() : mNPML(10) {
		for (int i=0; i<3; i++){
			for (int j=0; j<1; j++) mConds[i][j] = false;
		}
	}

	std::size_t & N() {return mNPML;};
	const std::size_t & N() const {return mNPML;};

	template <Dir d, Orientation o>
	bool & get() {return mConds[static_cast<char>(d)][static_cast<char>(o)];};
	template <Dir d, Orientation o>
	const bool & get() const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};


	bool & get(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};
	const bool & get(Dir d, Orientation o) const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};


	bool & operator()(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};	
	const bool & operator()(Dir d, Orientation o) const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};	


	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<PML>" << std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<N>" << mNPML << "</N>" << std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<X>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MIN>" << get<Dir::X, Orientation::MIN>() << " </MIN>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MAX>" << get<Dir::X, Orientation::MAX>() << " </MAX>" <<  std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "</X>" << std::endl ;
			
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Y>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MIN>" << get<Dir::Y, Orientation::MIN>() << " </MIN>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MAX>" << get<Dir::Y, Orientation::MAX>() << " </MAX>" <<  std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "</Y>" << std::endl ;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Z>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MIN>" << get<Dir::Z, Orientation::MIN>() << " </MIN>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MAX>" << get<Dir::Z, Orientation::MAX>() << " </MAX>" <<  std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "</Z>" << std::endl ;
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</PML>" << std::endl;
	}



	#ifdef TINYXML2_INCLUDED
	static PMLOptions readXML(std::string filename) {
		PMLOptions po;

		tinyxml2::XMLDocument doc;
		doc.LoadFile(filename.c_str());

		tinyxml2::XMLNode * n = doc.FirstChild();
		auto c = (n->FirstChild());
				

		while (c != nullptr){
			std::stringstream ss;

			if (!strcmp(c->Value(), "N")){
				ss << c->FirstChild()->Value();
				ss >> po.N();
			}
			else if(!strcmp(c->Value(), "X")){
				tinyxml2::XMLNode * mm = c->FirstChild();
				while (mm != nullptr){
					if (!strcmp(mm->Value(), "MIN")){
						ss << mm->FirstChild()->Value();
						ss >> po.get<Dir::X, Orientation::MIN>();
					}
					else if (!strcmp(mm->Value(), "MAX")){
						ss << mm->FirstChild()->Value();
						ss >> po.get<Dir::X, Orientation::MAX>();
					}
					mm = mm->NextSibling();
				}
			}
			else if(!strcmp(c->Value(), "Y")){
				tinyxml2::XMLNode * mm = c->FirstChild();
				while (mm != nullptr){
					if (!strcmp(mm->Value(), "MIN")){
						ss << mm->FirstChild()->Value();
						ss >> po.get<Dir::Y, Orientation::MIN>();
					}
					else if (!strcmp(mm->Value(), "MAX")){
						ss << mm->FirstChild()->Value();
						ss >> po.get<Dir::Y, Orientation::MAX>();
					}
					mm = mm->NextSibling();
				}
			}
			else if(!strcmp(c->Value(), "Z")){
				tinyxml2::XMLNode * mm = c->FirstChild();
				while (mm != nullptr){
					if (!strcmp(mm->Value(), "MIN")){
						ss << mm->FirstChild()->Value();
						ss >> po.get<Dir::Z, Orientation::MIN>();
					}
					else if (!strcmp(mm->Value(), "MAX")){
						ss << mm->FirstChild()->Value();
						ss >> po.get<Dir::Z, Orientation::MAX>();
					}
					mm = mm->NextSibling();
				}
			}

			c = (c->NextSibling());
		}

		return po;
	}
	#endif
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <typename Mode, FieldType ft, Dir d,
		  template <typename> class FieldPolicy = GetField,
		  template <typename,Dir> class DifferencePolicy = DefaultDifferenceOperator>
struct UpdatePML{
private:
	static_assert(std::is_base_of<EMMode, Mode>::value, "YeeUpdate needs a valid Mode");
	// static_assert(!std::is_same<d, Dir::NONE>::value, "PML Update needs a valid direction");

	double mDt, mDx;

	template <typename EMField>
	struct atomic_update{
		static constexpr Dir 			I = FieldDir<EMField>::value;
		static constexpr Dir 			J = d;
		static constexpr Dir 			K = MutuallyOrthogonal<I,J>::value;

		static constexpr FieldType CoeffField = (EMField::off[static_cast<char>(J)] < 0.0) ?  FieldType::Electric : FieldType::Magnetic;

		// perpendicular components get updated
		template <typename YeeCell, Dir di = d>
		static std::enable_if_t<di != I, void> 
		get(YeeCell && f, double dt, double dx){
			// modify the flux
			GetField<EMField>::get(f) += yu_details::Coeff<ft>::value*dt*GetPML::integrator<ft, I, J>(f)*GetPML::G<CoeffField, J>(f);
			// update the integrator
			GetPML::integrator<ft, I, J>(f) = 	GetPML::B<CoeffField, J>(f)*GetPML::integrator<ft, I, J>(f) + 
												GetPML::C<CoeffField, J>(f)*LeviCivita<I, J, K>::value/dx*DifferencePolicy<typename FieldComponent<yu_details::CurlType<ft>::value, K>::type, J>::get(f);
		}

		// parallel components do nothing
		template <typename YeeCell, Dir di = d>
		static std::enable_if_t<di == I, void> 
		get(YeeCell && f, double dt, double dx){};
	};


public:

	UpdatePML(double dt, double dx) : mDt(dt), mDx(dx) {};

	template <typename YeeCell>
	void update(YeeCell && f, double delt, double delx){
		// loop through all the fluxes and apply update to each
		Detail::for_each_tuple_type<
				 std::conditional_t<ft == FieldType::Electric, 
				 					typename FieldComponents<Mode>::electric_flux,
				 					typename FieldComponents<Mode>::magnetic_flux>, 
				 					atomic_update>(f, delt, delx);
	}

	template <typename YeeCell>
	void operator()(YeeCell && f, double delt, double delx){
		update(std::forward<YeeCell>(f), delt, delx);
	}

	template <typename YeeCell>
	void operator()(YeeCell && f){
		update(std::forward<YeeCell>(f), mDt, mDx);
	}
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
	: mM(3), mMa(1), msMax(0.8*(mM+1)/(imp0*dx)), mkMax(5), maMax(0.05) {}; //1.0, 0.05) {};

	PMLParameterModel(double m, double ma, double sMax, double kMax, double aMax)
	: mM(m), mMa(ma), msMax(sMax), mkMax(kMax), maMax(aMax) {};
	
	double K(double x) const {
		if (x<0.0) return 1.0;
		return (x > 1.0 ? 1.0 : 1.0+pow(x,mM)*(mkMax-1.0));
	};

	double S(double x) const {
		if (x<0.0) return 0.0;
		return (x > 1.0 ? 0.0 : pow(x,mM)*msMax);
	};

	double A(double x) const {
		if (x<0.0) return 0.0;
		return (x > 1.0 ? 0.0 : pow(1.0-x, mMa)*maMax);
	};

};








}// end namespace fdtd

#endif