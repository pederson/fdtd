#ifndef _DISPERSIVEMATERIALS_H
#define _DISPERSIVEMATERIALS_H

#include "FDTDConstants.hpp"

namespace fdtd{



struct VacuumPolarization{
	constexpr double Px() const {return 0;};
	constexpr double Py() const {return 0;};
	constexpr double Pz() const {return 0;};
};


struct VacuumMagnetization{
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





template <class FieldType>
inline void ConstantUpdate(FieldType & Out, const FieldType & In, double c){Out=In/c;}


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

























// struct Vacuum{};
// struct PEC{};
// struct PMC{};

// struct OwnedConstant{
// 	double val;
// 	constexpr OwnedConstant(): val(1.0) {};
// 	constexpr OwnedConstant(double c): val(c) {};
// 	constexpr double getConstant() const {return val;};
// 	void setConstant(double c) {val=c;};
// };

// template <class ParameterPolicy = OwnedConstant>
// struct ConstantUpdate : public ParameterPolicy{
// 	template <typename FieldType>
// 	void update(FieldType & Out, const FieldType & In, double dt){
// 		Out = In/getConstant();
// 	};
// };

// // template <typename FieldType,
// // 		  class ParameterPolicy,
// // 		  template <class> class UpdatePolicy>
// // struct ElementalDispersion : public UpdatePolicy<FieldType>, public ParameterPolicy{

// // };



// struct NullUpdate{
// 	template <typename FieldType>
// 	void update(FieldType & Out, const FieldType & In, double dt){};
// };



// template <class Pxx, class Pxy, class Pxz,
// 				  , class Pyy, class Pyz,
// 				  			 , class Pzz>
// struct GeneralLinearDispersion : public Pxx, public Pxy, public Pxz,
// 											   public Pyy, public Pyz,
// 											 			   public Pzz
// {

// };


// template <class Pxx, class Pyy, class Pzz> 
// struct GeneralDiagonalDispersion GeneralLinearDispersion<Pxx, void, void,
// 														  Pyy		, void,
// 																	, Pzz>;


// template <class Pxx> 
// struct GeneralIsotropicDispersion GeneralDiagonalDispersion<Pxx, Pxx, Pxx>;


// template <class ElementalDispersion>
// struct IsotropicDipsersion{
// 	ElementalDispersion			mDisp;

// 	template <typename FieldType>
// 	void update(double dt, double dx, std::pair<FieldType, FieldType>... fields){

// 	};
// };


// template <class Mode>
// struct VacuumMaterialModel{
// 	static Mode m;


// 	// implementation for 3D
// 	template <typename M>
// 	void updateE(double dt, double dx, EMFields & f, typename std::enable_if<std::is_same<M, ThreeD>::value>::type=m){
// 		f.Ex() = f.Dx()/eps0;
// 		f.Ey() = f.Dy()/eps0;
// 		f.Ez() = f.Dz()/eps0;
// 	};


// 	// implementation for 3D
// 	template <typename M>
// 	void updateH(double dt, double dx, EMFields & f, typename std::enable_if<std::is_same<M, ThreeD>::value>::type=m){
// 		f.Hx() = f.Bx()/mu0;
// 		f.Hy() = f.By()/mu0;
// 		f.Hz() = f.Bz()/mu0;
// 	}
// };

// template<class EMFields>
// class LinearMaterialUpdater{
// public:
// 	void updateE(EMFields & f, double dt, double dx){

// 	};

// protected:
// 	LinearDispersion		mDispE;
// 	LinearDispersion		mDispH;
// };


// class VacuumDispersion{};

// class PECDispersion{};
// class PMCDispersion{};


// template<> void YeeCellTM<VacuumDispersion, VacuumDispersion>::update_E(double dt){
// 	fields.E[2] = fields.D[2]/eps0;
// }
// template<> void YeeCellTM<VacuumDispersion, VacuumDispersion>::update_H(double dt){
// 	fields.H[0] = fields.B[0]/mu0;
// 	fields.H[1] = fields.B[1]/mu0;
// }

// template<> void YeeCellTM<PECDispersion, VacuumDispersion>::update_E(double dt){
// 	// fields.E[0] = 0.0;
// 	// fields.E[1] = 0.0;
// 	fields.E[2] = 0.0;
// }

// template<> void YeeCellTM<PECDispersion, VacuumDispersion>::update_H(double dt){
// 	fields.H[0] = fields.B[0];
// 	fields.H[1] = fields.B[1];
// }


// class ConstantDispersion{
// public:
// 	double val;

// 	ConstantDispersion() 
// 	: val(1.0) {};

// 	ConstantDispersion(double v)
// 	: val(v) {};

// 	void update(std::array<double, 3> & E, std::array<double, 3> & D, double c, double dt){
// 		for (auto i=0; i<3; i++) E[i] = D[i]/(c*val);
// 	}
// };



}// end namespace fdtd

#endif