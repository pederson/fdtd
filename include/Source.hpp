#ifndef _SOURCE_H
#define _SOURCE_H

#include <map>


#include "FDTDConstants.hpp"

namespace fdtd{


namespace Detail{

//***********************
	// an assignable vector requires the non-const begin() and end()
	// iterator accessors
	template<typename T, typename _ = void>
	struct has_method_freq_data : std::false_type {};

	template<typename T>
	struct has_method_freq_data<
	      T,
	      std::conditional_t<
	          false,
	          has_method_helper<
	              decltype(std::declval<T>().freq()),
	              decltype(std::declval<T>().min_freq()),
	              decltype(std::declval<T>().max_freq())
	              >,
	          void
	          >
	      > : public std::true_type {};
} // end namespace Detail





// wraps around the print_summary() functionality so that it can
// be stored for heterogeneous types
template <typename T>
struct FreqDataWrapper{
private:
  const T* mT;
public:
  typedef std::function<void(std::ostream &, unsigned int)> type;

  FreqDataWrapper(const T & t) : mT(&t) {
    static_assert(Detail::has_method_freq_data<T>::value, "Type must have frequency data with methods freq(), min_freq(), max_freq()");
  };
  FreqDataWrapper(const FreqDataWrapper & p) : mT(p.mT) {
    static_assert(Detail::has_method_freq_data<T>::value, "Type must have frequency data with methods freq(), min_freq(), max_freq()");
  };

  double freq() {return mT->freq();};
  double min_freq() {return mT->min_freq();};
  double max_freq() {return mT->max_freq();};

  double operator()(){return mT->freq();};
};




//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
namespace Detail{
	enum class Function : int{
		sin,
		ricker
	};
}
template <> struct NameArray<Detail::Function>{
static constexpr std::array<const char *, 2> value = {"Sinusoid",
													  "Ricker"};};
constexpr std::array<const char *, 2> NameArray<Detail::Function>::value;


//************************************************************
//************************************************************

// signal function must be in the list of possible ones
template <typename DataPolicy>
struct SignalFunction : public DataPolicy {
public:
	template <typename ... Args>
	constexpr SignalFunction(Args... args) : DataPolicy(args...) {};
	
	double operator()(double t){return DataPolicy::get(t, *this);};
};


//************************************************************
//************************************************************

struct SinData{
private:
	double mFreq, mPhase, mDelay;
public:
	SinData(double Freq, double Phase=0) : mFreq(Freq), mPhase(Phase), mDelay(10.0/Freq*fdtd::pi/2.0) {};
	typedef Detail::Function type;
	static constexpr type value = Detail::Function::sin;
	static constexpr double get(double t, SinData & s){return std::tanh(t/s.mDelay)*std::sin(s.mFreq*t+s.mPhase);};

	double min_freq() const {return mFreq;};
	double freq()	  const {return mFreq;};
	double max_freq() const {return mFreq;};

	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<Sinusoid>" << std::endl;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Freq>" << mFreq << "</Freq>" << std::endl;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Phase>" << mPhase << "</Phase>" << std::endl;

		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</Sinusoid>" << std::endl;
	}
};

using Sinusoid = SignalFunction<SinData>;




//************************************************************
//************************************************************

struct RickerData{
private:
	double mFreq, mDelay;
public:
	RickerData(double Freq) : mFreq(Freq), mDelay(2.0/Freq) {};
	typedef Detail::Function type;
	static constexpr type value = Detail::Function::ricker;
	static constexpr double get(double t, RickerData & s){
		double param = fdtd::pi*s.mFreq*(t-s.mDelay);
		return (1.0-2.0*param*param)*exp(-param*param);
	};

	double min_freq() const {return 0.0;};
	double freq()	  const {return mFreq;};
	double max_freq() const {return 3.0*mFreq;};


	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<Ricker>" << std::endl;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Freq>" << mFreq << "</Freq>" << std::endl;

		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</Ricker>" << std::endl;
	}
};

using Ricker = SignalFunction<RickerData>;


//************************************************************
//************************************************************


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************



template <typename SignalType, typename CellType, typename FieldGetter>
struct SoftSource{
private:

	// optional printing behavior
	static constexpr bool isPrintable = Detail::has_method_print_summary<SignalType>::value;

	SignalType 		mFunc;
	CellType * 		mLoc;

	
	SoftSource() {};
public:
	
	SoftSource(SignalType f, CellType & c)
	: mFunc(f), mLoc(&c) {};

	SignalType & signal() {return mFunc;};
	const SignalType & signal() const {return mFunc;};

	double freq() const {return mFunc.freq();};
	double min_freq() const {return mFunc.min_freq();};
	double max_freq() const {return mFunc.max_freq();};

	CellType * location() {return mLoc;};

	void apply(double t){
		FieldGetter::get(*mLoc) += mFunc(t);
	}

	// in order to make this into a functor
	void operator()(double t){return apply(t);};

	double at(double t) const {return mFunc(t);};

	template <bool enable = isPrintable>
	std::enable_if_t<!enable, void>
	print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<SoftSource>" << std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			std::cout << "<Location>" << std::endl;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			std::cout << "</Location>" << std::endl;
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</SoftSource>" << std::endl;
	}

	template <bool enable = isPrintable>
	std::enable_if_t<enable, void>
	print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<SoftSource>" << std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			std::cout << "<Signal>" << std::endl;

				mFunc.print_summary(os, ntabs+2);

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			std::cout << "</Signal>" << std::endl;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			std::cout << "<Location>" << std::endl;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			std::cout << "</Location>" << std::endl;
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</SoftSource>" << std::endl;
	}
};


template <typename F, typename SignalType, typename CellType>
SoftSource<SignalType, CellType, GetField<F>> make_soft_source(SignalType f, CellType & c){
	return SoftSource<SignalType, CellType, GetField<F>>(f, c);
}


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


#include "PML.hpp"
#include "YeeFields.hpp"
#include "YeeCell.hpp"
#include "YeeUpdates.hpp"
#include "DispersiveMaterials.hpp"
#include "DefaultInterfaces.hpp"

// TFSF Source
struct Aux1D{
private:
	typedef TEM							mode;
	typedef double 						scalar_type;
	typedef NewPML<mode, scalar_type> 	pml_type;
	typedef YeeFields<mode, scalar_type, std::array> 	field_type;
	typedef YeeCell<field_type, 
					VacuumPolarization<scalar_type>, 
					VacuumMagnetization<scalar_type>,
					pml_type
					>					yee_type;

	// make a neighbor type;
	template <typename T>
	struct Stencil1D : public T{
		typedef Stencil1D<T>		NeighborType;
		typedef NeighborType		GhostType;

		// ValueT										center;
		NeighborType * 		neighbMin;
		NeighborType * 		neighbMax;

		std::shared_ptr<GhostType> 		ghostMin;
		std::shared_ptr<GhostType>		ghostMax;
	public:
		Stencil1D() : T() {
			neighbMin = nullptr;
			neighbMax = nullptr;
		};

		// Stencil1D(const NeighborType & o) = default;

		Stencil1D(const T & c) : T(c){};

		Stencil1D & operator=(T & c){
			Stencil1D s(c);
			// std::swap(s.neighbMin, neighbMin);
			// std::swap(s.neighbMax, neighbMax);
			std::swap(static_cast<T &>(*this), static_cast<T &>(s));
			return *this;
		}

		// runtime set and get
		constexpr NeighborType & getNeighborMin(std::size_t d) const{
			return *neighbMin;
		};
		constexpr NeighborType & getNeighborMax(std::size_t d) const{
			return *neighbMax;
		};
		void setNeighborMin(NeighborType & val, std::size_t d){
			neighbMin = &val;
		};
		void setNeighborMax(NeighborType & val, std::size_t d){
			neighbMax = &val;
		};


		// runtime add ghost cells
		void addGhostMin(std::size_t d){
			ghostMin = std::make_unique<GhostType>();
			ghostMin->setNeighborMax(*this,d);
			setNeighborMin(*ghostMin,d);
		};
		void addGhostMax(std::size_t d){
			ghostMax = std::make_unique<GhostType>();
			ghostMax->setNeighborMin(*this,d);
			setNeighborMax(*ghostMax,d);
		};
	};


	typedef std::vector<Stencil1D<yee_type>> 		cell_container_type;
	typedef std::function<double(double)>			signal_type;

	cell_container_type 		mCells;
	static constexpr std::size_t 		mNPML = 15;
	double 						mT; // current time
	double 						mDx;
	signal_type 				mSignal;

public:
	template <typename Signal>
	Aux1D(std::size_t ncells, double dx, double dt, Signal sig) : mDx(dx), mT(0) {
		mCells.resize(ncells + 2*mNPML + 20);

		// connect neighbors
		for (std::size_t i=1; i<mCells.size()-1; i++){
			mCells[i].setNeighborMin(mCells[i-1], 0);
			mCells[i].setNeighborMax(mCells[i+1], 0);
		}
		mCells[0].addGhostMin(0);
		mCells[0].setNeighborMax(mCells[1], 0);
		mCells[mCells.size()-1].addGhostMax(0);
		mCells[mCells.size()-1].setNeighborMin(mCells[mCells.size()-2], 0);

		// setup PML
		PMLParameterModel p(dx);
		for (auto i=1; i<mNPML+2; i++){
			double x = static_cast<double>(mNPML-i+1)/static_cast<double>(mNPML);
			double xm = x - 0.5/static_cast<double>(mNPML);
			mCells[i].setPMLParameters<FieldType::Electric, Dir::X>(p.K(x), p.S(x), p.A(x), dt);
			mCells[i].setPMLParameters<FieldType::Magnetic, Dir::X>(p.K(xm), p.S(xm), p.A(xm), dt);

			mCells[mCells.size()-i].setPMLParameters<FieldType::Electric, Dir::X>(p.K(x), p.S(x), p.A(x), dt);
			mCells[mCells.size()-1-i].setPMLParameters<FieldType::Magnetic, Dir::X>(p.K(xm), p.S(xm), p.A(xm), dt);
		}

		// setup source
		mSignal = static_cast<signal_type>(sig);
	}


	// take a time-step
	void time_step(double dt){
		mT += dt;
		// std::cout << "time-stepping aux sim D: " << std::endl;

		// update D/E
		std::for_each(mCells.begin(), mCells.end(), YeeUpdate<mode, Dz>(dt, mDx));
		// std::cout << "time-stepping aux sim D PML: " << std::endl;
		std::for_each(mCells.begin(), mCells.begin()+mNPML+2, UpdatePML<mode, FieldType::Electric, Dir::X>(dt, mDx));
		std::for_each(mCells.end()-(mNPML+3), mCells.end(), UpdatePML<mode, FieldType::Electric, Dir::X>(dt, mDx));
		// std::cout << "time-stepping aux sim D materials: " << std::endl;
		std::for_each(mCells.begin(), mCells.end(), VacuumUpdate<mode, FieldType::Electric>());

		// std::cout << "time-stepping aux sim D signal: " << std::endl;
		// update the sources
		mCells[mNPML + 5].Ez() += mSignal(mT);

		// std::cout << "time-stepping aux sim B: " << std::endl;

		// update B/H
		std::for_each(mCells.begin(), mCells.end(), YeeUpdate<mode, By>(dt, mDx));
		std::for_each(mCells.begin(), mCells.begin()+mNPML+2, UpdatePML<mode, FieldType::Magnetic, Dir::X>(dt, mDx));
		std::for_each(mCells.end()-(mNPML+3), mCells.end(), UpdatePML<mode, FieldType::Magnetic, Dir::X>(dt, mDx));
		std::for_each(mCells.begin(), mCells.end(), VacuumUpdate<mode, FieldType::Magnetic>());

	}


	// interpolation of field values 
	// interpolate the E value to location "loc"
	// where loc is measured from the center of the 1D
	// domain, towards the direction of the wavevector
	double interpolateE(double loc){
		double nd = static_cast<double>(mCells.size())/2.0 + loc/mDx;
		std::size_t n = nd; // truncate to integer

		double frac = (nd - n)/mDx;

		return (1.0-frac)*mCells[n].Ez() + frac*mCells[n+1].Ez();
	}

	double interpolateH(double loc){
		double nd = static_cast<double>(mCells.size())/2.0 + loc/mDx - 0.5;
		std::size_t n = nd; // truncate to integer

		double frac = (nd - n)/mDx;

		return (1.0-frac)*mCells[n].Hy() + frac*mCells[n+1].Hy();
	}



};


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

template <typename Mode,
		  typename LocationFunctor,
		  typename IteratorType,
		  class TimePolicy 			= BD1, 
		  class PMLCoeffPolicy 		= GetPML>
struct TFSFPlane{
private:	
	static constexpr std::size_t 			ndim = Mode::dim;
	typedef std::array<double, ndim> 		PointType;

	PointType 								mWaveVector, mNormal;
	std::shared_ptr<Aux1D> 					mAuxSim;
	std::array<double, 3>					mPolarization, mPolarizationH;
	LocationFunctor 						mLocFunctor;
	IteratorType							mBegin, mEnd;

	template <typename EMField>
	struct field_update{
	private:
		static constexpr Dir 			I = FieldDir<EMField>::value;
		static constexpr FieldType 		FT = EMField::field_type;
		typedef GetField<EMField> 		FP;

		template <Dir d>
		struct dir_update{
			static constexpr Dir 			J = d;
			static constexpr Dir 			K = MutuallyOrthogonal<I,J>::value;
			static constexpr FieldType 		CurlType =  yu_details::CurlType<FT>::value;
			typedef typename FieldComponent<CurlType, K>::type 	CurlField;
			typedef std::conditional_t<CurlType == FieldType::Electric,
									   typename FieldComponents<Mode>::electric,
									   typename FieldComponents<Mode>::magnetic> 	tuple_type;

			template <typename YeeCell>
			static constexpr std::enable_if_t<Detail::tuple_contains_type<CurlField, tuple_type>::value, std::remove_reference_t<decltype(FP::get(std::declval<YeeCell>()))>>
			get(YeeCell && f, double dt, double dx, PointType & normal, std::array<double, 3> & pol, double srcval){
				// std::cout << "Levi-civita: " << LeviCivita<I, J, K>::value << std::endl;
				// std::cout << "PML: " << GetPML::F<FT,J>(f) << std::endl;
				// std::cout << "Normal: " << normal[static_cast<int>(J)] << std::endl;
				// std::cout << "srcval: " << srcval << std::endl;
				// std::cout << "polarization: " << pol[static_cast<int>(K)] << std::endl;
				// std::cout << "return value: " <<  LeviCivita<I, J, K>::value*GetPML::F<FT, J>(f)/(dx)* normal[static_cast<int>(J)] * srcval * pol[static_cast<int>(K)] << std::endl;
				return LeviCivita<I, J, K>::value*GetPML::F<FT, J>(f)/(dx)* normal[static_cast<int>(J)] * srcval * pol[static_cast<int>(K)];
			}

			template <typename YeeCell>
			static constexpr std::enable_if_t<!Detail::tuple_contains_type<CurlField, tuple_type>::value, std::remove_reference_t<decltype(FP::get(std::declval<YeeCell>()))>>
			get(YeeCell && f, double dt, double dx, PointType & normal, std::array<double, 3> & pol, double srcval){
				return 0.0;
			}
		};
	public:
		template <typename YeeCell>
		static void get(YeeCell && f, double dt, double dx, PointType & normal, std::array<double, 3> & pol, double srcval){
			static constexpr Dir d1 = (I == Dir::X ? Dir::Y : (I == Dir::Y ? Dir::Z : Dir::X));
			static constexpr Dir d2 = MutuallyOrthogonal<I, d1>::value;

			GetField<EMField>::get(f) +=  yu_details::Coeff<FT>::value*TimePolicy::curl_coeff*dt*
										  (dir_update<d1>::get(f, dt, dx, normal, pol, srcval) + 
										  	dir_update<d2>::get(f, dt, dx, normal, pol, srcval));
		}
	};
public:
	TFSFPlane(PointType wavevec, PointType normal, LocationFunctor lf, 
			  IteratorType beg, IteratorType end, std::shared_ptr<Aux1D> auxsim)
	: mWaveVector(wavevec), mNormal(normal), mLocFunctor(lf), mBegin(beg), mEnd(end), mAuxSim(auxsim)
	{
		// make sure normal is unit size
		double N = 0;
		for (auto d=0; d<ndim; d++) N += mNormal[d]*mNormal[d];
		N = std::sqrt(N);
		for (auto d=0; d<ndim; d++) mNormal[d] /= N;

		// normalize the wavevector direction
		double K = 0;
		for (auto d=0; d<ndim; d++) K += mWaveVector[d]*mWaveVector[d];
		K = std::sqrt(K);
		for (auto d=0; d<ndim; d++) mWaveVector[d] /= K;

		// set default polarization
		if (std::is_same<Mode, TE>::value){
			mPolarizationH = {0,0,1};
			mPolarization = {mWaveVector[1], -mWaveVector[0], 0};
		}
		else if (std::is_same<Mode, TM>::value){
			mPolarization = {0,0,1};
			mPolarizationH = {mWaveVector[1], -mWaveVector[0], 0};
		}
		else if (std::is_same<Mode, TEM>::value){
			mPolarization = {0,0,1};
			mPolarizationH = {0,1,0};
		}

		// std::cout << "Normal: " << mNormal[0] << ", " << mNormal[1] << std::endl;
		// std::cout << "WaveVector: " << mWaveVector[0] << ", " << mWaveVector[1] << std::endl;

		// std::cout << "Polarization: " << mPolarization[0] << ", " << mPolarization[1] << ", " << mPolarization[2] << std::endl;
		// std::cout << "PolarizationH: " << mPolarizationH[0] << ", " << mPolarizationH[1] << ", " << mPolarizationH[2] << std::endl;
		
	}


	void updateD(double dt, double dx){
		std::cout << "update D: " << std::endl;
		for (auto it=mBegin; it!=mEnd; it++){
			// get the location of the current iterator
			PointType pt = mLocFunctor(it);

			// std::cout << "point location: " << pt[0] << ", " << pt[1] << std::endl;

			// get the projection of this vector on the wavevector
			double proj = 0;
			for (auto d=0; d<ndim; d++) proj += mWaveVector[d]*pt[d];

			// std::cout << "projection: " << proj << std::endl;

			// get the source value at the iterator location
			double src = mAuxSim->interpolateH(proj);

			// update for each component that needs updating
			nested_for_each_tuple_type<field_update, typename FieldComponents<Mode>::electric_flux>(*it, dt, dx, mNormal, mPolarizationH, src);
		}

		// throw -1;
		mAuxSim->time_step(dt);
	}

	void updateB(double dt, double dx){
		for (auto it=mBegin; it!=mEnd; it++){
			// get the location of the current iterator
			PointType pt = mLocFunctor(it);

			// get the projection of this vector on the wavevector
			double proj = 0;
			for (auto d=0; d<ndim; d++) proj += mWaveVector[d]*pt[d];

			// get the source value at the iterator location
			double src = mAuxSim->interpolateE(proj);

			// update for each component that needs updating
			nested_for_each_tuple_type<field_update, typename FieldComponents<Mode>::magnetic_flux>(*it, dt, dx, mNormal, mPolarization, src);
		}
	}

};


template <typename Mode, Dir dd, Orientation o, typename LocFunctor, typename Iterator>
TFSFPlane<Mode, LocFunctor, Iterator> make_tfsf_plane(Iterator && begin, Iterator && end, 
													  LocFunctor lf, std::array<double, Mode::dim> K,
													  std::shared_ptr<Aux1D> sim){
	// construct the normal
	std::array<double, Mode::dim> normal;
	for (auto d=0; d<Mode::dim; d++) normal[d] = 0;
	normal[static_cast<int>(dd)] = (o == Orientation::MAX ? 1 : -1);

	return TFSFPlane<Mode, LocFunctor, Iterator>(K, normal, lf, begin, end, sim);
}


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


// A source of this type has :
//		- a defined functor void(double) which retrieves the
// 			value of the source, given the time
//
//		- a print_summary() function
//	
//		- frequency data in the form of a freq() function
struct Source{
private:
	typedef std::function<void(double)> 						ftype;
	typedef std::function<void(std::ostream &, unsigned int)> 	ptype;
	typedef std::function<double(void)> 						freqtype;

	ftype 		mSrc;
	ptype		mPrinter;
	freqtype 	mFreqData;
public:
	template <typename T>
	Source(T && t) : mSrc(t), mPrinter(PrintSummaryWrapper<T>(t)), mFreqData(FreqDataWrapper<T>(t)) {};

	void operator()(double t){return mSrc(t);};
	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const {return mPrinter(os, ntabs);};
	double freq() const {return mFreqData();};
};

template <typename T>
Source make_source(T && t){return Source(std::forward<T>(t));};





//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

template <typename Identifier>
struct SourceMap : public std::map<Identifier, Source>
{
	typedef std::map<Identifier, Source> MapType;
	using MapType::begin;
	using MapType::end;
	using MapType::cbegin;
	using MapType::cend;

	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<SourceMap>" << std::endl;
			for (auto it=cbegin(); it!=cend(); it++){

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				std::cout << "<Item>" << std::endl;


					for (auto i=0; i<ntabs+2; i++) os << "\t" ;
					std::cout << "<Identifier>" << it->first << "</Identifier>" << std::endl;

					it->second.print_summary(os, ntabs+2);

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				std::cout << "</Item>" << std::endl;
			}
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</SourceMap>" << std::endl;
	}



	#ifdef TINYXML2_INCLUDED

		
	// public:

	// 	template <typename Mode, bool forward = false>
	// 	static SourceMap readXML(tinyxml2::XMLNode * n, double dt){
	// 		SourceMap mp;

	// 		std::string id;
	// 		tinyxml2::XMLNode * mnode;

	// 		auto c = (n->FirstChild());

	// 		while (c != nullptr){
	// 			std::stringstream ss;

	// 			if(!strcmp(c->Value(), "Item")){
	// 				tinyxml2::XMLNode * mm = c->FirstChild();
					
	// 				while (mm != nullptr){
	// 					if(!strcmp(mm->Value(), "Identifier")){
	// 						id = mm->FirstChild()->Value();
	// 					}
	// 					else if(!strcmp(mm->Value(), "Source")){
	// 						mnode = mm;
	// 					}

	// 					mm = (mm->NextSibling());
	// 				}
	// 				// mp.insert(std::make_pair(id, MaterialPair<CellType>::template readXML<Mode, forward>(mnode, dt)));
	// 			}

	// 			c = (c->NextSibling());
	// 		}

	// 		return mp;
	// 	}

	// 	template <typename Mode, bool forward = false>
	// 	static SourceMap readXML(std::string filename, double dt) {
	// 		tinyxml2::XMLDocument doc;
	// 		doc.LoadFile(filename.c_str());

	// 		tinyxml2::XMLNode * n = doc.FirstChild();
	// 		return readXML<Mode, forward>(n, dt);
	// 	}

	#endif
};



}// end namespace fdtd

#endif