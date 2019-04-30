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