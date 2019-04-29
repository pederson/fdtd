#ifndef _SOURCE_H
#define _SOURCE_H

#include <map>


#include "FDTDConstants.hpp"

namespace fdtd{

//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
namespace Detail{
	enum class Function : int{
		sin,
		tanh
	};
}
template <> struct NameArray<Detail::Function>{
static constexpr std::array<const char *, 2> value = {"Sinusoid",
													  "Tanh"};};
constexpr std::array<const char *, 2> NameArray<Detail::Function>::value;


//************************************************************
//************************************************************

// composable function must be in the list of possible ones
template <typename Functor, typename Argument>
struct ComposableFunction{
private:
	Argument mArg;

	struct PrintArg{
		template <typename T = Argument>
		static std::enable_if_t<Detail::has_method_print_summary<T>::value && !std::is_void<T>::value, void>
		get(std::ostream & os, unsigned int ntabs, Argument * arg){
			arg->print_summary(os, ntabs);
		}

		template <typename T = Argument>
		static std::enable_if_t<!Detail::has_method_print_summary<T>::value || std::is_void<T>::value, void>
		get(std::ostream & os, unsigned int ntabs, Argument * arg){};
	};

	// overload this for every template specialization
	static double functor(double t){return Functor::get(t);};
public:
	constexpr ComposableFunction(Argument arg) : mArg(arg) {};
	
	template <typename T = Argument>
	std::enable_if_t<!std::is_void<T>::value, double>
	operator()(double t){return functor(mArg(t));};

	template <typename T = Argument>
	std::enable_if_t<std::is_void<T>::value, double>
	operator()(double t){return functor(t);};

	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<" << NameArray<typename Functor::type>::value[static_cast<int>(Functor::value)] << ">" << std::endl;
			PrintArg::get(os, ntabs+1, &mArg);
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<" << NameArray<typename Functor::type>::value[static_cast<int>(Functor::value)] << ">" << std::endl;

	}
};

//************************************************************
//************************************************************

struct SinFunctor{
	typedef Detail::Function type;
	static constexpr type value = Detail::Function::sin;
	static constexpr double get(double t){return std::sin(t);};
};

template <typename Argument>
using Sinusoid = ComposableFunction<SinFunctor, Argument>;

//************************************************************
//************************************************************

struct TanhFunctor{
	typedef Detail::Function type;
	static constexpr type value = Detail::Function::sin;
	static constexpr double get(double t){return std::tanh(t);};
};

template <typename Argument>
using Tanh = ComposableFunction<TanhFunctor, Argument>;


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









//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

template <typename Identifier>
struct SourceMap : public std::map<Identifier, std::function<void(double)>>
{
	typedef std::map<Identifier, std::function<void(double)>> MapType;
	using MapType::begin;
	using MapType::end;

	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<SourceMap>" << std::endl;
			for (auto it=begin(); it!=end(); it++){

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