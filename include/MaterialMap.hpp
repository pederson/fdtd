#ifndef _MATERIALMAP_H
#define _MATERIALMAP_H

#include <map>
#include <functional>

#include "FDTDConstants.hpp"
#include "DispersiveMaterials.hpp"

namespace fdtd{


// wraps around the print_summary() functionality so that it can
// be stored for heterogeneous types
template <typename T>
struct PrintSummaryWrapper{
private:
	T mT;
public:
	PrintSummaryWrapper(const T & t) : mT(t) {};
	PrintSummaryWrapper(const PrintSummaryWrapper & p) : mT(p.mT) {};

	template <typename StreamType>
	void operator()(StreamType & os, unsigned int ntabs=0){mT.print_summary(os, ntabs);};
};

// template <typename Derived>
// struct AddPrinter {
// public:
// 	// template <typename StreamType>
// 	void print_summary(std::ostream & os, unsigned int ntabs){

// 	}
// };

//*************************************************************************
// Usage example:
//			// construct
// 			material_map_type 		matmap;
// 			matmap["vacuum"] 		= make_material_pair<node_type>(ConstantUpdate<Mode>(1.0)							, ConstantUpdate<Mode>(1.0));
// 			matmap["dielectric"] 	= make_material_pair<node_type>(ConstantUpdate<Mode>(epsilon)						, ConstantUpdate<Mode>(mu));
//			matmap["conductor"] 	= make_material_pair<node_type>(ConductiveUpdate<Mode>(epsilon, coll_freq, delta_t)	, ConstantUpdate<Mode>(mu));
//			matmap["pec"] 			= make_material_pair<node_type>(ConductiveUpdate<Mode>(1.0, coll_freq_pec, delta_t)	, ConstantUpdate<Mode>(mu));
//
//			// use
//			CellType cell;
//			matmap["dielectric"].electric()(cell);
//**************************************************************************

// CellType is kept as a generic type so that it can accept
// both lvalue and rvalue references
template <typename CellType>
struct MaterialPair{
private:
	typedef std::function<void(std::add_rvalue_reference_t<CellType>)> function_type;
	function_type 		mElectric;
	function_type 		mMagnetic;

	typedef std::function<void(std::ostream &, unsigned int)> printer_type;
	printer_type 		mElectricPrinter;
	printer_type 		mMagneticPrinter;

public:
	MaterialPair(function_type e, function_type m)
	: mElectric(e), mMagnetic(m) {};

	MaterialPair(function_type e, function_type m, printer_type pe, printer_type pm)
	: mElectric(e), mMagnetic(m), mElectricPrinter(pe), mMagneticPrinter(pm) {};

	function_type & electric() {return mElectric;};
	const function_type & electric() const {return mElectric;};

	function_type & magnetic() {return mMagnetic;};
	const function_type & magnetic() const {return mMagnetic;};


	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<Electric>" << std::endl;
			mElectricPrinter(os, ntabs+1);
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</Electric>" << std::endl;

		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<Magnetic>" << std::endl;
			mMagneticPrinter(os, ntabs+1);
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</Magnetic>" << std::endl;
	}
};




// template <typename CellType>
// MaterialPair<CellType> make_material_pair(std::function<void(std::add_rvalue_reference_t<CellType>)> e, std::function<void(std::add_rvalue_reference_t<CellType>)> m){
// 	return MaterialPair<CellType>(e,m);
// };

template <typename CellType, typename Electric, typename Magnetic>
MaterialPair<CellType> make_material_pair(Electric && e, Magnetic && m){

	// static_assert(std::is_same<Electric, void>::value, "Debug");
	typedef std::function<void(std::add_rvalue_reference_t<CellType>)> ftype;
	typedef std::function<void(std::ostream &, unsigned int)> ptype;
	return MaterialPair<CellType>(static_cast<ftype>(e), static_cast<ftype>(m), 
								  static_cast<ptype>(PrintSummaryWrapper<Electric>(e)), static_cast<ptype>(PrintSummaryWrapper<Magnetic>(m)));
};

template <typename Identifier, typename CellType>
using MaterialMap = std::map<Identifier, MaterialPair<CellType>>;








// // define a print_summary() function
// template <typename Identifier, typename CellType>
// void MaterialMap<Identifier, CellType>::print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const {
	
// 	// for (auto i=0; i<ntabs; i++) os << "\t" ;
// 	// os << "<MaterialMap>" << std::endl;
// 	// 	for (auto it=begin(); it!=end(); it++){
// 	// 		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
// 	// 		std::cout << "<Identifier>" << it->first << "</Identifier>" << std::endl;
// 	// 		it->second.print_summary(os, ntabs+2);
// 	// 	}
// 	// for (auto i=0; i<ntabs; i++) os << "\t" ;
// 	// os << "</MaterialMap>" << std::endl;
// }

}// end namespace fdtd

#endif