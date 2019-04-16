#ifndef _MATERIALMAP_H
#define _MATERIALMAP_H

#include <map>
#include <functional>

#include "FDTDConstants.hpp"
#include "DispersiveMaterials.hpp"

namespace fdtd{



//*************************************************************************
// Usage example:
//			// construct
// 			material_map_type 		matmap;
// 			matmap["vacuum"] = make_material_pair<Mode, node_type>(ConstantUpdateE<Mode>(1.0), ConstantUpdateH<Mode>(1.0));
// 			matmap["dielectric"] = make_material_pair<Mode, node_type>(ConstantUpdateE<Mode>(epsilon), ConstantUpdateH<Mode>(mu));
//			matmap["conductor"] = make_material_pair<Mode, node_type>(ConductiveUpdateE<Mode>(epsilon, coll_freq, delta_t), ConstantUpdateH<Mode>(mu));
//			matmap["pec"] = make_material_pair<Mode, node_type>(ConductiveUpdateE<Mode>(1.0, coll_freq_pec, delta_t), ConstantUpdateH<Mode>(mu));
//
//			// use
//			CellType cell;
//			matmap["dielectric"].electric()(cell);
//**************************************************************************

// CellType is kept as a generic type so that it can accept
// both lvalue and rvalue references
template <typename Mode, typename CellType>
struct MaterialPair{
private:
	typedef std::function<void(std::add_rvalue_reference_t<CellType>)> function_type;
	function_type 		mElectric;
	function_type 		mMagnetic;

public:
	MaterialPair() 
	: mElectric(ConstantUpdate<Mode>(1.0))
	, mMagnetic(ConstantUpdate<Mode>(1.0)) {};

	MaterialPair(function_type e, function_type m)
	: mElectric(e), mMagnetic(m) {};

	function_type & electric() {return mElectric;};
	const function_type & electric() const {return mElectric;};

	function_type & magnetic() {return mMagnetic;};
	const function_type & magnetic() const {return mMagnetic;};
};


template <typename Mode, typename CellType>
MaterialPair<Mode, CellType> make_material_pair(std::function<void(std::add_rvalue_reference_t<CellType>)> e, std::function<void(std::add_rvalue_reference_t<CellType>)> m){
	return MaterialPair<Mode, CellType>(e,m);
};

template <typename Identifier, typename Mode, typename CellType>
using MaterialMap = std::map<Identifier, MaterialPair<Mode, CellType>>;





// remove the Mode template parameter, it's unnecessary. 
// Also allow for the possibility of more arguments to the function (e.g. delta_t, delta_x ...)
// template <typename CellType>
// struct MaterialPair{
// private:
// 	typedef std::function<void(std::add_rvalue_reference_t<CellType>)> function_type;
// 	function_type 		mElectric;
// 	function_type 		mMagnetic;

// public:
// 	MaterialPair() 
// 	: mElectric(ConstantUpdateE<fdtd::ThreeD>(1.0))
// 	, mMagnetic(ConstantUpdateH<fdtd::ThreeD>(1.0)) {};

// 	MaterialPair(function_type e, function_type m)
// 	: mElectric(e), mMagnetic(m) {};

// 	function_type & electric() {return mElectric;};
// 	const function_type & electric() const {return mElectric;};

// 	function_type & magnetic() {return mMagnetic;};
// 	const function_type & magnetic() const {return mMagnetic;};
// };


// template <typename CellType>
// MaterialPair<CellType> make_material_pair(std::function<void(std::add_rvalue_reference_t<CellType>)> e, std::function<void(std::add_rvalue_reference_t<CellType>)> m){
// 	return MaterialPair<CellType>(e,m);
// };

// template <typename Identifier, typename CellType>
// using MaterialMap = std::map<Identifier, MaterialPair<CellType>>;









}// end namespace fdtd

#endif