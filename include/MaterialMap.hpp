#ifndef _MATERIALMAP_H
#define _MATERIALMAP_H

#include <map>

#include "FDTDConstants.hpp"
#include "DispersiveMaterials.hpp"

namespace fdtd{

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
	: mElectric(ConstantUpdateE<Mode>(1.0))
	, mMagnetic(ConstantUpdateH<Mode>(1.0)) {};

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


}// end namespace fdtd

#endif