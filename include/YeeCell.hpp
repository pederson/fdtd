#ifndef _YEECELL_H
#define _YEECELL_H

// #include "FDTDConstants.hpp"
// #include "YeeFields.hpp"
// #include "YeeUpdates.hpp"

// #include <limits>
// #include <cmath>
#include <array>
// #include <memory>
// #include <iostream>

namespace fdtd{




template<class EMFields, std::size_t dim>
class StoredNeighbors{
public:
	static_assert(dim <= 3, "Neighbors can only have up to 3 dimensions");

	// run-time functions
	EMFields & getNeighborMax(std::size_t d){return *mNeighbHi[d];};
	EMFields & getNeighborMin(std::size_t d){return *mNeighbLo[d];};
	void setNeighborMax(std::size_t d, const EMFields & f){mNeighbHi[d]=&f;};
	void setNeighborMin(std::size_t d, const EMFields & f){mNeighbLo[d]=&f;};

	// compile-time functions
	template<std::size_t d>
	const EMFields & getNeighborMax(){return *std::get<d>(mNeighbHi);};
	template<std::size_t d>
	const EMFields & getNeighborMin(){return *std::get<d>(mNeighbLo);};
	template<std::size_t d>
	void setNeighborMax(const EMFields & f){return std::get<d>(mNeighbHi) = &f;};
	template<std::size_t d>
	void setNeighborMin(const EMFields & f){return std::get<d>(mNeighbLo) = &f;};

	// convenience accessors
	const EMFields & xMin() const {return std::get<0>(mNeighbLo);};
	const EMFields & xMax() const {return std::get<0>(mNeighbHi);};
	const EMFields & yMin() const {return std::get<1>(mNeighbLo);};
	const EMFields & yMax() const {return std::get<1>(mNeighbHi);};
	const EMFields & zMin() const {return std::get<2>(mNeighbLo);};
	const EMFields & zMax() const {return std::get<2>(mNeighbHi);};

protected:
	std::array<const EMFields *, dim>			mNeighbHi;
	std::array<const EMFields *, dim>			mNeighbLo;
};






template<class C> 
struct Neighbor1Typedef{typedef StoredNeighbors<C,1> type;};
template<class C>
using Neighbor1 = typename Neighbor1Typedef<C>::type;

template<class C> 
struct Neighbor2Typedef{typedef StoredNeighbors<C,2> type;};
template<class C>
using Neighbor2 = typename Neighbor2Typedef<C>::type;

template<class C> 
struct Neighbor3Typedef{typedef StoredNeighbors<C,3> type;};
template<class C>
using Neighbor3 = typename Neighbor3Typedef<C>::type;





/** @class YeeCell
 *	@brief Base YeeCell structure to allow implementation of different fields and other stuff
 *
 *	@tparam 	FieldPolicy 	defines the fields Ex(), By(), etc...
 *	@tparam 	PolarizationPolicy 	defines the polarizations Px(), PXOld(), Py(), etc...
 *	@tparam 	MagnetizationPolicy 	defines the magnetization Mx(), Myold(), etc...
 *	@tparam 	NeighborPolicy<FieldPolicy> 	defines the neighbor fields via getNeighborMax(d) and getNeighborMin(d)
 *	@tparam 	SourcePolicy 		optionally defines source data
 *
 *
 *	@details YeeCell storage class for storing E,D,B,H fields and optional
 *			polarization/magnetization state. In addition, more templates
 *			can be added for Sources, PML, TFSF, etc... but make sure to 
 *			add a default void. Fields are exposed by each of the policies
 *			via accessors (e.g. Dx(), Hy(), Pz(), My()...). Neighbors are 
 *			available via getNeighborMin(std::size_t dim) and getNeighborMax(std::size_t dim)
 *			member functions, but this can be customized through the NeighborPolicy class
 */
template <class FieldPolicy,
		  class PolarizationPolicy,
		  class MagnetizationPolicy,
		  template <class> class NeighborPolicy,
		  class SourcePolicy = void>
class YeeCell : public FieldPolicy,
				public PolarizationPolicy,
				public MagnetizationPolicy,
				public NeighborPolicy<FieldPolicy>,
				public SourcePolicy 
{
public:
	typedef FieldPolicy		 				FieldT;
	typedef PolarizationPolicy 				PolarizationT;
	typedef MagnetizationPolicy 			MagnetizationT;
	typedef NeighborPolicy<FieldPolicy>		NeighborT;
	typedef SourcePolicy 					SourceT;

};


}// end namespace fdtd

#endif