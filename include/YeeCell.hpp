#ifndef _YEECELL_H
#define _YEECELL_H

#include <array>

namespace fdtd{


struct NoSource{

};


// PML model has dispersion X_v = K_v*E + S_v/(jW + A_v)
// and convolutional terms I_{w,v} where w and v = {x,y,z}
struct NoPML{

	// PML parameters
	constexpr double pmlKx() const {return 1.0;};
	constexpr double pmlKy() const {return 1.0;};
	constexpr double pmlKz() const {return 1.0;};

	constexpr double pmlSx() const {return 0.0;};
	constexpr double pmlSy() const {return 0.0;};
	constexpr double pmlSz() const {return 0.0;};

	constexpr double pmlAx() const {return 0.0;};
	constexpr double pmlAy() const {return 0.0;};
	constexpr double pmlAz() const {return 0.0;};

	// derived PML parameters
	constexpr double pmlBx() const {return 1.0;};
	constexpr double pmlBy() const {return 1.0;};
	constexpr double pmlBz() const {return 1.0;};

	constexpr double pmlCx() const {return 0.0;};
	constexpr double pmlCy() const {return 0.0;};
	constexpr double pmlCz() const {return 0.0;};


	// convolution terms
	constexpr double pmlIxx() const {return 0.0;};
	constexpr double pmlIxy() const {return 0.0;};
	constexpr double pmlIxz() const {return 0.0;};
	constexpr double pmlIyx() const {return 0.0;};
	constexpr double pmlIyy() const {return 0.0;};
	constexpr double pmlIyz() const {return 0.0;};
	constexpr double pmlIzx() const {return 0.0;};
	constexpr double pmlIzy() const {return 0.0;};
	constexpr double pmlIzz() const {return 0.0;};
};


template<class EMFields, std::size_t dim>
class StoredNeighbors{
public:
	static_assert(dim <= 3, "Neighbors can only have up to 3 dimensions");

	// run-time functions
	const EMFields & getNeighborMax(std::size_t d) const {return *(mNeighbHi[d]);};
	const EMFields & getNeighborMin(std::size_t d) const {return *(mNeighbLo[d]);};
	void setNeighborMax(std::size_t d, const EMFields & f){mNeighbHi[d]=&f;};
	void setNeighborMin(std::size_t d, const EMFields & f){mNeighbLo[d]=&f;};

	// compile-time functions
	template<std::size_t d>
	const EMFields & getNeighborMax() const {return *std::get<d>(mNeighbHi);};
	template<std::size_t d>
	const EMFields & getNeighborMin() const {return *std::get<d>(mNeighbLo);};
	template<std::size_t d>
	void setNeighborMax(const EMFields & f){std::get<d>(mNeighbHi) = &f;};
	template<std::size_t d>
	void setNeighborMin(const EMFields & f){std::get<d>(mNeighbLo) = &f;};

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
		  class PMLPolicy = NoPML,
		  class SourcePolicy = NoSource>
class YeeCell : public FieldPolicy,
				public PolarizationPolicy,
				public MagnetizationPolicy,
				public PMLPolicy,
				public SourcePolicy 
{
public:
	typedef FieldPolicy		 				FieldT;
	typedef PolarizationPolicy 				PolarizationT;
	typedef MagnetizationPolicy 			MagnetizationT;
	// typedef NeighborPolicy<FieldPolicy>		NeighborT;
	typedef PMLPolicy						PMLT;
	typedef SourcePolicy 					SourceT;

};


}// end namespace fdtd

#endif