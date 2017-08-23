#ifndef _YEECELL_H
#define _YEECELL_H

#include <array>

namespace fdtd{


struct NoSource{};
struct NoSurfaceQuantities{};
struct NoExtra{};


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
		  class SourcePolicy = NoSource,
		  class SurfaceQuantityPolicy = NoSurfaceQuantities,
		  class ExtraFields = NoExtra>
class YeeCell : public FieldPolicy,
				public PolarizationPolicy,
				public MagnetizationPolicy,
				public PMLPolicy,
				public SourcePolicy,
				public SurfaceQuantityPolicy,
				public ExtraFields 
{
public:
	typedef FieldPolicy		 				FieldT;
	typedef PolarizationPolicy 				PolarizationT;
	typedef MagnetizationPolicy 			MagnetizationT;
	typedef PMLPolicy						PMLT;
	typedef SourcePolicy 					SourceT;
	typedef SurfaceQuantityPolicy 			SurfaceQuantityT;
	typedef ExtraFields						ExtraFieldT;
};


}// end namespace fdtd

#endif