#ifndef _BOUNDARYMAP_H
#define _BOUNDARYMAP_H

#include <map>
#include <functional>

#include "FDTDConstants.hpp"
#include "BoundaryUpdates.hpp"

namespace fdtd{



//*************************************************************************
// Usage example:
//			// construct
// 			boundary_map_type 		boundmap;
// 			boundmap["xmin"] 	= make_pec_boundary 		<TE, X, MIN>(xmin_begin(), xmin_end());
// 			boundmap["xmax"] 	= make_symmetric_boundary	<TE, X, MAX>(xmax_begin(), xmax_end());
//			boundmap["ymin"] 	= make_periodic_boundary	<TE, Y, MIN>(ymin_begin(), ymin_end(), ymax_functor());
//			boundmap["ymax"] 	= make_periodic_boundary	<TE, Y, MAX>(ymax_begin(), ymax_end(), ymax_functor());
//
//			// use
//			CellType cell;
//			boundmap["xmin"].updateE();
//**************************************************************************

template <typename Identifier>
using BoundaryMap = std::map<Identifier, BoundaryData>;


}// end namespace fdtd

#endif