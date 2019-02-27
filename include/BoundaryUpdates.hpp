#ifndef _BOUNDARYUPDATES_H
#define _BOUNDARYUPDATES_H

#include "FDTDConstants.hpp"
#include "YeeUpdates.hpp"
#include "DefaultInterfaces.hpp"

namespace fdtd{


enum class Boundary : char {
	NONE,
	Periodic,
	BlochPeriodic,
	PEC,
	PMC,
	Symmetric,
	Antisymmetric,
	Parallel
};


struct BoundaryOptions{
private:
	Boundary mConds[3][2];
public:

	BoundaryOptions() {
		for (int i=0; i<3; i++){
			for (int j=0; j<1; j++) mConds[i][j] = Boundary::PEC;
		}
	}

	template <Dir d, Orientation o>
	Boundary & get() {return mConds[static_cast<char>(d)][static_cast<char>(o)];};

	Boundary & get(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};

	Boundary & operator()(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};	


	// void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
	// 	for (auto i=0; i<ntabs; i++) os << "\t" ;
	// 	os << "<BoundaryOptions>" << std::endl;
	// 	for (auto i=0; i<ntabs+1; i++) os << "\t" ;
	// 	os << "<Center>" << m_center << "</Center>" << std::endl;
	// 	for (auto i=0; i<ntabs+1; i++) os << "\t" ;
	// 	os << "<Radius>" << m_radius << "</Radius>" << std::endl ;
	// 	for (auto i=0; i<ntabs; i++) os << "\t" ;
	// 	os << "</BoundaryOptions>" << std::endl;
	// }
};



//************************************************************
//************************************************************
//********************* PEC BOUNDARY *************************
//************************************************************
//************************************************************
//************************************************************


template <class Mode, Dir d, Orientation o, class IndexableObject>
struct UpdateBoundaryPECD{
private:
	IndexableObject * mIdxObj;

public:
	UpdateBoundaryPECD(IndexableObject & idx_obj): mIdxObj(&idx_obj) {};

	void operator()(){
		// do nothing
	};

	// static void update(){
	// };
};

template <class Mode, Dir d, Orientation o, class IndexableObject>
struct UpdateBoundaryPECB{
private:
	IndexableObject * mIdxObj;

public:
	UpdateBoundaryPECB(IndexableObject & idx_obj): mIdxObj(&idx_obj) {};

	void operator()(){
		// do nothing
	};

	// static void update(){
	// };
};

//************************************************************
//************************************************************
//********************* PMC BOUNDARY *************************
//************************************************************
//************************************************************
//************************************************************


// template <class Mode, Dir d, Orientation o, class IndexableObject>
// struct UpdateBoundaryPMCB{
// private:
// 	std::array<std::size_t, Mode::dim> mDims;
// 	IndexableObject * mIdxObj;

// 	template <Orientation T = o>
// 	typename std::enable_if<T==Orientation::MIN, std::size_t>::type
// 	fixed_index(){return 0;};

// 	template <Orientation T = o>
// 	typename std::enable_if<T==Orientation::MAX, std::size_t>::type
// 	fixed_index(){return mDims[static_cast<char>(d)]-1;};


// 	template <Orientation T = o>
// 	typename std::enable_if<T==Orientation::MAX, std::size_t>::type
// 	opp_index(){return 0;};

// 	template <Orientation T = o>
// 	typename std::enable_if<T==Orientation::MIN, std::size_t>::type
// 	opp_index(){return mDims[static_cast<char>(d)]-1;};

// 	template <typename CellType, Orientation T = o>
// 	typename std::enable_if<T == Orientation::MAX, CellType>::type &
// 	get_neighbor(CellType && t){return t.getNeighborMax(static_cast<char>(d));};

// 	template <typename CellType, Orientation T = o>
// 	typename std::enable_if<T == Orientation::MIN, CellType>::type &
// 	get_neighbor(CellType && t){return t.getNeighborMin(static_cast<char>(d));};

// 	template <class EMField, std::size_t... Idx>
// 	void apply_impl(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds, std::index_sequence<Idx...>){
// 		GetField<EMField>::get(idxobj(inds[Idx]...)) = 0;
// 	}

// 	template <class M = Mode>
// 	typename std::enable_if<std::is_same<M, ThreeD>::value, void>
// 	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
// 		apply_impl<Hx>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
// 		apply_impl<Hy>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
// 		apply_impl<Hz>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
// 	}

// 	template <class M = Mode>
// 	typename std::enable_if<std::is_same<M, TM>::value, void>
// 	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
// 		apply_impl<Hx>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
// 		apply_impl<Hy>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
// 	}

// 	template <class M = Mode>
// 	typename std::enable_if<std::is_same<M, TE>::value, void>
// 	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
// 		apply_impl<Hz>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
// 	}

// 	template <class M = Mode>
// 	typename std::enable_if<std::is_same<M, TEM>::value, void>
// 	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
// 		apply_impl<Hy>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
// 	}

// public:
// 	UpdateBoundaryPMCB(IndexableObject & idx_obj): mIdxObj(&idx_obj) {};

// 	void operator()(){
// 		// starting inds
// 		std::array<std::size_t, Mode::dim> idx, opp_idx;
// 		idx[static_cast<char>(d)] = fixed_index();
// 		opp_idx[static_cast<char>(d)] = opp_index();
		

// 		// which dimensions are the outer dimensions
// 		std::array<std::size_t, Mode::dim-1> odims;
// 		int dout = 0;
// 		for (auto i=0; i<Mode::dim; i++){
// 			if (i == static_cast<char>(d)) continue;
// 			odims[dout] = i;
// 			dout++;
// 		}

// 		// the sizes of the outer dimensions
// 		std::array<std::size_t, Mode::dim-1> odsizes;
// 		for (auto i=0; i<Mode::dim-1; i++){
// 			odsizes[i] = mDims[odims[i]];
// 		}

// 		std::array<std::size_t, Mode::dim-1> odidx;
// 		for (auto i=0; i<Mode::dim-1; i++) odidx[i] = 0;

// 		// need to loop over Mode::dim-1 dimensions
// 		// outer loop over indices of other dimensions
// 		while (odidx[0]<odsizes[0]){
// 			for (auto i=0; i<Mode::dim-1; i++) {
// 				idx[odims[i]] = odidx[i];
// 				opp_idx[odims[i]] = odidx[i];
// 			}

// 			apply_update(mIdxObj, idx, opp_idx);

// 			// advance odidx
// 			int dd = odidx.size()-1;
// 			while (dd > 0){
// 				odidx[dd]++;
// 				odidx[dd-1] += odidx[dd]/odims[dd];
// 				odidx[dd] = odidx[dd]%odims[dd];
// 				dd--;
// 			}
// 		}
// 	};

// 	// static void update(){
// 	// };
// };



//************************************************************
//************************************************************
//***************** PERIODIC BOUNDARY ************************
//************************************************************
//************************************************************
//************************************************************

template <class Mode, Dir d, Orientation o, class IndexableObject>
struct UpdateBoundaryPeriodicD{
private:
	std::array<std::size_t, Mode::dim> mDims;
	IndexableObject * mIdxObj;

	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MIN, std::size_t>::type
	fixed_index(){return 0;};

	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MAX, std::size_t>::type
	fixed_index(){return mDims[static_cast<char>(d)]-1;};


	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MAX, std::size_t>::type
	opp_index(){return 0;};

	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MIN, std::size_t>::type
	opp_index(){return mDims[static_cast<char>(d)]-1;};

	template <typename CellType, Orientation T = o>
	typename std::enable_if<T == Orientation::MAX, CellType>::type &
	get_neighbor(CellType && t){return t.getNeighborMax(static_cast<char>(d));};

	template <typename CellType, Orientation T = o>
	typename std::enable_if<T == Orientation::MIN, CellType>::type &
	get_neighbor(CellType && t){return t.getNeighborMin(static_cast<char>(d));};

	template <class EMField, std::size_t... Idx>
	void apply_impl(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds, std::index_sequence<Idx...>){
		GetField<EMField>::get(get_neighbor(idxobj(inds[Idx]...))) = GetField<EMField>::get(idxobj(pinds[Idx]...));
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, ThreeD>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		apply_impl<Ex>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
		apply_impl<Ey>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
		apply_impl<Ez>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, TE>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		apply_impl<Ex>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
		apply_impl<Ey>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, TM>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		apply_impl<Ez>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, TEM>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		apply_impl<Ez>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
	}

public:
	UpdateBoundaryPeriodicD(IndexableObject & idx_obj): mIdxObj(&idx_obj) {};

	void operator()(){
		// starting inds
		std::array<std::size_t, Mode::dim> idx, opp_idx;
		idx[static_cast<char>(d)] = fixed_index();
		opp_idx[static_cast<char>(d)] = opp_index();
		

		// which dimensions are the outer dimensions
		std::array<std::size_t, Mode::dim-1> odims;
		int dout = 0;
		for (auto i=0; i<Mode::dim; i++){
			if (i == static_cast<char>(d)) continue;
			odims[dout] = i;
			dout++;
		}

		// the sizes of the outer dimensions
		std::array<std::size_t, Mode::dim-1> odsizes;
		for (auto i=0; i<Mode::dim-1; i++){
			odsizes[i] = mDims[odims[i]];
		}

		std::array<std::size_t, Mode::dim-1> odidx;
		for (auto i=0; i<Mode::dim-1; i++) odidx[i] = 0;

		// need to loop over Mode::dim-1 dimensions
		// outer loop over indices of other dimensions
		while (odidx[0]<odsizes[0]){
			for (auto i=0; i<Mode::dim-1; i++) {
				idx[odims[i]] = odidx[i];
				opp_idx[odims[i]] = odidx[i];
			}

			apply_update(mIdxObj, idx, opp_idx);

			// advance odidx
			int dd = odidx.size()-1;
			while (dd > 0){
				odidx[dd]++;
				odidx[dd-1] += odidx[dd]/odims[dd];
				odidx[dd] = odidx[dd]%odims[dd];
				dd--;
			}
		}
	};

	// static void update(){
	// };
};





template <class Mode, Dir d, Orientation o, class IndexableObject>
struct UpdateBoundaryPeriodicB{
private:
	std::array<std::size_t, Mode::dim> mDims;
	IndexableObject * mIdxObj;

	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MIN, std::size_t>::type
	fixed_index(){return 0;};

	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MAX, std::size_t>::type
	fixed_index(){return mDims[static_cast<char>(d)]-1;};


	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MAX, std::size_t>::type
	opp_index(){return 0;};

	template <Orientation T = o>
	typename std::enable_if<T==Orientation::MIN, std::size_t>::type
	opp_index(){return mDims[static_cast<char>(d)]-1;};

	template <typename CellType, Orientation T = o>
	typename std::enable_if<T == Orientation::MAX, CellType>::type &
	get_neighbor(CellType && t){return t.getNeighborMax(static_cast<char>(d));};

	template <typename CellType, Orientation T = o>
	typename std::enable_if<T == Orientation::MIN, CellType>::type &
	get_neighbor(CellType && t){return t.getNeighborMin(static_cast<char>(d));};

	template <class EMField, std::size_t... Idx>
	void apply_impl(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds, std::index_sequence<Idx...>){
		GetField<EMField>::get(get_neighbor(idxobj(inds[Idx]...))) = GetField<EMField>::get(idxobj(pinds[Idx]...));
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, ThreeD>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		apply_impl<Hx>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
		apply_impl<Hy>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
		apply_impl<Hz>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, TE>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		apply_impl<Hz>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});		
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, TM>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		
		apply_impl<Hx>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
		apply_impl<Hy>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
	}

	template <class M = Mode>
	typename std::enable_if<std::is_same<M, TEM>::value, void>
	apply_update(IndexableObject & idxobj, std::array<std::size_t, Mode::dim> & inds, std::array<std::size_t, Mode::dim> & pinds){
		apply_impl<Hy>(idxobj, inds, pinds, std::make_index_sequence<Mode::dim>{});
	}

public:
	UpdateBoundaryPeriodicB(IndexableObject & idx_obj): mIdxObj(&idx_obj) {};

	void operator()(){
		// starting inds
		std::array<std::size_t, Mode::dim> idx, opp_idx;
		idx[static_cast<char>(d)] = fixed_index();
		opp_idx[static_cast<char>(d)] = opp_index();
		

		// which dimensions are the outer dimensions
		std::array<std::size_t, Mode::dim-1> odims;
		int dout = 0;
		for (auto i=0; i<Mode::dim; i++){
			if (i == static_cast<char>(d)) continue;
			odims[dout] = i;
			dout++;
		}

		// the sizes of the outer dimensions
		std::array<std::size_t, Mode::dim-1> odsizes;
		for (auto i=0; i<Mode::dim-1; i++){
			odsizes[i] = mDims[odims[i]];
		}

		std::array<std::size_t, Mode::dim-1> odidx;
		for (auto i=0; i<Mode::dim-1; i++) odidx[i] = 0;

		// need to loop over Mode::dim-1 dimensions
		// outer loop over indices of other dimensions
		while (odidx[0]<odsizes[0]){
			for (auto i=0; i<Mode::dim-1; i++) {
				idx[odims[i]] = odidx[i];
				opp_idx[odims[i]] = odidx[i];
			}

			apply_update(mIdxObj, idx, opp_idx);

			// advance odidx
			int dd = odidx.size()-1;
			while (dd > 0){
				odidx[dd]++;
				odidx[dd-1] += odidx[dd]/odims[dd];
				odidx[dd] = odidx[dd]%odims[dd];
				dd--;
			}
		}
	};

	// static void update(){
	// };
};





//************************************************************
//************************************************************
//************ BLOCH-PERIODIC BOUNDARY ***********************
//************************************************************
//************************************************************
//************************************************************

// Define Bloch-periodic difference operator


//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************

// Define parallel difference operator


}// end namespace fdtd

#endif