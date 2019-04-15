#ifndef _BOUNDARYUPDATES_H
#define _BOUNDARYUPDATES_H

#include "FDTDConstants.hpp"
#include "YeeUpdates.hpp"
#include "DefaultInterfaces.hpp"


#include <functional>

namespace fdtd{


// restrict the set of possible boundary
// conditions to these
enum class Boundary : char {
	Periodic,
	BlochPeriodic,
	PEC,
	PMC,
	Symmetric,
	Antisymmetric,
	Parallel,
	NONE
};
template <> struct NameArray<Boundary>{
  static constexpr std::array<const char *, 8> value = {"Periodic", 
													   "BlochPeriodic", 
													   "PEC",
													   "PMC",
													   "Symmetric",
													   "Antisymmetric",
													   "Parallel", 
													   "NONE"};};
constexpr std::array<const char *, 8> NameArray<Boundary>::value;




//***************************************************************

// these are all masks, which are compile-time constants
// which are applied non-uniformly to field components in order
// to enforce certain boundary conditions (PEC, PMC, symmetries)
template <typename EMField, Boundary b, Dir d, Orientation o>
struct Mask {
	static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
	static constexpr double value = 1;
};



//***************************************************************

// the decorator class is a runtime base class intended to 
// be used for applying the Bloch-periodic boundary condition
// phase factor. It is applied uniformly to all field components
struct DefaultDecorator{
	constexpr double decorate() const {return 1.0;};
};



//***************************************************************


// Single cell update from source to destination
// by applying a mask (compile-time) and a 
// decorator (runtime)
template <typename Mode, Boundary b, Dir d, Orientation o, 
		  typename DecoratorPolicy = DefaultDecorator>
struct SingleUpdate : public DecoratorPolicy {
private:
	static_assert(std::is_base_of<EMMode, Mode>::value, "Mode must be a valid EMMode");
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// adapter for the for_each_tuple_type.
	// atomic PEC update for a single field
	template <typename EMField>
	struct atomic_update{
		static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");

		template <typename SourceIterator, typename DestIterator>
		static void get(SourceIterator && sit, DestIterator && dit, DecoratorPolicy & dec){
			GetField<EMField>::get(*dit) = dec.decorate() * Mask<EMField, b, d, o>::value * GetField<EMField>::get(*sit);
		}
	};
	

public:

	template <typename SourceIterator, typename DestIterator>
	void update(SourceIterator sit, DestIterator dit){
		Detail::for_each_tuple_type<typename FieldComponents<Mode>::electric, atomic_update>(sit, dit, *this);
		Detail::for_each_tuple_type<typename FieldComponents<Mode>::magnetic, atomic_update>(sit, dit, *this);
	}

	// turn this into a functor to enable runtime polymorphism
	template <typename SourceIterator, typename DestIterator>
	void operator()(SourceIterator sit, DestIterator dit){update(std::forward<SourceIterator>(sit), std::forward<DestIterator>(dit));};

};

//***************************************************************


// an intermediate update struct that adapts the statically 
// polymorphic structs to a std::function in order to enable
// runtime polymorphism
//
// this is a generalized framework for this sort of behavior,
// but it might be more efficient to redefine the template 
// arguments depending on the situation. Either way, they all
// become a functor in the end
template <typename SourceIterator, typename DestIterator>
struct BoundaryUpdater{
private:
	SourceIterator 		mSrcBegin, mSrcEnd;
	DestIterator		mDestBegin, mDestEnd;
	std::function<void(SourceIterator, DestIterator)>	mFunct; // takes two iterators and applies the update

public:
	BoundaryUpdater(SourceIterator sBeg, SourceIterator sEnd,
					DestIterator   dBeg, DestIterator 	dEnd,
					std::function<void(SourceIterator, DestIterator)> f)
	: mSrcBegin(sBeg), mSrcEnd(sEnd), mDestBegin(dBeg), mDestEnd(dEnd)
	, mFunct(f)
	{};

	void update(){
		for(auto it = std::make_pair(mSrcBegin, mDestBegin);
			it.first != mSrcEnd ; 
			it.first++, it.second++){
			mFunct(it.first, it.second);
		}
	}

	// turn this into a functor to enable runtime polymorphism
	void operator()(void) {update();}; 
};


template <typename SourceIterator, typename DestIterator>
BoundaryUpdater<SourceIterator, DestIterator> make_boundary_updater(
					SourceIterator sBeg, SourceIterator sEnd,
					DestIterator   dBeg, DestIterator 	dEnd,
					std::function<void(SourceIterator, DestIterator)> f){
	return BoundaryUpdater<SourceIterator, DestIterator>(sBeg, sEnd, dBeg, dEnd, f);
}



//***************************************************************

// this is the runtime interface for boundary updates
struct BoundaryData{
private:
	Boundary 						mType;
	Dir 							mDir;
	Orientation 					mOr;

	std::function<void(void)> 		mUpdateE;
	std::function<void(void)>		mUpdateH;

public:

	BoundaryData(Boundary b, 
				 Dir 	  d,
				 Orientation o,
				 std::function<void(void)> updateE,
				 std::function<void(void)> updateH)
	: mType(b), mDir(d), mOr(o), mUpdateE(updateE), mUpdateH(updateH) 
	{};

	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<Boundary>" << std::endl;
		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << "<Type>" << NameArray<Boundary>::value[static_cast<int>(mType)] << "</Type>" << std::endl;
		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << "<Dir>" << NameArray<Dir>::value[static_cast<int>(mDir)] << "</Dir>" << std::endl ;
		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << "<Orientation>" << NameArray<Orientation>::value[static_cast<int>(mOr)] << "</Orientation>" << std::endl ;
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</Boundary>" << std::endl;
	}
};

//***************************************************************



// iterator that wraps another type of iterator
// the dereference operator applies a functor in order
// to access the correct value
//
// the dereference functor is expected to have the 
// static function ::get()
template <typename Iterator, typename DereferenceFunctor>
struct DerivedIterator{
private:
	Iterator mIt;
public:

	typedef DerivedIterator							self_type;
	typedef typename Iterator::difference_type		difference_type;
    typedef typename Iterator::value_type 			value_type;
    typedef typename Iterator::reference 			reference;
    typedef typename Iterator::pointer 				pointer;
    typedef typename Iterator::iterator_category 	iterator_category;

	// construction
	DerivedIterator(Iterator it) : mIt(it) {};

	// copy-construct
	DerivedIterator(const DerivedIterator & cit) : mIt(cit.mIt) {};

	// copy assignment
	DerivedIterator & operator=(const DerivedIterator & cit){
		DerivedIterator i(cit);
		std::swap(i,*this);
		return *this;
	}

	pointer operator->() const {return &DereferenceFunctor::get(mIt);};
	reference operator*() const {return DereferenceFunctor::get(mIt);};

	// increment operators
	self_type operator++(){
		mIt++;
		return *this;
	}
	self_type operator++(int blah){
		mIt++;
		return *this;
	}

	// decrement operators
	self_type operator--(){
		mIt--;
		return *this;
	}
	self_type operator--(int blah){
		mIt--;
		return *this;
	}

	// equivalence operators
	bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt;};
	bool operator==(const self_type & leaf) const {return mIt == leaf.mIt;};

	// // random access iterator operators
	// template <typename T = base_iterator>
	// typename std::enable_if<std::is_same<typename std::iterator_traits<T>::iterator_category, std::random_access_iterator_tag>::value, difference_type>::type 
	// operator-(const self_type & it) const {return mIt-it.mIt;};

};

//***************************************************************


template <Dir d>
struct GetNeighborMin{
	template <typename IteratorType>
	static decltype(auto) get(IteratorType && it){return (*it).getNeighborMin(d);};
};


template <typename Iterator, Dir d>
struct MinIteratorDef{
	typedef DerivedIterator<Iterator, GetNeighborMin<d>> type;
};
template <typename Iterator, Dir d>
using MinIterator = typename MinIteratorDef<Iterator, d>::type;



//***************************************************************

typedef std::map<std::string, BoundaryData>	BoundaryMap;



//************************************************************
//************************************************************
//********************* PEC BOUNDARY *************************
//************************************************************
//************************************************************
//************************************************************


// PEC
template <typename EMField, Dir d, Orientation o>
struct Mask<EMField, Boundary::PEC, d, o>{
	static constexpr double value = (IsElectric<EMField>::value ? 
										(FieldDir<EMField>::value == d ?  1 : -1) :
										(FieldDir<EMField>::value == d ? -1 :  1) 
									);
};


// use the default decorator

// helper function to create this structure without passing the data types explicitly
template <typename Mode, Dir d, Orientation o, typename Iterator>
BoundaryData make_pec_boundary(Iterator begit, Iterator endit){
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// source iterator
	typedef std::conditional_t<o==Orientation::MIN, Iterator, MinIterator<Iterator, d>>			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, Iterator>			DestIterator;

	// single update
	typedef SingleUpdate<Mode, Boundary::PEC, d, o> 	SingleUpdateType;
	SingleUpdateType s;

	// function type
	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

	// whole-boundary update
	auto upE = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));
	auto upH = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));

	return BoundaryData(Boundary::PEC, d, o, upE, upH);
};


//************************************************************
//************************************************************
//********************* PMC BOUNDARY *************************
//************************************************************
//************************************************************
//************************************************************


// PMC
template <typename EMField, Dir d, Orientation o>
struct Mask<EMField, Boundary::PMC, d, o>{
	static constexpr double value = (IsElectric<EMField>::value ? 
										(FieldDir<EMField>::value == d ?  -1 : 1) :
										(FieldDir<EMField>::value == d ?  1 : -1) 
									);
};

// use the default decorator

// helper function to create this structure without passing the data types explicitly
template <typename Mode, Dir d, Orientation o, typename Iterator>
BoundaryData make_pmc_boundary(Iterator begit, Iterator endit){
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// source iterator
	typedef std::conditional_t<o==Orientation::MIN, Iterator, MinIterator<Iterator, d>>			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, Iterator>			DestIterator;

	// single update
	typedef SingleUpdate<Mode, Boundary::PMC, d, o> 	SingleUpdateType;
	SingleUpdateType s;

	// function type
	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

	// whole-boundary update
	auto upE = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));
	auto upH = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));

	return BoundaryData(Boundary::PMC, d, o, upE, upH);
};

//************************************************************
//************************************************************
//****************** SYMMETRY BOUNDARY ***********************
//************************************************************
//************************************************************
//************************************************************


// Symmetry
template <typename EMField, Dir d, Orientation o>
struct Mask<EMField, Boundary::Symmetric, d, o>{
	static constexpr double value = (IsElectric<EMField>::value ? 
										(FieldDir<EMField>::value == d ? -1 : 1) :
										(FieldDir<EMField>::value == d ? -1 : 1) 
									);
};


// use the default decorator

// helper function to create this structure without passing the data types explicitly
template <typename Mode, Dir d, Orientation o, typename Iterator>
BoundaryData make_symmetric_boundary(Iterator begit, Iterator endit){
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// source iterator
	typedef std::conditional_t<o==Orientation::MIN, Iterator, MinIterator<Iterator, d>>			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, Iterator>			DestIterator;

	// single update
	typedef SingleUpdate<Mode, Boundary::Symmetric, d, o> 	SingleUpdateType;
	SingleUpdateType s;

	// function type
	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

	// whole-boundary update
	auto upE = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));
	auto upH = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));

	return BoundaryData(Boundary::Symmetric, d, o, upE, upH);
};



//************************************************************
//************************************************************
//****************** ANTI-SYMMETRY BOUNDARY ******************
//************************************************************
//************************************************************
//************************************************************

// Anti-Symmetry
template <typename EMField, Dir d, Orientation o>
struct Mask<EMField, Boundary::Antisymmetric, d, o>{
	static constexpr double value = (IsElectric<EMField>::value ? 
										(FieldDir<EMField>::value == d ? 1 : -1) :
										(FieldDir<EMField>::value == d ? 1 : -1) 
									);
};


// use the default decorator

// helper function to create this structure without passing the data types explicitly
template <typename Mode, Dir d, Orientation o, typename Iterator>
BoundaryData make_antisymmetric_boundary(Iterator begit, Iterator endit){
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// source iterator
	typedef std::conditional_t<o==Orientation::MIN, Iterator, MinIterator<Iterator, d>>			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, Iterator>			DestIterator;

	// single update
	typedef SingleUpdate<Mode, Boundary::Antisymmetric, d, o> 	SingleUpdateType;
	SingleUpdateType s;

	// function type
	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

	// whole-boundary update
	auto upE = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));
	auto upH = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));

	return BoundaryData(Boundary::Antisymmetric, d, o, upE, upH);
};


//************************************************************
//************************************************************
//***************** PERIODIC BOUNDARY ************************
//************************************************************
//************************************************************
//************************************************************



// // an intermediate update struct that adapts the statically 
// // polymorphic structs to a std::function in order to enable
// // runtime polymorphism
// //
// // this is a generalized framework for this sort of behavior,
// // but it might be more efficient to redefine the template 
// // arguments depending on the situation. Either way, they all
// // become a functor in the end
// template <typename SourceIterator, typename DestIterator>
// struct BoundaryUpdater{
// private:
// 	SourceIterator 		mSrcBegin, mSrcEnd;
// 	DestIterator		mDestBegin, mDestEnd;
// 	std::function<void(SourceIterator, DestIterator)>	mFunct; // takes two iterators and applies the update

// public:
// 	BoundaryUpdater(SourceIterator sBeg, SourceIterator sEnd,
// 					DestIterator   dBeg, DestIterator 	dEnd,
// 					std::function<void(SourceIterator, DestIterator)> f)
// 	: mSrcBegin(sBeg), mSrcEnd(sEnd), mDestBegin(dBeg), mDestEnd(dEnd)
// 	, mFunct(f)
// 	{};

// 	void update(){
// 		for(auto it = std::make_pair(mSrcBegin, mDestBegin);
// 			it.first != mSrcEnd ; 
// 			it.first++, it.second++){
// 			mFunct(it.first, it.second);
// 		}
// 	}

// 	// turn this into a functor to enable runtime polymorphism
// 	void operator()(void) {update();}; 
// };


// template <typename SourceIterator, typename DestIterator>
// BoundaryUpdater<SourceIterator, DestIterator> make_boundary_updater(
// 					SourceIterator sBeg, SourceIterator sEnd,
// 					DestIterator   dBeg, DestIterator 	dEnd,
// 					std::function<void(SourceIterator, DestIterator)> f){
// 	return BoundaryUpdater<SourceIterator, DestIterator>(sBeg, sEnd, dBeg, dEnd, f);
// }



// // use the default mask

// // use the default decorator

// ************************ FINISH ME
// // helper function to create this structure without passing the data types explicitly
// template <typename Mode, Dir d, Orientation o, typename Iterator>
// BoundaryData make_periodic_boundary(Iterator begit, Iterator endit){
// 	static_assert(d != Dir::NONE, "Must have a valid direction!");
// 	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

// 	// source iterator
// 	typedef std::conditional_t<o==Orientation::MIN, Iterator, MinIterator<Iterator, d>>			SourceIterator;

// 	// dest iterator
// 	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, Iterator>			DestIterator;

// 	// single update
// 	typedef SingleUpdate<Mode, Boundary::Periodic, d, o> 	SingleUpdateType;
// 	SingleUpdateType s;

// 	// function type
// 	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

// 	// whole-boundary update
// 	auto upE = make_periodic_boundary_updater(SourceIterator(begit), SourceIterator(endit),
// 									 DestIterator(begit), DestIterator(endit),
// 									 static_cast<FunctorType>(s));
// 	auto upH = make_periodic_boundary_updater(SourceIterator(begit), SourceIterator(endit),
// 									 DestIterator(begit), DestIterator(endit),
// 									 static_cast<FunctorType>(s));

// 	return BoundaryData(Boundary::Periodic, d, o, upE, upH);
// };




// //************************************************************
// //************************************************************
// //************ BLOCH-PERIODIC BOUNDARY ***********************
// //************************************************************
// //************************************************************
// //************************************************************

// // Define Bloch-periodic difference operator


// //************************************************************
// //************************************************************
// //**************** PARALLEL BOUNDARY *************************
// //************************************************************
// //************************************************************
// //************************************************************

// // an intermediate update struct that adapts the statically 
// // polymorphic structs to a std::function in order to enable
// // runtime polymorphism
// //
// // this is a generalized framework for this sort of behavior,
// // but it might be more efficient to redefine the template 
// // arguments depending on the situation. Either way, they all
// // become a functor in the end
// template <typename ValueType, typename SourceIterator, typename DestIterator>
// struct BufferBoundaryUpdater{
// private:
// 	SourceIterator 		mSrcBegin, mSrcEnd;
// 	DestIterator		mDestBegin, mDestEnd;
// 	std::function<void(SourceIterator, DestIterator)>	mFunct; // takes two iterators and applies the update
// 	std::vector<ValueType> 		mBuffer;

// public:
// 	BufferBoundaryUpdater(SourceIterator sBeg, SourceIterator sEnd,
// 					DestIterator   dBeg, DestIterator 	dEnd,
// 					std::function<void(SourceIterator, DestIterator)> f)
// 	: mSrcBegin(sBeg), mSrcEnd(sEnd), mDestBegin(dBeg), mDestEnd(dEnd)
// 	, mFunct(f)
// 	{};

// 	void update(){
// 		for(auto it = std::make_pair(mSrcBegin, mDestBegin);
// 			it.first != mSrcEnd ; 
// 			it.first++, it.second++){
// 			mFunct(it.first, it.second);
// 		}
// 	}

// 	// turn this into a functor to enable runtime polymorphism
// 	void operator()(void) {update();}; 
// };


// template <typename SourceIterator, typename DestIterator>
// BufferBoundaryUpdater<SourceIterator, DestIterator> make_boundary_updater(
// 					SourceIterator sBeg, SourceIterator sEnd,
// 					DestIterator   dBeg, DestIterator 	dEnd,
// 					std::function<void(SourceIterator, DestIterator)> f){
// 	return BufferBoundaryUpdater<SourceIterator, DestIterator>(sBeg, sEnd, dBeg, dEnd, f);
// }


}// end namespace fdtd

#endif