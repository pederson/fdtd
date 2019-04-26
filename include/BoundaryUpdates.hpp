#ifndef _BOUNDARYUPDATES_H
#define _BOUNDARYUPDATES_H

#include "FDTDConstants.hpp"
#include "YeeUpdates.hpp"
#include "DefaultInterfaces.hpp"


#include <functional>
#include <complex>
#include <vector>

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
static constexpr std::array<const char *, 8> value = {"Periodic", 		// translational symmetry
													   "BlochPeriodic", 	// translational symmetry
													   "PEC",				// reflection symmetry
													   "PMC",				// reflection symmetry
													   "Symmetric",			// reflection symmetry
													   "Antisymmetric",		// reflection symmetry
													   "Parallel", 			
													   "NONE"};};
constexpr std::array<const char *, 8> NameArray<Boundary>::value;


//***************************************************************


struct BoundaryOptions{
private:
	Boundary mConds[3][2];

public:

	BoundaryOptions() {
		for (int i=0; i<3; i++){
			for (int j=0; j<2; j++){
				mConds[i][j] = Boundary::PEC;
			} 
		}
	}

	template <Dir d, Orientation o>
	Boundary & get() {return mConds[static_cast<char>(d)][static_cast<char>(o)];};
	template <Dir d, Orientation o>
	const Boundary & get() const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};


	Boundary & get(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};
	const Boundary & get(Dir d, Orientation o) const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};


	Boundary & operator()(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};	
	const Boundary & operator()(Dir d, Orientation o) const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};	


	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<BoundaryConditions>" << std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<X>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MIN>" << NameArray<Boundary>::value[static_cast<char>(get<Dir::X, Orientation::MIN>())] << "</MIN>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MAX>" << NameArray<Boundary>::value[static_cast<char>(get<Dir::X, Orientation::MAX>())] << "</MAX>" <<  std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "</X>" << std::endl ;
			
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Y>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MIN>" << NameArray<Boundary>::value[static_cast<char>(get<Dir::Y, Orientation::MIN>())] << "</MIN>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MAX>" << NameArray<Boundary>::value[static_cast<char>(get<Dir::Y, Orientation::MAX>())] << "</MAX>" <<  std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "</Y>" << std::endl ;

			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Z>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MIN>" << NameArray<Boundary>::value[static_cast<char>(get<Dir::Z, Orientation::MIN>())] << "</MIN>" <<  std::endl;
				for (auto i=0; i<ntabs+2; i++) os << "\t" ;
				os << "<MAX>" << NameArray<Boundary>::value[static_cast<char>(get<Dir::Z, Orientation::MAX>())] << " </MAX>" <<  std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "</Z>" << std::endl ;
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</BoundaryConditions>" << std::endl;
	}



	#ifdef TINYXML2_INCLUDED
	static BoundaryOptions readXML(std::string filename) {
		BoundaryOptions bo;

		tinyxml2::XMLDocument doc;
		doc.LoadFile(filename.c_str());

		tinyxml2::XMLNode * n = doc.FirstChild();
		auto c = (n->FirstChild());
				

		while (c != nullptr){
			std::stringstream ss;

			if(!strcmp(c->Value(), "X")){
				tinyxml2::XMLNode * mm = c->FirstChild();
				while (mm != nullptr){
					if (!strcmp(mm->Value(), "MIN")){
						bo.get<Dir::X, Orientation::MIN>() = MapNameTo<Boundary>(mm->FirstChild()->Value());
					}
					else if (!strcmp(mm->Value(), "MAX")){
						bo.get<Dir::X, Orientation::MAX>() = MapNameTo<Boundary>(mm->FirstChild()->Value());
					}
					mm = mm->NextSibling();
				}
			}
			else if(!strcmp(c->Value(), "Y")){
				tinyxml2::XMLNode * mm = c->FirstChild();
				while (mm != nullptr){
					if (!strcmp(mm->Value(), "MIN")){
						bo.get<Dir::Y, Orientation::MIN>() = MapNameTo<Boundary>(mm->FirstChild()->Value());
					}
					else if (!strcmp(mm->Value(), "MAX")){
						bo.get<Dir::Y, Orientation::MAX>() = MapNameTo<Boundary>(mm->FirstChild()->Value());
					}
					mm = mm->NextSibling();
				}
			}
			else if(!strcmp(c->Value(), "Z")){
				tinyxml2::XMLNode * mm = c->FirstChild();
				while (mm != nullptr){
					if (!strcmp(mm->Value(), "MIN")){
						bo.get<Dir::Z, Orientation::MIN>() = MapNameTo<Boundary>(mm->FirstChild()->Value());
					}
					else if (!strcmp(mm->Value(), "MAX")){
						bo.get<Dir::Z, Orientation::MAX>() = MapNameTo<Boundary>(mm->FirstChild()->Value());
					}
					mm = mm->NextSibling();
				}
			}

			c = (c->NextSibling());
		}

		return bo;
	}
	#endif
};





//***************************************************************

// these are all masks, which are compile-time constants
// which are applied non-uniformly to field components in order
// to enforce certain boundary conditions (PEC, PMC, symmetries)
template <typename EMField, Boundary b, Dir d, Orientation o>
struct Mask {
	static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
	static constexpr double value = 1;
};
// declare here for the linker if needed at runtime
template <typename EMField, Boundary b, Dir d, Orientation o> constexpr double Mask<EMField, b, d, o>::value;



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
			// std::cout << "******** " << EMField::name << " **********" << std::endl;
			// std::cout << "before: " << GetField<EMField>::get(*dit) << std::endl;
			// std::cout << "\t decorator: " << dec.decorate() << std::endl;
			// std::cout << "\t mask: " <<  Mask<EMField, b, d, o>::value << std::endl;
			// std::cout << "\t source: " << GetField<EMField>::get(*sit) << std::endl;
			GetField<EMField>::get(*dit) = dec.decorate() * Mask<EMField, b, d, o>::value * GetField<EMField>::get(*sit);
			// std::cout << "after: " << GetField<EMField>::get(*dit) << std::endl;
		}
	};
	

public:

	SingleUpdate() {};
	template <typename... Args>
	SingleUpdate(Args... init) : DecoratorPolicy(init...) {};

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

	void updateE() {mUpdateE();};
	void updateH() {mUpdateH();};

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
// static function ::get() which passes the base iterator
// and the iterator itself. The DereferenceFunctor
// can store data, which can be useful. 
template <typename Iterator, typename DereferenceFunctor>
struct DerivedIterator : public DereferenceFunctor{
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
	// pass all extra constructor arguments on to the base class
	template <typename... Args>
	DerivedIterator(Iterator it, Args... init) : DereferenceFunctor(init...), mIt(it) {};

	// copy-construct
	DerivedIterator(const DerivedIterator & cit) : DereferenceFunctor(cit), mIt(cit.mIt) {};

	// copy assignment
	DerivedIterator & operator=(const DerivedIterator & cit){
		DerivedIterator i(cit);
		std::swap(i,*this);
		return *this;
	}

	pointer operator->() {return &DereferenceFunctor::get(mIt, *this);};
	reference operator*() {return DereferenceFunctor::get(mIt, *this);};

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
	static decltype(auto) get(IteratorType && it, const GetNeighborMin & g){return (*it).getNeighborMin(static_cast<char>(d));};
};


template <typename Iterator, Dir d>
struct MinIteratorDef{
	typedef DerivedIterator<Iterator, GetNeighborMin<d>> type;
};
template <typename Iterator, Dir d>
using MinIterator = typename MinIteratorDef<Iterator, d>::type;


template <Dir d>
struct GetNeighborMax{
	template <typename IteratorType>
	static decltype(auto) get(IteratorType && it, const GetNeighborMax & g){return (*it).getNeighborMax(static_cast<char>(d));};
};


template <typename Iterator, Dir d>
struct MaxIteratorDef{
	typedef DerivedIterator<Iterator, GetNeighborMax<d>> type;
};
template <typename Iterator, Dir d>
using MaxIterator = typename MaxIteratorDef<Iterator, d>::type;


//***************************************************************

// typedef std::map<std::string, BoundaryData>	BoundaryMap;



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
	typedef Iterator 			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

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


// // helper function to create this structure without passing the data types explicitly
// template <typename Mode, typename Iterator>
// BoundaryData make_pec_boundary(Dir d, Orientation o, Iterator begit, Iterator endit){
// 	assert(d != Dir::NONE, "Must have a valid direction!");
// 	assert(o != Orientation::NONE, "Must have a valid orientation!");

// 	// source iterator
// 	typedef Iterator 			SourceIterator;

// 	// dest iterator
// 	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

// 	// single update
// 	typedef SingleUpdate<Mode, Boundary::PEC, d, o> 	SingleUpdateType;
// 	SingleUpdateType s;

// 	// function type
// 	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

// 	// whole-boundary update
// 	auto upE = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
// 									 DestIterator(begit), DestIterator(endit),
// 									 static_cast<FunctorType>(s));
// 	auto upH = make_boundary_updater(SourceIterator(begit), SourceIterator(endit),
// 									 DestIterator(begit), DestIterator(endit),
// 									 static_cast<FunctorType>(s));

// 	return BoundaryData(Boundary::PEC, d, o, upE, upH);
// };


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
	typedef Iterator 			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

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
	typedef Iterator 			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

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
	typedef Iterator 			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

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



template <typename Functor>
struct GetNeighborPeriodic : public Functor{
	// Functor fctor;

	GetNeighborPeriodic() {};
	GetNeighborPeriodic(Functor & f) : Functor(f) {};

	template <typename IteratorType>
	static decltype(auto) get(IteratorType && it, GetNeighborPeriodic & g) {return g(it);};
};


template <typename Iterator, typename Functor>
struct PeriodicIteratorDef{
	typedef DerivedIterator<Iterator, GetNeighborPeriodic<Functor>> type;
};
template <typename Iterator, typename Functor>
using PeriodicIterator = typename PeriodicIteratorDef<Iterator, Functor>::type;


// use the default mask

// use the default decorator


// helper function to create this structure without passing the data types explicitly
// the Functor passed in here takes an iterator and returns a reference to the periodic neighbor
template <typename Mode, Dir d, Orientation o, typename Iterator, typename Functor>
BoundaryData make_periodic_boundary(Iterator begit, Iterator endit, Functor f){
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// source iterator
	typedef PeriodicIterator<Iterator, Functor> 	SourceIterator;
	// typedef std::conditional_t<o==Orientation::MIN, Iterator, PeriodicIterator<Iterator, Functor>>			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

	// single update
	typedef SingleUpdate<Mode, Boundary::Periodic, d, o> 	SingleUpdateType;
	SingleUpdateType s;

	// function type
	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

	// whole-boundary update
	auto upE = make_boundary_updater(SourceIterator(begit, f), SourceIterator(endit, f),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));
	auto upH = make_boundary_updater(SourceIterator(begit, f), SourceIterator(endit, f),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));

	return BoundaryData(Boundary::Periodic, d, o, upE, upH);
};




// //************************************************************
// //************************************************************
// //************ BLOCH-PERIODIC BOUNDARY ***********************
// //************************************************************
// //************************************************************
// //************************************************************

// use the default mask

// use a Bloch-periodic decorator that multiplies by the bloch-factor specified by the user
struct BlochDecorator{
	std::complex<double> mBlochFactor;

	BlochDecorator(std::complex<double> val) : mBlochFactor(val) {};

	constexpr std::complex<double> decorate() const {return mBlochFactor;};
};


// helper function to create this structure without passing the data types explicitly
// the Functor passed in here takes an iterator and returns a reference to the periodic neighbor
template <typename Mode, Dir d, Orientation o, typename Iterator, typename Functor>
BoundaryData make_bloch_periodic_boundary(Iterator begit, Iterator endit, Functor f, std::complex<double> BlochFactor){
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// source iterator
	typedef PeriodicIterator<Iterator, Functor> 	SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

	// single update
	typedef SingleUpdate<Mode, Boundary::BlochPeriodic, d, o, BlochDecorator> 	SingleUpdateType;
	SingleUpdateType s(BlochFactor);

	// function type
	typedef std::function<void(SourceIterator, DestIterator)> FunctorType;

	// whole-boundary update
	auto upE = make_boundary_updater(SourceIterator(begit, f), SourceIterator(endit, f),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));
	auto upH = make_boundary_updater(SourceIterator(begit, f), SourceIterator(endit, f),
									 DestIterator(begit), DestIterator(endit),
									 static_cast<FunctorType>(s));

	return BoundaryData(Boundary::BlochPeriodic, d, o, upE, upH);
};


// //************************************************************
// //************************************************************
// //**************** PARALLEL BOUNDARY *************************
// //************************************************************
// //************************************************************
// //************************************************************


#ifdef MPI_VERSION


// use the default mask

// use the default decorator



template <typename ValueType>
struct SendRecvBuffer{
	typedef std::vector<ValueType> 		VectorType;
	VectorType 		mSendBuffer;
	VectorType 		mRecvBuffer;

	void resize(std::size_t n) {mSendBuffer.resize(n); mRecvBuffer.resize(n);};
	VectorType & send() 	{return mSendBuffer;};
	VectorType & recv() 	{return mRecvBuffer;};
};

// an intermediate update struct made specifically
// for the parallel boundary updates
template <typename Mode, FieldType ft, typename SourceIterator, typename DestIterator, typename BufferPolicy>
struct BufferBoundaryUpdater : public BufferPolicy{
private:
	SourceIterator 		mSrcBegin, mSrcEnd;
	DestIterator		mDestBegin, mDestEnd;
	// std::function<void(SourceIterator, DestIterator)>	mFunct; // takes two iterators and applies the update
	int 				mNeighborProcRank;


	template <typename EMField>
	struct single_field_update{
		static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");


		static void buffer(BufferBoundaryUpdater & b){
			// std::cout << "attempting to buffer field " << EMField::name << std::endl;
			// first copy to the send buffer
			for(auto it = std::make_pair(b.mSrcBegin, b.send().begin());
				it.first != b.mSrcEnd ; 
				it.first++, it.second++){
				*(it.second) = GetField<EMField>::get(*it.first);
			}
			// std::cout << "successfully buffered field " << EMField::name << std::endl;
		}

		static void communicate(BufferBoundaryUpdater & b){
			// void * bs = &b.send().front();
			// void * rs = &b.recv().front();
			// std::cout << "attempting to communicate field " << EMField::name << std::endl;
			// std::cout << " to neighbor rank: " << b.mNeighborProcRank << std::endl;
			// std::cout << " send size: " << b.send().size() << std::endl;
			// std::cout << " recv size: " << b.recv().size() << std::endl;
			// communicate the possibly non-native datatypes stored in buffer
			// across neighbors



			// need to define an MPI datatype for non-native datatypes
			MPI_Status status;
			MPI_Sendrecv(&b.send().front(), b.send().size(), MPI_DOUBLE, b.mNeighborProcRank, 0,
						 &b.recv().front(), b.recv().size(), MPI_DOUBLE, b.mNeighborProcRank, 0,
						 MPI_COMM_WORLD, &status);
			// std::cout << "successfully communicated field " << EMField::name << std::endl;
		}

		static void update(BufferBoundaryUpdater & b){
			// std::cout << "attempting to update field " << EMField::name << std::endl;

			// now extract from the receive buffer
			for(auto it = std::make_pair(b.mDestBegin, b.recv().begin());
				it.first != b.mDestEnd; 
				it.first++, it.second++){
				GetField<EMField>::get(*it.first) = *(it.second);
			}
			// std::cout << "successfully updated field " << EMField::name << std::endl;
		}

		static void get(BufferBoundaryUpdater & b){
			int rnk;
			MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
			// std::cout << " on processor rank " << rnk << " with neighbor rank " << b.mNeighborProcRank << std::endl;
			// std::cout << " which has size: " << b.recv().size() << std::endl;
			buffer(b); update(b); communicate(b); update(b);
			
		}
	};

public:
	BufferBoundaryUpdater(SourceIterator sBeg, SourceIterator 	sEnd,
						  DestIterator   dBeg, DestIterator 	dEnd,
						  int neighb_rank)
	: mSrcBegin(sBeg), mSrcEnd(sEnd)
	, mDestBegin(dBeg), mDestEnd(dEnd)
	, mNeighborProcRank(neighb_rank)
	{
		std::size_t ct = 0;
		for (auto it=dBeg; it!=dEnd; it++) ct++;
		BufferPolicy::resize(ct);
	};


	// turn this into a functor to enable runtime polymorphism
	template <FieldType f = ft>
	std::enable_if_t<f == FieldType::Electric, void> operator()(void) {
		// std::cout << "attempting to communicate fields " << std::endl; ;
		Detail::for_each_tuple_type<typename FieldComponents<Mode>::electric, single_field_update>(*this);
	}; 

	template <FieldType f = ft>
	std::enable_if_t<f == FieldType::Magnetic, void> operator()(void) {
		Detail::for_each_tuple_type<typename FieldComponents<Mode>::magnetic, single_field_update>(*this);
	}; 
};

// std::remove_reference_t<decltype(*std::declval<SourceIterator>())>
template <typename Mode, FieldType ft, typename SourceIterator, typename DestIterator>
BufferBoundaryUpdater<Mode, ft, SourceIterator, DestIterator, SendRecvBuffer<double>>
make_boundary_buffer_updater(SourceIterator sBeg, SourceIterator sEnd,
					DestIterator   dBeg, DestIterator 	dEnd,
					int neighb_rank){
	return BufferBoundaryUpdater<Mode, ft, SourceIterator, DestIterator, SendRecvBuffer<double>>(sBeg, sEnd, dBeg, dEnd, neighb_rank);
}




// helper function to create this structure without passing the data types explicitly
template <typename Mode, Dir d, Orientation o, typename Iterator>
BoundaryData make_parallel_boundary(Iterator begit, Iterator endit, int nproc){
	static_assert(d != Dir::NONE, "Must have a valid direction!");
	static_assert(o != Orientation::NONE, "Must have a valid orientation!");

	// source iterator
	typedef Iterator 	SourceIterator;
	// typedef std::conditional_t<o==Orientation::MIN, Iterator, PeriodicIterator<Iterator, Functor>>			SourceIterator;

	// dest iterator
	typedef std::conditional_t<o==Orientation::MIN, MinIterator<Iterator, d>, MaxIterator<Iterator, d>>			DestIterator;

	// // single update
	// typedef SingleUpdate<Mode, Boundary::Parallel, d, o> 	SingleUpdateType;
	// SingleUpdateType s;

	// whole-boundary update
	auto upE = make_boundary_buffer_updater<Mode, FieldType::Electric>(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 nproc);
	auto upH = make_boundary_buffer_updater<Mode, FieldType::Magnetic>(SourceIterator(begit), SourceIterator(endit),
									 DestIterator(begit), DestIterator(endit),
									 nproc);

	return BoundaryData(Boundary::Parallel, d, o, upE, upH);
};

#endif

}// end namespace fdtd

#endif