#ifndef _SOURCE_H
#define _SOURCE_H

#include <map>

#include "FDTDConstants.hpp"

namespace fdtd{


template <typename CellType, typename FieldGetter>
struct Source{
private:
	typedef std::function<double(double)> function_type;
	function_type 		mFunc;
	CellType * 			mLoc;

public:
	Source() {};

	Source(function_type f, CellType & c)
	: mFunc(f), mLoc(&c) {};

	function_type & function() {return mFunc;};
	const function_type & function() const {return mFunc;};

	CellType * location() {return mLoc;};

	void apply(double t){
		FieldGetter::get(*mLoc) += mFunc(t);
	}

	// in order to make this into a functor
	void operator()(double t){return apply(t);};

	double at(double t) const {return mFunc(t);};

};


template <typename CellType, typename F>
Source<CellType, GetField<F>> make_source(std::function<double(double)> f, CellType & c){
	return Source<CellType, GetField<F>>(f, c);
}

}// end namespace fdtd

#endif