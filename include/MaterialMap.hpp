#ifndef _MATERIALMAP_H
#define _MATERIALMAP_H

#include <map>
#include <functional>
#include <memory>

#include "FDTDConstants.hpp"
#include "DispersiveMaterials.hpp"

namespace fdtd{



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
	
// template <typename Derived>
// struct AddPrinter : public Derived{
// public:
// 	// template <typename StreamType>
// 	void print_summary(std::ostream & os, unsigned int ntabs) const;
// };

//*************************************************************************
// Usage example:
//			// construct
// 			material_map_type 		matmap;
// 			matmap["vacuum"] 		= make_material_pair<node_type>(ConstantUpdate<Mode>(1.0)							, ConstantUpdate<Mode>(1.0));
// 			matmap["dielectric"] 	= make_material_pair<node_type>(ConstantUpdate<Mode>(epsilon)						, ConstantUpdate<Mode>(mu));
//			matmap["conductor"] 	= make_material_pair<node_type>(ConductiveUpdate<Mode>(epsilon, coll_freq, delta_t)	, ConstantUpdate<Mode>(mu));
//			matmap["pec"] 			= make_material_pair<node_type>(ConductiveUpdate<Mode>(1.0, coll_freq_pec, delta_t)	, ConstantUpdate<Mode>(mu));
//
//			// use
//			CellType cell;
//			matmap["dielectric"].electric()(cell);
//**************************************************************************

// CellType is kept as a generic type so that it can accept
// both lvalue and rvalue references
template <typename CellType>
struct MaterialPair{
private:
	typedef std::function<void(std::add_lvalue_reference_t<CellType>)> function_type;
	function_type 		mElectric;
	function_type 		mMagnetic;

	typedef std::function<void(std::ostream &, unsigned int)> printer_type;
	printer_type 		mElectricPrinter;
	printer_type 		mMagneticPrinter;

	MaterialPair(){};

public:
	MaterialPair(function_type e, function_type m)
	: mElectric(e), mMagnetic(m) {};

	MaterialPair(function_type e, function_type m, printer_type pe, printer_type pm)
	: mElectric(e), mMagnetic(m), mElectricPrinter(pe), mMagneticPrinter(pm) {};

	function_type & electric() {return mElectric;};
	const function_type & electric() const {return mElectric;};

	function_type & magnetic() {return mMagnetic;};
	const function_type & magnetic() const {return mMagnetic;};


	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		std::cout << "<MaterialPair>" << std::endl;

		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << "<Electric>" << std::endl;
			mElectricPrinter(os, ntabs+2);
		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << "</Electric>" << std::endl;

		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << "<Magnetic>" << std::endl;
			mMagneticPrinter(os, ntabs+2);
		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << "</Magnetic>" << std::endl;

		for (auto i=0; i<ntabs; i++) os << "\t" ;
		std::cout << "</MaterialPair>" << std::endl;
	}


#ifdef TINYXML2_INCLUDED

private:
	template <typename Mode, bool forward>
	struct XMLReader{
		template <typename E, typename M>
		struct atomic_build{
			typedef std::function<void(std::add_lvalue_reference_t<CellType>)> ftype;
			typedef std::function<void(std::ostream &, unsigned int)> ptype;
			
			static void get(double dt, Detail::Dispersion de, Detail::Dispersion dm, tinyxml2::XMLNode * e, tinyxml2::XMLNode * m, MaterialPair & mp){
				if (de == E::value && dm == M::value){
					std::cout << "MaterialPair builder" << std::endl;
					// auto mp = make_material_pair(BuildDispersionXML<E>::get(e), BuildDispersionXML<M>::get(m));
					std::cout << NameArray<typename E::type>::value[static_cast<int>(de)] << " ----- " << NameArray<typename M::type>::value[static_cast<int>(dm)] << std::endl;

					typedef decltype(BuildDispersionXML<E, Mode, FieldType::Electric, forward>::get(e, dt)) EType;
					typedef decltype(BuildDispersionXML<M, Mode, FieldType::Magnetic, forward>::get(m, dt)) MType;
					
					EType eval = BuildDispersionXML<E, Mode, FieldType::Electric, forward>::get(e, dt);
					MType mval = BuildDispersionXML<M, Mode, FieldType::Magnetic, forward>::get(m, dt);

					// eval.print_summary();
					// mval.print_summary();

					mp = MaterialPair(static_cast<ftype>(eval), 
									static_cast<ftype>(mval),
									static_cast<ptype>(PrintSummaryWrapper<EType>(eval)),
									static_cast<ptype>(PrintSummaryWrapper<MType>(mval)));

					// mp->print_summary();
					// return mp;
					// mp.print_summary();
				}
			}
		};
	};

	
public:

	template <typename Mode, bool forward = false>
	static MaterialPair readXML(tinyxml2::XMLNode * n, double dt){
		typedef XMLReader<Mode, forward> ReaderType;

		// MaterialPair mp;
		std::string estring, mstring;
		tinyxml2::XMLNode * enode, * mnode;

		auto c = (n->FirstChild());

		while (c != nullptr){
			std::stringstream ss;

			if(!strcmp(c->Value(), "Electric")){
				enode = c->FirstChild();
				estring = c->FirstChild()->Value();
				// auto emat = ParseXMLMaterials<Mode, ft, forward>::get(c->FirstChild()); 
			}
			else if(!strcmp(c->Value(), "Magnetic")){
				mnode = c->FirstChild();
				mstring = c->FirstChild()->Value();
			}

			c = (c->NextSibling());
		}

		// std::cout << "XML MATERIAL PAIR: " << estring << ", " << mstring << std::endl;
		MaterialPair mp;
		fdtd::nested_for_each_tuple_type<ReaderType::template atomic_build, DispersionTuple, DispersionTuple>(dt, MapNameTo<Detail::Dispersion>(estring), 
																						 MapNameTo<Detail::Dispersion>(mstring),
																						 enode, mnode, mp);
		// mp.print_summary();
		return mp;
	}

	template <typename Mode, bool forward = false>
	static MaterialPair readXML(std::string filename, double dt) {
		tinyxml2::XMLDocument doc;
		doc.LoadFile(filename.c_str());

		tinyxml2::XMLNode * n = doc.FirstChild();
		return readXML<Mode, forward>(n, dt);
	}

#endif

};



template <typename CellType, typename Electric, typename Magnetic>
MaterialPair<CellType> make_material_pair(Electric && e, Magnetic && m){

	// static_assert(std::is_same<Electric, void>::value, "Debug");
	typedef std::function<void(std::add_rvalue_reference_t<CellType>)> ftype;
	typedef std::function<void(std::ostream &, unsigned int)> ptype;
	return MaterialPair<CellType>(static_cast<ftype>(e), static_cast<ftype>(m), 
								  static_cast<ptype>(PrintSummaryWrapper<Electric>(e)), static_cast<ptype>(PrintSummaryWrapper<Magnetic>(m)));
};



//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************
//************************************************************


template <typename Identifier, typename CellType>
struct MaterialMap : public std::map<Identifier, MaterialPair<CellType>>
{
	typedef std::map<Identifier, MaterialPair<CellType>> MapType;
	using MapType::begin;
	using MapType::end;

	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const
	{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<MaterialMap>" << std::endl;
			for (auto it=begin(); it!=end(); it++){

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				std::cout << "<Item>" << std::endl;


					for (auto i=0; i<ntabs+2; i++) os << "\t" ;
					std::cout << "<Identifier>" << it->first << "</Identifier>" << std::endl;

					it->second.print_summary(os, ntabs+2);

				for (auto i=0; i<ntabs+1; i++) os << "\t" ;
				std::cout << "</Item>" << std::endl;
			}
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</MaterialMap>" << std::endl;
	}



	#ifdef TINYXML2_INCLUDED

		
	public:

		template <typename Mode, bool forward = false>
		static MaterialMap readXML(tinyxml2::XMLNode * n, double dt){
			MaterialMap mp;

			std::string id;
			tinyxml2::XMLNode * mnode;

			auto c = (n->FirstChild());

			while (c != nullptr){
				std::stringstream ss;

				if(!strcmp(c->Value(), "Item")){
					tinyxml2::XMLNode * mm = c->FirstChild();
					
					while (mm != nullptr){
						if(!strcmp(mm->Value(), "Identifier")){
							id = mm->FirstChild()->Value();
						}
						else if(!strcmp(mm->Value(), "MaterialPair")){
							mnode = mm;
						}

						mm = (mm->NextSibling());
					}
					mp.insert(std::make_pair(id, MaterialPair<CellType>::template readXML<Mode, forward>(mnode, dt)));
				}

				c = (c->NextSibling());
			}

			return mp;
		}

		template <typename Mode, bool forward = false>
		static MaterialMap readXML(std::string filename, double dt) {
			tinyxml2::XMLDocument doc;
			doc.LoadFile(filename.c_str());

			tinyxml2::XMLNode * n = doc.FirstChild();
			return readXML<Mode, forward>(n, dt);
		}

	#endif
};


}// end namespace fdtd

#endif