#ifndef _DATACAPTURE_H
#define _DATACAPTURE_H


#include "FDTDConstants.hpp"


namespace fdtd{


// Data Capture Types
// 
// Volume Data - n_per_period, fields
// 
// Time History - location, fields
//
// Fluxes - 
//
//
//


template <typename ... Fields>
struct VolumeData{
private:
	double mOutputPeriod;

	template <typename T>
	struct atomic_print{
		static void get(std::ostream & os, unsigned int ntabs){
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<" << T::name << ">" << "</" << T::name << ">" << std::endl;
		}
	};
public:
	VolumeData(double period) : mOutputPeriod(period) {};


	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<VolumeData>" << std::endl;
			for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << "<Period>" << mOutputPeriod << "</Period>" << std::endl;

			nested_for_each_tuple_type<atomic_print, std::tuple<Fields...>>(os, ntabs+1);

		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</VolumeData>" << std::endl;
	}
};


template <typename ... Fields>
VolumeData<Fields...> make_volume_data(double period){return VolumeData<Fields...>(period);};


// template <std::tuple<typename ... Fields>>
// VolumeData<Fields...> make_volume_data(double period){return VolumeData<Fields...>(period);};


struct DataCapture{
private:
	// bool mConds[3][2];
	// std::size_t mNPML;

public:

	DataCapture() {};

	// std::size_t & N() {return mNPML;};
	// const std::size_t & N() const {return mNPML;};

	// template <Dir d, Orientation o>
	// bool & get() {return mConds[static_cast<char>(d)][static_cast<char>(o)];};
	// template <Dir d, Orientation o>
	// const bool & get() const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};


	// bool & get(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};
	// const bool & get(Dir d, Orientation o) const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};


	// bool & operator()(Dir d, Orientation o) {return mConds[static_cast<char>(d)][static_cast<char>(o)];};	
	// const bool & operator()(Dir d, Orientation o) const {return mConds[static_cast<char>(d)][static_cast<char>(o)];};	


	void print_summary(std::ostream & os = std::cout, unsigned int ntabs=0) const{
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "<DataCapture>" << std::endl;
			// for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			// os << "<N>" << mNDataCapture << "</N>" << std::endl;
			// for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			// os << "<X>" <<  std::endl;
			// 	for (auto i=0; i<ntabs+2; i++) os << "\t" ;
			// 	os << "<MIN>" << get<Dir::X, Orientation::MIN>() << " </MIN>" <<  std::endl;
			// 	for (auto i=0; i<ntabs+2; i++) os << "\t" ;
			// 	os << "<MAX>" << get<Dir::X, Orientation::MAX>() << " </MAX>" <<  std::endl;
			// for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			// os << "</X>" << std::endl ;
			
			// for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			// os << "<Y>" <<  std::endl;
			// 	for (auto i=0; i<ntabs+2; i++) os << "\t" ;
			// 	os << "<MIN>" << get<Dir::Y, Orientation::MIN>() << " </MIN>" <<  std::endl;
			// 	for (auto i=0; i<ntabs+2; i++) os << "\t" ;
			// 	os << "<MAX>" << get<Dir::Y, Orientation::MAX>() << " </MAX>" <<  std::endl;
			// for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			// os << "</Y>" << std::endl ;

			// for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			// os << "<Z>" <<  std::endl;
			// 	for (auto i=0; i<ntabs+2; i++) os << "\t" ;
			// 	os << "<MIN>" << get<Dir::Z, Orientation::MIN>() << " </MIN>" <<  std::endl;
			// 	for (auto i=0; i<ntabs+2; i++) os << "\t" ;
			// 	os << "<MAX>" << get<Dir::Z, Orientation::MAX>() << " </MAX>" <<  std::endl;
			// for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			// os << "</Z>" << std::endl ;
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		os << "</DataCapture>" << std::endl;
	}



	// #ifdef TINYXML2_INCLUDED
	// static PMLOptions readXML(std::string filename) {
	// 	PMLOptions po;

	// 	tinyxml2::XMLDocument doc;
	// 	doc.LoadFile(filename.c_str());

	// 	tinyxml2::XMLNode * n = doc.FirstChild();
	// 	auto c = (n->FirstChild());
				

	// 	while (c != nullptr){
	// 		std::stringstream ss;

	// 		if (!strcmp(c->Value(), "N")){
	// 			ss << c->FirstChild()->Value();
	// 			ss >> po.N();
	// 		}
	// 		else if(!strcmp(c->Value(), "X")){
	// 			tinyxml2::XMLNode * mm = c->FirstChild();
	// 			while (mm != nullptr){
	// 				if (!strcmp(mm->Value(), "MIN")){
	// 					ss << mm->FirstChild()->Value();
	// 					ss >> po.get<Dir::X, Orientation::MIN>();
	// 				}
	// 				else if (!strcmp(mm->Value(), "MAX")){
	// 					ss << mm->FirstChild()->Value();
	// 					ss >> po.get<Dir::X, Orientation::MAX>();
	// 				}
	// 				mm = mm->NextSibling();
	// 			}
	// 		}
	// 		else if(!strcmp(c->Value(), "Y")){
	// 			tinyxml2::XMLNode * mm = c->FirstChild();
	// 			while (mm != nullptr){
	// 				if (!strcmp(mm->Value(), "MIN")){
	// 					ss << mm->FirstChild()->Value();
	// 					ss >> po.get<Dir::Y, Orientation::MIN>();
	// 				}
	// 				else if (!strcmp(mm->Value(), "MAX")){
	// 					ss << mm->FirstChild()->Value();
	// 					ss >> po.get<Dir::Y, Orientation::MAX>();
	// 				}
	// 				mm = mm->NextSibling();
	// 			}
	// 		}
	// 		else if(!strcmp(c->Value(), "Z")){
	// 			tinyxml2::XMLNode * mm = c->FirstChild();
	// 			while (mm != nullptr){
	// 				if (!strcmp(mm->Value(), "MIN")){
	// 					ss << mm->FirstChild()->Value();
	// 					ss >> po.get<Dir::Z, Orientation::MIN>();
	// 				}
	// 				else if (!strcmp(mm->Value(), "MAX")){
	// 					ss << mm->FirstChild()->Value();
	// 					ss >> po.get<Dir::Z, Orientation::MAX>();
	// 				}
	// 				mm = mm->NextSibling();
	// 			}
	// 		}

	// 		c = (c->NextSibling());
	// 	}

	// 	return po;
	// }
	// #endif
};

} // end namespace fdtd
#endif