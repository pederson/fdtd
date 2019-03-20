// compile with:
// 			clang++-3.9 -std=c++14 -I../ testConstants.cpp -o testConstants


#include <fdtd.h>

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>


using namespace fdtd;


template <typename Field>
struct PrintName{
	static void get(){
		std::cout << Field::name ;
		std::cout << " offset: " ;
		std::cout << "(";
		std::cout << Field::off[0] << ", ";
		std::cout << Field::off[1] << ", ";
		std::cout << Field::off[2] ;
		std::cout << ")";
		std::cout << " isMagnetic?: " << (Field::field_type == FieldType::Magnetic) ;
		std::cout << " Neighbor side: " << (Field::neighb_side == Orientation::MAX ? 1 : -1) ;
		std::cout << " direction: " << static_cast<int>(Field::direction) << std::endl;
		std::cout << std::endl;
	}
};


int main(int argc, char * argv[]){

	std::cout << "3D: " << std::endl;
	Detail::for_each_tuple_type<FieldComponents<ThreeD>::electric, PrintName>();
	Detail::for_each_tuple_type<FieldComponents<ThreeD>::magnetic, PrintName>();
	std::cout << "*************" << std::endl;

	std::cout << "TE: " << std::endl;
	Detail::for_each_tuple_type<FieldComponents<TE>::electric, PrintName>();
	Detail::for_each_tuple_type<FieldComponents<TE>::magnetic, PrintName>();
	std::cout << "*************" << std::endl;

	std::cout << "TM: " << std::endl;
	Detail::for_each_tuple_type<FieldComponents<TM>::electric, PrintName>();
	Detail::for_each_tuple_type<FieldComponents<TM>::magnetic, PrintName>();
	std::cout << "*************" << std::endl;

	std::cout << "TEM: " << std::endl;
	Detail::for_each_tuple_type<FieldComponents<TEM>::electric, PrintName>();
	Detail::for_each_tuple_type<FieldComponents<TEM>::magnetic, PrintName>();
	std::cout << "*************" << std::endl;


	return 0;
}