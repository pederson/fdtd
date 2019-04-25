#ifndef _FDTDCONSTANTS_H
#define _FDTDCONSTANTS_H

#include <limits>
#include <cmath>
#include <array>
#include <memory>
#include <iostream>

namespace fdtd{


// declare a stored private member variable
// of given type
// and expose public accessors functions with
// given name
#define FDTD_DECLARE_MEMBER(type, name)     \
  private:                                  \
    type  m ## name;                        \
  public:                                   \
    const type & name() const {return m ## name;};  \
    type &     name()     {return m ## name;};



// struct for naming things... easier to output generalized
// stuff to streams this way
template <typename T> struct NameArray{
  static_assert(std::is_enum<T>::value, "NameArray can only be defined for enum classes!");
};


// function to map strings to an enum which has a defined NameArray
template <typename T>
T MapNameTo(std::string name){
  int idx;
  bool found = false;
  for (auto i=0; i!=NameArray<T>::value.size(); i++){
    if (!strcmp(name.c_str(), NameArray<T>::value[i])){
      idx = i;
      found = true;
    }
  }

  if (!found){
    std::cerr << "ERROR: Could not MapNameTo, because name \"" << name << "\" was not found" << std::endl;
    throw -1;
  }

  return static_cast<T>(idx);
}


namespace Detail
{
    double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
    {
        return curr == prev
            ? curr
            : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
    }
}

namespace Detail{

  template <typename... Args>
  struct type_list {};



  // this does the implementation
  template<typename TupleType, int I, template<typename> class Fctor>
  struct for_each_tuple_type_struct {
    template <typename... Args>
    static void get(Args && ... args) {
      //Call get() of Fctor
      Fctor<std::tuple_element_t<I,TupleType>>::get(std::forward<Args>(args)...);
      
      //Continue Loop
      for_each_tuple_type_struct<TupleType, I - 1, Fctor>::get(std::forward<Args>(args)...);
    }
  };

  template<typename TupleType, template<typename> class Fctor>
  struct for_each_tuple_type_struct<TupleType, 0, Fctor> {
    template <typename... Args>
    static void get(Args && ... args) {
      // //Call get() of Fctor
      Fctor<std::tuple_element_t<0,TupleType>>::get(std::forward<Args>(args)...);
    }
  };


  // Fctor is a templated class with a single template parameter. 
  // it is expected to contain a static function get() that accepts the same 
  // number of arguments as are passed into this function
  //
  // Fctor<T>::get(Args...) is called for each type T in the tuple 
  template<typename TupleType, template<typename> class Fctor, typename... Args>
  void for_each_tuple_type(Args && ... args) {
    for_each_tuple_type_struct<TupleType, std::tuple_size<TupleType>::value-1, Fctor>::get(std::forward<Args>(args)...);
  };




  template <typename T, typename... Ts>
  struct FirstType{
    typedef T type;
  };

  template <typename T>
  struct FirstType<T>{
    typedef T type;
  };


  // this does the implementation
  template<typename AggregateTuple, int I, template<typename...> class Fctor, typename TupleType1, typename... TupleTypes>
  struct nested_for_each_tuple_type_struct {
    template <typename... Args>
    static void get(Args && ... args) {
      // set the Ith type
      typedef std::tuple_element_t<I, TupleType1>   Type1;
      typedef std::tuple<Type1>                   SingletonTuple1;
      typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

      // loop through inner loop
      nested_for_each_tuple_type_struct<NewAggregate, 
                                        std::tuple_size<typename FirstType<TupleTypes...>::type>::value-1, 
                                        Fctor,
                                        TupleTypes...>::get(std::forward<Args>(args)...);

      // on to the next iterate
      nested_for_each_tuple_type_struct<AggregateTuple, 
                                        I - 1, 
                                        Fctor, 
                                        TupleType1, 
                                        TupleTypes...>::get(std::forward<Args>(args)...);
    }
  };

  // base case
  template<typename AggregateTuple, template<typename...> class Fctor, typename TupleType1, typename... TupleTypes>
  struct nested_for_each_tuple_type_struct<AggregateTuple, 0, Fctor, TupleType1, TupleTypes...> {
    template <typename... Args>
    static void get(Args && ... args) {
      // set the Ith type
      typedef std::tuple_element_t<0, TupleType1>   Type1;
      typedef std::tuple<Type1>                   SingletonTuple1;
      typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

      // loop through inner loop
      nested_for_each_tuple_type_struct<NewAggregate, 
                                        std::tuple_size<typename FirstType<TupleTypes...>::type>::value-1, 
                                        Fctor,
                                        TupleTypes...>::get(std::forward<Args>(args)...);
    }
  };

  // this does the implementation
  template<typename AggregateTuple, int I, template<typename...> class Fctor, typename TupleType>
  struct nested_for_each_tuple_type_struct<AggregateTuple, I, Fctor, TupleType>{
  private:
    template <typename Tuple>
    struct ApplyFunctor{
    };

    // specialize
    template <typename... TupleArgs>
    struct ApplyFunctor<std::tuple<TupleArgs...>>{
      template <typename... Args>
      static void get(Args && ... args){
        Fctor<TupleArgs...>::get(std::forward<Args>(args)...);
      }
    };
  public:
    template <typename... Args>
    static void get(Args && ... args) {

      // set the Ith type
      typedef std::tuple_element_t<I, TupleType>   Type1;
      typedef std::tuple<Type1>                   SingletonTuple1;
      typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

      //Call get() of Fctor
      ApplyFunctor<NewAggregate>::get(std::forward<Args>(args)...);
      
      //Continue to next iteration
      nested_for_each_tuple_type_struct<AggregateTuple, I - 1, Fctor, TupleType>::get(std::forward<Args>(args)...);
    }
  };


  // this does the implementation
  template<typename AggregateTuple, template<typename...> class Fctor, typename TupleType>
  struct nested_for_each_tuple_type_struct<AggregateTuple, 0, Fctor, TupleType>{
  private:
    template <typename Tuple>
    struct ApplyFunctor{
    };

    // specialize
    template <typename... TupleArgs>
    struct ApplyFunctor<std::tuple<TupleArgs...>>{
      template <typename... Args>
      static void get(Args && ... args){
        Fctor<TupleArgs...>::get(std::forward<Args>(args)...);
      }
    };
  public:
    template <typename... Args>
    static void get(Args && ... args) {

      // set the Ith type
      typedef std::tuple_element_t<0, TupleType>   Type1;
      typedef std::tuple<Type1>                   SingletonTuple1;
      typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

      //Call get() of Fctor
      ApplyFunctor<NewAggregate>::get(std::forward<Args>(args)...);
    }
  };



  // template<typename AggregateTuple, typename TupleType, template<typename> class Fctor>
  // struct nested_for_each_tuple_type_struct<TupleType, 0, Fctor> {
  //   template <typename... Args>
  //   static void get(Args && ... args) {
  //     // //Call get() of Fctor
  //     Fctor<std::tuple_element_t<0,TupleType>>::get(std::forward<Args>(args)...);
  //   }
  // };


  // Fctor is a templated class with a single template parameter. 
  // it is expected to contain a static function get() that accepts the same 
  // number of arguments as are passed into this function
  //
  // Fctor<T>::get(Args...) is called for each type T in the tuple 
  template<template<typename...> class Fctor, typename TupleType1, typename... TupleTypes, typename... Args>
  void nested_for_each_tuple_type(Args && ... args) {
    nested_for_each_tuple_type_struct<std::tuple<>, std::tuple_size<TupleType1>::value-1, Fctor, TupleType1, TupleTypes...>::get(std::forward<Args>(args)...);
  };

  

  // check if a tuple has a certain type T in its template parameters
  template <typename T, typename Tuple>
  struct has_type;

  template <typename T>
  struct has_type<T, std::tuple<>> : std::false_type {};

  template <typename T, typename U, typename... Ts>
  struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {};

  template <typename T, typename... Ts>
  struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

  template <typename T, typename Tuple>
  using tuple_contains_type = typename has_type<T, Tuple>::type;
}

/*
* Constexpr version of the square root
* Return value:
*   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
*   - Otherwise, returns NaN
*/
double constexpr sqrt(double x)
{
    return x >= 0 && x < std::numeric_limits<double>::infinity()
        ? Detail::sqrtNewtonRaphson(x, x, 0)
        : std::numeric_limits<double>::quiet_NaN();
}


static constexpr double pi = 3.14159265358979323846264338327950288; 
static constexpr double eps0 = 8.854e-12; /**< permittivity of free space */
static constexpr double mu0 = pi*4.0e-7; /**< permeability of free space */
static constexpr double c0 = 2.99792458e+8;	/**< speed of light in a vacuum */
static constexpr double imp0 = sqrt(mu0/eps0); /**< free space impedance */


struct EMMode{};
struct TE : public EMMode{static const std::size_t dim=2;
                          static const std::size_t numE=2;
                          static const std::size_t numH=1;};
struct TM : public EMMode{static const std::size_t dim=2;
                          static const std::size_t numE=1;
                          static const std::size_t numH=2;};
struct TEM : public EMMode{static const std::size_t dim=1;
                           static const std::size_t numE=1;
                           static const std::size_t numH=1;};
struct ThreeD : public EMMode{static const std::size_t dim=3;
                              static const std::size_t numE=3;
                              static const std::size_t numH=3;};





enum class FieldType : char{
  Electric = 0,
  Magnetic,
  NONE
};
template <> struct NameArray<FieldType>{
  static constexpr std::array<const char *, 3> value = {"Electric", "Magnetic", "NONE"};
};
constexpr std::array<const char *, 3> NameArray<FieldType>::value;



enum class Dir : char{
  X = 0,
  Y,
  Z,
  NONE
};
template <> struct NameArray<Dir>{
  static constexpr std::array<const char *, 4> value = {"X", "Y", "Z", "NONE"};
};
constexpr std::array<const char *, 4> NameArray<Dir>::value;



template <Dir I, Dir J>
struct MutuallyOrthogonal{static constexpr Dir value = Dir::NONE;};

template <> struct MutuallyOrthogonal<Dir::X, Dir::Y>{static constexpr Dir value = Dir::Z;};
template <> struct MutuallyOrthogonal<Dir::Y, Dir::X>{static constexpr Dir value = Dir::Z;};
template <> struct MutuallyOrthogonal<Dir::X, Dir::Z>{static constexpr Dir value = Dir::Y;};
template <> struct MutuallyOrthogonal<Dir::Z, Dir::X>{static constexpr Dir value = Dir::Y;};
template <> struct MutuallyOrthogonal<Dir::Z, Dir::Y>{static constexpr Dir value = Dir::X;};
template <> struct MutuallyOrthogonal<Dir::Y, Dir::Z>{static constexpr Dir value = Dir::X;};


template <Dir Idx1, Dir Idx2, Dir... Idxs>
struct CrossAll{
  static constexpr Dir value = MutuallyOrthogonal<Idx1, MutuallyOrthogonal<Idx2, Idxs...>::value>::value;
};

template <Dir Idx1, Dir Idx2>
struct CrossAll<Idx1, Idx2>{
  static constexpr Dir value = MutuallyOrthogonal<Idx1, Idx2>::value;
};

//************************************************************
//************************************************************


// recursively check that each type is distinct from all of its following types
template <typename T, typename T1, typename... TypeList>
struct IsListed{
  static constexpr bool value = std::is_same<T, T1>::value || IsListed<T, TypeList...>::value;
};

template <typename T, typename T1>
struct IsListed<T, T1>{
  static constexpr bool value = std::is_same<T, T1>::value;
};



// recursively check that each direction is distinct from all of its following directions
template <Dir T, Dir T1, Dir... TypeList>
struct DirIsListed{
  static constexpr bool value = T==T1 || DirIsListed<T, TypeList...>::value;
};

template <Dir T, Dir T1>
struct DirIsListed<T, T1>{
  static constexpr bool value = T==T1;
};



// recursively check that each direction is distinct from all of its following directions
template <Dir T, Dir T1, Dir... TypeList>
struct ContainsDuplicates{
  static constexpr bool value = DirIsListed<T,T1,TypeList...>::value || ContainsDuplicates<T1, TypeList...>::value;
};

template <Dir T, Dir T1>
struct ContainsDuplicates<T, T1>{
  static constexpr bool value = DirIsListed<T,T1>::value;
};


// // generalized Levi-Civita
// template <Dir... Idx>
// struct GeneralizedLeviCivita{
// private:

// public:
//  // check that each direction passed in is distinct... if not then value = 0
//   static constexpr int value = (ContainsDuplicates<Idx...>::value ? 0 : 1);
// };


// Define components of the Levi-Civita tensor
template <Dir I, Dir J, Dir K>
struct LeviCivita{
  // Default, return 0
  static constexpr int value = 0.0;
  static constexpr decltype(auto) get(){return 0.0;};
};

template <> struct LeviCivita<Dir::X, Dir::Y, Dir::Z>{
  static constexpr int value = 1.0;
  static constexpr decltype(auto) get(){return 1.0;};
};

template <> struct LeviCivita<Dir::Y, Dir::Z, Dir::X>{
  static constexpr int value = 1.0;
  static constexpr decltype(auto) get(){return 1.0;};
};

template <> struct LeviCivita<Dir::Z, Dir::X, Dir::Y>{
  static constexpr int value = 1.0;
  static constexpr decltype(auto) get(){return 1.0;};
};

template <> struct LeviCivita<Dir::Z, Dir::Y, Dir::X>{
  static constexpr int value = -1.0;
  static constexpr decltype(auto) get(){return -1.0;};
};

template <> struct LeviCivita<Dir::X, Dir::Z, Dir::Y>{
  static constexpr int value = -1.0;
  static constexpr decltype(auto) get(){return -1.0;};
};

template <> struct LeviCivita<Dir::Y, Dir::X, Dir::Z>{
  static constexpr int value = -1.0;
  static constexpr decltype(auto) get(){return -1.0;};
};


//************************************************************
//************************************************************
//************************************************************
//************************************************************


enum class Orientation : char{
  MIN = 0,
  MAX,
  NONE
};
template <> struct NameArray<Orientation>{
  static constexpr std::array<const char *, 3> value = {"MIN", "MAX", "NONE"};
};
constexpr std::array<const char *, 3> NameArray<Orientation>::value;


struct Field{typedef std::array<float, 3> offset_type;};

  struct EType : public Field{
    static const Orientation neighb_side = Orientation::MAX;
    static const FieldType field_type = FieldType::Electric;
  };
    struct Ex : public EType{
      static constexpr offset_type    off = {0.0, -0.5, -0.5};
      static constexpr const char *   name = "Ex";
      static constexpr Dir            direction = Dir::X;
  };
    struct Ey : public EType{
      static constexpr offset_type    off = {-0.5, 0.0, -0.5};
      static constexpr const char *   name = "Ey";
      static constexpr Dir            direction = Dir::Y;
  };
    struct Ez : public EType{
      static constexpr offset_type    off = {-0.5, -0.5, 0.0};
      static constexpr const char *   name = "Ez";
      static constexpr Dir            direction = Dir::Z;
  };
    struct Dx : public EType{
      static constexpr offset_type    off = Ex::off;
      static constexpr const char *   name = "Dx";
      static constexpr Dir            direction = Dir::X;
    };
    struct Dy : public EType{
      static constexpr offset_type    off = Ey::off;
      static constexpr const char *   name = "Dy";
      static constexpr Dir            direction = Dir::Y;
    };
    struct Dz : public EType{
      static constexpr offset_type    off = Ez::off;
      static constexpr const char *   name = "Dz";
      static constexpr Dir            direction = Dir::Z;
    };

  struct HType : public Field{
    static const Orientation neighb_side = Orientation::MIN;
    static const FieldType field_type = FieldType::Magnetic;
  };
    struct Hx : public HType{
      static constexpr offset_type    off = {-0.5, 0.0, 0.0};
      static constexpr const char *   name = "Hx";
      static constexpr Dir            direction = Dir::X;
  };
    struct Hy : public HType{
      static constexpr offset_type    off = {0.0, -0.5, 0.0};
      static constexpr const char *   name = "Hy";
      static constexpr Dir            direction = Dir::Y;
  };
    struct Hz : public HType{
      static constexpr offset_type    off = {0.0, 0.0, -0.5};
      static constexpr const char *   name = "Hz";
      static constexpr Dir            direction = Dir::Z;
  };
    struct Bx : public HType{
      static constexpr offset_type    off = Hx::off;
      static constexpr const char *   name = "Bx";
      static constexpr Dir            direction = Dir::X;
    };
    struct By : public HType{
      static constexpr offset_type    off = Hy::off;
      static constexpr const char *   name = "By";
      static constexpr Dir            direction = Dir::Y;
    };
    struct Bz : public HType{
      static constexpr offset_type    off = Hz::off;
      static constexpr const char *   name = "Bz";
      static constexpr Dir            direction = Dir::Z;
    };


// have to put this declaration here b/c c++14 and lower
// have some issues with static constexpr...
// it is fixed in c++17 but that isn't implemented here yet
constexpr Field::offset_type Ex::off;
constexpr Field::offset_type Ey::off;
constexpr Field::offset_type Ez::off;

constexpr Field::offset_type Dx::off;
constexpr Field::offset_type Dy::off;
constexpr Field::offset_type Dz::off;

constexpr Field::offset_type Hx::off;
constexpr Field::offset_type Hy::off;
constexpr Field::offset_type Hz::off;

constexpr Field::offset_type Bx::off;
constexpr Field::offset_type By::off;
constexpr Field::offset_type Bz::off;


template <FieldType ft, Dir d>
struct FieldComponent{};

template <> struct FieldComponent<FieldType::Electric, Dir::X>{typedef Ex type;};
template <> struct FieldComponent<FieldType::Electric, Dir::Y>{typedef Ey type;};
template <> struct FieldComponent<FieldType::Electric, Dir::Z>{typedef Ez type;};

template <> struct FieldComponent<FieldType::Magnetic, Dir::X>{typedef Hx type;};
template <> struct FieldComponent<FieldType::Magnetic, Dir::Y>{typedef Hy type;};
template <> struct FieldComponent<FieldType::Magnetic, Dir::Z>{typedef Hz type;};


template <FieldType ft, Dir d>
struct FluxComponent{};

template <> struct FluxComponent<FieldType::Electric, Dir::X>{typedef Dx type;};
template <> struct FluxComponent<FieldType::Electric, Dir::Y>{typedef Dy type;};
template <> struct FluxComponent<FieldType::Electric, Dir::Z>{typedef Dz type;};

template <> struct FluxComponent<FieldType::Magnetic, Dir::X>{typedef Bx type;};
template <> struct FluxComponent<FieldType::Magnetic, Dir::Y>{typedef By type;};
template <> struct FluxComponent<FieldType::Magnetic, Dir::Z>{typedef Bz type;};



// helper structs to easily access field data
template <typename EMField>
struct FieldDir{
  static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
  static constexpr Dir value = EMField::direction;
};

template <typename EMField>
struct IsElectric{
  static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
  static constexpr bool value = (EMField::field_type == FieldType::Electric);
};

template <typename EMField>
struct IsMagnetic{
  static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
  static constexpr bool value = (EMField::field_type == FieldType::Magnetic);
};






// define field components for different modes
template <typename Mode>
struct FieldComponents{
  static_assert(std::is_base_of<EMMode, Mode>::value, "Mode must be a valid EMMode");
};

template <> struct FieldComponents<ThreeD>{
  typedef std::tuple<Ex, Ey, Ez>          electric;
  typedef std::tuple<Hx, Hy, Hz>          magnetic;

  typedef std::tuple<Dx, Dy, Dz>          electric_flux;
  typedef std::tuple<Bx, By, Bz>          magnetic_flux;
};

template <> struct FieldComponents<TE>{
  typedef std::tuple<Ex, Ey>   electric;
  typedef std::tuple<Hz>       magnetic;

  typedef std::tuple<Dx, Dy>   electric_flux;
  typedef std::tuple<Bz>       magnetic_flux;
};

template <> struct FieldComponents<TM>{
  typedef std::tuple<Ez>       electric;
  typedef std::tuple<Hx, Hy>   magnetic;

  typedef std::tuple<Dz>       electric_flux;
  typedef std::tuple<Bx, By>   magnetic_flux;
};

template <> struct FieldComponents<TEM>{
  typedef std::tuple<Ez>       electric;
  typedef std::tuple<Hy>       magnetic;

  typedef std::tuple<Dz>       electric_flux;
  typedef std::tuple<By>       magnetic_flux;
};





}// end namespace fdtd

#endif