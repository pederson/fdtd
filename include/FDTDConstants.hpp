#ifndef _FDTDCONSTANTS_H
#define _FDTDCONSTANTS_H

#include <limits>
#include <cmath>
#include <array>
#include <memory>
#include <iostream>

namespace fdtd{


namespace Detail
{
    double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
    {
        return curr == prev
            ? curr
            : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
    }
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
  NONE,
  Electric,
  Magnetic
};

enum class Dir : char{
    NONE,
  X = 0,
  Y,
  Z
};

template <Dir I, Dir J>
struct MutuallyOrthogonal{};

template <> struct MutuallyOrthogonal<Dir::X, Dir::Y>{static constexpr Dir value = Dir::Z;};
template <> struct MutuallyOrthogonal<Dir::Y, Dir::X>{static constexpr Dir value = Dir::Z;};
template <> struct MutuallyOrthogonal<Dir::X, Dir::Z>{static constexpr Dir value = Dir::Y;};
template <> struct MutuallyOrthogonal<Dir::Z, Dir::X>{static constexpr Dir value = Dir::Y;};
template <> struct MutuallyOrthogonal<Dir::Z, Dir::Y>{static constexpr Dir value = Dir::X;};
template <> struct MutuallyOrthogonal<Dir::Y, Dir::Z>{static constexpr Dir value = Dir::X;};


enum class Orientation : char{
    NONE,
  MIN = 0,
  MAX
};

struct Field{typedef std::array<float, 3> offset_type;};

  struct EType : public Field{
    static const Orientation neighb_side = Orientation::MAX;
    static const FieldType field_type = FieldType::Electric;
  };
    struct Ex : public EType{static constexpr offset_type off = {0.0, -0.5, -0.5};};
    struct Ey : public EType{static constexpr offset_type off = {-0.5, 0.0, -0.5};};
    struct Ez : public EType{static constexpr offset_type off = {-0.5, -0.5, 0.0};};
    struct Dx : public EType{static constexpr offset_type off = Ex::off;};
    struct Dy : public EType{static constexpr offset_type off = Ey::off;};
    struct Dz : public EType{static constexpr offset_type off = Ez::off;};

  struct HType : public Field{
    static const Orientation neighb_side = Orientation::MIN;
    static const FieldType field_type = FieldType::Magnetic;
  };
    struct Hx : public HType{static constexpr offset_type off = {-0.5, 0.0, 0.0};};
    struct Hy : public HType{static constexpr offset_type off = {0.0, -0.5, 0.0};};
    struct Hz : public HType{static constexpr offset_type off = {0.0, 0.0, -0.5};};
    struct Bx : public HType{static constexpr offset_type off = Hx::off;};
    struct By : public HType{static constexpr offset_type off = Hy::off;};
    struct Bz : public HType{static constexpr offset_type off = Hz::off;};


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



template <typename EMField>
struct FieldDir{
  static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
};

template <> struct FieldDir<Dx>{static constexpr Dir value = Dir::X;};
template <> struct FieldDir<Dy>{static constexpr Dir value = Dir::Y;};
template <> struct FieldDir<Dz>{static constexpr Dir value = Dir::Z;};

template <> struct FieldDir<Ex>{static constexpr Dir value = Dir::X;};
template <> struct FieldDir<Ey>{static constexpr Dir value = Dir::Y;};
template <> struct FieldDir<Ez>{static constexpr Dir value = Dir::Z;};

template <> struct FieldDir<Bx>{static constexpr Dir value = Dir::X;};
template <> struct FieldDir<By>{static constexpr Dir value = Dir::Y;};
template <> struct FieldDir<Bz>{static constexpr Dir value = Dir::Z;};

template <> struct FieldDir<Hx>{static constexpr Dir value = Dir::X;};
template <> struct FieldDir<Hy>{static constexpr Dir value = Dir::Y;};
template <> struct FieldDir<Hz>{static constexpr Dir value = Dir::Z;};



template <typename EMField>
struct IsElectric{
  static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
  static constexpr bool value = false;
};

template <> struct IsElectric<Dx>{static constexpr bool value = true;};
template <> struct IsElectric<Dy>{static constexpr bool value = true;};
template <> struct IsElectric<Dz>{static constexpr bool value = true;};

template <> struct IsElectric<Ex>{static constexpr bool value = true;};
template <> struct IsElectric<Ey>{static constexpr bool value = true;};
template <> struct IsElectric<Ez>{static constexpr bool value = true;};



template <typename EMField>
struct IsMagnetic{
  static_assert(std::is_base_of<Field, EMField>::value, "Field must be a valid EMField");
  static constexpr bool value = false;
};

template <> struct IsMagnetic<Bx>{static constexpr bool value = true;};
template <> struct IsMagnetic<By>{static constexpr bool value = true;};
template <> struct IsMagnetic<Bz>{static constexpr bool value = true;};

template <> struct IsMagnetic<Hx>{static constexpr bool value = true;};
template <> struct IsMagnetic<Hy>{static constexpr bool value = true;};
template <> struct IsMagnetic<Hz>{static constexpr bool value = true;};

}// end namespace fdtd

#endif