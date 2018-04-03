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

enum class Orientation : char{
    NONE,
  MIN = 0,
  MAX
};

struct Field{};

  struct EType : public Field{
    static const Orientation neighb_side = Orientation::MAX;
    static const FieldType field_type = FieldType::Electric;
  };
    struct Ex : public EType{};
    struct Ey : public EType{};
    struct Ez : public EType{};
    struct Dx : public EType{};
    struct Dy : public EType{};
    struct Dz : public EType{};

  struct HType : public Field{
    static const Orientation neighb_side = Orientation::MIN;
    static const FieldType field_type = FieldType::Magnetic;
  };
    struct Hx : public HType{};
    struct Hy : public HType{};
    struct Hz : public HType{};
    struct Bx : public HType{};
    struct By : public HType{};
    struct Bz : public HType{};







}// end namespace fdtd

#endif