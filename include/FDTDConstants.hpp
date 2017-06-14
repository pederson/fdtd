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
struct TE : public EMMode{static const std::size_t dim=2;};
struct TM : public EMMode{static const std::size_t dim=2;};
struct TEM : public EMMode{static const std::size_t dim=1;};
struct ThreeD : public EMMode{static const std::size_t dim=3;};



}// end namespace fdtd

#endif