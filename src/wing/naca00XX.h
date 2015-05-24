
/**
 * @file naca00XX.h
 * @brief Add description here
 */
#pragma once

namespace UVLM {
namespace wing {

/**
 * @brief Equation for a symmetrical 4-digit NACA airfoil
 * http://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_symmetrical_4-digit_NACA_airfoil
 * @param c chord length
 * @param xx number
 */
double NACA00XX(double x, double c, int xx);

}  // namespace wing
}  // namespace UVLM
