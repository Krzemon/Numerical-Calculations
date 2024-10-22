#ifndef CALCULATIONS_H
#define CALCULATIONS_H

/**
 * @brief Generates a random number in the range [0, 1].
 * 
 * @return Generated random number.
 */
double uniform();

/**
 * @brief Trial wave function Psi_T(r, a, c).
 * 
 * @param r Argument of the function.
 * @param a Parameter of the function.
 * @param c Parameter of the function.
 * @return Value of Psi_T(r, a, c).
 */
double Psi_T(double r, double a, double c);

/**
 * @brief Local energy Eps_loc(r, a, c).
 * 
 * @param r Argument of the function.
 * @param a Parameter of the function.
 * @param c Parameter of the function.
 * @return Local energy Eps_loc(r, a, c).
 */
double Eps_loc(double r, double a, double c);

/**
 * @brief Computes the moment of the local energy Eps_loc(r, a, c).
 * 
 * @param a Parameter of the function.
 * @param c Parameter of the function.
 * @param N Number of iterations.
 * @param moment Moment computed for Eps_loc(r, a, c).
 * @param delta_r Increment in r for each iteration.
 * @return Computed moment of the local energy.
 */
double Eps_loc_mom(double a, double c, int N, int moment, double delta_r = 0.1);

#endif // CALCULATIONS_H
