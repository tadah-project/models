#ifndef CUTOFFS_H
#define CUTOFFS_H

#include <string>
#include <tadah/core/config.h>
#include <tadah/core/registry.h>

/** Base class to be inherited by all cutoffs */
class Cut_Base {
    public:
        virtual ~Cut_Base();
        virtual std::string label()=0;
        virtual double calc(double r)=0;
        virtual double get_rcut()=0;
        virtual void set_rcut(const double r)=0;
        virtual double get_rcut_sq()=0;
        virtual double calc_prime(double r)=0;
        void test_rcut(const double r);
};

/**
 * @class Cut_Dummy
 * @brief Represents a basic cutoff function with a sharp transition at the cutoff radius.
 *
 * The `Cut_Dummy` class implements a simple cutoff function defined by:
 * \f[
 * f_c(r) =
 * \begin{cases}
 * 1, & \text{if } r \leq r_{\text{cut}} \\
 * 0, & \text{if } r > r_{\text{cut}}
 * \end{cases}
 * \f]
 *
 * and its derivative:
 * \f[
 * f_c'(r) = 0
 * \f]
 *
 * where:
 * - \( r \) is the radial distance.
 * - \( r_{\text{cut}} \) is the cutoff radius.
 *
 * **Characteristics:**
 * - The function value is constant (1) within the cutoff radius and zero beyond it.
 * - The derivative of the function is zero everywhere except at \( r = r_{\text{cut}} \), where it is undefined due to the discontinuity.
 *
 * **Note:** Since the function is discontinuous at \( r = r_{\text{cut}} \), its derivative is not defined at that point. In practical computations, the derivative function `calc_prime` returns `0.0` for all \( r \).
 */
class Cut_Dummy : public Cut_Base {
    private:
        std::string lab = "Cut_Dummy";
        double rcut, rcut_sq, rcut_inv;
    public:
        Cut_Dummy();
        Cut_Dummy(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};

/**
 * @class Cut_Cos
 * @brief Cosine cutoff function.
 *
 * The `Cut_Cos` class implements a smooth cosine cutoff function defined by:
 * \f[
 * f_c(r) =
 * \begin{cases}
 *  \dfrac{1}{2}\left[ \cos\left( \dfrac{\pi r}{r_c} \right) + 1 \right], & \text{if } r \leq r_c \\
 *  0, & \text{if } r > r_c
 * \end{cases}
 * \f]
 *
 * This function smoothly transitions from 1 to 0 over the interval \( 0 \leq r \leq r_c \).
 * It is commonly used in molecular simulations to smoothly truncate interactions without introducing discontinuities in the potential energy or its derivatives.
 *
 * **Reference:**
 *
 * Behler, J., & Parrinello, M. (2007). *Generalized neural-network representation of high-dimensional potential-energy surfaces*. Physical Review Letters, 98(14), 146401. [DOI:10.1103/PhysRevLett.98.146401](https://doi.org/10.1103/PhysRevLett.98.146401)
 */
class Cut_Cos : public Cut_Base {
    private:
        std::string lab = "Cut_Cos";
        double rcut, rcut_sq, rcut_inv;
    public:
        Cut_Cos();
        Cut_Cos(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};

/**
 * @class Cut_Tanh
 * @brief Hyperbolic tangent cutoff function.
 *
 * The `Cut_Tanh` class implements a smooth cutoff function using the hyperbolic tangent, defined by:
 * \f[
 * f_c(r) =
 * \begin{cases}
 *  \tanh^3\left( 1 - \dfrac{r}{r_c} \right), & \text{if } r \leq r_c \\
 *  0, & \text{if } r > r_c
 * \end{cases}
 * \f]
 *
 * This function smoothly transitions from 1 to 0 over the interval \( 0 \leq r \leq r_c \), with the transition shape controlled by the cubic power of the hyperbolic tangent. It is particularly useful in constructing atom-centered symmetry functions for high-dimensional neural network potentials.
 *
 * **Reference:**
 *
 * Behler, J. (2011). Atom-centered symmetry functions for constructing high-dimensional neural network potentials. *Journal of Chemical Physics*, 134(7), 074106. [DOI:10.1063/1.3553717](https://doi.org/10.1063/1.3553717)
 */
class Cut_Tanh : public Cut_Base {
    private:
        std::string lab = "Cut_Tanh";
        double rcut, rcut_sq, rcut_inv;
    public:
        Cut_Tanh();
        Cut_Tanh(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};

/**
 * @class Cut_Poly2
 * @brief Polynomial cutoff function ensuring smoothness up to the second derivative.
 *
 * The `Cut_Poly2` class implements a polynomial cutoff function defined as:

 * \f[
 * f_c(r) =
 * \begin{cases}
 * 1, & \text{if } r \leq r_{\text{inner}} \\
 * \left( r - r_c + 1 \right )^3 \left[ r \left( 15 - 6 r \right ) - 10 \right ] + 1, & \text{if } r_{\text{inner}} < r \leq r_c \\
 * 0, & \text{if } r > r_c
 * \end{cases}
 * \f]

 * where:
 * - \( r \) is the radial distance.
 * - \( r_c \) is the cutoff radius.
 * - \( r_{\text{inner}} = r_c - 1 \) is the inner cutoff radius.

 * This function provides a smooth transition from 1 to 0 between \( r_{\text{inner}} \) and \( r_c \), ensuring continuity and smoothness in the potential and its derivatives up to the second order.

 * **Reference:**

 * Singraber, A., Behler, J., & Dellago, C. (2019). *Library-Based LAMMPS Implementation of High-Dimensional Neural Network Potentials*. Journal of Chemical Theory and Computation, 15(1), 182–190. [DOI:10.1021/acs.jctc.8b00770](https://doi.org/10.1021/acs.jctc.8b00770)
 */
class Cut_Poly2 : public Cut_Base {
    private:
        std::string lab = "Cut_Poly2";
        double rcut, rcut_sq, rcut_inv;
        double rcut_inner;
    public:
        Cut_Poly2();
        Cut_Poly2(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};
/**
 * @class Cut_Poly1
 * @brief Represents a 5th order polynomial cutoff function with zero first and second derivatives at the cutoff radius.
 *
 * This class implements a smooth cutoff function defined by:
 * \f[
 * f(r) =
 * \begin{cases}
 * 1 - 10 \left( \dfrac{r}{r_{\text{cut}}} \right)^3 + 15 \left( \dfrac{r}{r_{\text{cut}}} \right)^4 - 6 \left( \dfrac{r}{r_{\text{cut}}} \right)^5, & \text{for } r \leq r_{\text{cut}} \\
 * 0, & \text{for } r > r_{\text{cut}}
 * \end{cases}
 * \f]
 * and its derivative:
 * \f[
 * f'(r) =
 * \begin{cases}
 * \left( -30 \left( \dfrac{r}{r_{\text{cut}}} \right)^2 + 60 \left( \dfrac{r}{r_{\text{cut}}} \right)^3 - 30 \left( \dfrac{r}{r_{\text{cut}}} \right)^4 \right) \dfrac{1}{r_{\text{cut}}}, & \text{for } r \leq r_{\text{cut}} \\
 * 0, & \text{for } r > r_{\text{cut}}
 * \end{cases}
 * \f]
 * where:
 * - \( r \) is the radial distance.
 * - \( r_{\text{cut}} \) is the cutoff radius.
 *
 * The function smoothly transitions from 1 to 0 within the cutoff radius \( r_{\text{cut}} \), ensuring that the function and its first and second derivatives are zero at \( r \geq r_{\text{cut}} \).
 *
 */
class Cut_Poly1 : public Cut_Base {
    private:
        std::string lab = "Cut_Poly1";
        double rcut, rcut_sq, rcut_inv;
    public:
        Cut_Poly1();
        Cut_Poly1(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};
class Cut_Cos_S : public Cut_Base {
    private:
        std::string lab = "Cut_Cos_S";
        double rcut, rcut_sq, rcut_inv;
        double rcut_inner;
    public:
        Cut_Cos_S();
        Cut_Cos_S(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};
/**
 * @class Cut_PT
 * @brief Represents a smooth cutoff function using hyperbolic tangent smoothing.
 *
 * This class implements a smooth cutoff function defined by:
 * \f[
 * f(r) = \frac{1}{2} + \frac{1}{2} \tanh\left( \frac{1}{2} \left( \frac{b}{r - a} + \frac{b}{r - R} \right) \right)
 * \f]
 * and its derivative:
 * \f[
 * f'(r) = -\frac{b}{4} \left( \frac{1}{\cosh\left( \dfrac{b}{2} \left( \dfrac{1}{r - a} + \dfrac{1}{r - R} \right) \right)} \right)^2 \left( \dfrac{1}{(r - a)^2} + \dfrac{1}{(r - R)^2} \right)
 * \f]
 * where:
 * - \( r \) is the radial distance.
 * - \( a \) is the inner cutoff radius (`rcut_inner`).
 * - \( R \) is the outer cutoff radius (`rcut`).
 * - \( b \) is the smoothing parameter.
 *
 * The function smoothly transitions from 0 to 1 between the inner cutoff radius \( a \) and the outer cutoff radius \( R \).
 *
 */
class Cut_PT : public Cut_Base {
    private:
        std::string lab="Cut_PT";
        double rcut;
        double rcut_sq;
        double rcut_inv;
        double rcut_inner;

    public:
        Cut_PT();
        Cut_PT(double rcut_value);
        std::string label();
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};

//template<> inline Registry<Cut_Base,double>::Map Registry<Cut_Base,double>::registry{};
#endif
