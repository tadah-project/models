#ifndef DM_func_H
#define DM_func_H

#include <tadah/core/config.h>
#include <tadah/models/descriptors/d_base.h>

#include <cstdlib>

class F_Base {
    public:
        virtual double calc_F(double rho, size_t c)=0;
        virtual double calc_dF(double rho, size_t c)=0;
};

/**
 * @class F_RLR
 * @brief Implements an embedding function of the form: \( s \cdot x \cdot \log(c \cdot x) \)
 *
 * This class supports embedding functions characterized by two main parameters:
 * - **SEMBFUNC**: Controls the depth of the embedding function.
 * - **CEMBFUNC**: Determines the x-intercept, with the x-intercept at \( 1/c \).
 * 
 * The number of keys for these parameters must match the entries in the mEAD descriptor.
 *
 * @note Ensure the configuration keys match the descriptor requirements.
 */
class F_RLR:  public F_Base {
    private:
        Config *config=nullptr;
        v_type sgrid;   // controls depth
        v_type cgrid;   // control x-intercep at 1/c
    public:
        F_RLR(Config &conf);
        F_RLR();
        double calc_F(double rho,size_t c);
        double calc_dF(double rho,size_t c);
        //~F_RLR() {
            // Do not delete config!
        //}
};
#endif
