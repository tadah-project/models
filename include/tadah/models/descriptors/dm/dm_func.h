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

/** Embedding function of the form:  sx*log(cx)
 *
 * REQUIRED KEYS:
 *
 * SEMBFUNC, CEMBFUNC
 * control depth and x-intercept of an embedding function respectively.
 * The number of keys must match those in the mEAD descriptor.
 *
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
