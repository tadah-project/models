#ifndef BASIS_FUNCTIONS_H
#define BASIS_FUNCTIONS_H

#include <tadah/core/core_types.h>
#include <tadah/models/functions/function_base.h>

#include <iostream>

/** \brief Base class to be used by all basis functions */
struct BF_Base: public virtual Function_Base {
    virtual ~BF_Base();
    virtual double epredict(const t_type &, const aed_type2& ) const=0;
    virtual double fpredict(const t_type &, const fd_type &,
            const aed_type2& , const size_t ) const=0;
    virtual force_type fpredict(const t_type &, const fd_type &,
            const aed_type2& ) const=0;
};
#endif
