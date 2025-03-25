#ifndef DC_SELECTOR_H
#define DC_SELECTOR_H

#include <tadah/core/context.h>
#include <tadah/models/descriptors/d_all.h>
#include <tadah/models/cut_all.h>
#include <tadah/core/registry.h>

namespace tadah {
namespace models {

class DC_Selector {
    public:
        Cut_Base *c2b=nullptr;
        Cut_Base *cmb=nullptr;
        D2_Base *d2b=nullptr;
        DM_Base *dmb=nullptr;

        //tadah::core::Context context;

        DC_Selector ()
        {};

        // Copy constructor
        DC_Selector(const DC_Selector& )
        {
            std::cout << "Copy constructor called\n";
        }
         // Const Assignment operator
         // Instead of deep copying
         // we use a trick where we just copy
         // the context file and run init()
         // as in a costructor
         // DC_Selector& operator=(const DC_Selector& dc)
         // {
         //     //context=dc.context;
         //     //init();
         //     return *this;
         // }
        // Assignment operator
        DC_Selector& operator=(DC_Selector& dc)
        {
            //std::swap(context,dc.context);
            std::swap(c2b,dc.c2b);
            std::swap(cmb,dc.cmb);
            std::swap(d2b,dc.d2b);
            std::swap(dmb,dc.dmb);
            return *this;
        }

        DC_Selector (tadah::core::Context &c)
            //context(c)
    {
        init(c);
    };
        void init(tadah::core::Context &context) {
            double rcutzero = 0.0; // for dummies

            size_t bias=0;
            if (context.get<bool>("BIAS"))
                bias++;

            if (context.get<bool>("INIT2B")) {
                double rcut2b = context.get<double>("RCUT2BMAX");
                c2b = tadah::core::factory<Cut_Base,double>( context.get<std::string>("RCTYPE2B"), rcut2b );
                d2b = tadah::core::factory<D2_Base,tadah::core::Context&>( context.get<std::string>("TYPE2B"), context );
                d2b->set_fidx(bias);
            }
            else {
                c2b = tadah::core::factory<Cut_Base,double>( "Cut_Dummy", rcutzero );
                d2b = tadah::core::factory<D2_Base,tadah::core::Context&>( "D2_Dummy", context );
                d2b->set_fidx(bias);
            }

            if (context.get<bool>("INITMB")) {
                double rcutmb = context.get<double>("RCUTMBMAX");
                cmb = tadah::core::factory<Cut_Base,double>( context.get<std::string>("RCTYPEMB"), rcutmb );
                dmb = tadah::core::factory<DM_Base,tadah::core::Context&>( context.get<std::string>("TYPEMB"), context );
                dmb->set_fidx(bias + d2b->size());
            }
            else {
                cmb = tadah::core::factory<Cut_Base,double>( "Cut_Dummy", rcutzero );
                dmb = tadah::core::factory<DM_Base,tadah::core::Context&>( "DM_Dummy", context );
                dmb->set_fidx(bias + d2b->size());
            }
        }
        ~DC_Selector()
        {
            if (d2b)
                delete d2b;
            if (c2b)
                delete c2b;
            if (cmb)
                delete cmb;
            if (dmb)

                delete dmb;
        }
};
} // namespace models
} // namespace tadah
#endif
