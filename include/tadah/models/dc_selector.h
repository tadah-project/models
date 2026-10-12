#ifndef DC_SELECTOR_H
#define DC_SELECTOR_H

#include <tadah/core/context.h>
#include <tadah/core/intercept_layout.h>
#include <tadah/models/descriptors/d_all.h>
#include <tadah/models/cut_all.h>
#include <tadah/models/descriptors/dm/embf0.h>
#include <tadah/core/registry.h>

namespace tadah {
  namespace models {

    class DC_Selector {
      public:
        int bias = 0;
        size_t dsize =0;  // size of the descriptor (bias + 2b + mb)
        size_t rsize =0;  // size of the density array (rhoi + rhoip)
        Cut_Base *c2b=nullptr;
        Cut_Base *cmb=nullptr;
        D2_Base *d2b=nullptr;
        DM_Base *dmb=nullptr;

        DC_Selector ()
        {};

        // Copy constructor
        DC_Selector(const DC_Selector& ) {}

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

        explicit DC_Selector(const tadah::core::Context& c) { init(c); }

        void init(const tadah::core::Context &context) {
          double rcutzero = 0.0; // for dummies

          // EMBF0 / RHOFORM without a DM_REAM block have no effect: one
          // stderr notice per key and process, wherever a model is built.
          report_ream_only_keys(context);

          bias = static_cast<int>(tadah::core::n_intercept_columns(context));

          if (context.get<bool>("INIT2B")) {
            double rcut2b = context.get<double>("RCUT2BMAX");
            c2b = tadah::core::factory<Cut_Base,double&>( context.get<std::string>("RCTYPE2B"), rcut2b );
            d2b = tadah::core::factory<D2_Base,const tadah::core::Context&>( context.get<std::string>("TYPE2B"), context );
            d2b->set_fidx(bias);
            d2b->warn_if_discontinuous_noncompact();
          }
          else {
            c2b = tadah::core::factory<Cut_Base,double&>( "Cut_Dummy", rcutzero );
            d2b = tadah::core::factory<D2_Base,const tadah::core::Context&>( "D2_Dummy", context );
          }

          if (context.get<bool>("INITMB")) {
            double rcutmb = context.get<double>("RCUTMBMAX");
            cmb = tadah::core::factory<Cut_Base,double&>( context.get<std::string>("RCTYPEMB"), rcutmb );
            dmb = tadah::core::factory<DM_Base,const tadah::core::Context&>( context.get<std::string>("TYPEMB"), context );
            dmb->set_fidx(bias + d2b->size());
            dmb->warn_if_discontinuous_noncompact();
          }
          else {
            cmb = tadah::core::factory<Cut_Base,double&>( "Cut_Dummy", rcutzero );
            dmb = tadah::core::factory<DM_Base,const tadah::core::Context&>( "DM_Dummy", context );
          }

          // Set the size of the descriptor and density arrays
          dsize = bias + d2b->size() + dmb->size();
          rsize = dmb->rhoi_size()+ dmb->rhoip_size();
        }

        size_t descriptor_size() const { return dsize; }
        size_t density_size() const { return rsize; }

        /** Update descriptor parameters in-place from context.
         *  Returns true if descriptor size changed (Trainer must be rebuilt).
         *  The descriptors only (D_Base::reinit()): tadah hpo does not call
         *  it, it builds a Trainer for every evaluation. After a size change
         *  the many-body block starts after the two-body block, as init()
         *  places it (dmb's first index bias + d2b->size()).
         *
         *  The objects stay the same (a DescriptorsCalc, such as a
         *  Trainer's, holds references to them), so the update is in
         *  place and its exception guarantee is the basic one: nothing
         *  leaks, but a reinit that throws (a key it cannot read, a
         *  refused grid or WATOM) may leave d2b, its cutoff or a part of a
         *  descriptor holding the new values and the rest the old ones.
         *  Build a new selector (a new Trainer) from a context that reads
         *  instead of using one whose reinit threw. */
        bool reinit(const tadah::core::Context& ctx) {
          bool changed = false;
          if (d2b) {
            changed |= d2b->reinit(ctx);
            d2b->update_weights(ctx);
            if (c2b && ctx.get<bool>("INIT2B"))
              c2b->set_rcut(ctx.get<double>("RCUT2BMAX"));
          }
          if (dmb) {
            changed |= dmb->reinit(ctx);
            dmb->update_weights(ctx);
            if (cmb && ctx.get<bool>("INITMB"))
              cmb->set_rcut(ctx.get<double>("RCUTMBMAX"));
          }
          if (changed) {
            // The many-body values follow the two-body ones (init()).
            if (dmb) dmb->set_fidx(static_cast<size_t>(bias) + d2b->size());
            dsize = static_cast<size_t>(bias) + d2b->size() + dmb->size();
            rsize = dmb->rhoi_size() + dmb->rhoip_size();
          }
          return changed;
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
