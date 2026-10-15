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

        explicit DC_Selector(const tadah::core::Context& c) {
          try {init(c);}catch(...){release();throw;}
        }
        explicit DC_Selector(tadah::core::Context& c) {
          try {init(c);}catch(...){release();throw;}
        }

        void init(tadah::core::Context &context) {
          resolve_embed_blocks(context);
          init(static_cast<const tadah::core::Context &>(context));
        }

        void init(const tadah::core::Context &input) {
          const bool component_blocks = !input.blocks("DESCRIPTOR").empty();
          tadah::core::Context resolved;
          if (component_blocks) { resolved=input; resolve_embed_blocks(resolved); }
          const auto &context=component_blocks?resolved:input;
          double rcutzero = 0.0; // for dummies

          // EMBF0 / RHOFORM without a DM_REAM block have no effect: one
          // stderr notice per key and process, wherever a model is built.
          report_ream_only_keys(context);

          bias = static_cast<int>(tadah::core::n_intercept_columns(context));

          if (context.get<bool>("INIT2B")) {
            double rcut2b = context.get<double>("RCUT2BMAX");
            c2b = tadah::core::factory<Cut_Base,double&>( context.get<std::string>("RCTYPE2B"), rcut2b );
            d2b = tadah::core::factory<D2_Base,const tadah::core::Context&>( context.get<std::string>("TYPE2B"), context );
            d2b->warn_if_discontinuous_noncompact();
          }
          else {
            c2b = tadah::core::factory<Cut_Base,double&>( "Cut_Dummy", rcutzero );
            d2b = tadah::core::factory<D2_Base,const tadah::core::Context&>( "D2_Dummy", context );
          }

          if (context.get<bool>("INITMB")) {
            double rcutmb = context.get<double>("RCUTMBMAX");
            if (component_blocks) {
              auto children=embed_contexts(context);
              cmb = new Cut_Dummy(rcutmb);
              if(children.size()==1)
                dmb=tadah::core::factory<DM_Base,const tadah::core::Context&>(children[0].get<std::string>("TYPEMB"),children[0]);
              else dmb=new DM_mJoin(context);
            } else {
              cmb = tadah::core::factory<Cut_Base,double&>( context.get<std::string>("RCTYPEMB"), rcutmb );
              dmb = tadah::core::factory<DM_Base,const tadah::core::Context&>( context.get<std::string>("TYPEMB"), context );
            }
            dmb->warn_if_discontinuous_noncompact();
          }
          else {
            cmb = tadah::core::factory<Cut_Base,double&>( "Cut_Dummy", rcutzero );
            dmb = tadah::core::factory<DM_Base,const tadah::core::Context&>( "DM_Dummy", context );
          }

          // The blocks' first indices, the size of the descriptor and
          // density arrays
          place();
        }

        size_t descriptor_size() const { return dsize; }
        size_t density_size() const { return rsize; }

        /** Update descriptor parameters in-place from context.
         *  Returns true if the descriptor or the density array changed size
         *  or a descriptor reports a change (a Trainer must be rebuilt).
         *  The descriptors only (D_Base::reinit()): tadah hpo does not call
         *  it, it builds a Trainer for every evaluation. Afterwards the
         *  selector is laid out as init() lays out one built from @p ctx:
         *  the intercept columns of its BIAS and ATOM, the two-body block
         *  from the last of them, the many-body block after it, the density
         *  array with rho_i from 0 and F'(rho_i) right after it, and dsize
         *  and rsize of those sizes.
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
          if (d2b && dmb) {
            const size_t old_dsize = dsize, old_rsize = rsize;
            bias = static_cast<int>(tadah::core::n_intercept_columns(ctx));
            place();
            changed |= dsize != old_dsize || rsize != old_rsize;
          }
          return changed;
        }

        ~DC_Selector()
        {
          release();
        }

      private:
        void release() {
          if (d2b)
            delete d2b;
          if (c2b)
            delete c2b;
          if (cmb)
            delete cmb;
          if (dmb)
            delete dmb;
          d2b=nullptr;c2b=nullptr;cmb=nullptr;dmb=nullptr;
        }

      private:
        /* The layout of init() and reinit(), for `bias` intercept
         * columns: the two-body block from column bias, the many-body
         * block after it, the density array with rho_i from 0 and
         * F'(rho_i) right after it (where each descriptor's constructor
         * puts them, and where pair_style tadah exchanges them); dsize and
         * rsize of those sizes. */
        void place() {
          d2b->set_fidx(static_cast<size_t>(bias));
          dmb->set_fidx(static_cast<size_t>(bias) + d2b->size());
          dmb->set_rfidx(0);
          dmb->set_drfidx(dmb->rhoi_size());
          dsize = static_cast<size_t>(bias) + d2b->size() + dmb->size();
          rsize = dmb->rhoi_size() + dmb->rhoip_size();
        }

    };
  } // namespace models
} // namespace tadah
#endif
