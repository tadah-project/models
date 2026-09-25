#ifndef TADAH_MODELS_EMBF0_H
#define TADAH_MODELS_EMBF0_H

#include <tadah/core/context.h>
#include <optional>
#include <string>

namespace tadah {
namespace models {

/** \brief DM_REAM embedding reference, selected by the EMBF0 key.
 *
 *  - Zero (default): the zero-referenced basis g~_m(rho) = g_m(rho) - g_m(0),
 *    so F(0) = 0 exactly for any weights.
 *  - Free: the raw basis g_m of tadah 1.5.0 and earlier.
 *
 *  Key-handling contract: DM_REAM is the only reader that queries EMBF0
 *  through Context::get (embf0_from_context), which materialises the
 *  blueprint default ZERO into the model Context so every DM_REAM potential is
 *  written with an explicit EMBF0 line. Every other reader uses
 *  stored_embf0() / context_has_ream(), which never insert defaults. */
enum class EmbF0Mode { Zero, Free };

/** "ZERO" | "FREE". */
const char *embf0_name(EmbF0Mode m);

/** Case-insensitive parse; throws std::runtime_error for anything other than
 *  ZERO or FREE. */
EmbF0Mode parse_embf0(const std::string &value);

/** Mode read through Context::get: an absent key materialises the blueprint
 *  default ZERO into @p ctx. Reserved for DM_REAM. */
EmbF0Mode embf0_from_context(const tadah::core::Context &ctx);

/** Mode stored in @p ctx, or nullopt when the key is absent. Never inserts
 *  a default; throws for an invalid stored value. */
std::optional<EmbF0Mode> stored_embf0(const tadah::core::Context &ctx);

/** True when @p ctx stores INITMB true and a TYPEMB token DM_REAM (directly or
 *  inside a DM_mJoin token list). Never inserts defaults. */
bool context_has_ream(const tadah::core::Context &ctx);

/** True when the stored MODEL is linear in the descriptor (BF_Linear or
 *  Kern_Linear); false when MODEL is absent or non-linear. Never inserts. */
bool model_is_linear(const tadah::core::Context &ctx);

/** Potential writers: a context with a DM_REAM block gets an explicit,
 *  canonical EMBF0 line (ZERO when absent); a context without one has a
 *  stray EMBF0 removed (WARNING once per process), so potentials without a
 *  DM_REAM block never carry the key. */
void normalise_embf0_key(tadah::core::Context &c);

/** Engines that merge a task context over a loaded potential (predict,
 *  analysis): call with the potential as read and the task context BEFORE
 *  merging. Different explicit values -> std::runtime_error (EMBF0 is a
 *  property of the potential); a task value for a potential without the key
 *  -> INFO (the value applies); a task value for a potential without a
 *  DM_REAM block -> WARNING (ignored). */
void check_embf0_precedence(const tadah::core::Context &pot,
                            const tadah::core::Context &task,
                            const std::string &pot_name);

}  // namespace models
}  // namespace tadah

#endif
