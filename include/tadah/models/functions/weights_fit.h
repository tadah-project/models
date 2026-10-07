#ifndef TADAH_MODELS_WEIGHTS_FIT_H
#define TADAH_MODELS_WEIGHTS_FIT_H

#include <tadah/core/context.h>
#include <tadah/models/functions/function_base.h>

#include <cstddef>

namespace tadah {
namespace models {

/** \brief Sizes of the trained arrays a model predicts with. */
struct TrainedSizes {
  std::size_t weights = 0;   ///< WEIGHTS
  std::size_t nmean = 0;     ///< NMEAN
  std::size_t nstdev = 0;    ///< NSTDEV
};

/** \brief The sizes of the WEIGHTS, NMEAN and NSTDEV keys stored in @p c
 *  (0 for an absent key). */
TrainedSizes trained_sizes(const tadah::core::Context &c);

/** \brief Refuses a potential whose trained arrays do not fit the
 *  descriptor it predicts from.
 *
 *  The predictors' dot products run over the weights, so weights written
 *  for another descriptor layout give silently wrong energies and forces.
 *  Required, for a descriptor of @p dsize entries (DSIZE) and the basis
 *  function or kernel @p f:
 *    - @p s.weights == f.weights_count(dsize);
 *    - a kernel with a basis: basis vectors of dsize entries;
 *    - NORM true in @p c: dsize values of NMEAN and of NSTDEV.
 *
 *  The message names WEIGHTS (NMEAN, NSTDEV, BASIS), DSIZE, BIAS and
 *  ATOM. When BIAS is true with N >= 2 ATOM species and the sizes are
 *  those of a descriptor with one intercept column instead of N (DSIZE -
 *  N + 1: how tadah before 1.6.0 trained BIAS true), it says so and that
 *  the potential has to be retrained.
 *
 *  Called where trained weights are used to predict (MLIP M_BLR / M_KRR
 *  predictions, pair_style tadah's pair_coeff), never where a model is
 *  built for training: a training configuration may carry stale WEIGHTS,
 *  NMEAN and NSTDEV, which the fit replaces.
 *
 *  @throws std::runtime_error */
void require_weights_fit(const tadah::core::Context &c, std::size_t dsize,
                         const Function_Base &f, const TrainedSizes &s);

/** \brief require_weights_fit with the sizes of the WEIGHTS, NMEAN and
 *  NSTDEV keys of @p c (a potential as it is loaded). */
void require_weights_fit(const tadah::core::Context &c, std::size_t dsize,
                         const Function_Base &f);

} // namespace models
} // namespace tadah

#endif
