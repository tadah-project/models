#ifndef M_CORE_H
#define M_CORE_H

#include <tadah/core/core_types.h>

/**
 * @class M_Core
 * @brief Abstract base class for model cores.
 *
 * Provides common functionality and interface for model implementations, including training status
 * and weight management.
 */
class M_Core {
public:
    int verbose = 0; ///< Verbose level for logging.
    bool trained = false; ///< Indicates if the model has been trained.
    t_type weights; ///< Weights vector for the model.

  int is_verbose() {return verbose;}

    /**
     * @brief Virtual destructor for polymorphic deletion.
     */
    virtual ~M_Core() {}

    /**
     * @brief Checks if the model has been trained.
     *
     * @return True if the model is trained, otherwise false.
     */
    bool is_trained() const {
        return trained;
    }

    /**
     * @brief Retrieves the weights of the model.
     *
     * @return Constant reference to the weights vector.
     */
    const t_type &get_weights() const {
        return weights;
    }

    /**
     * @brief Sets the model weights.
     *
     * @param w New weights vector to be set.
     */
    void set_weights(const t_type w) {
        weights = w;
    }

    /**
     * @brief Pure virtual function for making predictions.
     *
     * Must be implemented by derived classes.
     *
     * @param v Input vector for prediction.
     * @return Predicted value.
     */
    virtual double predict(const aed_type2 &v) const = 0;

    /**
     * @brief Pure virtual function to get weights' uncertainty.
     *
     * Must be implemented by derived classes.
     *
     * @return Vector of uncertainties for the weights.
     */
    virtual t_type get_weights_uncertainty() const = 0;
};

#endif // M_CORE_H

