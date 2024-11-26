#ifndef M_TRAIN_H
#define M_TRAIN_H

#include <tadah/core/core_types.h>

/**
 * @class M_Train
 * @brief Abstract base class for training models.
 *
 * Defines the interface for training models with a design matrix and target values.
 */
class M_Train {
public:
    /**
     * @brief Virtual destructor for polymorphic deletion.
     */
    virtual ~M_Train() {}

    /**
     * @brief Pure virtual function to train the model.
     *
     * Must be implemented by derived classes.
     *
     * @param Phi Design matrix containing input features.
     * @param T Target vector for training.
     */
    virtual void train(phi_type &Phi, t_type &T) = 0;
};

#endif // M_TRAIN_H
