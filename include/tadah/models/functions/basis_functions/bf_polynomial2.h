#ifndef BF_POLYNOMIAL2_H
#define BF_POLYNOMIAL2_H
#include <tadah/models/functions/basis_functions/bf_base.h>

/**
 * @class BF_Polynomial2
 * @brief Polynomial order 2 basis function with coefficient c=0.
 *
 * This structure provides methods for predicting energy and forces using a second-order polynomial basis function.
 */
struct BF_Polynomial2 : public virtual BF_Base {

    /**
     * @brief Default constructor.
     */
    BF_Polynomial2();

    /**
     * @brief Constructor with configuration.
     * @param c Configuration object reference.
     */
    BF_Polynomial2(const Config &c);

    std::string label = "BF_Polynomial2"; ///< Label identifying the basis function.

    /**
     * @brief Retrieve the label of the basis function.
     * @return String label of the basis function.
     */
    std::string get_label() const;

    /**
     * @brief Predicts energy using the polynomial basis function.
     * @param weights Weights vector.
     * @param aed Atomic descriptors.
     * @return Predicted energy.
     */
    double epredict(const t_type &weights, const aed_type2 &aed) const;

    /**
     * @brief Predicts force for a specific dimension using the polynomial basis function.
     * @param weights Weights vector.
     * @param fdij Force descriptor between atoms.
     * @param aedi Atomic descriptors for an individual atom.
     * @param k Dimension index.
     * @return Predicted force component.
     */
    double fpredict(const t_type &weights, const fd_type &fdij,
                    const aed_type2 &aedi, const size_t k) const;

    /**
     * @brief Predicts forces using the polynomial basis function.
     * @param weights Weights vector.
     * @param fdij Force descriptor between atoms.
     * @param aedi Atomic descriptors for an individual atom.
     * @return Predicted forces.
     */
    force_type fpredict(const t_type &weights, const fd_type &fdij,
                        const aed_type2 &aedi) const;
};
#endif
