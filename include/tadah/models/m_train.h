#ifndef M_TRAIN_H
#define M_TRAIN_H

#include <tadah/core/core_types.h>
#include <tadah/core/utils.h>
#include <tadah/models/memory/imodels_workspace_manager.h>
#include <tadah/models/memory/models_workspace_manager.h>

namespace tadah {
namespace models {

/**
 * @class M_Train
 * @brief Abstract base class for training models.
 *
 * Defines the interface for training models with a design matrix and target values.
 */
class M_Train {
protected:
  bool own_workspace_manager_ = false;
  tadah::models::memory::IModelsWorkspaceManager* workspaceManager_ = nullptr;

public:
  M_Train():
    own_workspace_manager_(true),
    workspaceManager_(new tadah::models::memory::ModelsWorkspaceManager())
  {}
  M_Train(tadah::models::memory::IModelsWorkspaceManager &workspaceManager):
    own_workspace_manager_(false),
    workspaceManager_(&workspaceManager)
  {}
  /**
     * @brief Virtual destructor for polymorphic deletion.
     */
  virtual ~M_Train() {
    if (workspaceManager_ != nullptr && own_workspace_manager_) {
      delete workspaceManager_;
      workspaceManager_ = nullptr;
      own_workspace_manager_ = false;
    }
  }

  /**
     * @brief Pure virtual function to train the model.
     *
     * Must be implemented by derived classes.
     *
     * @param Phi Design matrix containing input features.
     * @param T Target vector for training.
     */
  virtual void train(tadah::core::phi_type &Phi, tadah::core::t_type &T) = 0;

  /**
 * Swaps columns and fits a reduced matrix to obtain residuals.
 *
 * This function processes a full matrix Phi by swapping columns specified by indices,
 * creating a new pointer for the fixed matrix Phi_f, and then use w_f weights
 * to obtain prediction T_f = Phi_f w_f. The function then copmutes the residual vector
 * T_r = T - T_f. On exit, the Phi matrix is rearranged such that the indices.size()
 * rightmost columns of the original matrix are the same as Phi_f.
 * The function return a mapping between columns of original Phi matrix and rearranged.
 *
 * @param Phi         The full matrix to be reduced and fitted to the adjusted target
 *                    vector.
 * @param T           The full target vector to which the reduced matrix Phi will
 *                    be fitted.
 * @param indices     A vector of indices specifying which columns in the
 *                    original matrix are fixed and should be considered during
 *                    the regression.
 * @param w_f         A vector of fixed weights corresponding to the specified
 *                    indices. These weights are used to adjust the target
 *                    vector, resulting in residuals.
 */
  template <typename M>
  std::vector<size_t> prep_train_with_residuals(
    tadah::core::Matrix_Base<M> &Phi, tadah::core::t_type &T, std::vector<size_t> &indices, const tadah::core::t_type w_f, tadah::core::t_type &T_r) {

    if (w_f.size() != indices.size()) {
      throw std::runtime_error("Size of w_f must be equal to size of indices.");
    }

    std::sort(indices.begin(), indices.end());

    if (indices.back() >= Phi.cols()) {
      throw std::runtime_error("Largest index is greater than the number of columns in Phi.");
    }

    std::vector<size_t> move_map = move_columns_in_place(Phi,indices);

    //   indices.size() = 1
    //   indices[0] 1
    //   1 4 7        1 7    4
    //   2 5 8   ->   2 8    5
    //   3 6 9        3 9    6
    //    Phi        Phi_o   Phi_f
    //    move_map = {0,2,1}

    // the fixed columns are to the right
    tadah::core::MatrixView Phi_f(Phi, Phi.cols()-indices.size());

    tadah::core::t_type T_f = Phi_f * w_f;

    T_r = T-T_f;

    return move_map;
  }

  void setWorkspaceManager(tadah::models::memory::IModelsWorkspaceManager  *workspaceManager) {
    if (workspaceManager_ != nullptr && own_workspace_manager_) {
      delete workspaceManager_;
    }
    workspaceManager_ = workspaceManager;
  }


};

} // namespace models
} // namespace tadah
#endif // M_TRAIN_H
