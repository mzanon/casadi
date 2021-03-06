/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "socp_solver_internal.hpp"
#include "../matrix/matrix_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"
#include <numeric>
#include <functional>

INPUTSCHEME(SOCPInput)
OUTPUTSCHEME(SOCPOutput)

using namespace std;
namespace casadi {

  // Constructor
  SocpSolverInternal::SocpSolverInternal(const std::vector<Sparsity> &st) : st_(st) {
    addOption("ni", OT_INTEGERVECTOR, GenericType(),
              "Provide the size of each SOC constraint. Must sum up to N.");
    addOption("print_problem", OT_BOOLEAN, false, "Print out problem statement for debugging.");

    input_.scheme = SCHEME_SOCPInput;
    output_.scheme = SCHEME_SOCPOutput;

  }

  void SocpSolverInternal::init() {
    // Call the init method of the base class
    FunctionInternal::init();

    ni_ = getOption("ni");
    print_problem_ = getOption("print_problem");

    m_ = ni_.size();

    const Sparsity& A = st_[SOCP_STRUCT_A];
    const Sparsity& G = st_[SOCP_STRUCT_G];
    const Sparsity& E = st_[SOCP_STRUCT_E];

    N_ = std::accumulate(ni_.begin(), ni_.end(), 0);
    casadi_assert_message(N_==G.size2(),
                          "SocpSolverInternal: Supplied G sparsity: number of cols ("
                          << G.size2()
                          <<  ")  must match sum of vector provided with option 'ni' ("
                          << N_ << ").");
    casadi_assert_message(m_==E.size2(),
                          "SocpSolverInternal: Supplied E sparsity: number of cols ("
                          << E.size2()
                          <<  ")  must match number of cone (2-norm) constraints ("
                          << m_ << ").");

    nc_ = A.size1();
    n_ = A.size2();

    casadi_assert_message(n_==G.size1(),
       "SocpSolverInternal: Supplied G sparsity: number of rows ("
        << G.size1()
        <<  ") must match number of decision variables (cols of A): " << n_ << ".");
    casadi_assert_message(n_==E.size1(),
       "SocpSolverInternal: Supplied E sparsity: number of rows ("
        << E.size1()
        <<  ") must match number of decision variables (cols of A): " << n_ << ".");

    // Input arguments
    setNumInputs(SOCP_SOLVER_NUM_IN);
    input(SOCP_SOLVER_G) = DMatrix::zeros(G);
    input(SOCP_SOLVER_H) = DMatrix::zeros(N_, 1);
    input(SOCP_SOLVER_E) = DMatrix::zeros(E);
    input(SOCP_SOLVER_F) = DMatrix::zeros(m_, 1);
    input(SOCP_SOLVER_A) = DMatrix::zeros(A);
    input(SOCP_SOLVER_C) = DMatrix::zeros(n_);
    input(SOCP_SOLVER_LBX) = -DMatrix::inf(n_);
    input(SOCP_SOLVER_UBX) = DMatrix::inf(n_);
    input(SOCP_SOLVER_LBA) = -DMatrix::inf(nc_);
    input(SOCP_SOLVER_UBA) = DMatrix::inf(nc_);

    // Output arguments
    setNumOutputs(SOCP_SOLVER_NUM_OUT);
    output(SOCP_SOLVER_X) = DMatrix::zeros(n_, 1);
    output(SOCP_SOLVER_COST) = 0.0;
    output(SOCP_SOLVER_DUAL_COST) = 0.0;
    output(SOCP_SOLVER_LAM_X) = DMatrix::zeros(n_, 1);
    output(SOCP_SOLVER_LAM_A) = DMatrix::zeros(nc_, 1);
    output(SOCP_SOLVER_LAM_CONE) = DMatrix::zeros(m_, 1);
  }

  SocpSolverInternal::~SocpSolverInternal() {
  }

  void SocpSolverInternal::evaluate() {
    throw CasadiException("SocpSolverInternal::evaluate: Not implemented");
  }

  void SocpSolverInternal::solve() {
    throw CasadiException("SocpSolverInternal::solve: Not implemented");
  }

  void SocpSolverInternal::printProblem(std::ostream &stream) const {
    stream << "SOCP Problem statement -- start" << std::endl;
    stream << "ni: "<< ni_ << std::endl;
    stream << "g: "<< std::endl;  input(SOCP_SOLVER_G).printDense(stream);
    stream << "h: "<< std::endl;  input(SOCP_SOLVER_H).printDense(stream);
    stream << "e: "<< input(SOCP_SOLVER_E) << std::endl;
    stream << "f: "<< std::endl;  input(SOCP_SOLVER_F).printDense(stream);
    stream << "c: " << input(SOCP_SOLVER_C) << std::endl;
    stream << "a: " << input(SOCP_SOLVER_A) << std::endl;
    stream << "lba: " << input(SOCP_SOLVER_LBA) << std::endl;
    stream << "uba: " << input(SOCP_SOLVER_UBA) << std::endl;
    stream << "lbx: " << input(SOCP_SOLVER_LBX) << std::endl;
    stream << "ubx: " << input(SOCP_SOLVER_UBX) << std::endl;

    stream << "SOCP Problem statement -- end" << std::endl;
  }


  void SocpSolverInternal::checkInputs() const {
    for (int i=0;i<input(SOCP_SOLVER_LBX).nnz();++i) {
      casadi_assert_message(input(SOCP_SOLVER_LBX).at(i)<=input(SOCP_SOLVER_UBX).at(i),
                            "LBX[i] <= UBX[i] was violated for i=" << i
                            << ". Got LBX[i]=" << input(SOCP_SOLVER_LBX).at(i)
                            << " and UBX[i]=" << input(SOCP_SOLVER_UBX).at(i));
    }
    for (int i=0;i<input(SOCP_SOLVER_LBA).nnz();++i) {
      casadi_assert_message(input(SOCP_SOLVER_LBA).at(i)<=input(SOCP_SOLVER_UBA).at(i),
                            "LBA[i] <= UBA[i] was violated for i=" << i
                            << ". Got LBA[i]=" << input(SOCP_SOLVER_LBA).at(i)
                            << " and UBA[i]=" << input(SOCP_SOLVER_UBA).at(i));
    }
  }

  void SocpSolverInternal::convertToDualSocp() {

    /*****************************************************************
     * Define aliases                                                *
     *****************************************************************/
    DMatrix& primal_C           = input(SOCP_SOLVER_C);
    DMatrix& primal_G           = input(SOCP_SOLVER_G);
    DMatrix& primal_H           = input(SOCP_SOLVER_H);
    DMatrix& primal_E           = input(SOCP_SOLVER_E);
    DMatrix& primal_F           = input(SOCP_SOLVER_F);
    DMatrix& primal_A           = input(SOCP_SOLVER_A);
    DMatrix& primal_LBA         = input(SOCP_SOLVER_LBA);
    DMatrix& primal_UBA         = input(SOCP_SOLVER_UBA);
    DMatrix& primal_LBX         = input(SOCP_SOLVER_LBX);
    DMatrix& primal_UBX         = input(SOCP_SOLVER_UBX);
    const Sparsity& primal_A_sparse   = primal_A.sparsity();
    const Sparsity& primal_G_sparse   = primal_G.sparsity();

    /*****************************************************************
     * Categorize number of optimization variables for dual problem  *
     *****************************************************************/

    // Empty list of indices
    primal_idx_lba_.resize(0);
    primal_idx_uba_.resize(0);
    primal_idx_lbx_.resize(0);
    primal_idx_ubx_.resize(0);

    // Loop over linear equality constraints
    for (int i=0;i<nc_;++i) {
      if (primal_LBA.at(i) != -std::numeric_limits<double>::infinity()) {
        primal_idx_lba_.push_back(i);
      }
      if (primal_UBA.at(i) != std::numeric_limits<double>::infinity()) {
        primal_idx_uba_.push_back(i);
      }
    }

    // Loop over simple bounds
    for (int i=0;i<n_;++i) {
      if (primal_LBX.at(i) != -std::numeric_limits<double>::infinity()) {
        primal_idx_lbx_.push_back(i);
      }
      if (primal_UBX.at(i) != std::numeric_limits<double>::infinity()) {
        primal_idx_ubx_.push_back(i);
      }
    }

    /*****************************************************************
     * Set up dual problem                                           *
     *****************************************************************/

    // c-vector: [fi' -hi' -LBA' UBA' -LBX' UBX']'
    dual_c_.resize(m_+N_);
    int start_idx = 0;
    int end_idx = m_;
    std::transform(primal_F.data().begin(), primal_F.data().end(), dual_c_.begin()+start_idx,
      std::negate<double>());
    start_idx = end_idx;
    end_idx += N_;
    std::copy(primal_H.data().begin(), primal_H.data().end(), dual_c_.begin()+start_idx);
    for (int k=0;k<primal_idx_lba_.size();++k) {
      int i = primal_idx_lba_[k];
      dual_c_.push_back(primal_LBA.data()[i]);
    }
    for (int k=0;k<primal_idx_uba_.size();++k) {
      int i = primal_idx_uba_[k];
      dual_c_.push_back(-primal_UBA.data()[i]);
    }
    for (int k=0;k<primal_idx_lbx_.size();++k) {
      int i = primal_idx_lbx_[k];
      dual_c_.push_back(primal_LBX.data()[i]);
    }
    for (int k=0;k<primal_idx_ubx_.size();++k) {
      int i = primal_idx_ubx_[k];
      dual_c_.push_back(-primal_UBX.data()[i]);
    }

    // A-matrix: [-ei Gi -Alba' Auba' -Ilbx Iubx]
    int begin_colind;
    int end_colind;
    dual_A_data_.resize(primal_E.size()+primal_G.size());
    dual_A_row_.resize(primal_E.size()+primal_G.size());
    dual_A_colind_.resize(m_+N_+1);

    // TODO(jgillis): replace T()-call by casadi_trans()
    DMatrix primal_A_T = primal_A.T();
    const Sparsity& primal_A_T_sparse = primal_A_T.sparsity();
    start_idx = 0;
    end_idx = primal_E.size();
    std::copy(primal_E.data().begin(), primal_E.data().end(), dual_A_data_.begin()+start_idx);
    std::transform(
      primal_E.data().begin(), primal_E.data().end(),
      dual_A_data_.begin()+start_idx, std::negate<double>());
    start_idx = end_idx;
    end_idx += primal_G.size();
    std::copy(primal_G.data().begin(), primal_G.data().end(), dual_A_data_.begin()+start_idx);
    for (int k=0;k<primal_idx_lba_.size();++k) {
      int i = primal_idx_lba_[k];
      begin_colind = primal_A_T_sparse.colind(i);
      end_colind = primal_A_T_sparse.colind(i+1);
      for (int ii=begin_colind; ii<end_colind;++ii) {
        dual_A_data_.push_back(-primal_A_T.data()[ii]);
        dual_A_row_.push_back(primal_A_T_sparse.getRow()[ii]);
      }
      dual_A_colind_.push_back(dual_A_colind_.back()+end_colind-begin_colind);
    }
    for (int k=0;k<primal_idx_uba_.size();++k) {
      int i = primal_idx_uba_[k];
      begin_colind = primal_A_T_sparse.colind(i);
      end_colind = primal_A_T_sparse.colind(i+1);
      for (int ii=begin_colind; ii<end_colind;++ii) {
        dual_A_data_.push_back(primal_A_T.data()[ii]);
        dual_A_row_.push_back(primal_A_T_sparse.getRow()[ii]);
      }
      dual_A_colind_.push_back(dual_A_colind_.back()+end_colind-begin_colind);
    }
    for (int k=0;k<primal_idx_lbx_.size();++k) {
      int i = primal_idx_lbx_[k];
      dual_A_data_.push_back(-1);
      dual_A_row_.push_back(i);
      dual_A_colind_.push_back(dual_A_colind_.back()+1);
    }
    for (int k=0;k<primal_idx_ubx_.size();++k) {
      int i = primal_idx_ubx_[k];
      dual_A_data_.push_back(1);
      dual_A_row_.push_back(i);
      dual_A_colind_.push_back(dual_A_colind_.back()+1);
    }

    // b-vector
    std::copy(primal_C.data().begin(), primal_C.data().end(), dual_b_.begin());

  }

  std::map<std::string, SocpSolverInternal::Plugin> SocpSolverInternal::solvers_;

  const std::string SocpSolverInternal::infix_ = "socpsolver";

} // namespace casadi

