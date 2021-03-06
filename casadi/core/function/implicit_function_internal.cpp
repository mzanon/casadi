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


#include "implicit_function_internal.hpp"
#include "mx_function.hpp"
#include "../mx/mx_node.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include <iterator>

#include "../casadi_options.hpp"
#include "../profiling.hpp"

using namespace std;
namespace casadi {

  ImplicitFunctionInternal::ImplicitFunctionInternal(const Function& f) : f_(f) {
    addOption("linear_solver",            OT_STRING, "csparse",
              "User-defined linear solver class. Needed for sensitivities.");
    addOption("linear_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the linear solver.");
    addOption("constraints",              OT_INTEGERVECTOR, GenericType(),
              "Constrain the unknowns. 0 (default): no constraint on ui, "
              "1: ui >= 0.0, -1: ui <= 0.0, 2: ui > 0.0, -2: ui < 0.0.");
    addOption("implicit_input",           OT_INTEGER,      0,
              "Index of the input that corresponds to the actual root-finding");
    addOption("implicit_output",          OT_INTEGER,      0,
              "Index of the output that corresponds to the actual root-finding");
    addOption("jacobian_function",        OT_FUNCTION,  GenericType(),
              "Function object for calculating the Jacobian (autogenerated by default)");
    addOption("linear_solver_function",   OT_FUNCTION,  GenericType(),
              "Function object for solving the linearized problem (autogenerated by default)");
  }

  ImplicitFunctionInternal::~ImplicitFunctionInternal() {
  }

  void ImplicitFunctionInternal::deepCopyMembers(std::map<SharedObjectNode*,
                                                 SharedObject>& already_copied) {
    FunctionInternal::deepCopyMembers(already_copied);
    f_ = deepcopy(f_, already_copied);
    jac_ = deepcopy(jac_, already_copied);
    linsol_ = deepcopy(linsol_, already_copied);
  }

  void ImplicitFunctionInternal::init() {

    // Initialize the residual function
    f_.init(false);

    // Which input/output correspond to the root-finding problem?
    iin_ = getOption("implicit_input");
    iout_ = getOption("implicit_output");

    // Get the number of equations and check consistency
    casadi_assert_message(iin_>=0 && iin_<f_.getNumInputs() && f_.getNumInputs()>0,
                          "Implicit input not in range");
    casadi_assert_message(iout_>=0 && iout_<f_.getNumOutputs() && f_.getNumOutputs()>0,
                          "Implicit output not in range");
    casadi_assert_message(f_.output(iout_).isDense() && f_.output(iout_).isVector(),
                          "Residual must be a dense vector");
    casadi_assert_message(f_.input(iin_).isDense() && f_.input(iin_).isVector(),
                          "Unknown must be a dense vector");
    n_ = f_.output(iout_).nnz();
    casadi_assert_message(n_ == f_.input(iin_).nnz(),
                          "Dimension mismatch. Input size is "
                          << f_.input(iin_).nnz()
                          << ", while output size is "
                          << f_.output(iout_).nnz());

    // Allocate inputs
    setNumInputs(f_.getNumInputs());
    for (int i=0; i<getNumInputs(); ++i) {
      input(i) = f_.input(i);
    }

    // Allocate output
    setNumOutputs(f_.getNumOutputs());
    for (int i=0; i<getNumOutputs(); ++i) {
      output(i) = f_.output(i);
    }

    // Same input and output schemes
    setInputScheme(f_.getInputScheme());
    setOutputScheme(f_.getOutputScheme());

    // Call the base class initializer
    FunctionInternal::init();

    // Get the Jacobian function object, if any
    if (hasSetOption("jacobian_function")) {
      jac_ = getOption("jacobian_function");
    }

    // Generate Jacobian if not provided
    if (jac_.isNull()) jac_ = f_.jacobian(iin_, iout_);
    jac_.init(false);

    // Check for structural singularity in the Jacobian
    casadi_assert_message(
      !jac_.output().sparsity().isSingular(),
      "ImplicitFunctionInternal::init: singularity - the jacobian is structurally rank-deficient. "
      "sprank(J)=" << sprank(jac_.output()) << " (instead of "<< jac_.output().size1() << ")");

    // Get the linear solver function object, if any
    if (hasSetOption("linear_solver_function")) {
      Function linsol = getOption("linear_solver_function");
      linsol_ = shared_cast<LinearSolver>(linsol);
    }

    // Get the linear solver creator function
    if (linsol_.isNull()) {
      if (hasSetOption("linear_solver")) {
        std::string linear_solver_name = getOption("linear_solver");

        // Allocate an NLP solver
        linsol_ = LinearSolver(linear_solver_name, jac_.output().sparsity(), 1);

        // Pass options
        if (hasSetOption("linear_solver_options")) {
          const Dictionary& linear_solver_options = getOption("linear_solver_options");
          linsol_.setOption(linear_solver_options);
        }

        // Initialize
        linsol_.init();
      }
    } else {
      // Initialize the linear solver, if provided
      linsol_.init(false);
      casadi_assert(linsol_.input().sparsity()==jac_.output().sparsity());
    }

    // No factorization yet;
    fact_up_to_date_ = false;

    // Constraints
    if (hasSetOption("constraints")) u_c_ = getOption("constraints");

    casadi_assert_message(u_c_.size()==n_ || u_c_.empty(),
                          "Constraint vector if supplied, must be of length n, but got "
                          << u_c_.size() << " and n = " << n_);
  }

  void ImplicitFunctionInternal::evaluate() {

    // Mark factorization as out-of-date. TODO: make this conditional
    fact_up_to_date_ = false;

    // Get initial guess
    output(iout_).set(input(iin_));

    // Solve the nonlinear system of equations
    solveNonLinear();
  }

  Function ImplicitFunctionInternal::getDerForward(int nfwd) {
    // Symbolic expression for the input
    vector<MX> arg = symbolicInput();
    arg[iin_] = MX::sym(arg[iin_].getName() + "_guess",
                        Sparsity(arg[iin_].shape()));
    vector<MX> res = symbolicOutput();
    vector<vector<MX> > fseed = symbolicFwdSeed(nfwd, arg), fsens;
    callForward(arg, res, fseed, fsens, false, false);

    // Construct return function
    arg.insert(arg.end(), res.begin(), res.end());
    for (int d=0; d<nfwd; ++d) arg.insert(arg.end(), fseed[d].begin(), fseed[d].end());
    res.clear();
    for (int d=0; d<nfwd; ++d) res.insert(res.end(), fsens[d].begin(), fsens[d].end());
    return MXFunction(arg, res);
  }

  Function ImplicitFunctionInternal::getDerReverse(int nadj) {
    // Symbolic expression for the input
    vector<MX> arg = symbolicInput();
    arg[iin_] = MX::sym(arg[iin_].getName() + "_guess",
                        Sparsity(arg[iin_].shape()));
    vector<MX> res = symbolicOutput();
    vector<vector<MX> > aseed = symbolicAdjSeed(nadj, res), asens;
    callReverse(arg, res, aseed, asens, false, false);

    // Construct return function
    arg.insert(arg.end(), res.begin(), res.end());
    for (int d=0; d<nadj; ++d) arg.insert(arg.end(), aseed[d].begin(), aseed[d].end());
    res.clear();
    for (int d=0; d<nadj; ++d) res.insert(res.end(), asens[d].begin(), asens[d].end());
    return MXFunction(arg, res);
  }

  void ImplicitFunctionInternal::spEvaluate(bool fwd) {

    // Initialize the callback for sparsity propagation
    f_.spInit(fwd);

    if (fwd) {

      // Pass inputs to function
      f_.input(iin_).setZeroBV();
      for (int i=0; i<getNumInputs(); ++i) {
        if (i!=iin_) f_.input(i).setBV(input(i));
      }

      // Propagate dependencies through the function
      f_.spEvaluate(true);

      // "Solve" in order to propagate to z
      output(iout_).setZeroBV();
      linsol_.spSolve(output(iout_), f_.output(iout_), false);

      // Propagate to auxiliary outputs
      if (getNumOutputs()>1) {
        f_.input(iin_).setBV(output(iout_));
        f_.spEvaluate(true);
        for (int i=0; i<getNumOutputs(); ++i) {
          if (i!=iout_) output(i).setBV(f_.output(i));
        }
      }

    } else {

      // Propagate dependencies from auxiliary outputs
      if (getNumOutputs()>1) {
        f_.output(iout_).setZeroBV();
        for (int i=0; i<getNumOutputs(); ++i) {
          if (i!=iout_) f_.output(i).setBV(output(i));
        }
        f_.spEvaluate(false);
        for (int i=0; i<getNumInputs(); ++i) {
          input(i).setBV(f_.input(i));
        }
      } else {
        for (int i=0; i<getNumInputs(); ++i) {
          input(i).setZeroBV();
        }
      }

      // Add dependency on implicitly defined variable
      input(iin_).borBV(output(iout_));

      // "Solve" in order to get seed
      f_.output(iout_).setZeroBV();
      linsol_.spSolve(f_.output(iout_), input(iin_), true);

      // Propagate dependencies through the function
      f_.spEvaluate(false);

      // Collect influence on inputs
      for (int i=0; i<getNumInputs(); ++i) {
        if (i!=iin_) input(i).borBV(f_.input(i));
      }

      // No dependency on the initial guess
      input(iin_).setZeroBV();
    }
  }

  std::map<std::string, ImplicitFunctionInternal::Plugin> ImplicitFunctionInternal::solvers_;

  const std::string ImplicitFunctionInternal::infix_ = "implicitfunction";

  void ImplicitFunctionInternal::
  callForward(const std::vector<MX>& arg, const std::vector<MX>& res,
          const std::vector<std::vector<MX> >& fseed,
          std::vector<std::vector<MX> >& fsens,
          bool always_inline, bool never_inline) {
    // Number of directional derivatives
    int nfwd = fseed.size();
    fsens.resize(nfwd);

    // Quick return if no seeds
    if (nfwd==0) return;

    // Propagate through f_
    vector<MX> f_arg(arg);
    f_arg.at(iin_) = res.at(iout_);
    vector<MX> f_res(res);
    f_res.at(iout_) = MX(input(iin_).shape()); // zero residual
    std::vector<std::vector<MX> > f_fseed(fseed);
    for (int d=0; d<nfwd; ++d) {
      f_fseed[d].at(iin_) = MX(input(iin_).shape()); // ignore seeds for guess
    }
    f_.callForward(f_arg, f_res, f_fseed, fsens, always_inline, never_inline);

    // Get expression of Jacobian
    MX J = jac_(f_arg).front();

    // Solve for all the forward derivatives at once
    vector<MX> rhs(nfwd);
    for (int d=0; d<nfwd; ++d) rhs[d] = vec(fsens[d][iout_]);
    rhs = horzsplit(J->getSolve(-horzcat(rhs), false, linsol_));
    for (int d=0; d<nfwd; ++d) fsens[d][iout_] = reshape(rhs[d], input(iin_).shape());

    // Propagate to auxiliary outputs
    int num_out = getNumOutputs();
    if (num_out>1) {
      for (int d=0; d<nfwd; ++d) f_fseed[d][iin_] = fsens[d][iout_];
      f_.callForward(f_arg, f_res, f_fseed, fsens, always_inline, never_inline);
      for (int d=0; d<nfwd; ++d) fsens[d][iout_] = f_fseed[d][iin_]; // Otherwise overwritten
    }
  }

  void ImplicitFunctionInternal::
  callReverse(const std::vector<MX>& arg, const std::vector<MX>& res,
          const std::vector<std::vector<MX> >& aseed,
          std::vector<std::vector<MX> >& asens,
          bool always_inline, bool never_inline) {

    // Number of directional derivatives
    int nadj = aseed.size();
    asens.resize(nadj);

    // Quick return if no seeds
    if (nadj==0) return;

    // Get expression of Jacobian
    vector<MX> f_arg(arg);
    f_arg[iin_] = res.at(iout_);
    MX J = jac_(f_arg).front();

    // Get adjoint seeds for calling f
    int num_out = getNumOutputs();
    int num_in = getNumInputs();
    vector<MX> f_res(res);
    f_res[iout_] = MX(input(iin_).shape()); // zero residual
    vector<vector<MX> > f_aseed(nadj);
    for (int d=0; d<nadj; ++d) {
      f_aseed[d].resize(num_out);
      for (int i=0; i<num_out; ++i) f_aseed[d][i] = i==iout_ ? f_res[iout_] : aseed[d][i];
    }

    // Propagate dependencies from auxiliary outputs
    vector<MX> rhs(nadj);
    vector<vector<MX> > asens_aux;
    if (num_out>1) {
      f_.callReverse(f_arg, f_res, f_aseed, asens_aux, always_inline, never_inline);
      for (int d=0; d<nadj; ++d) rhs[d] = vec(asens_aux[d][iin_] + aseed[d][iout_]);
    } else {
      for (int d=0; d<nadj; ++d) rhs[d] = vec(aseed[d][iout_]);
    }

    // Solve for all the adjoint seeds at once
    rhs = horzsplit(J->getSolve(-horzcat(rhs), true, linsol_));
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<num_out; ++i) {
        if (i==iout_) {
          f_aseed[d][i] = reshape(rhs[d], output(i).shape());
        } else {
          // Avoid counting the auxiliary seeds twice
          f_aseed[d][i] = MX(output(i).shape());
        }
      }
    }

    // Propagate through f_
    f_.callReverse(f_arg, f_res, f_aseed, asens, always_inline, never_inline);

    // No dependency on guess
    for (int d=0; d<nadj; ++d) {
      asens[d][iin_] = MX(input(iin_).shape());
    }

    // Add contribution from auxiliary outputs
    if (num_out>1) {
      for (int d=0; d<nadj; ++d) {
        for (int i=0; i<num_in; ++i) if (i!=iin_) asens[d][i] += asens_aux[d][i];
      }
    }
  }

} // namespace casadi
