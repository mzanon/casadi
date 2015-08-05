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


#ifndef CASADI_HSL_HPP
#define CASADI_HSL_HPP

#include "casadi/core/function/linear_solver_internal.hpp"
#include <casadi/interfaces/hsl/casadi_linearsolver_hsl_export.h>

/** \defgroup plugin_LinearSolver_hsl
*
* This class solves the linear system <tt>A.x=b</tt> by making 
* use of HSL
*/

/** \pluginsection{LinearSolver,hsl} */

/// \cond INTERNAL
namespace casadi {

  /** \brief  \pluginbrief{LinearSolver,hsl}
   *
   @copydoc LinearSolver_doc
   @copydoc plugin_LinearSolver_hsl
   *
   */
  class CASADI_LINEARSOLVER_HSL_EXPORT Hsl : public LinearSolverInternal {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    Hsl(const Sparsity& sparsity, int nrhs);

    // Clone
    virtual Hsl* clone() const;

    /** \brief  Create a new LinearSolver */
    static LinearSolverInternal* creator(const Sparsity& sp, int nrhs)
    { return new Hsl(sp, nrhs);}

    // Destructor
    virtual ~Hsl();

    // Initialize the solver
    virtual void init();

    // Prepare the solution of the linear system
    virtual void prepare();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    /// A documentation string
    static const std::string meta_doc;

  protected:


  };

} // namespace casadi

/// \endcond
#endif // CASADI_HSL_HPP
