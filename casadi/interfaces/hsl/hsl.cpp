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


#include "hsl.hpp"
#include "../../core/std_vector_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINEARSOLVER_HSL_EXPORT
  casadi_register_linearsolver_hsl(LinearSolverInternal::Plugin* plugin) {
    plugin->creator = Hsl::creator;
    plugin->name = "hsl";
    plugin->doc = Hsl::meta_doc.c_str();;
    plugin->version = 22;
    return 0;
  }

  extern "C"
  void CASADI_LINEARSOLVER_HSL_EXPORT casadi_load_linearsolver_hsl() {
    LinearSolverInternal::registerPlugin(casadi_register_linearsolver_hsl);
  }

  Hsl::Hsl(const Sparsity& sparsity, int nrhs) :
      LinearSolverInternal(sparsity, nrhs) {
  }

  Hsl::~Hsl() {
  }

  void Hsl::init() {
    // Call the base class initializer
    LinearSolverInternal::init();

  }

  void Hsl::prepare() {
    prepared_ = false;

  }

  void Hsl::solve(double* x, int nrhs, bool transpose) {

  }

  Hsl* Hsl::clone() const {
    return new Hsl(*this);
  }

} // namespace casadi
