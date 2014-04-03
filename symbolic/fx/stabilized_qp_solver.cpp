/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "stabilized_qp_solver_internal.hpp"

using namespace std;
namespace CasADi{


StabilizedQPSolver::StabilizedQPSolver() {
}

StabilizedQPSolverInternal* StabilizedQPSolver::operator->(){
  return static_cast<StabilizedQPSolverInternal*>(FX::operator->());
}

const StabilizedQPSolverInternal* StabilizedQPSolver::operator->() const{
    return static_cast<const StabilizedQPSolverInternal*>(FX::operator->());
}

bool StabilizedQPSolver::checkNode() const{
  return dynamic_cast<const StabilizedQPSolverInternal*>(get())!=0;
}

void StabilizedQPSolver::setLPOptions() {
  (*this)->setLPOptions();
}

void StabilizedQPSolver::generateNativeCode(const std::string &filename) const {
  std::ofstream file;
  file.open (filename.c_str());
  (*this)->generateNativeCode(file);
}

} // namespace CasADi

  

