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

#ifndef OPTIMICA_OCP_HPP
#define OPTIMICA_OCP_HPP

#include "casadi/printable_object.hpp"
#include "variable.hpp"

namespace CasADi{
  namespace OptimalControl{
    
    /// Tree structure for storing variables
    class VariableTree{
      public:
        /// Access a sub-collection by name
        VariableTree& subByName(const std::string& name, bool allocate=false);

        /// Access a sub-collection by index
        VariableTree& subByIndex(int ind, bool allocate=false);
        
        /// Get all variables
        void getAll(std::vector<Variable>& v, bool skip_dependent=false) const;
    
        /// Print node
        #ifndef SWIG
        void print(std::ostream &stream, int indent=0) const;
        #endif // SWIG

        /// Variable
        Variable var_;
        
        /// Children nodes
        std::vector<VariableTree> children_;
        
        /// Names of children
        std::map<std::string,int> name_part_;
    };
    
/** Symbolic, object oriented representation of an optimal control problem (OCP) */
class OCP : public PrintableObject{
  public:
    /// OCP
    OCP();

#ifndef SWIG
    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;
    
    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;
#endif
    /// Sort variables according to type
    void sortType();

    /// Divide the state into its differential and algebraic components
    void sortState();

    /// Eliminate dependent equations
    void eliminateDependent();
    
    /// Sort the variables and equations according to BLT
    void sortBLT();
    
    /// Try to make explicit by symbolically solving for xdot (experimental, only small systems)
    void makeExplicit();

    /// Replace all state derivatives by algebraic variables with the same name
    void makeSemiExplicit();
    
    /// Add a binding equation
    void addExplicitEquation(const SX& var, const SX& bind_eq);
    
    /// Create a new, scaled OCP
    void scale();
    
    /// Access the variables in a class hierarchy -- public data member
    VariableTree variables_;
    
    /// Time
    SX t_;
    
    /// Differential states
    std::vector<Variable> xd_;

    /// Algebraic states
    std::vector<Variable> xa_;
    
    /// States
    std::vector<Variable> x_;

    /// Controls
    std::vector<Variable> u_;
    
    /// Free parameters
    std::vector<Variable> p_;

    /// Explicit equations
    std::vector<SX> explicit_lhs_, explicit_rhs_;
    
    /// Dynamic equations
    std::vector<SX> dynamic_eq_;
    
    /// Initial equations
    std::vector<SX> initial_eq_;

    /// Constraint function with upper and lower bounds
    std::vector<SX> cfcn, cfcn_lb, cfcn_ub;

    /// Mayer objective terms
    std::vector<SX> mterm;

    /// Mayer time time points TODO: Remove this when WITH_TIMEDVARIABLE is default
    std::vector<double> mtp;
    
    /// Lagrange objective terms
    std::vector<SX> lterm;
    
    /// Initial time
    double t0;
    
    /// Initial time is free
    bool t0_free;
    
    /// Final time
    double tf;
    
    /// Final time is free
    bool tf_free;
    
    /// Is scaled?
    bool is_scaled_;

    /// BLT blocks
    std::vector<int> rowblock_;  // block k is rows r[k] to r[k+1]-1
    std::vector<int> colblock_;  // block k is cols s[k] to s[k+1]-1

    /// BLT sorted?
    bool blt_sorted_;
    
};

#ifdef SWIG
%extend OCP{
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif

  } // namespace OptimalControl
} // namespace CasADi

#endif // OPTIMICA_OCP_HPP


