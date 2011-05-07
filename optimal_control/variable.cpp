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

#include "variable_internal.hpp"
#include "../casadi/casadi_exception.hpp"

using namespace std;
namespace CasADi{
namespace OptimalControl{
  
Variable::Variable(){
}

Variable::Variable(const string& name, bool create_expression){
  assignNode(new VariableInternal(name));
  if(create_expression){
    setExpression(SX(name));
  }
}


Variable::~Variable(){
}

VariableInternal* Variable::operator->(){
  return (VariableInternal*)(SharedObject::operator->());
}

const VariableInternal* Variable::operator->() const{
  return (const VariableInternal*)(SharedObject::operator->());
}
  
// Variable::operator SX() const{
//   return var();  
// }

SX Variable::der(bool allocate) const{
  return (*this)->der(allocate);  
}

SX Variable::der(bool allocate){
  return (*this)->der(allocate);  
}

SX Variable::var() const{
  return (*this)->var();
}

const string& Variable::getName() const{
  return (*this)->name_;
}

void Variable::setName(const string& name){
  (*this)->name_ = name;
}

Variability Variable::getVariability() const{
  return (*this)->variability_;
}

void Variable::setVariability(Variability variability){
  (*this)->variability_ = variability;
}

Causality Variable::getCausality() const{
  return (*this)->causality_;
}

void Variable::setCausality(Causality causality){
  (*this)->causality_ = causality;
}
    
Alias Variable::getAlias() const{
  return (*this)->alias_;
}

void Variable::setAlias(Alias alias){
  (*this)->alias_ = alias;
}
    
const string& Variable::getDescription() const{
  return (*this)->description_;
}

void Variable::setDescription(const string& description){
  (*this)->description_ = description;
}
    
int Variable::getValueReference() const{
  return (*this)->valueReference_;
}

void Variable::setValueReference(int valueReference){
  (*this)->valueReference_ = valueReference;
}
    
double Variable::getMin() const{
  return (*this)->min_;
}

void Variable::setMin(double min){
  (*this)->min_ = min;
}
    
double Variable::getMax() const{
  return (*this)->max_;
}

void Variable::setMax(double max){
  (*this)->max_ = max;
}
    
double Variable::getNominal() const{
  return (*this)->nominal_;
}

void Variable::setNominal(double nominal){
  (*this)->nominal_ = nominal;
}
    
double Variable::getStart() const{
  return (*this)->start_;
}

void Variable::setStart(double start){
  (*this)->start_ = start;
}
    
const string& Variable::getUnit() const{
  return (*this)->unit_;
}

void Variable::setUnit(const string& unit){
  (*this)->unit_ = unit;
}
    
const string& Variable::getDisplayUnit() const{
  return (*this)->displayUnit_;
}

void Variable::setDisplayUnit(const string& displayUnit){
  (*this)->displayUnit_ = displayUnit;
}

void Variable::setExpression(const SX& sx){
  (*this)->sx_ = sx;
}

void Variable::setDerivative(const SX& dx){
  (*this)->dx_ = dx;
}

bool Variable::checkNode() const{
  return dynamic_cast<const VariableInternal*>(get());
}

SX Variable::atTime(double t, bool allocate) const{
  return (*this)->atTime(t,allocate);
}

SX Variable::atTime(double t, bool allocate){
  return (*this)->atTime(t,allocate);
}

int Variable::index() const{
  return (*this)->index_;
}

void Variable::setIndex(int ind){
  (*this)->index_ = ind;
}
    
bool Variable::isDifferential() const{
  return !(*this)->dx_.isNan();
}

void Variable::setDependent(int dep){
  (*this)->dependent_ = dep;
}
    
bool Variable::getDependent() const{
  return (*this)->dependent_;
}

SX Variable::highest() const{
  if(isDifferential()){
    return der();
  } else {
    return var();
  }
}


} // namespace OptimalControl
} // namespace CasADi

