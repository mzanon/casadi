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


#include "mapsum_internal.hpp"

namespace casadi {
  using namespace std;

  MapSum::MapSum() {
  }

  MapSum::MapSum(const std::string& name, const Function& f,
                 int n,
                 const std::vector<bool> &repeat_in,
                 const Dict& opts) {
    assignNode(new MapSumInternal(f, n, repeat_in));
    setOption("name", name);
    setOption(opts);
    init();
  }

  MapSumInternal* MapSum::operator->() {
    return static_cast<MapSumInternal*>(Function::operator->());
  }

  const MapSumInternal* MapSum::operator->() const {
    return static_cast<const MapSumInternal*>(Function::operator->());
  }

  bool MapSum::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const MapSumInternal*>(ptr)!=0;
  }

} // namespace casadi
