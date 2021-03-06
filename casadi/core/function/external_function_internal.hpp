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


#ifndef CASADI_EXTERNAL_FUNCTION_INTERNAL_HPP
#define CASADI_EXTERNAL_FUNCTION_INTERNAL_HPP

#include "external_function.hpp"
#include "function_internal.hpp"

#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#define NOMINMAX
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL

/// \cond INTERNAL

namespace casadi {

  class CASADI_EXPORT ExternalFunctionInternal : public FunctionInternal {
    friend class ExternalFunction;
  public:

    /** \brief  constructor */
    explicit ExternalFunctionInternal(const std::string& bin_name);

    /** \brief  clone function */
    virtual ExternalFunctionInternal* clone() const;

    /** \brief  Destructor */
    virtual ~ExternalFunctionInternal();

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const cpv_double& arg, const pv_double& res, int* itmp, double* rtmp);

    /** \brief  Work vector sizes */
    virtual void nTmp(size_t& ni, size_t& nr);

    /** \brief  Initialize */
    virtual void init();

  protected:

    ///@{
    /** \brief  Function pointer types */
    typedef int (*evalPtr)(const double* const* arg, double* const* res, int* iw, double* w);
    typedef int (*initPtr)(int *n_in, int *n_out);
    typedef int (*getSparsityPtr)(int n_in, int *n_row, int *n_col, int **colind, int **row);
    typedef int (*nworkPtr)(int *ni, int *nr);
    ///@}

    /** \brief  Name of binary */
    std::string bin_name_;

    /** \brief  Function pointers */
    evalPtr eval_;

#if defined(WITH_DL) && defined(_WIN32) // also for 64-bit
    typedef HINSTANCE handle_t;
#else
    typedef void* handle_t;
#endif

    /** \brief  handle to the dll */
    handle_t handle_;

    /** \brief Work vector sizes */
    size_t ni_, nr_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_EXTERNAL_FUNCTION_INTERNAL_HPP
