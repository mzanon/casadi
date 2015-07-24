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
#include "mx_function.hpp"

using namespace std;

namespace casadi {

  MapSumInternal::MapSumInternal(const Function& f, int n, const std::vector<bool> &repeat_in)
    : f_(f), n_(n), repeat_in_(repeat_in) {

    casadi_assert_message(repeat_in_.size()==f.nIn(),
      "MapSumInternal expected repeat_in of size " << f.nIn() <<
      ", but got " << repeat_in_.size() << " instead.");

    // Give a name
    setOption("name", "unnamed_mapsum");
  }

  MapSumInternal::~MapSumInternal() {
  }

  void MapSumInternal::init() {
    int num_in = f_.nIn(), num_out = f_.nOut();

    // Initialize the functions, get input and output sparsities
    // Input and output sparsities

    ibuf_.resize(num_in);
    obuf_.resize(num_out);

    for (int i=0;i<num_in;++i) {
      if (repeat_in_[i]) {
        input(i) = DMatrix::zeros(repmat(f_.input(i).sparsity(), 1, n_));
      } else {
        input(i) = DMatrix::zeros(f_.input(i).sparsity());
      }
    }

    for (int i=0;i<num_out;++i) {
      output(i) = DMatrix::zeros(f_.output(i).sparsity());
    }

    step_in_.resize(num_in, 0);
    step_out_.resize(num_out, 0);

    for (int i=0;i<num_in;++i) {
        if (repeat_in_[i])
          step_in_[i] = f_.input(i).nnz();
    }

    for (int i=0;i<num_out;++i) {
      step_out_[i] = f_.output(i).nnz();
    }

    // Call the initialization method of the base class
    FunctionInternal::init();

    // Allocate some space to evaluate each function to.
    nnz_out_ = 0;
    for (int i=0;i<num_out;++i) {
      nnz_out_+= step_out_[i];
    }

    alloc_w(f_.sz_w() + nnz_out_);
    alloc_iw(f_.sz_iw());
    alloc_arg(2*f_.sz_arg());
    alloc_res(2*f_.sz_res());

  }

  void MapSumInternal::evalD(const double** arg, double** res, int* iw, double* w) {
    int num_in = f_.nIn(), num_out = f_.nOut();

    const double** arg1 = arg+f_.sz_arg();
    double** sum = res;

    // Clear the accumulators
    for (int k=0;k<num_out;++k) {
      if (sum[k]!=0) std::fill(sum[k], sum[k]+step_out_[k], 0);
    }

    double** res1 = res+f_.sz_res();

    for (int i=0; i<n_; ++i) {

      double* temp_res = w+f_.sz_w();
      if (temp_res!=0) std::fill(temp_res, temp_res+nnz_out_, 0);

      for (int j=0; j<num_in; ++j) {
        arg1[j] = arg[j]+i*step_in_[j];
      }

      // Make the function outputs end up in temp_res
      for (int j=0; j<num_out; ++j) {
        res1[j] = (res[j]==0)? 0: temp_res;
        temp_res+= step_out_[j];
      }
      f_->evalD(arg1, res1, iw, w);

      // Sum results
      for (int k=0;k<num_out;++k) {
        if (res1[k] && sum[k])
          std::transform(res1[k], res1[k]+step_out_[k], sum[k], sum[k], std::plus<double>());
      }
    }
  }

  Function MapSumInternal
  ::getDerForward(const std::string& name, int nfwd, const Dict& opts) {

    Function df = f_.derForward(nfwd);

    std::vector<bool> repeat_in;
    repeat_in.insert(repeat_in.begin(), repeat_in_.begin(), repeat_in_.end());
    repeat_in.insert(repeat_in.end(), f_.nOut(), false);
    for (int i=0;i<nfwd;++i) {
      repeat_in.insert(repeat_in.end(), repeat_in_.begin(), repeat_in_.end());
    }

    MapSum ret("forward", df, n_, repeat_in, opts);
    ret.init();

    return ret;

  }

  Function MapSumInternal
  ::getDerReverse(const std::string& name, int nadj, const Dict& opts) {

    Function df = f_.derReverse(nadj);

    std::vector<MX> sym_in = symbolicInput();
    std::vector<MX> sym_out = symbolicOutput();
    std::vector<MX> ret_in;
    ret_in.insert(ret_in.end(), sym_in.begin(), sym_in.end());
    ret_in.insert(ret_in.end(), sym_out.begin(), sym_out.end());
    for (int i=0;i<nadj;++i) {
       std::vector<MX> sym_out = symbolicOutput();
       ret_in.insert(ret_in.end(), sym_out.begin(), sym_out.end());
    }

    std::vector<MX> ret_out = df.map(ret_in);

    // Outputs corresponding to repeat_in == true -> split and sum
    for (int j=0;j<f_.nIn();++j) {
      if (!repeat_in_[j]) {
        std::vector<MX> sym = std::vector<MX>(1, sym_in[j]);
        MXFunction id = MXFunction("identity", sym, sym);
        for (int i=0;i<nadj;++i) {
          ret_out[i*f_.nIn()+j] = id.mapsum(std::vector<MX>(1, ret_out[i*f_.nIn()+j]))[0];
        }
      }
    }

    MXFunction ret = MXFunction("reverse", ret_in, ret_out);

    return ret;

  }

  void MapSumInternal::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapSumInternal::generateBody(CodeGenerator& g) const {

    int num_in = f_.nIn(), num_out = f_.nOut();

    g.body << "  const real_t** arg1 = arg+" << f_.sz_arg() << ";" << endl
           << "  real_t** sum = res;" << endl;

    // Clear the accumulators
    for (int k=0;k<num_out;++k) {
      g.body << "  if (sum[" << k << "]!=0) " << g.fill_n(STRING("sum[" << k << "]"), step_out_[k], "0") << endl;
    };

    g.body << "  real_t** res1 = res+"  << f_.sz_res() <<  ";" << endl;

    g.body << "  int i;" << endl;
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;

    g.body << "    real_t* temp_res = w+"  << f_.sz_w() <<  ";" << endl
           << "    if (temp_res!=0)" << g.fill_n("temp_res", nnz_out_, "0") << endl;

    for (int j=0; j<num_in; ++j) {
      g.body << "    arg1[" << j << "] = arg[" << j << "]+i*" << step_in_[j] << ";" << endl;
    }
    for (int j=0; j<num_out; ++j) {
      g.body << "    res1[" << j << "] = (res[" << j << "]==0)? 0: temp_res;" << endl
             << "    temp_res+= " << step_out_[j] << ";" << endl;
    }

    g.body << "    " << g.call(f_, "arg1", "res1", "iw", "w") << ";" << endl;

    g.addAuxiliary(CodeGenerator::AUX_AXPY);
    // Sum results
    for (int k=0; k<num_out; ++k) {
      g.body << "    if (res1[" << k << "] && sum[" << k << "])" << endl
             << "       axpy(" << step_out_[k] << ",1,res1["<< k << "],1,sum[" << k << "],1);" << endl;
    }
    g.body << "  }" << std::endl;
  }

  inline string name(const Function& f) {
    if (f.isNull()) {
      return "NULL";
    } else {
      return f.getOption("name");
    }
  }

  void MapSumInternal::print(ostream &stream) const {
    stream << "MapSum(" << name(f_) << ", " << n_ << ")";
  }

} // namespace casadi
