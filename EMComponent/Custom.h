//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_USE_CXX11

#ifndef _DEBUG
#define ARMA_NO_DEBUG
#endif

#include <memory>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <assert.h>
#include <armadillo>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <numeric>
#include <cmath>
#include <chrono>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <deque>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>



#define     Qcout               std::cout
#define     Qcin                std::cin
#define     Qcerr               std::cerr

#define     HEADING             "============================================================\n"
#define     TRAILING            "------------------------------\n"
#define     LEVEL1              std::setw(30) << "  -> "
#define     LEVEL2              std::setw(30) << "    => "
#define     TIME_LOG(message)   tool::elapsedTime(time_logger_, "Time of " message)
#define     LOG(func, message)  do{ \
                                    Qcout << std::setw(30) << message ": "; \
                                    func; \
                                    Qcout << std::endl; \
                                }while(0)

typedef std::fstream            Qfstream;
typedef std::ifstream           Qifstream;
typedef std::ofstream           Qofstream;
typedef std::stringstream       Qsstream;
typedef std::ostream            Qostream;
typedef std::string             Qstring;
typedef float                   value_t;
typedef std::complex<value_t>   Complex;
typedef arma::Mat<value_t>      Qmat;
typedef arma::Col<value_t>      Qvec;
typedef arma::Mat<Complex>      Qcx_mat; 

typedef arma::subview<value_t>  Qsubview;
typedef arma::subview<Complex>  Qcx_subview;
typedef arma::Col<Complex>      Qcx_vec;
typedef arma::SpMat<Complex>    Qsp_cx_mat;
typedef arma::SpMat<value_t>    Qsp_mat;
typedef arma::Mat<size_t>       Qumat;

typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> CVector;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> CMatrix;
//typedef Eigen::Matrix<component::VectorC3, Eigen::Dynamic, Eigen::Dynamic> C3Matrix;
typedef Eigen::Matrix<value_t, Eigen::Dynamic, 1> VaVector;
typedef Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> VaMatrix;
typedef Eigen::SparseMatrix<Complex, Eigen::RowMajor> SpCMat_Row;
typedef Eigen::SparseMatrix<Complex> SpCMat;
typedef Eigen::Triplet<Complex> LocSpE;
//typedef Eigen::Matrix<VectorR3, Eigen::Dynamic, Eigen::Dynamic> VR3Matrix;
//typedef Eigen::Matrix<VectorC3, Eigen::Dynamic, Eigen::Dynamic> VC3Matrix;
//typedef Eigen::Matrix<VectorC3, Eigen::Dynamic, 1> VC3Vector;
typedef std::function<CVector(const CVector &)> MulVct;
//  desktop: wt
namespace wt {

//  sparse matrix template
template <typename T>
class QSpMat {
public:
    QSpMat() { }
    void set_size(size_t r, size_t c)
    {
        const size_t num = r * c;
        loc_.set_size(2, num);
        elems_.set_size(num);
    }
    void set_value(const size_t r, const size_t c, const size_t i, const T& val)
    {
        assert(i < elems_.n_elem);
        loc_(0, i) = r;
        loc_(1, i) = c;
        elems_(i) = val;
    }
    T& operator()(const size_t i)
    {
        assert(i < elems_.n_elem);
        return elems_(i);
    }
    const T& operator()(const size_t i) const
    {
        assert(i < elems_.n_elem);
        return elems_(i);
    }
    size_t row(const size_t i) const
    {
        assert(i < elems_.n_elem);
        return loc_(0, i);
    }
    size_t col(const size_t i) const
    {
        assert(i < elems_.n_elem);
        return loc_(1, i);
    }
    size_t size() const
    {
        return elems_.n_elem;
    }
private:
    Qumat           loc_;
    arma::Col<T>    elems_;
};

//  column-major ordering matrix template
template <typename T>
class QMat {
public:
    QMat() : row_(0), col_(0), vect_() { }
    void set_size(const size_t r, const size_t c)
    {
        row_ = r;
        col_ = c;
        vect_.resize(r * c);
    }
    void set_size(const size_t r, const size_t c, const T& val)
    {
        row_ = r;
        col_ = c;
        vect_.resize(r * c, val);
    }
    void reset(const T& val)
    {
        for (auto& elem : vect_)
            elem = val;
    }
    T& operator()(const size_t ir, const size_t ic)
    {
        assert(ic * row_ + ir < vect_.size());
        return vect_[ic * row_ + ir];
    }
    const T& operator()(const size_t ir, const size_t ic) const
    {
        assert(ic * row_ + ir < vect_.size());
        return vect_[ic * row_ + ir];
    }
    T& operator()(const size_t ith)
    {
        assert(ith < vect_.size());
        return vect_[ith];
    }
    const T& operator() (const size_t ith) const
    {
        assert(ith < vect_.size());
        return vect_[ith];
    }
    QMat& operator+=(const QMat& other)
    {
        assert(row_ == other.row_ && col_ == other.col_);
        std::transform(vect_.begin(), vect_.end(), other.vect_.begin(), 
                       vect_.begin(), std::plus<T>());
        return *this;
    }
    size_t row() const { return row_; }
    size_t col() const { return col_; }
    size_t size() const { return vect_.size(); };
    void clear() { vect_.clear(); }
private:
    size_t          row_, col_;
    std::vector<T>  vect_;
};

} // namespace wt

namespace policy {

enum BasisFunction { HALF_RWG = 0, FULL_RWG };
enum DDM { LESS_TIME = 2, LESS_MEMORY };
enum ResultType { BiRCS = 0, MonoRCS, Rad };

} // namespace category