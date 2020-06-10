//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"

namespace component {

class ACA {
public:
    struct Args {
        const std::vector<int> *fptr;
        const std::vector<int> *sptr;
    };
    using Functor = std::function<Complex (int, int)>;

    ACA();
    ~ACA();
    ACA(const ACA&) = delete;
    ACA& operator=(const ACA&) = delete;

public:
    void approximate(const Args& args, value_t threshold, Functor caller);
    Qcx_mat UMatrix() const;
    Qcx_mat VMatrix() const;
    size_t row() const { return row_; }
    size_t col() const { return col_; }
    size_t rank() const { return rank_; }
private:
    size_t searchJk(size_t row, const std::set<size_t>& jks) const;
    size_t searchIk(size_t col, const std::set<size_t>& iks) const;
    void setVk(size_t ik, Complex max_val);
    void setUk(size_t jk);
    void setRRow(size_t ik, int f, const std::vector<int>& srcs);
    void setRCol(size_t jk, int s, const std::vector<int>& flds);
    value_t midEpslion() const;
private:
    size_t row_, col_, rank_;
    Functor caller_;
    Qcx_mat U_, V_, R_;
};

} // namespace component