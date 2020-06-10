#include "ACA.h"

using namespace component;

ACA::ACA()
{
}


ACA::~ACA()
{
}

void ACA::approximate(const Args & args, value_t threshold, Functor caller)
{
    auto& flds = *(args.fptr);
    auto& srcs = *(args.sptr);
    row_ = flds.size();
    col_ = srcs.size();
    rank_ = 0;
    caller_ = caller;
    auto min_rc = row_ < col_ ? row_ : col_;
    R_.zeros(row_, col_);
    U_.zeros(row_, min_rc);
    V_.zeros(min_rc, col_);

    std::set<size_t> iks, jks;

    size_t ik = 0, jk = 0;
    iks.insert(ik);
    for (size_t c = 0; c < col_; ++c)
        R_(ik, c) = caller_(flds[ik], srcs[c]);

    jk = searchJk(ik, jks);
    jks.insert(jk);
    setVk(ik, R_(ik, jk));
    for (size_t r = 0; r < row_; ++r)
        R_(r, jk) = caller_(flds[r], srcs[jk]);
    setUk(jk);
    
    auto vk_norm = arma::norm(V_.row(rank_));
    auto uk_norm = arma::norm(U_.col(rank_));
    auto epslion = uk_norm * vk_norm;
    auto fbs_imp = epslion;
    auto fbs_imp2 = fbs_imp * fbs_imp;
    ++rank_;
    while (rank_ < min_rc && fbs_imp * threshold < epslion)
    {
        ik = searchIk(jk, iks);
        iks.insert(ik);
        setRRow(ik, flds[ik], srcs);
        jk = searchJk(ik, jks);
        jks.insert(jk);
        setVk(ik, R_(ik, jk));
        setRCol(jk, srcs[jk], flds);
        setUk(jk);

        uk_norm = arma::norm(U_.col(rank_));
        vk_norm = arma::norm(V_.row(rank_));
        epslion = uk_norm * vk_norm;
        fbs_imp2 += 2 * midEpslion() + epslion * epslion;
        fbs_imp = sqrt(fbs_imp2);
        ++rank_;
    }
}

Qcx_mat ACA::UMatrix() const
{
    return U_.submat(0, 0, arma::size(row_, rank_));
}

Qcx_mat ACA::VMatrix() const
{
    return V_.submat(0, 0, arma::size(rank_, col_));
}

size_t ACA::searchJk(size_t row, const std::set<size_t>& jks) const
{
    size_t id = 0;
    while (id < col_ && jks.find(id) != jks.end())
        ++id;

    auto cur_max = norm(R_(row, id));
    for (size_t i = id + 1; i < col_; ++i)
    {
        if (cur_max < norm(R_(row, i)))
        {
            if (jks.find(i) != jks.end())
                continue;
            cur_max = norm(R_(row, i));
            id = i;
        }
    }
    return id;
}

size_t ACA::searchIk(size_t col, const std::set<size_t>& iks) const
{
    size_t id = 0;
    while (id < row_ && iks.find(id) != iks.end())
        ++id;
    auto cur_max = norm(R_(id, col));
    for (size_t i = id + 1; i < row_; ++i)
    {
        if (cur_max < norm(R_(i, col)))
        {
            if (iks.find(i) != iks.end())
                continue;
            cur_max = norm(R_(i, col));
            id = i;
        }
    }
    return id;
}

void ACA::setVk(size_t ik, Complex max_val)
{
    V_.row(rank_) = R_.row(ik) / max_val;
}

void ACA::setUk(size_t jk)
{
    U_.col(rank_) = R_.col(jk);
}

void ACA::setRRow(size_t ik, int f, const std::vector<int>& srcs)
{
    for (size_t c = 0; c < col_; ++c)
    {
        Complex imped = caller_(f, srcs[c]);
        Complex uv(0, 0);
        for (size_t rk = 0; rk < rank_; ++rk)
            uv += U_(ik, rk) * V_(rk, c);
        R_(ik, c) = imped - uv;
    }
}

void ACA::setRCol(size_t jk, int s, const std::vector<int>& flds)
{
    for (size_t r = 0; r < row_; ++r)
    {
        Complex imped = caller_(flds[r], s);
        Complex vu(0, 0);
        for (size_t rk = 0; rk < rank_; ++rk)
            vu += U_(r, rk) * V_(rk, jk);
        R_(r, jk) = imped - vu;
    }
}

value_t ACA::midEpslion() const
{
    value_t val = 0;
    for (size_t rk = 0; rk < rank_; ++rk)
    {
        Complex uu = arma::dot(U_.col(rk), U_.col(rank_));
        Complex vv = arma::dot(V_.row(rk), V_.row(rank_));
        val += norm(uu) * norm(vv);
    }
    return val;
}
