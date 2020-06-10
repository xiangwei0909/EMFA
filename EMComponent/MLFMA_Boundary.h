//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "VectorR3.h"

namespace component {

class Mesh;

class MLFMA_Boundary {
public:
    MLFMA_Boundary();
    ~MLFMA_Boundary();

    void            init(std::shared_ptr<Mesh> pmesh, value_t lambda, value_t delta);
    const VectorR3& getMinPos() const { return min_pos_; }
    const VectorR3& getMaxPos() const { return max_pos_; }
    VectorR3        getCenter() const { return (min_pos_ + max_pos_) / 2; }
    int             getLevel() const { return level_; }
    value_t         getMaxLength() const { return max_pos_.x - min_pos_.x; }
    void            reportInfo(Qostream& strm) const;

private:
    int         level_;
    value_t     delta_;
    value_t     lambda_;
    VectorR3    min_pos_;
    VectorR3    max_pos_;
};


} // namespace component