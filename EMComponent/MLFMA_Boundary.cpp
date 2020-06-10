#include "MLFMA_Boundary.h"
#include "Mesh.h"

using namespace component;

MLFMA_Boundary::MLFMA_Boundary()
{
}


MLFMA_Boundary::~MLFMA_Boundary()
{
}

void MLFMA_Boundary::init(std::shared_ptr<Mesh> pmesh, value_t lambda, value_t delta)
{
    delta_ = delta * lambda;
    lambda_ = lambda;
    pmesh->getBoundary(min_pos_, max_pos_);
    auto center = (min_pos_ + max_pos_) / 2;

    value_t xlength = max_pos_.x - min_pos_.x;
    value_t ylength = max_pos_.y - min_pos_.y;
    value_t zlength = max_pos_.z - min_pos_.z;

    value_t max_length = std::max(std::max(xlength, ylength), zlength) + 0.1f * lambda;
    value_t half = max_length / 2;

    min_pos_ = center - VectorR3(half, half, half);
    max_pos_ = center + VectorR3(half, half, half);

    level_ = 0;
    for (; max_length > delta_; ++level_)
        max_length /= 2;
}



void MLFMA_Boundary::reportInfo(Qostream & strm) const
{
    value_t length = max_pos_.x - min_pos_.x;
    int box_num = 1;
    for (int i = 0; i < level_; ++i, box_num <<= 3)
        length /= 2;
    int old_flag = strm.flags();
    strm << std::fixed << std::setprecision(2);
    strm << HEADING "MLFMA Boundary Information:\n" TRAILING
        << LEVEL1 "Minimum position: " << '(' << min_pos_.x << ',' << min_pos_.y << ',' << min_pos_.z << ")\n"
        << LEVEL1 "Maximum position: " << '(' << max_pos_.x << ',' << max_pos_.y << ',' << max_pos_.z << ")\n"
        << LEVEL1 "MLFMA layers: " << level_ + 1 << '\n'
        << LEVEL1 "Distance between boxes: " << length << " (" << length / lambda_  << "*wavelength)\n"
        << LEVEL1 "Total box number: " << box_num << '\n';
    strm << std::flush;
    strm.flags(old_flag);
}
