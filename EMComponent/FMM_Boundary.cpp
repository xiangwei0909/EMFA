#include "FMM_Boundary.h"
#include "Mesh.h"

using namespace component;

FMM_Boundary::FMM_Boundary()
{
}


FMM_Boundary::~FMM_Boundary()
{
}

void FMM_Boundary::init(std::shared_ptr<Mesh> pmesh, value_t lambda, value_t delta)
{
    pmesh->getBoundary(min_pos_, max_pos_);
    delta *= lambda;

    auto box = max_pos_ - min_pos_;
    value_t max_length = (box.x < box.y) ? (box.y < box.z ? box.z : box.y) : (box.x < box.z ? box.z : box.x);
    int num = static_cast<int>(max_length / delta + 0.5f);
    delta_ = max_length / num;

    value_t offset = 0.1f * delta_;
    min_pos_ = min_pos_ - VectorR3(offset, offset, offset);
    max_pos_ = max_pos_ + VectorR3(offset, offset, offset);
    box = max_pos_ - min_pos_;

    x_box_ = static_cast<int>((box.x + 0.9f * delta_) / delta_);
    y_box_ = static_cast<int>((box.y + 0.9f * delta_) / delta_);
    z_box_ = static_cast<int>((box.z + 0.9f * delta_) / delta_);
}

int FMM_Boundary::inWhichBox(const VectorR3 & mid) const
{
    int x = static_cast<int>((mid.x - min_pos_.x) / delta_);
    int y = static_cast<int>((mid.y - min_pos_.y) / delta_);
    int z = static_cast<int>((mid.z - min_pos_.z) / delta_);
    return (x * y_box_ + y) * z_box_ + z;
}

VectorR3 FMM_Boundary::getCenter(int idx) const
{
    int z = idx % z_box_;
    idx /= z_box_;
    int y = idx % y_box_;
    int x = idx / y_box_;

    return VectorR3((x + 0.5f) * delta_, (y + 0.5f) * delta_, (z + 0.5f) * delta_);
}

void FMM_Boundary::reportInfo(Qostream & strm) const
{
    auto oldFlag = strm.flags();
    strm << std::left << std::fixed << std::setprecision(2);
    strm << HEADING "FMM Boundary Information:\n" TRAILING
     << LEVEL1 "Minimum position: " << '(' << min_pos_.x << ',' << min_pos_.y << ',' << min_pos_.z << ")\n"
     << LEVEL1 "Maximum position: " << '(' << max_pos_.x << ',' << max_pos_.y << ',' << max_pos_.z << ")\n"
     << LEVEL1 "Distance between boxes: " << delta_ << "m\n"
     << LEVEL1 "X|Y|Z box number: " << x_box_ << '|' << y_box_ << '|' << z_box_ << '\n';
    strm << std::flush;
    strm.flags(oldFlag);
}