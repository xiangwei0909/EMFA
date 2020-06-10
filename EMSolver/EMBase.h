//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once

namespace mom {

class EMBase {
public:
    EMBase() { }
    virtual ~EMBase() { }

    EMBase(const EMBase&) = delete;
    EMBase& operator=(const EMBase&) = delete;
};

} // namespace mom