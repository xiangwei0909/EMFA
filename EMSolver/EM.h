//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "EMBase.h"
#include "MyFWD.h"
#include "Custom.h"
#include "VectorR3.h"
#include "fieldStruct.h"
#include "Result.h"
#include "ThreadPool.h"

#define SEGMENT(message) do{ Qcout << HEADING message " Information\n" TRAILING; }while(0)

namespace mom {

using component::VectorR3;

class EM : public EMBase {
    RESULT_DECLARATION;
public:
    EM();
    ~EM();

public:
    virtual void init(component::ConfigLoader* ploader);
    virtual void solve();
    virtual void output();
    virtual void clear();
    virtual void reportInfo(Qostream& strm) const;
private:
    void reportTimeInfo(Qostream& strm);
protected:
    policy::ResultType rtype_;
    Qstring dir_;
	Qstring FEKOcur_path;
    Qstring folder_name_;

    Qofstream logger_;
    std::ostringstream time_logger_;

    component::Incidence incidence_;
    component::Scattering scattering_;
    component::Radiation radiation_;
	component::Nearfield nearfield_;

    VectorR3 inc_e_, inc_h_, inc_k_;
	VectorR3 r_inc_k_, r_inc_e_;
    VectorR3 sca_ev_, sca_eh_, sca_hv_, sca_hh_, sca_k_;

	value_t lamda;
};

} // namespace mom