#include "EM.h"
#include "ConfigLoader.h"
#include "tools.h"
#include "Mathdef.h"

using namespace component;
using namespace mom;

EM::EM()
{
    tool::initTime();
    logger_ << std::left;
    time_logger_ << std::left;
}

EM::~EM()
{
    tool::totalTime(time_logger_);

    reportTimeInfo(Qcout);
    reportTimeInfo(logger_);
}

void EM::init(component::ConfigLoader * ploader)
{
    rtype_ = ploader->getResultType();
    incidence_ = ploader->getInc();

    value_t pol = incidence_.pole * DegreesToRadians;
    value_t theta = incidence_.theta * DegreesToRadians;
    value_t phi = incidence_.phi * DegreesToRadians;

    inc_k_ = VectorR3(-sin(theta) * cos(phi), -sin(theta) * sin(phi), -cos(theta));
    inc_e_ = VectorR3(-cos(pol)*cos(theta)*cos(phi)-sin(pol)*sin(phi), -cos(pol)*cos(theta)*sin(phi)+sin(pol)*cos(phi), cos(pol)*sin(theta));

	//reflection parameters
	r_inc_k_ = VectorR3(-sin(theta) * cos(phi), -sin(theta) * sin(phi), cos(theta));
	r_inc_e_ = VectorR3(-cos(theta)*cos(phi - PI - pol), -cos(theta)*sin(phi - PI - pol), sin(theta));
	///////////////////////////////////

    inc_h_ = inc_k_ * inc_e_;
	lamda = cc / incidence_.freq;
    if(rtype_ == policy::ResultType::Rad)
        radiation_ = ploader->getRad();
    else
        scattering_ = ploader->getSca();
	nearfield_ = ploader->getNearfield();
    auto mesh_file = ploader->getMeshPath();
    auto pos_suffix = mesh_file.find_last_of('.');
    if (pos_suffix == Qstring::npos)
        throw ConfigInvalid("mesh path parameter cannot match");

    // support two separators '/' and '\'
    auto pos1 = mesh_file.find_last_of('/');
    auto pos2 = mesh_file.find_last_of('\\');
    if (pos1 == Qstring::npos && pos2 == Qstring::npos)
    {
        dir_ = ".";
        folder_name_ = mesh_file.substr(0, pos_suffix);
    }
    else
    {
        auto pos_dir = (pos1 != Qstring::npos && pos2 != Qstring::npos) ? std::max(pos1, pos2) : std::min(pos1, pos2);
        dir_ = mesh_file.substr(0, pos_dir);
        folder_name_ = mesh_file.substr(pos_dir + 1, pos_suffix - pos_dir - 1);
    }
}

void EM::solve()
{
}

void EM::output()
{
}

void EM::clear()
{
}

void EM::reportInfo(Qostream & strm) const
{
}

void EM::reportTimeInfo(Qostream & strm)
{
    auto old_state = strm.flags();
    strm << HEADING "Time Information\n" TRAILING + time_logger_.str() << std::flush;
    strm.flags(old_state);
}


