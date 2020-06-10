#include "ConfigLoader.h"
#include <cctype>


#define SET_FLAG(key)   do{ auto& result = dict_.at(key); result.has_reported = true; }while(0)
#define REPORT(report)  do{ report(Qcout); report(log_); }while(0)
#define CONFIG_LOG(key, report)    do{ SET_FLAG(key); REPORT(report); }while(0)

using std::setw;
using namespace component;

ConfigLoader::ConfigLoader()
{
    log_ << std::left << HEADING "Config Information\n" TRAILING;
    Qcout << std::left << HEADING "Config Information\n" TRAILING;
}


ConfigLoader::~ConfigLoader()
{
}

void ConfigLoader::load(const Qstring & filepath)
{
    Qifstream config(filepath, std::ios::in);
    if (config.fail())
        throw ConfigInvalid("fail opening config file(" + filepath + ')');
    auto max_len = std::numeric_limits<std::streamsize>::max();

    for (auto ch = config.peek(); std::isspace(ch) || ch == '#'; ch = config.peek())
    {
        if (ch == '#')
            config.ignore(max_len, '\n');
        else
            config.get();
    }
    Qstring key, value;
    while (config >> key >> value)
    {
        dict_[key] = Value(false, value);
        for (auto ch = config.peek(); std::isspace(ch) || ch == '#'; ch = config.peek())
        {
            if (ch == '#')
                config.ignore(max_len, '\n');
            else
                config.get();
        }
    }
    config.close();
}

Qstring ConfigLoader::getSolverType()
{
    const auto& val = dict_.at("EMSolver:");
    if (!val.has_reported)
        CONFIG_LOG("EMSolver:", reportSolverType);
    return val.value;
}

Incidence ConfigLoader::getInc()
{
    const auto& freq = dict_.at("frequency:");
	const auto& freq_from = dict_.at("frequency_from:");
	const auto& freq_to = dict_.at("frequency_to:");
	const auto& freq_delta = dict_.at("frequency_delta:");
	const auto& freq_inc = dict_.at("frequency_inc:");
    const auto& theta = dict_.at("theta:");
    const auto& phi = dict_.at("phi:");
    const auto& pole = dict_.at("pole:");
    if (!freq.has_reported)
    {
        SET_FLAG("frequency:");
        SET_FLAG("theta:");
        SET_FLAG("phi:");
        SET_FLAG("pole:");
        REPORT(reportInc);
    }
    Incidence inc;
    try {
        inc.freq = std::stof(freq.value);
		inc.freq_from = std::stof(freq_from.value);
		inc.freq_to = std::stof(freq_to.value);
		inc.freq_delta = std::stof(freq_delta.value);
		inc.freq_inc = std::stof(freq_inc.value);
        inc.theta = std::stof(theta.value);
        inc.phi = std::stof(phi.value);
        inc.pole = std::stof(pole.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing incidence information");
    }
    return inc;
}

policy::ResultType ConfigLoader::getResultType()
{
    const auto& val = dict_.at("result_type:");
    if (!val.has_reported)
        CONFIG_LOG("result_type:", reportResultType);
    if (val.value == "BistaticRCS")
        return policy::BiRCS;
    if (val.value == "MonostaticRCS")
        return policy::MonoRCS;
    if (val.value == "Radiation")
        return policy::Rad;
    throw ConfigInvalid("type parameter only supports three type values (BistaticRCS, MonostaticRCS and Radiation)");
}

Scattering ConfigLoader::getSca()
{
    const auto& theta_from = dict_.at("theta_from:");
    const auto& theta_to = dict_.at("theta_to:");
    const auto& theta_delta = dict_.at("theta_delta:");
    const auto& phi_from = dict_.at("phi_from:");
    const auto& phi_to = dict_.at("phi_to:");
    const auto& phi_delta = dict_.at("phi_delta:");
    if (!theta_from.has_reported)
    {
        SET_FLAG("theta_from:");
        SET_FLAG("theta_to:");
        SET_FLAG("theta_delta:");
        SET_FLAG("phi_from:");
        SET_FLAG("phi_to:");
        SET_FLAG("phi_delta:");
        REPORT(reportSca);
    }
    Scattering sca;
    try {
        sca.theta_from = std::stof(theta_from.value);
        sca.theta_to = std::stof(theta_to.value);
        sca.theta_delta = std::stof(theta_delta.value);
        sca.phi_from = std::stof(phi_from.value);
        sca.phi_to = std::stof(phi_to.value);
        sca.phi_delta = std::stof(phi_delta.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing scattering information");
    }
 
   return sca;
}

Radiation ConfigLoader::getRad()
{
    const auto& e_plane = dict_.at("E_Plane:");
    const auto& h_plane = dict_.at("H_Plane:");
    const auto& hh_plane = dict_.at("HH_Plane:");
    if (!e_plane.has_reported)
    {
        SET_FLAG("E_Plane:");
        SET_FLAG("H_Plane:");
        SET_FLAG("HH_Plane:");
        REPORT(reportRad);
    }
    Radiation rad;
    rad.e_plane = e_plane.value == "true" ? true : false;
    rad.h_plane = h_plane.value == "true" ? true : false;
    rad.hh_plane = hh_plane.value == "true" ? true : false;
    return rad;
}

Nearfield ConfigLoader::getNearfield()
{
	Nearfield NF;

	const auto& para1 = dict_.at("Origin_point:");
	const Qstring &_ori = para1.value;
	auto pos1 = _ori.find_first_of('(');
	auto pos2 = _ori.find_first_of(',');
	auto pos3 = _ori.find_last_of(',');
	auto pos4 = _ori.find_last_of(')');
	value_t x, y, z;
	try {
		x = std::stof(_ori.substr(pos1 + 1, pos2 - pos1 - 1));
		y = std::stof(_ori.substr(pos2 + 1, pos3 - pos2 - 1));
		z = std::stof(_ori.substr(pos3 + 1, pos4 - pos3 - 1));
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Origin point");
	}
	VectorR3 origin(x, y, z);
	NF.origin_point = origin;

	const auto& para2 = dict_.at("End_point:");
	const Qstring &_end = para2.value;
	pos1 = _end.find_first_of('(');
	pos2 = _end.find_first_of(',');
	pos3 = _end.find_last_of(',');
	pos4 = _end.find_last_of(')');
	try {
		x = std::stof(_end.substr(pos1 + 1, pos2 - pos1 - 1));
		y = std::stof(_end.substr(pos2 + 1, pos3 - pos2 - 1));
		z = std::stof(_end.substr(pos3 + 1, pos4 - pos3 - 1));
	}
	catch (...) {
		throw ConfigInvalid("fail parsing end point");
	}
	VectorR3 end(x, y, z);
	NF.end_point = end;

	const auto& para3 = dict_.at("Sampling_x:");
	const auto& para4 = dict_.at("Sampling_y:");
	const auto& para5 = dict_.at("Sampling_z:");
	try {
		NF.sampling_x = std::stoi(para3.value);
		NF.sampling_y = std::stoi(para4.value);
		NF.sampling_z = std::stoi(para5.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing sampling number");
	}

	if (!para1.has_reported)
	{
		REPORT(reportNearfield);
	}

	return NF;
}

MultipleIncidence component::ConfigLoader::getMultiInc()
{
	const auto& PW_theta_from = dict_.at("PW_theta_from:");
	const auto& PW_theta_delta = dict_.at("PW_theta_delta:");
	const auto& PW_theta_to = dict_.at("PW_theta_to:");
	const auto& PW_phi_from = dict_.at("PW_phi_from:");
	const auto& PW_phi_delta = dict_.at("PW_phi_delta:");
	const auto& PW_phi_to = dict_.at("PW_phi_to:");
	const auto& polarization = dict_.at("polarization:");

	if (!PW_theta_from.has_reported)
	{
		SET_FLAG("PW_theta_from:");
		SET_FLAG("PW_theta_delta:");
		SET_FLAG("PW_theta_to:");
		SET_FLAG("PW_phi_from:");
		SET_FLAG("PW_phi_delta:");
		SET_FLAG("PW_phi_to:");
		SET_FLAG("polarization:");
		REPORT(reportMultiInc);
	}

	MultipleIncidence MultiInc;
	try {
		MultiInc.PW_theta_from = std::stof(PW_theta_from.value);
		MultiInc.PW_theta_delta = std::stof(PW_theta_delta.value);
		MultiInc.PW_theta_to = std::stof(PW_theta_to.value);
		MultiInc.PW_phi_from = std::stof(PW_phi_from.value);
		MultiInc.PW_phi_delta = std::stof(PW_phi_delta.value);
		MultiInc.PW_phi_to = std::stof(PW_phi_to.value);
		MultiInc.polarization = polarization.value;
	}
	catch (...) {
		throw ConfigInvalid("fail parsing multiple PWs information");
	}

	return MultiInc;
}

Qstring ConfigLoader::getMeshPath()
{
    const auto& mesh_path = dict_.at("mesh_path:");
    if (!mesh_path.has_reported)
        CONFIG_LOG("mesh_path:", reportMeshPath);
    return mesh_path.value;
}

Qstring ConfigLoader::getSurroundPath()
{
	const auto& mesh_path = dict_.at("surround_path:");
	if (!mesh_path.has_reported)
		CONFIG_LOG("surround_path:", reportMeshPath);
	return mesh_path.value;
}

Qstring ConfigLoader::getFEKOcurPath()
{
	const auto& mesh_path = dict_.at("FEKOcur_path:");
	if (!mesh_path.has_reported)
		CONFIG_LOG("FEKOcur_path:", reportMeshPath);
	return mesh_path.value;
}

value_t ConfigLoader::getAlpha()
{
    const auto& alpha = dict_.at("alpha:");
    if (!alpha.has_reported)
        CONFIG_LOG("alpha:", reportAlpha);
    value_t res = 0.0f;
    try {
        res = std::stof(alpha.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing combined coefficient");
    }
    return res;
}

value_t ConfigLoader::getACABoxLength()
{
    const auto& para = dict_.at("aca_box_length:");
    if (!para.has_reported)
        CONFIG_LOG("aca_box_length:", reportACABoxLength);
    value_t box_len = 0.0f;
    try {
        box_len = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing aca box length");
    }
    return box_len;
}

value_t ConfigLoader::getACAThreshold()
{
    const auto& threshold = dict_.at("aca_threshold:");
    if (!threshold.has_reported)
        CONFIG_LOG("aca_threshold:", reportACAThreshold);
    value_t aca_threshold = 0.0f;
    try {
        aca_threshold = std::stof(threshold.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing aca threshold");
    }
    return aca_threshold;
}

bool ConfigLoader::getEDMSwitch()
{
    const auto& para = dict_.at("edm_acceleration:");
    if (!para.has_reported)
        CONFIG_LOG("edm_acceleration:", reportEDMSwitch);
    if (para.value == "true")
        return true;
    if (para.value == "false")
        return false;
    throw ConfigInvalid("edm acceleration option only supports two values (true or false)");
}

int ConfigLoader::getCBFMaxBoxNum()
{
    const auto& box_num = dict_.at("max_box_num:");
    if (!box_num.has_reported)
        CONFIG_LOG("max_box_num:", reportCBFMaxBoxNum);
    int max_box = 0;
    try {
        max_box = std::stoi(box_num.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing CBF max box number");
    }
    return max_box;
}

int ConfigLoader::getCBFPolarization()
{
    const auto& polar = dict_.at("polarization:");
    if (!polar.has_reported)
        CONFIG_LOG("polarization:", reportCBFPolarization);
    if (polar.value == "1")
        return 1;
    if (polar.value == "2")
        return 2;
    throw ConfigInvalid("CBF polarization only has two type (1 or 2)");
}

int ConfigLoader::getCBFSampleSpacing()
{
    const auto& spacing = dict_.at("sample_spacing:");
    if (!spacing.has_reported)
        CONFIG_LOG("sample_spacing:", reportCBFSampleSpacing);
    int sample = 0;
    try {
        sample = std::stoi(spacing.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing CBF sample spacing");
    }
    return sample;
}

value_t ConfigLoader::getEpsilon1()
{
    const auto& para = dict_.at("epsilon1:");
    if (!para.has_reported)
        CONFIG_LOG("epsilon1:", reportEpsilon1);
    value_t eps1 = 0.0f;
    try {
        eps1 = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing epsilon1");
    }
    return eps1;
}

value_t ConfigLoader::getEpsilon2()
{
    const auto& para = dict_.at("epsilon2:");
    if (!para.has_reported)
        CONFIG_LOG("epsilon2:", reportEpsilon2);
    value_t eps2 = 0.0f;
    try {
        eps2 = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing epsilon2");
    }
    return eps2;
}

value_t ConfigLoader::getMu1()
{
    const auto& para = dict_.at("mu1:");
    if (!para.has_reported)
        CONFIG_LOG("mu1:", reportMu1);
    value_t mu1 = 0.0f;
    try {
        mu1 = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing mu1");
    }
    return mu1;
}

value_t ConfigLoader::getMu2()
{
    const auto& para = dict_.at("mu2:");
    if (!para.has_reported)
        CONFIG_LOG("mu2:", reportMu2);
    value_t mu2 = 0.0f;
    try {
        mu2 = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing mu2");
    }
    return mu2;
}

value_t ConfigLoader::getAIMXSpacing()
{
    const auto& spacing = dict_.at("x_spacing:");
    if (!spacing.has_reported)
        CONFIG_LOG("x_spacing:", reportAIMXSpacing);
    value_t x_spacing = 0.0f;
    try {
        x_spacing = std::stof(spacing.value);
    }
    catch (...)
    {
        throw ConfigInvalid("fail parsing AIM x axis spacing");
    }
    return x_spacing;
}

value_t ConfigLoader::getAIMYSpacing()
{
    const auto& spacing = dict_.at("y_spacing:");
    if (!spacing.has_reported)
        CONFIG_LOG("y_spacing:", reportAIMYSpacing);
    value_t y_spacing = 0.0f;
    try {
        y_spacing = std::stof(spacing.value);
    }
    catch (...)
    {
        throw ConfigInvalid("fail parsing AIM y axis spacing");
    }
    return y_spacing;
}

value_t ConfigLoader::getAIMZSpacing()
{
    const auto& spacing = dict_.at("z_spacing:");
    if (!spacing.has_reported)
        CONFIG_LOG("z_spacing:", reportAIMZSpacing);
    value_t z_spacing = 0.0f;
    try {
        z_spacing = std::stof(spacing.value);
    }
    catch (...)
    {
        throw ConfigInvalid("fail parsing AIM z axis spacing");
    }
    return z_spacing;
}

value_t ConfigLoader::getNearThreshold()
{
    const auto& threshold = dict_.at("near_threshold:");
    if (!threshold.has_reported)
        CONFIG_LOG("near_threshold:", reportNearThreshold);
    value_t near_threshold = 0.0f;
    try {
        near_threshold = std::stof(threshold.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing AIM near threshold");
    }
    return near_threshold;
}

value_t ConfigLoader::getIterationThreshold()
{
    const auto& threshold = dict_.at("iteration_threshold:");
    if (!threshold.has_reported)
        CONFIG_LOG("iteration_threshold:", reportIterationThreshold);
    value_t iteration_threshold = 0.0f;
    try {
        iteration_threshold = std::stof(threshold.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing iteration threshold");
    }
    return iteration_threshold;
}

int ConfigLoader::getMaxIterationNum()
{
    const auto& num = dict_.at("max_iteration_number:");
    if (!num.has_reported)
        CONFIG_LOG("max_iteration_number:", reportMaxIterationNum);
    int iter_num = 0;
    try {
        iter_num = std::stoi(num.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing maximum iteration number");
    }
    return iter_num;
}

value_t ConfigLoader::getFMMBoxLength()
{
    const auto& para = dict_.at("fmm_box_length:");
    if (!para.has_reported)
        CONFIG_LOG("fmm_box_length:", reportFMMBoxLength);
    value_t box_len = 0.0f;
    try {
        box_len = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing FMM box length");
    }
    return box_len;
}

value_t ConfigLoader::getMLFMABoxThreshold()
{
    const auto& para = dict_.at("mlfma_box_threshold:");
    if (!para.has_reported)
        CONFIG_LOG("mlfma_box_threshold:", reportMLFMABoxThreshold);
    value_t box_threshold = 0.0f;
    try {
        box_threshold = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing MLFMA box threshold");
    }
    return box_threshold;
}

policy::BasisFunction ConfigLoader::getIEDGBasisFunctionType()
{
    const auto& type = dict_.at("basis_function_type:");
    if (!type.has_reported)
        CONFIG_LOG("basis_function_type:", reportIEDGBasisFunctionType);
    if (type.value == "full_RWG")
        return policy::FULL_RWG;
    if (type.value == "half_RWG")
        return policy::HALF_RWG;
    throw ConfigInvalid("IEDG basis function only supports two types (half_RWG or full_RWG)");
}

value_t ConfigLoader::getStabilizationFactor()
{
    const auto& para = dict_.at("stabilization_factor:");
    if (!para.has_reported)
        CONFIG_LOG("stabilization_factor:", reportStabilizationFactor);
    value_t factor = 0.0f;
    try {
        factor = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing IEDG stabilization factor");
    }
    return factor;
}

value_t ConfigLoader::getAverageSize()
{
    const auto& para = dict_.at("average_size:");
    if (!para.has_reported)
        CONFIG_LOG("average_size:", reportAverageSize);
    value_t average_size = 0.0f;
    try {
        average_size = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing average size");
    }
    return average_size;
}

policy::DDM ConfigLoader::getDDMPolicy()
{
    const auto& para = dict_.at("ddm_policy:");
    if (!para.has_reported)
        CONFIG_LOG("ddm_policy:", reportDDMPolicy);
    if (para.value == "less_time")
        return policy::DDM::LESS_TIME;
    if (para.value == "less_memory")
        return policy::DDM::LESS_MEMORY;
    throw ConfigInvalid("DDM policy only supports two options (less_time or less_memory)");
}

bool ConfigLoader::getPreconditioningSwitch()
{
    const auto& para = dict_.at("preconditioning:");
    if (!para.has_reported)
        CONFIG_LOG("preconditioning:", reportPreconditioningSwitch);
    if (para.value == "true")
        return true;
    if (para.value == "false")
        return false;
    throw ConfigInvalid("Precondition switch only supports two options (true or false)");
}

value_t ConfigLoader::getRowThreshold()
{
    const auto& para = dict_.at("row_threshold:");
    if (!para.has_reported)
        CONFIG_LOG("row_threshold:", reportRowThreshold);
    value_t row_threshold = 0.0f;
    try {
        row_threshold = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing row threshold");
    }
    if (row_threshold < 0.87f)
        throw ConfigInvalid("row threshold parameter cannot less than 0.87");
    return row_threshold;
}

value_t ConfigLoader::getColThreshold()
{
    const auto& para = dict_.at("col_threshold:");
    if (!para.has_reported)
        CONFIG_LOG("col_threshold:", reportColThreshold);
    value_t col_threshold = 0.0f;
    try {
        col_threshold = std::stof(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing column threshold");
    }
    if (col_threshold < 0.87f)
        throw ConfigInvalid("column threshold parameter cannot less than 0.87");
    return col_threshold;
}

size_t ConfigLoader::getThreadNumber()
{
    const auto& para = dict_.at("thread_number:");
    if (!para.has_reported)
        CONFIG_LOG("thread_number:", reportThreadNumber);
    size_t threads = 0;
    try {
        threads = std::stoull(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing thread number");
    }
    return threads;
}

size_t ConfigLoader::getTaskFactor()
{
    const auto& para = dict_.at("task_factor:");
    if (!para.has_reported)
        CONFIG_LOG("task_factor:", reportTaskFactor);
    size_t factor = 0;
    try {
        factor = std::stoull(para.value);
    }
    catch (...) {
        throw ConfigInvalid("fail parsing task factor");
    }
    return factor;
}

int ConfigLoader::getIsfast()
{
	const auto& para = dict_.at("Isfast:");
	if (!para.has_reported)
		CONFIG_LOG("Isfast:", reportIsfast);
	int _Isfast = 0;
	try {
		_Isfast = std::stoi(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Dx");
	}
	return _Isfast;
}

value_t ConfigLoader::getDx()
{
	const auto& para = dict_.at("Dx:");
	if (!para.has_reported)
		CONFIG_LOG("Dx:", reportDx);
	value_t Dx = 0;
	try {
		Dx = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Dx");
	}
	return Dx;
}

value_t ConfigLoader::getDy()
{
	const auto& para = dict_.at("Dy:");
	if (!para.has_reported)
		CONFIG_LOG("Dy:", reportDy);
	value_t Dy = 0;
	try {
		Dy = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Dy");
	}
	return Dy;
}

value_t ConfigLoader::getD_Angle()
{
	const auto& para = dict_.at("D_Angle:");
	if (!para.has_reported)
		CONFIG_LOG("D_Angle:", reportD_Angle);
	value_t D_Angle = 0;
	try {
		D_Angle = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing D_Angle");
	}
	return D_Angle;
}

value_t ConfigLoader::getPhase_0()
{
	const auto& para = dict_.at("Phase_0:");
	if (!para.has_reported)
		CONFIG_LOG("Phase_0:", reportPhase_0);
	value_t Phase_0 = 0;
	try {
		Phase_0 = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Phase_0");
	}
	return Phase_0;
}

value_t ConfigLoader::getPhase_x()
{
	const auto& para = dict_.at("Phase_x:");
	if (!para.has_reported)
		CONFIG_LOG("Phase_x:", reportPhase_x);
	value_t Phase_x = 0;
	try {
		Phase_x = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Phase_x");
	}
	return Phase_x;
}

value_t ConfigLoader::getPhase_y()
{
	const auto& para = dict_.at("Phase_y:");
	if (!para.has_reported)
		CONFIG_LOG("Phase_y:", reportPhase_y);
	value_t Phase_y = 0;
	try {
		Phase_y = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Phase_y");
	}
	return Phase_y;
}

value_t ConfigLoader::getA_Angle()
{
	const auto& para = dict_.at("A_Angle:");
	if (!para.has_reported)
		CONFIG_LOG("A_Angle:", reportA_Angle);
	value_t A_Angle = 0;
	try {
		A_Angle = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing A_Angle");
	}
	return A_Angle;
}

int ConfigLoader::getArray_x()
{
	const auto& para = dict_.at("Array_x:");
	if (!para.has_reported)
		CONFIG_LOG("Array_x:", reportArray_x);
	int Array_x = 0;
	try {
		Array_x = std::stoi(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Array_x");
	}
	return Array_x;
}

int ConfigLoader::getArray_y()
{
	const auto& para = dict_.at("Array_y:");
	if (!para.has_reported)
		CONFIG_LOG("Array_y:", reportArray_y);
	int Array_y = 0;
	try {
		Array_y = std::stoi(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Array_y");
	}
	return Array_y;
}

int ConfigLoader::gett_sum()
{
	const auto& para = dict_.at("t_sum:");
	if (!para.has_reported)
		CONFIG_LOG("t_sum:", reportt_sum);
	int t_sum = 0;
	try {
		t_sum = std::stoi(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing t_sum");
	}
	return t_sum;
}

value_t ConfigLoader::getEps1()
{
	const auto& para = dict_.at("eps1:");
	if (!para.has_reported)
		CONFIG_LOG("eps1:", reportEps1);
	value_t eps1 = 0.0f;
	try {
		eps1 = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing eps1");
	}
	return eps1;
}

value_t ConfigLoader::getMu()
{
	const auto& para = dict_.at("mu:");
	if (!para.has_reported)
		CONFIG_LOG("mu:", reportMu);
	value_t mu1 = 0.0f;
	try {
		mu1 = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing mu");
	}
	return mu1;
}

Complex ConfigLoader::getSigma()
{
	const auto& para = dict_.at("Sigma:");
	if (!para.has_reported)
		CONFIG_LOG("Sigma:", reportSigma);
	const Qstring &_sigma = para.value;
	auto pos1 = _sigma.find_first_of('(');
	auto pos2 = _sigma.find_first_of(',');
	auto pos3 = _sigma.find_first_of(')');
	value_t r, i;
	try {
		r = std::stof(_sigma.substr(pos1+1,pos2-pos1-1));
		i = std::stof(_sigma.substr(pos2+1, pos3 - pos2 - 1));
	}
	catch (...) {
		throw ConfigInvalid("fail parsing sigma");
	}
	Complex sigma(r, i);
	return sigma;
}

value_t ConfigLoader::getScale_alongx()
{
	const auto& para = dict_.at("Scale_alongx:");
	if (!para.has_reported)
		CONFIG_LOG("Scale_alongx:", reportScale_alongx);
	value_t scale_alongx = 0.0f;
	try {
		scale_alongx = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Scale_alongx");
	}
	return scale_alongx;
}

value_t ConfigLoader::getScale_alongy()
{
	const auto& para = dict_.at("Scale_alongy:");
	if (!para.has_reported)
		CONFIG_LOG("Scale_alongy:", reportScale_alongy);
	value_t scale_alongy = 0.0f;
	try {
		scale_alongy = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Scale_alongy");
	}
	return scale_alongy;
}

value_t ConfigLoader::getRotate_alongx()
{
	const auto& para = dict_.at("Rotate_alongx:");
	if (!para.has_reported)
		CONFIG_LOG("Rotate_alongx:", reportRotate_alongx);
	value_t Rotate_alongx = 0.0f;
	try {
		Rotate_alongx = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Rotate_alongx");
	}
	return Rotate_alongx;
}

value_t ConfigLoader::getRotate_alongy()
{
	const auto& para = dict_.at("Rotate_alongy:");
	if (!para.has_reported)
		CONFIG_LOG("Rotate_alongy:", reportRotate_alongy);
	value_t Rotate_alongy = 0.0f;
	try {
		Rotate_alongy = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Rotate_alongy");
	}
	return Rotate_alongy;
}

int ConfigLoader::getisContinuous()
{
	const auto& para = dict_.at("isContinuous:");
	if (!para.has_reported)
		CONFIG_LOG("isContinuous:", reportisContinuous);
	int _isCon = 0;
	try {
		_isCon = std::stoi(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing isContinuous");
	}
	return _isCon;
}

value_t ConfigLoader::getd_x()
{
	const auto& para = dict_.at("d_x:");
	if (!para.has_reported)
		CONFIG_LOG("d_x:", reportd_x);
	value_t _d_x;
	try {
		_d_x = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing d_x");
	}
	return _d_x;
}

value_t ConfigLoader::getd_y()
{
	const auto& para = dict_.at("d_y:");
	if (!para.has_reported)
		CONFIG_LOG("d_y:", reportd_y);
	value_t _d_y;
	try {
		_d_y = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing d_y");
	}
	return _d_y;
}

value_t ConfigLoader::getd_z()
{
	const auto& para = dict_.at("d_z:");
	if (!para.has_reported)
		CONFIG_LOG("d_z:", reportd_z);
	value_t _d_z;
	try {
		_d_z = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing d_z");
	}
	return _d_z;
}

value_t ConfigLoader::getThickness_sheet()
{
	const auto& para = dict_.at("Thickness_sheet:");
	if (!para.has_reported)
		CONFIG_LOG("Thickness_sheet:", reportThickness_sheet);
	value_t _Thickness_sheet;
	try{
		_Thickness_sheet = std::stof(para.value);
	}
	catch (...) {
		throw ConfigInvalid("fail parsing Thickness_sheet");
	}
	return _Thickness_sheet;
}
////////////////////////////////////////////////////////////////////////
void ConfigLoader::Debug() const
{
    Qcout << "Config Debug Information\n";
    for (auto cur = dict_.begin(); cur != dict_.end(); ++cur)
        Qcout << cur->first << " " << cur->second.value << '\n';
    Qcout.flush();
}
////////////////////////////////////////////////////////////////////////

void ConfigLoader::reportSolverType(Qostream & strm) const
{
    strm << LEVEL1 "EM Solver: " << dict_.at("EMSolver:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportInc(Qostream & strm) const
{
    strm << LEVEL1 "Incidence: " << '\n'
        << LEVEL2 "frequency: " << dict_.at("frequency:").value << " (Hz)\n"
        << LEVEL2 "theta: " << dict_.at("theta:").value << " (degree)\n"
        << LEVEL2 "phi: " << dict_.at("phi:").value << " (degree)\n"
        << LEVEL2 "pole: " << dict_.at("pole:").value << " (degree)\n";
    strm.flush();
}

void ConfigLoader::reportResultType(Qostream & strm) const
{
    strm << LEVEL1 "Result type: " << dict_.at("result_type:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportSca(Qostream & strm) const
{
    strm << LEVEL1 "Scattering: " << '\n'
        << LEVEL2 "theta from: " << dict_.at("theta_from:").value << " (degree)\n"
        << LEVEL2 "theta to: " << dict_.at("theta_to:").value << " (degree)\n"
        << LEVEL2 "theta delta: " << dict_.at("theta_delta:").value << " (degree)\n"
        << LEVEL2 "phi from: " << dict_.at("phi_from:").value << " (degree)\n"
        << LEVEL2 "phi to: " << dict_.at("phi_to:").value << " (degree)\n"
        << LEVEL2 "phi delta: " << dict_.at("phi_delta:").value << " (degree)\n";
    strm.flush();
}

void ConfigLoader::reportRad(Qostream & strm) const
{
    strm << LEVEL1 "Radiation: " << '\n'
        << LEVEL2 "XZ plane: " << dict_.at("E_Plane:").value << '\n'
        << LEVEL2 "YZ plane: " << dict_.at("H_Plane:").value << '\n'
        << LEVEL2 "XY plane: " << dict_.at("HH_Plane:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportNearfield(Qostream &strm) const
{
	strm << LEVEL1 "Near field:" << '\n'
		<< LEVEL2 "Origin point:" << dict_.at("Origin_point:").value << '\n'
		<< LEVEL2 "End point:" << dict_.at("End_point:").value << '\n'
		<< LEVEL2 "Sampling x:" << dict_.at("Sampling_x:").value << '\n'
		<< LEVEL2 "Sampling y:" << dict_.at("Sampling_y:").value << '\n'
		<< LEVEL2 "Sampling z:" << dict_.at("Sampling_z:").value << '\n';
	strm.flush();
}

void component::ConfigLoader::reportMultiInc(Qostream & strm) const
{
	strm << LEVEL1 "Multiple Incidence Setting: " << '\n'
		<< LEVEL2 "PW_theta_from: " << dict_.at("PW_theta_from:").value << '\n'
		<< LEVEL2 "PW_theta_delta: " << dict_.at("PW_theta_delta:").value << '\n'
		<< LEVEL2 "PW_theta_to: " << dict_.at("PW_theta_to:").value << '\n'
		<< LEVEL2 "PW_phi_from: " << dict_.at("PW_phi_from:").value << '\n'
		<< LEVEL2 "PW_phi_delta: " << dict_.at("PW_phi_delta:").value << '\n'
		<< LEVEL2 "PW_phi_to: " << dict_.at("PW_phi_to:").value << '\n'
		<< LEVEL2 "polarization: " << dict_.at("polarization:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportMeshPath(Qostream & strm) const
{
    strm << LEVEL1 "Mesh path: " << dict_.at("mesh_path:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportFEKOcurPath(Qostream & strm) const
{
	strm << LEVEL1 "FEKOcur path: " << dict_.at("FEKOcur_path:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportSurroundPath(Qostream & strm) const
{
	strm << LEVEL1 "Surround path: " << dict_.at("surround_path:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportAlpha(Qostream & strm) const
{
    strm << LEVEL1 "Combined coefficient: " << dict_.at("alpha:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportACABoxLength(Qostream & strm) const
{
    strm << LEVEL1 "ACA box length: " << dict_.at("aca_box_length:").value << " (wavelength)\n";
    strm.flush();
}

void ConfigLoader::reportACAThreshold(Qostream & strm) const
{
    strm << LEVEL1 "ACA threshold: " << dict_.at("aca_threshold:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportEDMSwitch(Qostream & strm) const
{
    strm << LEVEL1 "EDM acceleration: " << dict_.at("edm_acceleration:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportCBFMaxBoxNum(Qostream & strm) const
{
    strm << LEVEL1 "CBF max box number: " << dict_.at("max_box_num:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportCBFPolarization(Qostream & strm) const
{
    strm << LEVEL1 "CBF polarization: " << (dict_.at("polarization:").value == "1" ? "one" : "two") << '\n';
    strm.flush();
}

void ConfigLoader::reportCBFSampleSpacing(Qostream & strm) const
{
    strm << LEVEL1 "CBF sample spacing: " << dict_.at("sample_spacing:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportEpsilon1(Qostream & strm) const
{
    strm << LEVEL1 "PMCHW Epsilon1: " << dict_.at("epsilon1:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportEpsilon2(Qostream & strm) const
{
    strm << LEVEL1 "PMCHW Epsilon2: " << dict_.at("epsilon2:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportMu1(Qostream & strm) const
{
    strm << LEVEL1 "PMCHW Mu1: " << dict_.at("mu1:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportMu2(Qostream & strm) const
{
    strm << LEVEL1 "PMCHW Mu2: " << dict_.at("mu2:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportAIMXSpacing(Qostream & strm) const
{
    strm << LEVEL1 "AIM x axis spacing: " << dict_.at("x_spacing:").value << " (wavelength)\n";
    strm.flush();
}

void ConfigLoader::reportAIMYSpacing(Qostream & strm) const
{
    strm << LEVEL1 "AIM y axis spacing: " << dict_.at("y_spacing:").value << " (wavelength)\n";
    strm.flush();
}

void ConfigLoader::reportAIMZSpacing(Qostream & strm) const
{
    strm << LEVEL1 "AIM z axis spacing: " << dict_.at("z_spacing:").value << " (wavelength)\n";
    strm.flush();
}

void ConfigLoader::reportNearThreshold(Qostream & strm) const
{
    strm << LEVEL1 "AIM near threshold: " << dict_.at("near_threshold:").value << " (wavelength)\n";
    strm.flush();
}

void ConfigLoader::reportIterationThreshold(Qostream & strm) const
{
    strm << LEVEL1 "Iteration threshold: " << dict_.at("iteration_threshold:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportMaxIterationNum(Qostream & strm) const
{
    strm << LEVEL1 "Max iteration number: " << dict_.at("max_iteration_number:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportFMMBoxLength(Qostream & strm) const
{
    strm << LEVEL1 "FMM box length: " << dict_.at("fmm_box_length:").value << " (wavelength)\n";
    strm.flush();
}

void ConfigLoader::reportMLFMABoxThreshold(Qostream & strm) const
{
    strm << LEVEL1 "MLFMA box threshold: " << dict_.at("mlfma_box_threshold:").value << " (wavelength)\n";
    strm.flush();
}

void ConfigLoader::reportIEDGBasisFunctionType(Qostream & strm) const
{
    strm << LEVEL1 "IEDG basis function type: " << dict_.at("basis_function_type:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportStabilizationFactor(Qostream & strm) const
{
    strm << LEVEL1 "Stabilization factor: " << dict_.at("stabilization_factor:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportDDMPolicy(Qostream & strm) const
{
    strm << LEVEL1 "DDM policy: " << dict_.at("ddm_policy:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportPreconditioningSwitch(Qostream & strm) const
{
    strm << LEVEL1 "Preconditioning: " << dict_.at("preconditioning:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportAverageSize(Qostream & strm) const
{
    strm << LEVEL1 "Average size: " << dict_.at("average_size:").value << " (m)\n";
    strm.flush();
}

void ConfigLoader::reportRowThreshold(Qostream & strm) const
{
    strm << LEVEL1 "Row threshold: " << dict_.at("row_threshold:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportColThreshold(Qostream & strm) const
{
    strm << LEVEL1 "Column threshold: " << dict_.at("col_threshold:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportThreadNumber(Qostream & strm) const
{
    strm << LEVEL1 "Thread number: " << dict_.at("thread_number:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportTaskFactor(Qostream & strm) const
{
    strm << LEVEL1 "Task factor: " << dict_.at("task_factor:").value << '\n';
    strm.flush();
}

void ConfigLoader::reportIsfast(Qostream & strm) const
{
	strm << LEVEL1 "Isfast: " << dict_.at("Isfast:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportDx(Qostream& strm) const
{
	strm << LEVEL1 "Dx: " << dict_.at("Dx:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportDy(Qostream& strm) const
{
	strm << LEVEL1 "Dy: " << dict_.at("Dy:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportD_Angle(Qostream& strm) const
{
	strm << LEVEL1 "D_Angle: " << dict_.at("D_Angle:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportA_Angle(Qostream& strm) const
{
	strm << LEVEL1 "A_Angle: " << dict_.at("A_Angle:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportPhase_0(Qostream& strm) const
{
	strm << LEVEL1 "Phase_0: " << dict_.at("Phase_0:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportPhase_x(Qostream& strm) const
{
	strm << LEVEL1 "Phase_x: " << dict_.at("Phase_x:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportPhase_y(Qostream& strm) const
{
	strm << LEVEL1 "Phase_y: " << dict_.at("Phase_y:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportArray_x(Qostream& strm) const
{
	strm << LEVEL1 "Array_x: " << dict_.at("Array_x:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportArray_y(Qostream& strm) const
{
	strm << LEVEL1 "Array_y: " << dict_.at("Array_y:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportScale_alongx(Qostream& strm) const
{
	strm << LEVEL1 "Scale_alongx: " << dict_.at("Scale_alongx:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportScale_alongy(Qostream& strm) const
{
	strm << LEVEL1 "Scale_alongy: " << dict_.at("Scale_alongy:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportRotate_alongx(Qostream& strm) const
{
	strm << LEVEL1 "Rotate_alongx: " << dict_.at("Rotate_alongx:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportRotate_alongy(Qostream& strm) const
{
	strm << LEVEL1 "Rotate_alongy: " << dict_.at("Rotate_alongy:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportt_sum(Qostream& strm) const
{
	strm << LEVEL1 "t_sum: " << dict_.at("t_sum:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportEps1(Qostream & strm) const
{
	strm << LEVEL1 "VIE Eps1: " << dict_.at("eps1:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportMu(Qostream & strm) const
{
	strm << LEVEL1 "VIE Mu1: " << dict_.at("mu:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportSigma(Qostream& strm) const
{
	strm << LEVEL1 "FSPGF_EFIE Sigma: " << dict_.at("Sigma:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportisContinuous(Qostream& strm) const
{
	strm << LEVEL1 "FSPGF_EFIE isContinuous: " << dict_.at("isContinuous:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportd_x(Qostream& strm) const
{
	strm << LEVEL1 "PGF_VIE d_x: " << dict_.at("d_x:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportd_y(Qostream& strm) const
{
	strm << LEVEL1 "PGF_VIE d_y: " << dict_.at("d_y:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportd_z(Qostream& strm) const
{
	strm << LEVEL1 "PGF_VIE d_z: " << dict_.at("d_z:").value << '\n';
	strm.flush();
}

void ConfigLoader::reportThickness_sheet(Qostream& strm) const
{
	strm << LEVEL1 "TDS Thickness_sheet: " << dict_.at("Thickness_sheet:").value << '\n';
	strm.flush();
}