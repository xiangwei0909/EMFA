//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "VectorR3.h"
namespace component {

struct Incidence {
    value_t freq;
    value_t theta;
    value_t phi;
    value_t pole;

	value_t freq_from;
	value_t freq_to;
	value_t freq_delta;
	value_t freq_inc;
};

struct Scattering {
    value_t theta_from, theta_to, theta_delta;
    value_t phi_from, phi_to, phi_delta;
};

struct Radiation {
    bool e_plane;   //  phi = 0     theta 0 - 360
    bool h_plane;   //  phi = 90    theta 0 - 360
    bool hh_plane;  //  phi 0 - 360     theta = 90
};
struct MultipleIncidence {
	value_t PW_theta_from;
	value_t PW_theta_delta;
	value_t PW_theta_to;
	value_t	PW_phi_from;
	value_t PW_phi_delta;
	value_t PW_phi_to;
	Qstring polarization;
	int		PW_num;
};

struct Nearfield {
	VectorR3 origin_point;
	VectorR3 end_point;
	int sampling_x;
	int sampling_y;
	int sampling_z;
};
} // namespace component