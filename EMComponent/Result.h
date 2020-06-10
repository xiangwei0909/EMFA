//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "fieldStruct.h"
#include "VectorR3.h"
#include "VectorC3.h"
#include "tools.h"
#include "MyFWD.h"

#define RESULT_DECLARATION friend class component::Result

namespace component {

struct FieldData {
    value_t atheta, aphi;
    Complex ftheta, fphi;
};

struct CurrentData {
    VectorR3 v1, v2, v3;
	int n[3];
    value_t magn1, magn2, magn3, magnc;
};

struct NearFieldData {
	NearFieldData(VectorR3 _v, VectorC3 _cur) :v(_v), cur(_cur)
	{
	}
	VectorR3 v;
	VectorC3 cur;
	VectorC3 sca_inc;
};
Qostream& operator<<(Qostream& stream, const FieldData& data);
Qostream& operator<<(Qostream& stream, const CurrentData& data);

class Result {
public:
    Result();
    ~Result();
    
    Result(const Result&) = delete;
    Result& operator=(const Result&) = delete;

    template<typename ClassType, typename Functor>
    void getBistaticRCS(ClassType* ptr,  Functor caller)
    {
        auto dir = ptr->dir_;
        auto sca = ptr->scattering_;
        Qofstream rcs_data;
        VectorR3 sca_k;

        int num = static_cast<int>((sca.phi_to - sca.phi_from) / sca.phi_delta + 1);
        TOOL int total_progress = static_cast<int>((sca.theta_to - sca.theta_from) / sca.theta_delta) + 1;
        TOOL total_progress = total_progress * num + 360;
        TOOL int cur_progress = 0;
        TOOL tool::BarAndPercent bar_perc;
		//**************************************************
		//rcs_data.open((dir + "/bistatic_rcs_all") + ".dat");
		//prefixRCSData_ALL(rcs_data);
		//**************************************************
        for (int i = 0; i < num; ++i)
        {
            value_t phiDegree = sca.phi_from + i * sca.phi_delta;
            value_t phi = phiDegree * DegreesToRadians;
            value_t cosp = cos(phi), sinp = sin(phi);
            auto str = std::to_string(phiDegree);
            rcs_data.open(dir + "/bistatic_rcs_phi_" + str.substr(0, str.find('.')) + ".dat");
            prefixRCSData(rcs_data);
            for (value_t thetaDegree = sca.theta_from; thetaDegree <= sca.theta_to; thetaDegree += sca.theta_delta)
            {
                value_t theta = thetaDegree * DegreesToRadians;
                sca_k.x = sin(theta) * cosp;
                sca_k.y = sin(theta) * sinp;
                sca_k.z = cos(theta);
				//VectorR3 _theta(cos(theta)*cosp, cos(theta)*sinp, -sin(theta));
				//VectorR3 _phi(-sinp, cosp, 0.0f);
                rcs_data << std::setw(10) << thetaDegree << '\t' << std::setw(10) << (ptr->*caller)(sca_k) << '\n';
				//rcs_data << std::setw(10) << thetaDegree << '\t' << std::setw(10) << phiDegree << '\t' << std::setw(10) << (ptr->*caller)(sca_k) << '\n';
                TOOL bar_perc(++cur_progress, total_progress);
            }
            rcs_data.flush();
            rcs_data.close();
        }
		//*****************
		//rcs_data.flush();
		//rcs_data.close();
		//*****************


        rcs_data.open(dir + "/bistatic_rcs_theta_90.dat");
        prefixRCSData(rcs_data);
        for (value_t phiDegree = 0; phiDegree < 360; phiDegree += 1)
        {
            value_t phi = phiDegree * DegreesToRadians;
            sca_k = VectorR3(cos(phi), sin(phi), 0);
            rcs_data << std::setw(10) << phiDegree << '\t' << std::setw(10) << (ptr->*caller)(sca_k) << '\n';
            TOOL bar_perc(++cur_progress, total_progress);
        }
        rcs_data.flush();
        rcs_data.close();
    }

    template<typename ClassType, typename Functor>
    void getMonostaticRCS(const ClassType* ptr, Functor caller) const
    {
        auto dir = ptr->dir_;
        auto sca = ptr->scattering_;
        Qofstream rcs_data;
        VectorR3 sca_k;

        int phi_num = static_cast<int>((sca.phi_to - sca.phi_from) / sca.phi_delta) + 1;
        int theta_num = static_cast<int>((sca.theta_to - sca.theta_from) / sca.theta_delta) + 1;
        size_t cur_inc = 0, total_inc = phi_num * theta_num + 360;

        TOOL size_t cur_progress = 0, total_progress = total_inc;
        TOOL tool::BarAndPercent bar_perc;

        for (int i = 0; i < phi_num; ++i)
        {
            value_t phiDegree = sca.phi_from + i * sca.phi_delta;
            value_t phi = phiDegree * DegreesToRadians;
            value_t cosp = cos(phi), sinp = sin(phi);
            auto str = std::to_string(phiDegree);
            rcs_data.open(dir + "/monostatic_rcs_phi_" + str.substr(0, str.find('.')) + ".dat");
            prefixRCSData(rcs_data);
            for (value_t thetaDegree = sca.theta_from; thetaDegree <= sca.theta_to; thetaDegree += sca.theta_delta)
            {
                value_t theta = thetaDegree * DegreesToRadians;
                sca_k.x = sin(theta) * cosp;
                sca_k.y = sin(theta) * sinp;
                sca_k.z = cos(theta);
                rcs_data << std::setw(10) << thetaDegree << '\t' << std::setw(10) << (ptr->*caller)(cur_inc++, sca_k) << '\n';
                TOOL bar_perc(++cur_progress, total_progress);
            }
            rcs_data.flush();
            rcs_data.close();
        }
        rcs_data.open(dir + "/monostatic_rcs_theta_90.dat");
        prefixRCSData(rcs_data);
        for (value_t phiDegree = 0; phiDegree < 360; phiDegree += 1)
        {
            value_t phi = phiDegree * DegreesToRadians;
            sca_k = VectorR3(cos(phi), sin(phi), 0);
            rcs_data << std::setw(10) << phiDegree << '\t' << std::setw(10) << (ptr->*caller)(cur_inc++, sca_k) << '\n';
            TOOL bar_perc(++cur_progress, total_progress);
        }
        if (cur_inc != total_inc)
            throw std::runtime_error("(getMonostaticRCS) the number of incident wave is not equal to the number of scattering direction");

        rcs_data.flush();
        rcs_data.close();
    }

    template<typename ClassType, typename Functor>
    void getRadiationField(const ClassType* ptr, Functor caller) const
    {
        auto dir = ptr->dir_;
        auto filename = ptr->folder_name_;
        auto freq = ptr->incidence_.freq;
        auto rad = ptr->radiation_;
        Qofstream rad_data;
        VectorR3 rad_k, rad_ev, rad_eh;
        FieldData field;

        TOOL int total_progress = (rad.e_plane ? 361 : 0) + (rad.h_plane ? 361 : 0) + (rad.hh_plane ? 361 : 0);
        TOOL int cur_progress = 0;
        TOOL tool::BarAndPercent bar_perc;
        if (rad.e_plane)
        {
            rad_data.open(dir + "/radiation_ve.ffe");
            prefixFFEData(rad_data, filename, freq);
            field.aphi = 0;
            for (value_t thetaDegree = 0; thetaDegree <= 360; thetaDegree += 1)
            {
                field.atheta = thetaDegree;
                value_t theta = thetaDegree * DegreesToRadians;
                value_t phi0 = 0 * DegreesToRadians;
                rad_k = VectorR3(sin(theta) * cos(phi0), sin(theta) * sin(phi0), cos(theta));
                rad_ev = VectorR3(cos(theta) * cos(phi0), cos(theta) * sin(phi0), -sin(theta));
                rad_eh = VectorR3(-sin(phi0), -cos(phi0), 0);
                (ptr->*caller)(&field, rad_k, rad_ev, rad_eh);
                rad_data << field;
                TOOL bar_perc(++cur_progress, total_progress);
            }
            rad_data.flush();
            rad_data.close();
        }
        if (rad.h_plane)
        {
            rad_data.open(dir + "/radiation_vh.ffe");
            prefixFFEData(rad_data, filename, freq);
            field.aphi = 90;
            for (value_t thetaDegree = 0; thetaDegree <= 360; thetaDegree += 1)
            {
                field.atheta = thetaDegree;
                value_t theta = thetaDegree * DegreesToRadians;
                value_t phi90 = 90 * DegreesToRadians;
                rad_k = VectorR3(sin(theta) * cos(phi90), sin(theta) * sin(phi90), cos(theta));
                rad_ev = VectorR3(cos(theta) * cos(phi90), cos(theta) * sin(phi90), -sin(theta));
                rad_eh = VectorR3(-sin(phi90), -cos(phi90), 0);
                (ptr->*caller)(&field, rad_k, rad_ev, rad_eh);
                rad_data << field;
                TOOL bar_perc(++cur_progress, total_progress);
            }
            rad_data.flush();
            rad_data.close();
        }
        if (rad.hh_plane)
        {
            rad_data.open(dir + "/radiation_hh.ffe");
            prefixFFEData(rad_data, filename, freq);
            field.atheta = 90;
            for (value_t phiDegree = 0; phiDegree <= 360; phiDegree += 1)
            {
                field.aphi = phiDegree;
                value_t theta90 = 90 * DegreesToRadians;
                value_t phi = phiDegree * DegreesToRadians;
                rad_k = VectorR3(sin(theta90) * cos(phi), sin(theta90) * sin(phi), cos(theta90));
                rad_ev = VectorR3(cos(theta90) * cos(phi), cos(theta90) * sin(phi), -sin(theta90));
                rad_eh = VectorR3(-sin(phi), -cos(phi), 0);
                (ptr->*caller)(&field, rad_k, rad_ev, rad_eh);
                rad_data << field;
                TOOL bar_perc(++cur_progress, total_progress);
            }
            rad_data.flush();
            rad_data.close();
        }
    }

    template<typename ClassType, typename Functor>
    void getCurrentDistribution(const ClassType* ptr, Functor caller) const
    {
        auto dir = ptr->dir_;
        std::vector<CurrentData> data;
        (ptr->*caller)(&data);

        size_t data_num = data.size();
        Qofstream cur_data;
        cur_data.open(dir + "/current.sc");
        cur_data << "triangles: " << data_num << '\n';
        prefixSCData(cur_data);

        TOOL size_t cur_progress = 0, total_progress = data.size();
        TOOL tool::BarAndPercent bar_perc;
        for (size_t t = 0; t < data_num; ++t)
        {
            cur_data << std::setw(10) << t + 1 << data[t] << '\n';
            TOOL bar_perc(++cur_progress, total_progress);
        }
        cur_data.flush();
        cur_data.close();
	}

	template<typename ClassType, typename Functor>
	void getNearField(const ClassType* ptr, Functor caller) const
	{
		auto dir = ptr->dir_;
		std::vector<NearFieldData> data;
		(ptr->*caller)(&data);

		struct tm t;
		time_t now;
		localtime_s(&t, &now);

		Qofstream NEField_efe, NEField_tot;
		NEField_efe.open(dir + "/nearfield_sca.efe");
		NEField_efe << "##File Type: Electric near field" << '\n'
			<< "##File Format: 3" << '\n'
			<< "##Source: model" << '\n'
			<< "##Date: " << 1900 + t.tm_year << "-" << 1 + t.tm_mon << "-" << t.tm_mday << " "
			<< t.tm_hour << ":" << t.tm_min << ":" << t.tm_sec << '\n'
			<< "**File exported by POSTFEKO version 7.0-238289(x64)" << '\n'
			<< '\n'
			<< "#Request Name: Near Field" << '\n'
			<< "#Frequency: " << ptr->incidence_.freq << '\n'
			<< "#Coordinate System: Cartesian" << '\n'
			<< "#No. of X Samples: " << ptr->nearfield_.sampling_x << '\n'
			<< "#No. of Y Samples: " << ptr->nearfield_.sampling_y << '\n'
			<< "#No. of Z Samples: " << ptr->nearfield_.sampling_z << '\n'
			<< "#Result Type: Electric Field Values" << '\n'
			<< "#No. of Header Lines: 1" << '\n'
			<< std::left << setw(18) << "#"
			<< std::left << setw(18) << "\"X\"" << std::left << setw(18) << "\"Y\"" << std::left << setw(18) << "\"Z\""
			<< std::left << setw(18) << "\"Re(Ex)\"" << std::left << setw(18) << "\"Im(Ex)\""
			<< std::left << setw(18) << "\"Re(Ey)\"" << std::left << setw(18) << "\"Im(Ey)\""
			<< std::left << setw(18) << "\"Re(Ez)\"" << std::left << setw(18) << "\"Im(Ez)\""
			<< '\n';

		NEField_tot.open(dir + "/nearfield_tot.efe");
		NEField_tot << "##File Type: Electric near field" << '\n'
			<< "##File Format: 3" << '\n'
			<< "##Source: model" << '\n'
			<< "##Date: " << 1900 + t.tm_year << "-" << 1 + t.tm_mon << "-" << t.tm_mday << " "
			<< t.tm_hour << ":" << t.tm_min << ":" << t.tm_sec << '\n'
			<< "**File exported by POSTFEKO version 7.0-238289(x64)" << '\n'
			<< '\n'
			<< "#Request Name: Near Field" << '\n'
			<< "#Frequency: " << ptr->incidence_.freq << '\n'
			<< "#Coordinate System: Cartesian" << '\n'
			<< "#No. of X Samples: " << ptr->nearfield_.sampling_x << '\n'
			<< "#No. of Y Samples: " << ptr->nearfield_.sampling_y << '\n'
			<< "#No. of Z Samples: " << ptr->nearfield_.sampling_z << '\n'
			<< "#Result Type: Electric Field Values" << '\n'
			<< "#No. of Header Lines: 1" << '\n'
			<< std::left << setw(18) << "#"
			<< std::left << setw(18) << "\"X\"" << std::left << setw(18) << "\"Y\"" << std::left << setw(18) << "\"Z\""
			<< std::left << setw(18) << "\"Re(Ex)\"" << std::left << setw(18) << "\"Im(Ex)\""
			<< std::left << setw(18) << "\"Re(Ey)\"" << std::left << setw(18) << "\"Im(Ey)\""
			<< std::left << setw(18) << "\"Re(Ez)\"" << std::left << setw(18) << "\"Im(Ez)\""
			<< '\n';

		size_t data_num = data.size();
		TOOL tool::BarAndPercent bar_perc;
		for (int t = 0; t < data_num; t++)
		{
			TOOL bar_perc(t+1, data_num);

			NEField_efe << std::left << std::showpos << std::fixed << std::scientific << setw(18) << " "
				<< std::left << setw(18) << data[t].v.x
				<< std::left << setw(18) << data[t].v.y
				<< std::left << setw(18) << data[t].v.z
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].cur.x.real()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].cur.x.imag()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].cur.y.real()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].cur.y.imag()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].cur.z.real()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].cur.z.imag()
				<< '\n';

			NEField_tot << std::left << std::showpos << std::fixed << std::scientific << setw(18) << " "
				<< std::left << setw(18) << data[t].v.x
				<< std::left << setw(18) << data[t].v.y
				<< std::left << setw(18) << data[t].v.z
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].sca_inc.x.real()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].sca_inc.x.imag()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].sca_inc.y.real()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].sca_inc.y.imag()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].sca_inc.z.real()
				<< std::left << std::showpos << std::fixed << std::scientific << setw(18) << data[t].sca_inc.z.imag()
				<< '\n';
		}
		NEField_efe.flush();
		NEField_efe.close();

		NEField_tot.flush();
		NEField_tot.close();
	}
private:
    void prefixFFEData(Qostream &strm, const Qstring& name, value_t freq) const;
    void prefixRCSData(Qostream& strm) const;
	void prefixRCSData_ALL(Qostream& strm) const;
    void prefixSCData(Qostream& strm) const;
};

}