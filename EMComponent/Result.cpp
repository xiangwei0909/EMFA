#include "Result.h"

using namespace component;
using std::setw;

Result::Result()
{
}


Result::~Result()
{
}

void component::Result::prefixFFEData(Qostream & strm, const Qstring& name, value_t freq) const
{
    strm << std::scientific << std::left << std::setprecision(8);
    strm << "##File Type: Far field\n##File Format: 3\n##Source: " << name
        << "\n##Date:\n** File exported by POSTFEKO version 7.0-238289 (x64)\n\n#Request Name: farField\n"
        << "#Frequency: " << freq << "\n#Coordinate System: Spherical\n#No. of Theta Samples: 361\n"
        "#No. of Phi Samples: 1\n#Result Type: Gain\n#No. of Header Lines: 1\n"
        "#                 \"Theta\"           \"Phi\"             \"Re(Etheta)\"      \"Im(Etheta)\"      "
        "\"Re(Ephi)\"        \"Im(Ephi)\"        \"Gain(Theta)\"     \"Gain(Phi)\"       \"Gain(Total)\"\n";
}

void component::Result::prefixRCSData(Qostream & strm) const
{
    strm << std::right << std::fixed << std::setprecision(4);
    strm << std::setw(10) << "Theta\t" << std::setw(10) << "RCS(dBm^2)\n";
}

void component::Result::prefixRCSData_ALL(Qostream& strm) const
{
	strm << std::right << std::fixed << std::setprecision(4);
	strm << std::setw(10) << "Theta\t" << std::setw(10) << "Phi\t" << std::setw(10) << "RCS(dBm^2)"<<'\n';
}

void component::Result::prefixSCData(Qostream & strm) const
{
    strm << std::fixed << std::setprecision(8) << std::right;
    strm << std::setw(10) << "id"
        << std::setw(18) << "v1.x" << std::setw(18) << "v1.y" << std::setw(18) << "v1.z"
        << std::setw(18) << "v2.x" << std::setw(18) << "v2.y" << std::setw(18) << "v2.z"
        << std::setw(18) << "v3.x" << std::setw(18) << "v3.y" << std::setw(18) << "v3.z"
        << std::setw(18) << "magn1" << std::setw(18) << "magn2" << std::setw(18) << "magn3"
        << std::setw(18) << "magn center"
        << '\n';
}

Qostream & component::operator<<(Qostream & stream, const FieldData & data)
{
    stream << setw(18) << ' '
        << setw(18) << data.atheta << setw(18) << data.aphi 
        << setw(18) << data.ftheta.real() << setw(18) << data.ftheta.imag()
        << setw(18) << data.fphi.real() << setw(18) << data.fphi.imag()
        << setw(18) << 0.0f << setw(18) << 0.0f << setw(18) << 0.0f << '\n';
    return stream;
}

Qostream & component::operator<<(Qostream & stream, const CurrentData & data)
{
    stream << setw(18) << data.v1.x << setw(18) << data.v1.y << setw(18) << data.v1.z
        << setw(18) << data.v2.x << setw(18) << data.v2.y << setw(18) << data.v2.z
        << setw(18) << data.v3.x << setw(18) << data.v3.y << setw(18) << data.v3.z
        << setw(18) << data.magn1 << setw(18) << data.magn2 << setw(18) << data.magn3
        << setw(18) << data.magnc;
    return stream;
}
