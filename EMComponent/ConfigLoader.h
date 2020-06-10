//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "fieldStruct.h"
#include "ErrorHandle.h"

namespace component {

class ConfigLoader {
    struct Value {
        Value() = default;
        Value(bool parsed, const Qstring& val) : has_reported(parsed), value(val) { }
        bool has_reported;
        Qstring value;
    };
public:
    ConfigLoader();
    ~ConfigLoader();

    ConfigLoader(const ConfigLoader&) = delete;
    ConfigLoader& operator=(const ConfigLoader&) = delete;
public:
    void load(const Qstring& filepath);

    Qstring getSolverType();
    Incidence getInc();
    policy::ResultType getResultType();
    Scattering getSca();
    Radiation getRad();
	Nearfield getNearfield();
	MultipleIncidence getMultiInc();
    Qstring getMeshPath();
	Qstring getSurroundPath();
	Qstring getFEKOcurPath();
    value_t getAlpha();
    value_t getACABoxLength();
    value_t getACAThreshold();
    bool getEDMSwitch();
    int getCBFMaxBoxNum();
    int getCBFPolarization();
    int getCBFSampleSpacing();
    value_t getEpsilon1();
    value_t getEpsilon2();
    value_t getMu1();
    value_t getMu2();
    value_t getAIMXSpacing();
    value_t getAIMYSpacing();
    value_t getAIMZSpacing();
    value_t getNearThreshold();
    value_t getIterationThreshold();
    int getMaxIterationNum();
    value_t getFMMBoxLength();
    value_t getMLFMABoxThreshold();
    policy::BasisFunction getIEDGBasisFunctionType();
    value_t getStabilizationFactor();
    value_t getAverageSize();
    policy::DDM getDDMPolicy();
    bool getPreconditioningSwitch();
    value_t getRowThreshold();
    value_t getColThreshold();
    size_t getThreadNumber();
    size_t getTaskFactor();
	int getIsfast();
	value_t getDx();
	value_t getDy();
	value_t getD_Angle();
	value_t getA_Angle();
	value_t getPhase_0();
	value_t getPhase_x();
	value_t getPhase_y();
	int getArray_x();
	int getArray_y();
	int getisContinuous();
	value_t getScale_alongx();
	value_t getScale_alongy();
	value_t getRotate_alongx();
	value_t getRotate_alongy();
	int gett_sum();
	value_t getEps1();
	value_t getMu();
	Complex getSigma();
	value_t getd_x();
	value_t getd_y();
	value_t getd_z();
	value_t getThickness_sheet();
    Qstring getReport() const { return log_.str(); }
    void Debug() const;

private:
    void reportSolverType(Qostream& strm) const;
    void reportInc(Qostream& strm) const;
    void reportResultType(Qostream& strm) const;
    void reportSca(Qostream& strm) const;
    void reportRad(Qostream& strm) const;
	void reportNearfield(Qostream& strm) const;
	void reportMultiInc(Qostream& strm) const;
    void reportMeshPath(Qostream& strm) const;
	void reportFEKOcurPath(Qostream& strm) const;
	void reportSurroundPath(Qostream& strm) const;
    void reportAlpha(Qostream & strm) const;
    void reportACABoxLength(Qostream& strm) const;
    void reportACAThreshold(Qostream& strm) const;
    void reportEDMSwitch(Qostream& strm) const;
    void reportCBFMaxBoxNum(Qostream& strm) const;
    void reportCBFPolarization(Qostream& strm) const;
    void reportCBFSampleSpacing(Qostream& strm) const;
    void reportEpsilon1(Qostream& strm) const;
    void reportEpsilon2(Qostream& strm) const;
    void reportMu1(Qostream& strm) const;
    void reportMu2(Qostream& strm) const;
    void reportAIMXSpacing(Qostream& strm) const;
    void reportAIMYSpacing(Qostream& strm) const;
    void reportAIMZSpacing(Qostream& strm) const;
    void reportNearThreshold(Qostream& strm) const;
    void reportIterationThreshold(Qostream& strm) const;
    void reportMaxIterationNum(Qostream& strm) const;
    void reportFMMBoxLength(Qostream& strm) const;
    void reportMLFMABoxThreshold(Qostream& strm) const;
    void reportIEDGBasisFunctionType(Qostream& strm) const;
    void reportStabilizationFactor(Qostream& strm) const;
    void reportDDMPolicy(Qostream& strm) const;
    void reportPreconditioningSwitch(Qostream& strm) const;
    void reportAverageSize(Qostream& strm) const;
    void reportRowThreshold(Qostream& strm) const;
    void reportColThreshold(Qostream& strm) const;
    void reportThreadNumber(Qostream& strm) const;
    void reportTaskFactor(Qostream& strm) const;
	void reportIsfast(Qostream& strm) const;
	void reportDx(Qostream& strm) const;
	void reportDy(Qostream& strm) const;
	void reportD_Angle(Qostream& strm) const;
	void reportA_Angle(Qostream& strm) const;
	void reportPhase_0(Qostream& strm) const;
	void reportPhase_x(Qostream& strm) const;
	void reportPhase_y(Qostream& strm) const;
	void reportArray_x(Qostream& strm) const;
	void reportArray_y(Qostream& strm) const;
	void reportScale_alongx(Qostream& strm) const;
	void reportScale_alongy(Qostream& strm) const;
	void reportRotate_alongx(Qostream& strm) const;
	void reportRotate_alongy(Qostream& strm) const;
	void reportt_sum(Qostream& strm) const;
	void reportEps1(Qostream& strm) const;
	void reportMu(Qostream& strm) const;
	void reportSigma(Qostream& strm) const;
	void reportisContinuous(Qostream& strm) const;
	void reportd_x(Qostream& strm) const;
	void reportd_y(Qostream& strm) const;
	void reportd_z(Qostream& strm) const;
	void reportThickness_sheet(Qostream & strm) const;
private:
    Qsstream log_;
    std::map<Qstring, Value> dict_;
};

}