//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "MyDebugNew.h"
#include "ThreadPool.h"
#include "tools.h"

#include "iml.h"


#include "EFIE.h"
#include "MFIE.h"
#include "CFIE.h"
//#include "ACA_CFIE.h"
//#include "CBFM.h"
//#include "EICBFM.h"
//#include "ACA_EICBFM.h"
//#include "PMCHW.h"
//#include "CBFM_PMCHW.h"
//#include "EICBFM_PMCHW.h"
//#include "EFIE_PMCHW.h"
//#include "CBFM_EFIE_PMCHW.h"
//#include "IE_NDDM_CFIE.h"
//#include "NDDM_PMCHW.h"
//#include "ACA_DDM.h"
//#include "AIM.h"
//#include "PO.h"
#include "FMM.h"
#include "MLFMA.h"
//#include "IEDG.h"
//#include "IEDG_MLFMA.h"
//#include "IP_IE_DDM.h"
#include "FastWT.h"

namespace test {

void Test();

namespace qmat {
void OperatorAdd();
}

namespace threadpool {
void TestThreadPool();
}
} // namespace test