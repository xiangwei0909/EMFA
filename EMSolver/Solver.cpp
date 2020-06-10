#include "Solver.h"
#include "MyDebugNew.h"
#include "EFIE.h"
#include "MFIE.h"
#include "CFIE.h"
#include "FMM.h"
#include "MLFMA.h"
#include "FastWT.h"
#include "FSPGF_EFIE.h"
#include "FSPGF_EFIE_HS.h"
#include "EFIE_HS.h"
#include "VIE.h"
#include "VSIE.h"
#include "PGF_VIE.h"
#include "DGVIE.h"
#include "DGVSIE.h"
#include "PGF_VSIE.h"
#include "NEWSED.h"
#include "New_ASED_ID.h"
#include "New_ASED_PEC.h"
#include "New_ASED_EDM.h"
#include "EPA_EFIE.h"
#include "New_ASED_gra.h"
#include "SED_CM.h"
#include "SED_CM_FMM.h"
#include "SED_QP.h"
#include "MPIE.h"
#include "NSED_metal.h"
#include "SED_QP_MAIN.h"
//#include "NEWSED_SPIN.h"

#define INVOKE(solver_name) solver::invoke_solver<solver_name>(&loader)
#define SOLVER_OPTION(solver_string, solver_name) \
        if(solver_type == solver_string) \
            INVOKE(solver_name);
#define ADD_SOLVER_OPTION(solver_string, solver_name) \
        else if(solver_type == solver_string) \
            INVOKE(solver_name);
#define DEFAULT_SOLVER_OPTION() \
        else { throw ConfigInvalid("undefined EM solver type (" + solver_type + ')'); }

using namespace component;
using namespace mom;

void solver::loadSolver(const Qstring & conf_file)
{
    ConfigLoader loader;
    loader.load(conf_file);
    auto solver_type = loader.getSolverType();

	SOLVER_OPTION("EFIE", EFIE)
		ADD_SOLVER_OPTION("MFIE", MFIE)
		ADD_SOLVER_OPTION("CFIE", CFIE)
		//ADD_SOLVER_OPTION("ACA_CFIE", ACA_CFIE)
		//ADD_SOLVER_OPTION("CBFM", CBFM)
		//ADD_SOLVER_OPTION("EICBFM", EICBFM)
		//ADD_SOLVER_OPTION("ACA_EICBFM", ACA_EICBFM)
		//ADD_SOLVER_OPTION("PMCHW", PMCHW)
		//ADD_SOLVER_OPTION("EFIE_PMCHW", EFIE_PMCHW)
		//ADD_SOLVER_OPTION("CBFM_EFIE_PMCHW", CBFM_EFIE_PMCHW) // Error
		//ADD_SOLVER_OPTION("EICBFM_PMCHW", EICBFM_PMCHW)
		//ADD_SOLVER_OPTION("AIM", AIM)
		//ADD_SOLVER_OPTION("DDM", IE_NDDM_CFIE)
		//ADD_SOLVER_OPTION("DDM_PMCHW", NDDM_PMCHW) // Error
		//ADD_SOLVER_OPTION("ACA_DDM", ACA_DDM)
		//ADD_SOLVER_OPTION("PO", PO)
		ADD_SOLVER_OPTION("FMM", FMM)
		ADD_SOLVER_OPTION("MLFMA", MLFMA)
		//ADD_SOLVER_OPTION("IEDG", IEDG)
		//ADD_SOLVER_OPTION("IEDG_MLFMA", IEDG_MLFMA)
		//ADD_SOLVER_OPTION("IP_IE_DDM", IP_IE_DDM)
		ADD_SOLVER_OPTION("FastWT", FastWT)
		ADD_SOLVER_OPTION("FSPGF_EFIE", FSPGF_EFIE)
		ADD_SOLVER_OPTION("FSPGF_EFIE_HS", FSPGF_EFIE_HS)
		ADD_SOLVER_OPTION("EFIE_HS", EFIE_HS)
		ADD_SOLVER_OPTION("VIE", VIE)
		ADD_SOLVER_OPTION("VSIE", VSIE)
		ADD_SOLVER_OPTION("PGF_VIE", PGF_VIE)
		ADD_SOLVER_OPTION("DGVIE", DGVIE)
		ADD_SOLVER_OPTION("DGVSIE", DGVSIE)
		ADD_SOLVER_OPTION("PGF_VSIE", PGF_VSIE)
		ADD_SOLVER_OPTION("NEWSED",NEWSED)
		ADD_SOLVER_OPTION("New_ASED_ID", New_ASED_ID)
		ADD_SOLVER_OPTION("New_ASED_PEC",New_ASED_PEC)
		ADD_SOLVER_OPTION("New_ASED_EDM", New_ASED_EDM)
		ADD_SOLVER_OPTION("EPA_EFIE", EPA_EFIE)
		ADD_SOLVER_OPTION("New_ASED_gra", New_ASED_gra)
		ADD_SOLVER_OPTION("SED_CM", SED_CM)
		ADD_SOLVER_OPTION("SED_CM_FMM", SED_CM_FMM)
		ADD_SOLVER_OPTION("SED_QP", SED_QP)
		ADD_SOLVER_OPTION("NSED_metal", NSED_metal)
		ADD_SOLVER_OPTION("SED_QP_MAIN", SED_QP_MAIN)
		//ADD_SOLVER_OPTION("NEWSED_SPIN", NEWSED_SPIN)
    DEFAULT_SOLVER_OPTION()
}

void solver::loadSolver(EMSuit::EMAssistant *)
{

}
