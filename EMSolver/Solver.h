//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "Custom.h"
#include "ConfigLoader.h"

namespace EMSuit {

class EMAssistant;

}

namespace solver {

template <typename QSolver>
void invoke_solver(component::ConfigLoader* ploader)
{
    auto solver_sptr = std::make_unique<QSolver>();
    solver_sptr->init(ploader);
    solver_sptr->reportInfo(Qcout);
    solver_sptr->solve();
    solver_sptr->output();
    solver_sptr->clear();
}


void loadSolver(const Qstring& conf_file);

void loadSolver(EMSuit::EMAssistant*);

} // namespace solver

#endif // !_SOLVER_H_
