#include "Solver.h"
#include "Test.h"

//#define TEST
//#define SERVER 

/*********************** Main ***************************/
int main(int argc, char** argv)
{
#ifdef TEST
    test::Test();
#else
    Qstring config, state("OK");

    Qcout << HEADING "Input Path of Config File:";
    std::getline(Qcin, config);

    if (config.empty())
        config = "C:/Users/WT/Desktop/EFIE/CONFIG.yml";

    //solver::loadSolver(config);   /// Debug
    try {
        solver::loadSolver(config);
    }
    catch (const component::ConfigInvalid& e) {
        state = Qstring("[ConfigParameter Error]") + e.what();
    }
    catch (const component::FileError& e) {
        state = Qstring("[File Error]") + e.what();
    }
    catch (const std::invalid_argument& e) {
        state = Qstring("[InvalidArgument Error]") + e.what();
    }
    catch (const std::runtime_error& e) {
        state = Qstring("[Runtime Error]") + e.what();
    }
    catch (const std::exception& e) {
        state = Qstring("[Unrecognized Error]") + e.what();
    }
    Qcout << HEADING "Exit State: " << state << "\n" HEADING << std::endl;

#endif

#ifdef _DEBUG
    _CrtDumpMemoryLeaks();
#endif
	system("pause");
    return 0;
}