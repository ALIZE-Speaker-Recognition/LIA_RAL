/**
 * \file SpkDetServerMain.cpp
 * \author Christophe LEVY, alexandre PRETI
 * \version 1.0
 * \date september, 13th 2007
 *
 * \brief Server for Alize use
 *
**/

#include <iostream>
#include "../include/SpkDetServer.h"

using namespace alize;
using namespace std;

/*! \fn int main(int argc, char* argv[])
 *  \brief  main
 *
 *  \param[in]      argc        number of parameter in the command line
 *  \param[in]      argv        command line
 *
 *  \return EXIT_SUCCESS if no problem, otherwise EXIT_FAILLURE
 */
int main(int argc, char* argv[]) {
    try {
        CmdLine cmdLine(argc, argv);
        if (cmdLine.displayHelpRequired()) { 
            // --help
        }
        else {
            Config tmp;
            cmdLine.copyIntoConfig(tmp);
            Config config;
            if (tmp.existsParam("config")) 
                config.load(tmp.getParam("config"));
            cmdLine.copyIntoConfig(config);
            //config.setParam("featureFilesPath", "./");
            SpkDetServer(config);
        }   
    }
    catch (alize::Exception& e) {
        cout << e.toString() << endl;
    }
    return EXIT_SUCCESS;
}
