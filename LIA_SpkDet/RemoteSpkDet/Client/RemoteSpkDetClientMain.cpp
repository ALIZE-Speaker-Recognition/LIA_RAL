/**
 * \file RemoteSpkDetClientMain.cpp
 * \author Christophe LEVY, Alexandre PRETI, Teva Merlin
 * \version 1.5
 *
 * \brief Example of how to write a client for the RemoteSpkDet protocol
 *
**/

#include <iostream>
#include <cstdlib>

#include "RemoteSpkDetClient.h"

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
        RemoteSpkDetClient();
    }
    catch (exception& e) {
        cerr<<e.what()<<endl;
    }
    return EXIT_SUCCESS;
}
