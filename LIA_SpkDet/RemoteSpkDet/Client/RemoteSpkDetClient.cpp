// RemoteSpkDetClient.cpp
// This file is a part of LIA Software LIA_SpkDet, based on ALIZE toolkit 
// LIA_SpkDet  is a free, open tool for speaker recognition
// LIA_SpkDet is a development project initiated and funded by the LIA lab.
// See lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet
// for more information about ALIZE, see http://alize.univ-avignon.fr/
//
// Copyright (C) 2004
//  Laboratoire d'informatique d'Avignon [lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@univ-avignon.fr]
//      
// LIA_SpkDet is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// You should have received a copy of the GNU General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// The LIA team as well as the ALIZE project want to highlight the limits of voice authentication
// in a forensic context. 
// The following paper proposes a good overview of this point:
// [Bonastre J.F., Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-chagnolleau I.,
//  Person  Authentification by Voice: A Need of Caution,
//  Eurospeech 2003, Genova]
// The conclusion of the paper of the paper is proposed bellow:
// [Currently, it is not possible to completely determine whether the
//  similarity between two recordings is due to the speaker or to other
//  factors, especially when: (a) the speaker does not cooperate, (b) there
//  is no control over recording equipment, (c) recording conditions are not 
//  known, (d) one does not know whether the voice was disguised and, to a
//  lesser extent, (e) the linguistic content of the message is not
//  controlled. Caution and judgment must be exercised when applying speaker
//  recognition techniques, whether human or automatic, to account for these
//  uncontrolled factors. Under more constrained or calibrated situations,
//  or as an aid for investigative purposes, judicious application of these
//  techniques may be suitable, provided they are not considered as infallible.
//  At the present time, there is no scientific process that enables one to
//  uniquely characterize a person=92s voice or to identify with absolute
//  certainty an individual from his or her voice.]
//
// Contact Jean-Francois Bonastre (jean-francois.bonastre@univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet
// May 2007 - A.P.
// Updated May 2017 - TM

/**
 * \file RemoteSpkDetClient.cpp
 * \author Christophe LEVY, Alexandre PRETI, Teva Merlin
 * \version 1.5
 *
 * \brief Example of how to write a client for the RemoteSpkDet protocol
 *
**/

#include <cstdlib> 
#include <cstring> 
#include <arpa/inet.h>
#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <fstream>

#include "SpkDetServerConstants.h"

using namespace std;

/*! \fn uint8_t* paquetize(const uint8_t &uctype, const uint32_t &ulsize, const uint8_t *cdata)
 *  \brief  make a paquet containing the type followed the size of data followed by the data
 *
 *  \param[in]      uctype      the type of requested/replied fonctionnality
 *  \param[in]      ulsize      the data size
 *  \param[in]      cdata       the data
 *
 *  \return the paquet to send (allocated in the function)
 */
uint8_t* paquetize(const uint8_t &uctype, const uint32_t &ulsize, const uint8_t *cdata) {
    uint8_t *buff;
    uint32_t ullocSize = htonl(ulsize);
    
    buff = (uint8_t*) malloc(1+4+ulsize);
    if (!buff)
        return NULL;
    memmove((void*)buff, (const void*)&uctype, 1);
    memmove((void*)(buff+(1)), (const void*)&ullocSize, 4);
    memmove((void*)(buff+(1+4)), (const void*)cdata, ulsize);   
    return buff;
}

/*! \fn int send(const int isockfd, const uint8_t &uctype, const uint32_t &ulsize, const uint8_t *cdata)
 *  \brief  send a paquet (build by paquetize) to a socket
 *
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      uctype      the type of requested/replied fonctionnality
 *  \param[in]      ulsize      the data size
 *  \param[in]      cdata       the data
 *
 *  \return the number of byte sent
 */
int send(const int isockfd, const uint8_t &uctype, const uint32_t &ulsize, const uint8_t *cdata) {
    uint8_t *cbuff;
    ssize_t ssend;
    
    cbuff = paquetize(uctype, ulsize, cdata);
    if (!cbuff) {
        cerr<<"Paquetize error: paquet not allocated "<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }
    ssend = write(isockfd, cbuff, 1+4+ulsize);
    free(cbuff);
    if ((uint32_t)ssend != (1+4+ulsize)) {
        cerr<<"Send error: not enough data written "<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }
    return ssend;
}

/*! \fn void printCommandList()
 *  \brief  print the list of supported command
 *
 */
void printCommandList() {
    cout<<"G_QUIT 0"<<endl;
    cout<<"G_LIST 1"<<endl;
    cout<<"G_RESET 2"<<endl;
    cout<<"G_STATUS 3"<<endl;
    cout<<"G_SENDOPT 4"<<endl;
    cout<<"A_RESET 10"<<endl;
    cout<<"A_SAVE 11"<<endl;
    cout<<"A_LOAD 12"<<endl;
    cout<<"A_SEND 13"<<endl;
    cout<<"F_RESET 30"<<endl;
    cout<<"F_SAVE 31"<<endl;
    cout<<"F_LOAD 32"<<endl;
    cout<<"F_SEND 33"<<endl;
    cout<<"M_RESET 50"<<endl;   
    cout<<"M_SAVE 51"<<endl;
    cout<<"M_LOAD 52"<<endl;
    cout<<"M_WLOAD 53"<<endl;
    cout<<"M_DEL 54"<<endl;
    cout<<"M_ADAPT 55"<<endl;
    cout<<"M_TRAIN 56"<<endl;
    cout<<"I_DET 70"<<endl;
    cout<<"I_ID 71"<<endl;
    cout<<"I_DETCUM 72"<<endl;
    cout<<"I_IDCUM 73"<<endl;
    cout<<"I_DETCUMR 74"<<endl;
    cout<<"I_IDCUMR 75"<<endl;
}

/*! \fn bool G_Reset(const int isockfd)
 *  \brief  send a G_RESET message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if G_RESET ok in the server, otherwise false
 */
bool G_Reset(const int isockfd) {
    string configName;
    uint8_t creturn;
    
    cout<<"Enter config file name: ";
    cin>>configName;

    send(isockfd, G_RESET, configName.size()+1, (uint8_t*)configName.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Server not reset: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Server correctly reset"<<endl;
    return true;
}

/*! \fn bool G_Status(const int isockfd)
 *  \brief  send a G_STATUS message
 *
 *  \todo   receive all speaker id
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true
 */
bool G_Status(const int isockfd) {
    uint32_t lnbFeature, lnbSpk;
    string sUId;
    uint8_t c;
    
    send(isockfd, G_STATUS, 0, NULL);

    read(isockfd, &c, 1);
    if (c != RSD_NO_ERROR)
        return false;
    
    bzero(&lnbFeature,4);
    if (read(isockfd, &lnbFeature, 4)) {         //server's answer: nb feature
        lnbFeature = ntohl(lnbFeature);
        cout<<"nbFeature loaded: "<<lnbFeature<<endl; 
    }
    bzero(&lnbSpk,4);
    if (read(isockfd, &lnbSpk, 4)) {             //server's answer: nb spk
        lnbSpk = ntohl(lnbSpk);
        cout<<"nbSpk loaded: "<<lnbSpk<<endl; 
    }
    for (unsigned long lcptr=0 ; lcptr<lnbSpk ; lcptr++) {
        sUId="";
        read(isockfd, &c, 1);
        do {
            sUId += (unsigned char)c;
            read(isockfd, &c, 1);
        } while (c!=0);
        cout<<"\t"<<sUId<<endl;
    }
    read(isockfd, &c, 1);
    if(c==0)
        cout<<"No UBM loaded"<<endl;
    else
        cout<<"UBM loaded"<<endl;
    
    return true;
} 

/*! \fn bool G_SendOpt(const int isockfd)
 *  \brief  send a G_SENDOPT message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if option set to the wanted value, otherwise false
 */
bool G_SendOpt(const int isockfd) {
    string opt, optValue, concat;
    uint8_t creturn;
    
    cout<<"Enter option to change: ";
    cin>>opt;
    cout<<"Enter the new value for "<<opt<<": ";
    cin>>optValue;
    concat = opt + '\0' + optValue; 

    send(isockfd, G_SENDOPT, concat.size()+1, (uint8_t*)concat.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Option  "<<opt<<" unchanged !!!: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<opt<<" set to "<<optValue<<endl;
    return true;
}

/*! \fn bool A_Reset(const int isockfd)
 *  \brief  send a A_RESET message
 *
 *  \attention NOT IMPLEMENTED in the server
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if A_RESET ok in the server, otherwise false
 */
bool A_Reset (const int isockfd) {  
    uint8_t creturn;

    send(isockfd, A_RESET, 0, NULL);
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cerr<<"Audio server not reset: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Audio Server correctly reset"<<endl;
    return true;
} 


/*! \fn bool A_Save (const int isockfd)
 *  \brief  send a A_SAVE message
 *
 *  \attention NOT IMPLEMENTED in the server
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the audio buffer is saved in the server, otherwise false
 */
bool A_Save (const int isockfd) {   
    string fileName;
    uint8_t creturn;
    
    cout<<"Enter filename to save audio buffer: ";
    cin>>fileName;

    send(isockfd, A_SAVE, fileName.size()+1, (uint8_t*)fileName.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Audio buffer not saved: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Audio buffer saved"<<endl;
    return true;
}

/*! \fn bool A_Load (const int isockfd)
 *  \brief  send a A_LOAD message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the audio buffer is loaded in the server, otherwise false
 */
bool A_Load (const int isockfd) {   
    string fileName;
    uint8_t creturn;
    
    cout<<"Enter the name of a SERVER-SIDE file to load into the audio buffer: ";
    cin>>fileName;

    send(isockfd, A_LOAD, fileName.size()+1, (uint8_t*)fileName.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Audio buffer not loaded: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Audio buffer loaded"<<endl;
    return true;
}

/*! \fn bool A_Send (const int isockfd)
 *  \brief  send a A_SEND message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the audio buffer is sent to the server, otherwise false
 */
bool A_Send (const int isockfd) {   
    uint8_t creturn;
    size_t size=320;
    char *ctab = new char[size];
    string filename;
    FILE *fin;
    unsigned long ll=0, totalSent=0;
    
    cout<<"Enter the name of a LOCAL audio file to send to the server: ";
    cin>>filename;
    cout<<"Sending the file"<<endl;
    
    fin = fopen(filename.c_str(), "rb");
        
    creturn = RSD_NO_ERROR;
    while(((ll=fread(ctab, 1, size, fin)) > 0) && (creturn == RSD_NO_ERROR)) {
        send(isockfd, A_SEND, ll, (uint8_t*)ctab);
        totalSent += ll;
        read(isockfd, &creturn, 1);
    }
    fclose(fin);
	delete[] ctab;
    if(creturn != RSD_NO_ERROR) {
        cout<<"Server replied with an error after "<<totalSent<<" bytes: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }

    send(isockfd, A_SEND, 0, NULL);    // signal the server that the transfer is over
    read(isockfd, &creturn, 1);        // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Server replied with an error after the end of transfer ("<<totalSent<<" bytes): "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Audio buffer sent ("<<totalSent<<" bytes)"<<endl;
	
    return true;
}

/*! \fn bool F_Reset(const int isockfd)
 *  \brief  send a F_RESET message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if F_RESET ok in the server, otherwise false
 */
bool F_Reset (const int isockfd) {  
    uint8_t creturn;

    send(isockfd, F_RESET, 0, NULL);
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cerr<<"FeatureServer not reset: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"FeatureServer correctly reset"<<endl;
    return true;
} 


/*! \fn bool F_Save (const int isockfd)
 *  \brief  send a F_SAVE message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the audio buffer is saved in the server, otherwise false
 */
bool F_Save (const int isockfd) {   
    string fileName;
    uint8_t creturn;
    
    cout<<"Enter filename to save FeatureServer: ";
    cin>>fileName;

    send(isockfd, F_SAVE, fileName.size()+1, (uint8_t*)fileName.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"FeatureServer not saved: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"FeatureServer saved"<<endl;
    return true;
}

/*! \fn bool F_Load (const int isockfd)
 *  \brief  send a F_LOAD message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the audio buffer is loaded in the server, otherwise false
 */
bool F_Load (const int isockfd) {   
    string fileName;
    uint8_t creturn;
    
    cout<<"Enter basename of feature file to load: ";
    cin>>fileName;

    send(isockfd, F_LOAD, fileName.size()+1, (uint8_t*)fileName.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Features not loaded: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Features loaded"<<endl;
    return true;
}

/*! \fn bool F_Send (const int isockfd)
 *  \brief  send a F_SEND message
 *
 *  \attention FSEND not yet implemented
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the audio buffer is sent to the server, otherwise false
 */
bool F_Send (const int isockfd) {   
    uint8_t creturn;
    size_t size=320;
    char *ctab = new char[size];
    string filename;
    FILE *fin;
    unsigned long ll=0, totalSent=0;
    
    cout<<"Enter the filename for features (RAW format): ";
    cin>>filename;
    cout<<"Sending the file"<<endl;
    
    fin = fopen(filename.c_str(), "rb");
    
    creturn = RSD_NO_ERROR;
    while(((ll=fread(ctab, 1, size, fin)) > 0) && (creturn == RSD_NO_ERROR)) {
        send(isockfd, F_SEND, ll, (uint8_t*)ctab);
        totalSent += ll;
        read(isockfd, &creturn, 1);
    }
    fclose(fin);
	delete[] ctab;
    if(creturn != RSD_NO_ERROR) {
        cout<<"Server replied with an error after "<<totalSent<<" bytes: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    
    send(isockfd, F_SEND, 0, NULL);    // signal the server that the transfer is over
    read(isockfd, &creturn, 1);        // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Server replied with an error after the end of transfer ("<<totalSent<<" bytes): "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Audio buffer sent ("<<totalSent<<" bytes)"<<endl;
	
    return true;
}

/*! \fn bool M_Reset (const int isockfd)
 *  \brief  send a M_RESET message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the model server is reset, otherwise false
 */
bool M_Reset (const int isockfd) {  
    uint8_t creturn;
    
    send(isockfd, M_RESET, 0, NULL);
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"USER model removed of model server"<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"ModelServer reset"<<endl;
    return true;
}

/*! \fn bool M_Save (const int isockfd)
 *  \brief  send a M_SAVE message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the \a user_id model is saved, otherwise false
 */
bool M_Save (const int isockfd) {   
    string fileName, uId, concat;
    uint8_t creturn;
    
    cout<<"Enter user_id for the model to save: ";
    cin>>uId;
    cout<<"Enter filename to save model["<<uId<<"]: ";
    cin>>fileName;

    concat = uId + '\0' + fileName;
    
    send(isockfd, M_SAVE, concat.length()+1, (uint8_t*)concat.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Model["<<uId<<"] not saved: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Model["<<uId<<"] saved"<<endl;
    return true;
}

/*! \fn bool M_Load (const int isockfd)
 *  \brief  send a M_LOAD message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the model is loaded, otherwise false
 */
bool M_Load (const int isockfd) {   
    string fileName, uId, concat;
    uint8_t creturn;
    
    
    cout<<"Enter user_id for the model to load: ";
    cin>>uId;
    cout<<"Enter filename to load model["<<uId<<"]: ";
    cin>>fileName;

    concat = uId + '\0' + fileName;
    
    send(isockfd, M_LOAD, concat.length()+1, (uint8_t*)concat.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"Model["<<uId<<"] not loaded: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Model["<<uId<<"] loaded"<<endl;
    return true;
}

/*! \fn bool M_WLoad (const int isockfd)
 *  \brief  send a M_WLOAD message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the world model is loaded, otherwise false
 */
bool M_WLoad (const int isockfd) {  
    string fileName;
    uint8_t creturn;
    
    cout<<"Enter filename to load the world model: ";
    cin>>fileName;

    send(isockfd, M_WLOAD, fileName.length()+1, (uint8_t*)fileName.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"World model not loaded: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"World model  loaded"<<endl;
    return true;
}

/*! \fn bool M_Del (const int isockfd)
 *  \brief  send a M_DEL message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the \a user_id model is removed from model server, otherwise false
 */
bool M_Del (const int isockfd) {    
    string uId;
    uint8_t creturn;
    
    cout<<"Enter the user_id of the model to remove: ";
    cin>>uId;

    send(isockfd, M_DEL, uId.length()+1, (uint8_t*)uId.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"model["<<uId<<"] not removed: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"model["<<uId<<"] removed"<<endl;
    return true;
}

/*! \fn bool M_Adapt (const int isockfd)
 *  \brief  send a M_ADAPT message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if the \a user_id model is correctly adapt with feature in feature server, otherwise false
 */
bool M_Adapt (const int isockfd) {  
    string uId;
    uint8_t creturn;
    
    cout<<"Enter the user_id of the model to adapt: ";
    cin>>uId;

    send(isockfd, M_ADAPT, uId.length()+1, (uint8_t*)uId.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"model["<<uId<<"] not adapted: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"model["<<uId<<"] adapted"<<endl;
    return true;
}

/*! \fn bool M_Train (const int isockfd)
 *  \brief  send a M_TRAIN message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if a model is correctly created with feature in feature server (user_id is requested), otherwise false
 */
bool M_Train (const int isockfd) {  
    string uId;
    uint8_t creturn;
    
    cout<<"Enter the user_id of the model to train: ";
    cin>>uId;

    send(isockfd, M_TRAIN, uId.length()+1, (uint8_t*)uId.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"model["<<uId<<"] not trained: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"model["<<uId<<"] trained"<<endl;
    return true;
}

/*! \fn bool I_Det (const int isockfd)
 *  \brief  send a I_DET message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_Det (const int isockfd) {    
    string uId;
    uint8_t creturn, cdecision;
    float fscore;
    
    cout<<"Enter the user_id for whom a score must be computed: ";
    cin>>uId;

    send(isockfd, I_DET, uId.length()+1, (uint8_t*)uId.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                          // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"I_Det error: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    read(isockfd, &fscore, sizeof(float));                              // score
    read(isockfd, &cdecision, 1);                            // decision
    cout<<"Speaker "<<uId<<" obtained a score of: "<<fscore<<endl;
    cout<<"Decision: "<<((cdecision==RSD_ACCEPT)?"match":"no match")<<endl;
    return true;
}

/*! \fn bool I_Id (const int isockfd)
 *  \brief  send a I_ID message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_Id (const int isockfd) { 
    string uId;
	uint8_t creturn, cdecision;
	float fscore;
    
    send(isockfd, I_ID, 0, NULL);
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"I_Id error: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
	read(isockfd, &fscore, sizeof(float));                   // score
	read(isockfd, &cdecision, 1);                            // decision
    
    string uid;
	char ctmp;
	do {
		read(isockfd, &ctmp, 1);
		uid += ctmp;
	} while (ctmp != '\0');
	if (cdecision == RSD_ACCEPT)
		cout << "Match found: "<<uid<<", with a score of "<<fscore<<endl;
	else
		cout << "No match. Closest speaker: "<<uid<<", with a score of "<<fscore<<endl;
    return true;
}

/*! \fn bool I_IdGetList (const int isockfd)
 *  \brief  send a I_ID message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_IdGetList (const int isockfd) {
	string uId;
	uint8_t creturn;
	uint32_t nbSpk;
	
	send(isockfd, I_ID, 0, NULL);
	creturn=0;
	read(isockfd, &creturn, 1);                         // server's answer
	if(creturn != RSD_NO_ERROR) {
		cout<<"I_Id error: "<<__FILE__<<" "<<__LINE__<<endl;
		return false;
	}
	read(isockfd, &nbSpk, 4);
	nbSpk = ntohl(nbSpk);
	cout<<"nbSpk="<<nbSpk<<endl;
	
	string *tsUid = new string[nbSpk];
	float *ftscore = new float[nbSpk];
	uint8_t *ctdecision = new uint8_t[nbSpk];
	
	for (unsigned long lcptr=0 ; lcptr<nbSpk ; lcptr++) {
		char ctmp;
		read(isockfd, &ctmp, 1);                     // user_ID
		while (ctmp != '\0') {
			tsUid[lcptr]+=ctmp;
			read(isockfd, &ctmp, 1);
		}
		read(isockfd, ftscore+lcptr, sizeof(float));                    // score
		read(isockfd, ctdecision+lcptr, 1);              // decision
		cerr<<"Speaker "<<tsUid[lcptr]<<" obtained a score of "<<ftscore[lcptr]<<". Decision: "<<((ctdecision[lcptr]==RSD_ACCEPT)?"match":"no match")<<endl;
	}
	
	delete[] tsUid;
	delete[] ftscore;
	delete[] ctdecision;
	return true;
}

/*! \fn bool I_DetCum (const int isockfd)
 *  \brief  send a I_DETCUM message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_DetCum (const int isockfd) { 
    string uId;
    uint8_t creturn, cdecision;
    float fscore;
    
    cout<<"Enter the user_id for whom a score must be computed: ";
    cin>>uId;

    send(isockfd, I_DETCUM, uId.length()+1, (uint8_t*)uId.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                          // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"I_Det error: "<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    read(isockfd, &fscore, sizeof(float));                              // score
    read(isockfd, &cdecision, 1);                            // decision
    cout<<"Speaker "<<uId<<" obtained a score of: "<<fscore<<endl;
    cout<<"Decision: "<<((cdecision==RSD_ACCEPT)?"match":"no match")<<endl;
    return true;
}

/*! \fn bool I_IdCum (const int isockfd)
 *  \brief  send a I_IDCUM message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_IdCum (const int isockfd) {  
	string uId;
	uint8_t creturn, cdecision;
	float fscore;
	
	send(isockfd, I_ID, 0, NULL);
	creturn=0;
	read(isockfd, &creturn, 1);                         // server's answer
	if(creturn != RSD_NO_ERROR) {
		cout<<"I_Id error: "<<__FILE__<<" "<<__LINE__<<endl;
		return false;
	}
	read(isockfd, &fscore, sizeof(float));                   // score
	read(isockfd, &cdecision, 1);                            // decision
	
	string uid;
	char ctmp;
	do {
		read(isockfd, &ctmp, 1);
		uid += ctmp;
	} while (ctmp != '\0');
	if (cdecision == RSD_ACCEPT)
		cout << "Match found: "<<uid<<", with a score of "<<fscore<<endl;
	else
		cout << "No match. Closest speaker: "<<uid<<", with a score of "<<fscore<<endl;
	return true;
}

/*! \fn bool I_IdCumGetList (const int isockfd)
 *  \brief  send a I_IDCUM message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_IdCumGetList (const int isockfd) {
	string uId;
	uint8_t creturn;
	uint32_t nbSpk;
	
	send(isockfd, I_IDCUM, 0, NULL);
	creturn=0;
	read(isockfd, &creturn, 1);                         // server's answer
	if(creturn != RSD_NO_ERROR) {
		cout<<"I_Id error: "<<__FILE__<<" "<<__LINE__<<endl;
		return false;
	}
	read(isockfd, &nbSpk, 4);
	nbSpk = ntohl(nbSpk);
	cout<<"nbSpk="<<nbSpk<<endl;
	
	string *tsUid = new string[nbSpk];
	float *ftscore = new float[nbSpk];
	uint8_t *ctdecision = new uint8_t[nbSpk];
	
	for (unsigned long lcptr=0 ; lcptr<nbSpk ; lcptr++) {
		char ctmp;
		read(isockfd, &ctmp, 1);                     // user_ID
		while (ctmp != '\0') {
			tsUid[lcptr]+=ctmp;
			read(isockfd, &ctmp, 1);
		}
		read(isockfd, ftscore+lcptr, sizeof(float));                    // score
		read(isockfd, ctdecision+lcptr, 1);              // decision
		cerr<<"Speaker "<<tsUid[lcptr]<<" obtained a score of "<<ftscore[lcptr]<<". Decision: "<<((ctdecision[lcptr]==RSD_ACCEPT)?"match":"no match")<<endl;
	}
	
	delete[] tsUid;
	delete[] ftscore;
	delete[] ctdecision;
	return true;
}



/*! \fn bool I_DetCumR (const int isockfd)
 *  \brief  send a I_DetCumR message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_DetCumR (const int isockfd) {
    uint8_t creturn;
    string uId;
    
    cout<<"Enter the user_id for whom the score must be reset: ";
    cin>>uId;

    send(isockfd, I_DETCUMR, uId.length()+1, (uint8_t*)uId.c_str());
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"I_DETCUMR problem"<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Accumulator reset for "<<uId<<" speaker"<<endl;
    return true;
}

/*! \fn bool I_IdCumR (const int isockfd)
 *  \brief  send a I_IDCUMR message
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if there is no problem during Identification processing, otherwise false
 */
bool I_IdCumR (const int isockfd) { 
    uint8_t creturn;
    
    send(isockfd, I_IDCUMR, 0, NULL);
    creturn=0;
    read(isockfd, &creturn, 1);                         // server's answer
    if(creturn != RSD_NO_ERROR) {
        cout<<"I_IDCUMR problem"<<__FILE__<<" "<<__LINE__<<endl;
        return false;
    }
    cout<<"Accumulator reset"<<endl;
    return true;
}



















/*! \fn int RemoteSpkDetClient()
 *  \brief  run the client
 *
 */
int RemoteSpkDetClient() {
    bool end=false;
    short port;
    string ip;
    int isockfd;
    struct sockaddr_in servaddr;
    int a;

    cout<<"CONNECT TO  SERVER: "<<endl;
/*  cout<<"\t IP: ";
    cin>>ip;
    cout<<"\t port: ";
    cin>>port;
*/
    ip="127.0.0.1";
    port=32114;
    if ( (isockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        cerr<<"Unable to open a socket"<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }

    bzero(&servaddr, sizeof(struct sockaddr_in));
    servaddr.sin_family = AF_INET;
    servaddr.sin_port = htons(port);
    
    if (inet_pton(AF_INET, ip.c_str(), &servaddr.sin_addr) <= 0){
        cerr<<"inet_pton error"<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }
    
    if (connect(isockfd, (struct sockaddr *) &servaddr, sizeof(struct sockaddr_in)) < 0) {
        cerr<<"connect error"<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }

    while(!end) {
        cout<<endl<<"Enter your action choice (\'1\' for list): ";
        cin>>a;
        
        switch(a) {
            case G_LIST:
                printCommandList();
                break;
            case G_RESET:
                cout<<"G_RESET"<<endl;
                G_Reset(isockfd);
                break;
            case G_STATUS:
                cout<<"G_STATUS"<<endl;
                G_Status(isockfd);
                break;
            case G_SENDOPT:
                cout<<"G_SENDOPT"<<endl;
                G_SendOpt(isockfd);
                break;
            case A_RESET:
                cout<<"A_RESET"<<endl;
                A_Reset(isockfd);
                break;
            case A_SAVE:
                cout<<"A_SAVE"<<endl;
                A_Save(isockfd);
                break;
            case A_LOAD:
                cout<<"A_LOAD"<<endl;
                A_Load(isockfd);
                break;
            case A_SEND:
                cout<<"A_SEND"<<endl;
                A_Send(isockfd);
                break;
            case F_RESET:
                cout<<"F_RESET"<<endl;
                F_Reset(isockfd);
                break;
            case F_SAVE:
                cout<<"F_SAVE"<<endl;
                F_Save(isockfd);
                break;
            case F_LOAD:
                cout<<"F_LOAD"<<endl;
                F_Load(isockfd);
                break;
            case F_SEND:
                cout<<"F_SEND"<<endl;
                F_Send(isockfd);
                break;
            case M_RESET:
                cout<<"M_RESET"<<endl;
                M_Reset(isockfd);
                break;
            case M_SAVE:
                cout<<"M_SAVE"<<endl;
                M_Save(isockfd);
                break;
            case M_LOAD:
                cout<<"M_LOAD"<<endl;
                M_Load(isockfd);
                break;
            case M_WLOAD:
                cout<<"M_WLOAD"<<endl;
                M_WLoad(isockfd);
                break;
            case M_DEL:
                cout<<"M_DEL"<<endl;
                M_Del(isockfd);
                break;
            case M_ADAPT:
                cout<<"M_ADAPT"<<endl;
                M_Adapt(isockfd);
                break;
            case M_TRAIN:
                cout<<"M_TRAIN"<<endl;
                M_Train(isockfd);
                break;
            case I_DET:
                cout<<"I_DET"<<endl;
                I_Det(isockfd);
                break;
            case I_ID:
                cout<<"I_ID"<<endl;
                I_Id(isockfd);
                break;
            case I_DETCUM:
                cout<<"I_DETCUM"<<endl;
                I_DetCum(isockfd);
                break;
            case I_IDCUM:
                cout<<"I_IDCUM"<<endl;
                I_IdCum(isockfd);
                break;
            case I_DETCUMR:
                cout<<"I_DETCUMR"<<endl;
                I_DetCumR(isockfd);
                break;
            case I_IDCUMR:
                cout<<"I_IDCUMR"<<endl;
                I_IdCumR(isockfd);
                break;
            case G_QUIT:
                cout<<"G_QUIT"<<endl;
                send(isockfd, (uint8_t) G_QUIT, (uint32_t) 0, (uint8_t*) NULL);
                end=true;
                break;
            default:
                cout<<"unrecognized command"<<endl;
        } 
    }
    close(isockfd);
    
    return EXIT_SUCCESS;
}
