// SpkDetServer.cpp
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
// First version 03 May 2005 - JFB
// Completed May 2007 - A.P.
// Updated May 2017 - TM


/**
 * \file SpkDetServer.cpp
 * \author Christophe LEVY, Alexandre PRETI, Teva MERLIN
 * \version 2
 *
 * \brief Server for Alize use
 *
**/


#include <arpa/inet.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <sys/errno.h>
#include <vector>

#include "SpkDetServerConstants.h"
#include "SpkDetServer.h"
#include "GeneralTools.h"
#include "SimpleSpkDetSystem.h"

#if defined(SPRO)
extern "C" {
#include "spro.h"
}
#endif 

using namespace alize;
using namespace std;


// global variables
SimpleSpkDetSystem* worker;
Config* _config;                        ///< configuration file


/*! \fn int read_command (const int &ifd, uint8_t *command, uint32_t *size, uint8_t **data)
 *  \brief  read a command received though the socket
 *
 *  \param[in]      ifd         socket from which to read the data
 *  \param[out]     command     command number
 *  \param[out]     size        size of the data in bytes
 *  \param[out]     data        data
 *
 *  \return number of bytes read
 */
int read_command(const int &ifd, uint8_t *command, uint32_t *size, uint8_t **data) {
    size_t ss=0;
    
    read(ifd, command, 1);
    read(ifd, size, 4);
    *size = ntohl(*size);
    if (*size) {
        (*data) = (uint8_t*) malloc(*size);
        ss=read(ifd, *data, *size);
    }
    if (ss<*size) {
        cerr<<"not enough data read: read only "<<ss<<" instead of "<<(*size)<<endl;
    }
    
    return ss;
}

/*! \fn bool G_Reset(const int &isockfd, String filename)
 *  \brief  reset all servers
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      filename        config filename
 *
 *  \return true if servers reseted, otherwise false
 */
bool G_Reset(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;
    try {
        _config->load(filename);
		delete worker;
		worker = new SimpleSpkDetSystem(*_config);
        write(isockfd, &cc, 1);
        cout<<"RESET completed"<<endl;
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cerr<<"G_Reset Exception"<<e.toString().c_str()<<endl;
        return false;
    }
    return true;
}
 
/*! \fn bool G_Status(const int &isockfd)
 *  \brief  send the servers status to the client
 *
 *  \todo   receive all speaker id
 *
 *  \param[in]      isockfd     communication socket
 *
 *  \return true if no exception throwing, otherwise false
 */
bool G_Status(const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;
    try{
        uint32_t featureCount, spkCount;
        uint8_t UBMLoaded;
		
        featureCount = htonl(worker->featureCount());
        spkCount = htonl(worker->speakerCount());
		UBMLoaded = worker->isUBMLoaded();
		
        write(isockfd, &cc, 1);
        if(!write(isockfd, &featureCount, 4))
            cerr<<"can't write to the socket"<<endl;
        else
            cerr<<"featureCount = "<<ntohl(featureCount)<<endl;
        if(!write(isockfd, &spkCount, 4))
            cerr<<"can't write to the socket"<<endl;
        else
            cerr<<"spkCount = "<<ntohl(spkCount)<<endl;
        cerr<<"UBMLoaded = "<<UBMLoaded<<endl;
		
		std::vector<String> spkIDs = worker->speakerIDs();
		for (int i=0; i<spkIDs.size(); i++) {
            String stmp = spkIDs[i];
            write(isockfd, stmp.c_str(), stmp.length()+1);
        }
        write(isockfd, &UBMLoaded, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"G_Status Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool G_SendOpt(const int &isockfd, String opt, String optValue)
 *  \brief  set an option of the configuration to optValue
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      opt         option name
 *  \param[in]      optValue    option value
 *
 *  \return true if no exception throwing, otherwise false
 */
bool G_SendOpt(const int &isockfd, String opt, String optValue) {
    uint8_t cc = RSD_NO_ERROR;
    try {
        _config->setParam(opt, optValue);
		worker->setOption(opt, optValue);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"G_SendOpt Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}



/*! \fn bool A_Reset(const int &isockfd)
 *  \brief  reset the audio server
 *
 *  \param[in]      isockfd     communication socket
 *
 *  \return true if no exception throwing, otherwise false
 */
bool A_Reset(const int &isockfd){
    uint8_t cc = RSD_NO_ERROR;
    try{
		worker->resetAudio();
        write(isockfd, &cc, 1);
    }
    catch (Exception& e){ 
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"A_Reset Exception"<< e.toString().c_str() << endl;
        return false;  
    }
    return true;
}

/*! \fn bool A_Save(const int &isockfd, String filename)
 *  \brief  TODO :: save the audio server
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      filename    filename to save audio server
 *
 *  \return true if no exception throwing, otherwise false
 */
bool A_Save(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->saveAudio(filename);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"A_Save Exception"<< e.toString().c_str() << endl;
        return false;
        }
    return true;
}

/*! \fn bool A_Load(const int &isockfd, String filename)
 *  \brief  load an audio file, parameterize it and add it to the feature server
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      filename    filename of the acoustic parameters to load
 *
 *  \return true if no exception throwing, otherwise false
 */
bool A_Load(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->addAudio(filename);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"A_Load Exception"<< e.toString().c_str() << endl;
        return false;
    }
	return true;
}

/*! \fn bool A_Send(const int &isockfd, uint32_t &size, unit8_t *data)
 *  \brief  Receive an audio signal from the client, parameterize it and add it to the feature server
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      size        size of first paquet (read until a paquet with a size at 0 is receive)
 *  \param[in]      data        audio data
 *
 *  \return true if no exception throwing, otherwise false
 */
bool A_Send(const int &isockfd, uint32_t &size, uint8_t *data) {
    uint8_t cc = RSD_NO_ERROR;
    if (size == 0) {
        write(isockfd, &cc, 1);
        return true;
    }
    
    try {
        uint8_t loccommand;
        ofstream fout;
        uint32_t bytesread=0, locsize;
        uint8_t *locdata;
        
        time_t t; time(&t);
        struct tm *tt= gmtime(&t);
        char tmpAudioFileName[40];
        
        bzero(tmpAudioFileName, 40);
        snprintf(tmpAudioFileName, 39, "./audio/%02d%02d%02d_%02d%02d%02d.audio", tt->tm_year%100, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
        fout.open(tmpAudioFileName, ofstream::binary);
        if (!fout.is_open()) {
            cerr<<"Failed to open audio file for writing: "<<tmpAudioFileName<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        fout.write((char*)data, size);
        bytesread += size;
        write(isockfd, &cc, 1);
        
        read_command(isockfd, &loccommand, &locsize, &locdata);
        while ((locsize!=0) && (loccommand==A_SEND)) {
            fout.write((char*)locdata, locsize);
            bytesread += size;
            write(isockfd, &cc, 1);
            free(locdata);
            locdata=NULL;
            read_command(isockfd, &loccommand, &locsize, &locdata);
        }
        fout.close();
		
		worker->addAudio(tmpAudioFileName);
		unlink(tmpAudioFileName);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"A_Send Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}



/*! \fn bool F_Reset (const int &isockfd)
 *  \brief  Reset the feature server
 *
 *  \param[in]      isockfd     communication socket
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Reset(const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->resetFeatures();
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"F_Reset Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool F_Save (const int &isockfd, String filename)
 *  \brief  load a feature file
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      filename        name of the feature file to save the feature server
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Save(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->saveFeatures(filename);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"F_Save Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool F_Load (const int &isockfd, String filename)
 *  \brief  load a feature file
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      filename        name of the feature file
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Load(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->addFeatures(filename);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"F_Load Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool F_Send(const int &isockfd, uint32_t &size, uint8_t *data)
 *  \brief  send a feature file to the server
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      size        size of first paquet (read until a paquet with a size at 0 is receive)
 *  \param[in]      data        feature data
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Send(const int &isockfd, uint32_t &size, uint8_t *data) {
    uint8_t cc = RSD_NO_ERROR;
    if (size == 0) {
        write(isockfd, &cc, 1);
        return true;
    }
    
    try {
        uint8_t loccommand;
        ofstream fout;
        uint32_t bytesread=0, locsize;
        uint8_t *locdata;
        
        const String sloadFeatureFileFormat(_config->getParam("loadFeatureFileFormat"));
        _config->setParam("loadFeatureFileFormat", "RAW");
        if (!(_config->existsParam("loadFeatureFileExtension"))) {
            cerr<<"loadFeatureFileExtension set to : \'\'"<<endl;
            _config->setParam("loadFeatureFileExtension", "");
        }
        
        time_t t; time(&t);
        struct tm *tt= gmtime(&t);
        char tmpFeatureFileName[40];
        
        bzero(tmpFeatureFileName, 40);
        snprintf(tmpFeatureFileName, 39, "./prm/%02d%02d%02d_%02d%02d%02d_tmp.prm", tt->tm_year%100, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
        fout.open(tmpFeatureFileName, ofstream::binary);
        if (!fout.is_open()) {
            cerr<<"Failed to open temporary feature file for writing: "<<tmpFeatureFileName<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        fout.write((char*)data, size);
        bytesread += size;
        write(isockfd, &cc, 1);
        
        read_command(isockfd, &loccommand, &locsize, &locdata);
        while ((locsize!=0) && (loccommand==F_SEND)) {
            fout.write((char*)locdata, locsize);
            bytesread += size;
            write(isockfd, &cc, 1);
            free(locdata);
            locdata=NULL;
            read_command(isockfd, &loccommand, &locsize, &locdata);
        }
        fout.close();
		
		worker->addFeatures(tmpFeatureFileName);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"F_Send Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool M_Reset (const int &isockfd)
 *  \brief  reset all USER mixtures not the world
 *
 *  \param[in]      isockfd     communication socket

*
 *  \return true if no exception throwing, otherwise false
 */
bool M_Reset(const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->removeAllSpeakers();
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"M_Reset Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool M_Save (const int &isockfd, String fileName, String uId)
 *  \brief  save the \a uId mixture to the \a filename file
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      fileName        filename to save the mixture
 *  \param[in]      uId         user_id of the mixture to save
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Save (const int &isockfd, String fileName, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->saveSpeakerModel(uId, fileName);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"M_Save Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool M_Load (const int &isockfd, String fileName, String uId)
 *  \brief  load the \a uId mixture from the \a filename file
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      fileName        filename to load the mixture
 *  \param[in]      uId         user_id of the mixture to load
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Load (const int &isockfd, String fileName, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->loadSpeakerModel(uId, fileName);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"M_Load Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool M_WLoad (const int &isockfd, String fileName)
 *  \brief  load the world model from the \a filename file
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      fileName        filename to load the mixture
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_WLoad (const int &isockfd, String fileName) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->loadBackgroundModel(fileName);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"M_WLoad Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool M_Del (const int &isockfd, String uId)
 *  \brief  delete the mixture with user_id to \a uId
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      uId         user_id of the model to delete
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Del (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->removeSpeaker(uId);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"M_DEL Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool M_Adapt (const int &isockfd, String uId)
 *  \brief  adapt the \a uId mixture with feature in memory (in featureServer)
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      uId         user_id of the model to adapt
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Adapt (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->adaptSpeakerModel(uId);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"M_ADAPT Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool M_Train (const int &isockfd, String uId)
 *  \brief  train a mixture with feature in memory (in featureServer) and aasigned it the user_id \a uId
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      uId         user_id of the model to train
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Train (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->createSpeakerModel(uId);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"M_TRAIN Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}

/*! \fn bool I_Det (const int &isockfd, String uId)
 *  \brief  Compute the score between the features in memory and a given user
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      uId         user_id of the model to check
 *
 *  \return true if no exception throwing, otherwise false
 */
bool I_Det (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
        uint8_t decision;
        float score;
		
		if (worker->verifySpeaker(uId, score))
			decision=RSD_ACCEPT;
		else
			decision=RSD_REJECT;
        write(isockfd, &cc, 1);
        write(isockfd, &score, sizeof(float));
        write(isockfd, &decision, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"I_DET Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}


/*! \fn bool I_Id (const int &isockfd)
 *  \brief  Compute the score between the features in memory and a given user
 *
 *  \param[in]      isockfd     communication socket
 *
 *  \return 
 */
bool I_Id (const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		uint8_t decision;
		float score;
		String bestMatchId;
		
		if (worker->identifySpeaker(bestMatchId, score))
			decision=RSD_ACCEPT;
		else
			decision=RSD_REJECT;
		write(isockfd, &cc, 1);
		write(isockfd, &score, sizeof(float));
		write(isockfd, &decision, 1);
		write(isockfd, bestMatchId.c_str(), bestMatchId.length()+1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"I_ID Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}


/*! \fn bool I_DetCum (const int &isockfd, String uId)
 *  \brief  Compute the score between the features in memory and a given user
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      uId         user_id of the model to check
 *
 *  \return true if no exception throwing, otherwise false
 */
bool I_DetCum (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		uint8_t decision;
		float score;
		
		if (worker->verifySpeaker(uId, score, true))
			decision=RSD_ACCEPT;
		else
			decision=RSD_REJECT;
		write(isockfd, &cc, 1);
		write(isockfd, &score, sizeof(float));
		write(isockfd, &decision, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"I_DET Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}


/*! \fn bool I_IdCum (const int &isockfd)
 *  \brief  Compute the score between the features in memory and a given user
 *
 *  \param[in]      isockfd     communication socket
 *
 *  \return 
 */
bool I_IdCum (const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		uint8_t decision;
		float score;
		String bestMatchId;
		
		if (worker->identifySpeaker(bestMatchId, score, true))
			decision=RSD_ACCEPT;
		else
			decision=RSD_REJECT;
		write(isockfd, &cc, 1);
		write(isockfd, &score, sizeof(float));
		write(isockfd, &decision, 1);
		write(isockfd, bestMatchId.c_str(), bestMatchId.length()+1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"I_ID Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;

}



/*! \fn bool I_DetCumR (const int &isockfd, String uId)
 *  \brief  Reset the score accumulator for a given user
 *
 *  \param[in]      isockfd     communication socket
 *  \param[in]      uId         user_id for which to reset the accumulator
 *
 *  \return true if no exception throwing, otherwise false
 */
bool I_DetCumR (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->resetAccumulatedScore(uId);
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"I_DETCUMR Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;
}


/*! \fn bool I_IdCumR (const int &isockfd)
 *  \brief  Reset score accumulator for all speakers
 *
 *  \param[in]      isockfd     communication socket
 *
 *  \return 
 */
bool I_IdCumR (const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;
    try {
		worker->resetAllAccumulatedScores();
        write(isockfd, &cc, 1);
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"I_IDCumR Exception"<< e.toString().c_str() << endl;
        return false;
    }
    return true;

}



/*! \fn void getFilename(const char* data, char** w1, char** w2)
 *  \brief  private
 */
void getFilename(const uint8_t* data, char** w1, char** w2) {
    (*w1) = (char*)data;
    (*w2) = (char*)data + strlen(*w1) + 1;
}
    
/*! \fn int SpkDetServer(Config &config)
 *  \brief  run the server
 *
 *  \param[in]      config      configuration
 *
 */
void SpkDetServer(Config &config) {
    _config=&config;
    short port= (short)config.getParam("SocketPort").toLong();
    int listenfd, connfd;
    struct sockaddr_in servaddr;
    uint8_t command;
    uint32_t size;
    uint8_t *data=NULL;

	worker = new SimpleSpkDetSystem(*_config);

    cout<<"LAUNCH SERVER on port: "<< port<<endl;
#if defined(SPRO)
    cout<<"Compiled with support for parameterization through SPro"<<endl;
#else
    cout<<"Compiled without SPro - no parameterization support"<<endl;
#endif

    if ( (listenfd=socket(AF_INET, SOCK_STREAM, 0)) < 0 ) {
        cerr<<"Unable to open a socket : "<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }
    bzero(&servaddr, sizeof(struct sockaddr_in));
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = htonl(INADDR_ANY);
    servaddr.sin_port = htons(port);
    
	if (::bind(listenfd, (struct sockaddr *) &servaddr, sizeof(struct sockaddr_in)) < 0) {
        cerr<<"Bind error: "<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }

    if (listen(listenfd, 5) < 0) {
        cerr<<"Listen error"<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }

    cout<<"WAITING FOR CLIENT "<<endl<<endl;
    if ( (connfd = accept(listenfd, (struct sockaddr *) NULL, NULL)) < 0) {
        cerr<<"Accept error: "<<__FILE__<<" "<<__LINE__<<endl;
        exit(EXIT_FAILURE);
    }

    while (1) {
        cerr<<endl<<"waiting for order "<<endl;
        command=0;
        size=0;
        read_command(connfd, &command, &size, &data);
        switch(command) {
            case G_LIST :
                break;
            case G_RESET :
                cout<<"received G_RESET"<<endl;
                G_Reset(connfd, String((char*)data));
                break;
            case G_STATUS :
                cout<<"received G_STATUS"<<endl;
                G_Status(connfd);
                break;
            case G_SENDOPT : {
                cout<<"received G_SENDOPT"<<endl;
                char *opt, *optValue;
                getFilename(data, &opt, &optValue);
                G_SendOpt(connfd, String(opt), String(optValue));
                break; }
            case A_RESET :
                cout<<"received A_RESET"<<endl;
                A_Reset(connfd);
                break;
            case A_SAVE :
                cout<<"received A_SAVE"<<endl;
                A_Save(connfd, String((char*)data));
                break;
            case A_LOAD :
                cout<<"received A_LOAD"<<endl;
                A_Load(connfd, String((char*)data));
                break;
            case A_SEND :
                cout<<"received A_SEND"<<endl;
                A_Send(connfd, size, data);
                break;
            case F_RESET :
                cout<<"received F_RESET"<<endl;
                F_Reset(connfd);
                break;
            case F_SAVE :
                cout<<"received F_SAVE"<<endl;
                F_Save(connfd, String((char*)data));
                break;
            case F_LOAD :
                cout<<"received F_LOAD"<<endl;
                F_Load(connfd, String((char*)data));
                break;
            case F_SEND :
                cout<<"received F_SEND"<<endl;
                F_Send(connfd, size, data);
                break;
            case M_RESET :
                cout<<"received M_RESET"<<endl;
                M_Reset(connfd);
                break;
            case M_SAVE : {
                cout<<"received M_SAVE"<<endl;
                char *uId, *fileName;
                getFilename(data, &uId, &fileName);
                M_Save(connfd, String(fileName), String(uId));
                break; }
            case M_LOAD : {
                cout<<"received M_LOAD"<<endl;
                char *uId, *fileName;
                getFilename(data, &uId, &fileName);
                M_Load(connfd, String(fileName), String(uId));
                break; }
            case M_WLOAD :
                cout<<"received M_WLOAD"<<endl;
                M_WLoad(connfd, String((char*)data));
                break;
            case M_DEL :
                cout<<"received M_DEL"<<endl;
                M_Del(connfd, String((char*)data));
                break;
            case M_ADAPT :
                cout<<"received M_ADAPT"<<endl;
                M_Adapt(connfd, String((char*)data));
                break;
            case M_TRAIN :
                cout<<"received M_TRAIN"<<endl;
                M_Train(connfd, String((char*)data));
                break;
            case I_DET :
                cout<<"received I_DET"<<endl;
                I_Det(connfd, String((char*)data));
                break;
            case I_ID :
                cout<<"received I_ID"<<endl;
                I_Id(connfd);
                break;
            case I_DETCUM :
                cout<<"received I_DETCUM"<<endl;
                I_DetCum(connfd, String((char*)data));
                break;
            case I_IDCUM :
                cout<<"received I_IDCUM"<<endl;
                I_IdCum(connfd);
                break;
            case I_DETCUMR :
                cout<<"received I_DETCUMR"<<endl;
                I_DetCumR(connfd, String((char*)data));
                break;
            case I_IDCUMR :
                cout<<"received I_IDCUMR"<<endl;
                I_IdCumR(connfd);
                break;
            case G_QUIT :
                cout<<"bye bye"<<endl;
                close(connfd);
                close(listenfd);
                delete worker ;
                exit(EXIT_SUCCESS);
            default :                                               // only with a client using a different version of the protocol
                cout<<"unrecognized command : "<<(int)command<<endl;
                close(connfd);
                close(listenfd);
                delete worker ;
                exit(EXIT_FAILURE);
        }
        if (data) {
            free(data); 
            data=NULL;
        }
    }
}

