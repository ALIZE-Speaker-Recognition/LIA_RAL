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
 * \version 1.5
 *
 * \brief Server for Alize use
 *
**/


#include <arpa/inet.h>
#include <unistd.h>
#include <sys/stat.h>

#include "SpkDetServerConstants.h"
#include "SpkDetServer.h"
#include "GeneralTools.h"

#if defined(SPRO)
extern "C" {
#include "spro.h"
}
#endif 

using namespace alize;
using namespace std;

typedef struct {
    String uId;
    float score;
    unsigned long nbFrame;
} Cum;


// global variables
FeatureServer* _fs;                     ///< feature server
MixtureServer* _ms;                     ///< mixture server
StatServer* _ss;                        ///< stat server
Config* _config;                        ///< configuration file
static XLine lstFeatureFile;            ///< list of feature files loaded in the feature server
Cum cumulScore[100];



#if defined(SPRO)
int SPRO_format = SPRO_SIG_PCM16_FORMAT;    ///< signal file format
float SPRO_Fs = 8000.0;                             ///< input signal sample rate
int SPRO_channel = 1;                               ///< channel to process
int SPRO_lswap = 0;                                     ///< change input sample byte order 
size_t SPRO_ibs = 10000000;                     ///< input buffer size
size_t SPRO_obs = 10000000;                     ///< output buffer size (in bytes)
float SPRO_emphco = 0.95;                       ///< pre-emphasis coefficient
float SPRO_fm_l = 20.0;                             ///< frame length in ms
float SPRO_fm_d = 10.0;                             ///< frame shift in ms
int SPRO_win = SPRO_HAMMING_WINDOW;         ///< weighting window
unsigned short SPRO_nfilters = 24;              ///< number of filters
float SPRO_alpha = 0.0;                             ///< frequency deformation parameter
int SPRO_usemel = 0;                                ///< use Usemel scale?
float SPRO_f_min = 0.0;                             ///< lower srateuency bound
float SPRO_f_max = 0.0;                             ///< higher srateuency bound
int SPRO_fftnpts = 512;                             ///< FFT length
unsigned short SPRO_numceps = 12;           ///< number of cepstral coefficients
int SPRO_lifter = 0;                                ///< liftering value
int SPRO_flag = 0;                                  ///< output stream description
unsigned long SPRO_winlen = 0;                  ///< length in frames of the CMS window
float SPRO_escale = 0.0;                            ///< energy scale factor
int SPRO_trace = 1;                                 ///< trace level

void initSpro() {
    cout << "SPro initialization" << endl;
    if (_config->existsParam("SPRO_format")) {
        if (_config->getParam("SPRO_format") == "SPRO_SIG_PCM16_FORMAT")
            SPRO_format = SPRO_SIG_PCM16_FORMAT;
        else if (_config->getParam("SPRO_format") == "SPRO_SIG_WAVE_FORMAT")
            SPRO_format = SPRO_SIG_WAVE_FORMAT;
        else if (_config->getParam("SPRO_format") == "SPRO_SIG_SPHERE_FORMAT")
#if defined(SPHERE)
            SPRO_format = SPRO_SIG_SPHERE_FORMAT;
#else
            cerr<<"SPRO_format set to Sphere in config file, but the server was not compiled with Sphere support."<<endl;
#endif //SPHERE
    }
    if (_config->existsParam("SPRO_Fs")) 
        SPRO_Fs = _config->getParam("SPRO_Fs").toDouble();
    if (_config->existsParam("SPRO_channel")) 
        SPRO_channel = _config->getParam("SPRO_channel").toLong();
    if (_config->existsParam("SPRO_lswap")) 
        SPRO_lswap = _config->getParam("SPRO_lswap").toLong();
    if (_config->existsParam("SPRO_ibs")) 
        SPRO_ibs =  _config->getParam("SPRO_ibs").toLong();
    if (_config->existsParam("SPRO_obs")) 
        SPRO_obs =  _config->getParam("SPRO_obs").toLong();
    if (_config->existsParam("SPRO_emphco"))
        SPRO_emphco = _config->getParam("SPRO_emphco").toDouble();
    if (_config->existsParam("SPRO_fm_l")) 
        SPRO_fm_l = _config->getParam("SPRO_fm_l").toDouble();
    if (_config->existsParam("SPRO_fm_d")) 
        SPRO_fm_d = _config->getParam("SPRO_fm_d").toDouble();
    if (_config->existsParam("SPRO_win")) {
        if (_config->getParam("SPRO_win").toString() == "SPRO_NULL_WINDOW")
            SPRO_win = SPRO_NULL_WINDOW;
        else if (_config->getParam("SPRO_win").toString() == "SPRO_HAMMING_WINDOW")
            SPRO_win = SPRO_HAMMING_WINDOW;
        else if (_config->getParam("SPRO_win").toString() == "SPRO_HANNING_WINDOW")
            SPRO_win = SPRO_HANNING_WINDOW;
        else if (_config->getParam("SPRO_win").toString() == "SPRO_BLACKMAN_WINDOW")
            SPRO_win = SPRO_BLACKMAN_WINDOW;
    }
    if (_config->existsParam("SPRO_nfilters")) 
        SPRO_nfilters =  _config->getParam("SPRO_nfilters").toULong();
    if (_config->existsParam("SPRO_alpha")) 
        SPRO_alpha =  _config->getParam("SPRO_alpha").toDouble();
    if (_config->existsParam("SPRO_usemel")) 
        SPRO_usemel =  _config->getParam("SPRO_usemel").toLong();
    if (_config->existsParam("SPRO_f_min")) 
        SPRO_f_min =  _config->getParam("SPRO_f_min").toDouble();
    if (_config->existsParam("SPRO_f_max")) 
        SPRO_f_max =  _config->getParam("SPRO_f_max").toDouble();
    if (_config->existsParam("SPRO_fftnpts")) 
        SPRO_fftnpts =  _config->getParam("SPRO_fftnpts").toLong();
    if (_config->existsParam("SPRO_numceps")) 
        SPRO_numceps =  _config->getParam("SPRO_numceps").toULong();
    if (_config->existsParam("SPRO_lifter")) 
        SPRO_lifter = _config->getParam("SPRO_lifter").toLong();
    if (_config->existsParam("SPRO_flag")) 
        SPRO_flag = _config->getParam("SPRO_flag").toLong();
    if (_config->existsParam("SPRO_winlen")) 
        SPRO_winlen = _config->getParam("SPRO_winlen").toULong();
    if (_config->existsParam("SPRO_escale")) 
        SPRO_escale = _config->getParam("SPRO_escale").toDouble();
    if (_config->existsParam("SPRO_trace"))
        SPRO_trace = _config->getParam("SPRO_trace").toLong();
    
    if (_config->existsParam("SPRO_add_energy"))
        SPRO_flag |= WITHE;
    if (_config->existsParam("SPRO_add_delta"))
        SPRO_flag |= WITHD;
    if (_config->existsParam("SPRO_add_acceleration"))
        SPRO_flag |= WITHA;
    
    if (_config->existsParam("SPRO_normalize"))
        SPRO_flag |= WITHR;
    if (_config->existsParam("SPRO_no_static_energy"))
        SPRO_flag |= WITHN;
    if (_config->existsParam("SPRO_cms"))
        SPRO_flag |= WITHZ;
}

/*
 * Main standard loop for filter bank analysis. COPY FROM SFBCEP.C (SPRO)
 */
int spro_cepstral_analysis(sigstream_t *is, spfstream_t *os) {
    unsigned short *idx;
    unsigned long l, d, j;
    float *w = NULL, *r = NULL;
    sample_t *buf;
    spsig_t *s;
    spf_t *e, *c;
    double energy;
    int status;
    /* ----- initialize some more stuff ----- */
    l = (unsigned long)(SPRO_fm_l * is->Fs / 1000.0);               /* frame length in samples */
    d = (unsigned long)(SPRO_fm_d * is->Fs / 1000.0);               /* frame shift in samples */
    if ((s = sig_alloc(l)) == NULL)                                     /* frame signal */
        return(SPRO_ALLOC_ERR);
    if (SPRO_win) {
        if ((buf = (sample_t *)malloc(l * sizeof(sample_t))) == NULL)   /* frame buffer */
            return(SPRO_ALLOC_ERR);
        if ((w = set_sig_win(l, SPRO_win)) == NULL) {
            free(buf); 
            sig_free(s);
            return(SPRO_ALLOC_ERR);
        }
    }
    else
        buf = s->s;
    if ((e = (spf_t *)malloc(SPRO_nfilters * sizeof(spf_t))) == NULL) {         /* filter-bank output */
        if (SPRO_win) 
            free(buf); 
        sig_free(s); 
        if (SPRO_win) 
            free(w);
        return(SPRO_ALLOC_ERR);
    }
    if ((c = (spf_t *)malloc((SPRO_numceps+1) * sizeof(spf_t))) == NULL) {
        if (SPRO_win) 
            free(buf); 
        sig_free(s); 
        free(e); 
        if (SPRO_win) 
            free(w);
        return(SPRO_ALLOC_ERR);
    }
    if (SPRO_lifter)
        if ((r = set_lifter(SPRO_lifter, SPRO_numceps)) == NULL) {
            if (SPRO_win) 
                free(buf); 
            sig_free(s); 
            free(e); 
            free(c); 
            if (SPRO_win) 
                free(w);
            return(SPRO_ALLOC_ERR);      
        }
    if (SPRO_usemel) {
        if ((idx = set_mel_idx(SPRO_nfilters, SPRO_f_min / is->Fs, SPRO_f_max / is->Fs, is->Fs)) == NULL) {
            if (SPRO_win) 
                free(buf); 
            sig_free(s); 
            free(e); 
            free(c); 
            if (SPRO_win) 
                free(w); 
            if (r) 
                free(r);
        return(SPRO_ALLOC_ERR);
        }
    }  
    else if ((idx = set_alpha_idx(SPRO_nfilters, SPRO_alpha, SPRO_f_min / is->Fs, SPRO_f_max / is->Fs)) == NULL) {
        if (SPRO_win) 
            free(buf); 
        sig_free(s); 
        free(e); 
        free(c); 
        if (SPRO_win) 
            free(w); 
        if (r) 
            free(r);
        return(SPRO_ALLOC_ERR);
    }
    /* ----- loop on each frame ----- */
    while (get_next_sig_frame(is, SPRO_channel, l, d, SPRO_emphco, buf)) {
        /* weight signal */
        if (SPRO_win)
            sig_weight(s, buf, w);
        /* apply the filter bank */
        if ((status = log_filter_bank(s, SPRO_nfilters, idx, e)) != 0) {
            if (SPRO_win) 
                free(buf); 
            sig_free(s); 
            free(e); 
            free(c); 
            if (w) 
                free(w); 
            if (r) 
                free(r); 
            free(idx);
            return(status);
        }
        /* DCT */
        if ((status = dct(e, c)) != 0) {
            if (SPRO_win) 
                free(buf); 
            sig_free(s); 
            free(e); 
            free(c); 
            if (w) 
                free(w); 
            if (r) 
                free(r); 
            free(idx);
            return(status);
        }
        /* liftering */
        if (SPRO_lifter)
            for (j = 0; j < SPRO_numceps; j++)
                *(c+j) *= *(r+j);
        /* energy */
        if (SPRO_flag & WITHE) {
            if ((energy = sig_normalize(s, 0)) < SPRO_ENERGY_FLOOR)
                energy = SPRO_ENERGY_FLOOR;
            *(c+SPRO_numceps) = (spf_t)(2.0 * log(energy));
        }
        /* write vector to stream */
        if (spf_stream_write(os, c, 1) != 1) { 
            if (SPRO_win) 
                free(buf); 
            sig_free(s); 
            free(e); 
            free(c); 
            if (w) 
                free(w); 
            if (r) 
                free(r); 
            free(idx);
            return(SPRO_FEATURE_WRITE_ERR);
        }
    }
    /* ----- reset memory ----- */
    if (SPRO_win) {
        free(buf); 
        free(w); 
    }
    sig_free(s); 
    free(e); 
    free(c); 
    if (r) 
        free(r);
    free(idx);
    
    return(0);
}
/*
 * process input file -> output file. COPY FROM SFBCEP.C (SPRO)
 */
int spro_process_audiofile(const char *ifn, char *ofn) {
    sigstream_t *is;
    spfstream_t *os;  
    float frate;
    unsigned short dim;
    int status;

    /* ----- show what was asked to do ----- */
    if (SPRO_trace) {
        fprintf(stdout, "%s --> %s\n", ifn, ofn);
        fflush(stdout);
    }
    /* ----- open input stream ----- */
    if ((is = sig_stream_open(ifn, SPRO_format, SPRO_Fs, SPRO_ibs, SPRO_lswap)) == NULL) {
        fprintf(stderr, "sfbcep error -- cannot open input stream %s\n", ifn);
        return(SPRO_STREAM_OPEN_ERR);
    }
    /* ----- open output stream ----- */
    frate = (unsigned long)(SPRO_fm_d * is->Fs / 1000.0) / is->Fs;      /* real frame period */
    dim = (SPRO_flag & WITHE) ? SPRO_numceps + 1 : SPRO_numceps;
    if ((os = spf_output_stream_open(ofn, dim, SPRO_flag & WITHE, SPRO_flag, 1.0 / frate, NULL, SPRO_obs)) == NULL) {
        fprintf(stderr, "sfbcep error -- cannot open output stream %s\n", ofn);
        sig_stream_close(is);
        return(SPRO_STREAM_OPEN_ERR);
    }
    if (SPRO_winlen)
        set_stream_seg_length(os, SPRO_winlen);
    if (SPRO_escale != 0.0)
        set_stream_energy_scale(os, SPRO_escale);
    /* ----- run cepstral analysis ----- */
    if ((status = spro_cepstral_analysis(is, os)) != 0) {
        fprintf(stderr, "sfbcep error -- error processing stream %s\n", ifn);
        sig_stream_close(is); spf_stream_close(os);
        return(status);  
    }
    /* ----- clean ----- */
    sig_stream_close(is);
    spf_stream_close(os);
    return(0);
}

#endif //SPRO

bool parameterize_audio(String audioFileName) {
#if defined(SPRO)       
    try {
        time_t t; time(&t);
        struct tm *tt= gmtime(&t);
        char tmpPrmFileBasename[40];
        char tmpPrmFileName[40];
        
        bzero(tmpPrmFileBasename, 40);
        snprintf(tmpPrmFileBasename, 39, "%02d%02d%02d_%02d%02d%02d", tt->tm_year%100, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
        bzero(tmpPrmFileName, 40);
        snprintf(tmpPrmFileName, 39, "./prm/%s.prm", tmpPrmFileBasename);

        /* ----- initialize necessary stuff ----- */
        if (fft_init(SPRO_fftnpts)) {
            cerr<<"SpkDetServer error -- cannot initialize FFT with "<<SPRO_fftnpts<<" points"<<endl;
            return false;
        }
        if (dct_init(SPRO_nfilters, SPRO_numceps)) {
            cerr<<"SpkDetServer error -- cannot initialize "<<SPRO_nfilters<<"x"<<SPRO_numceps<<" DCT kernel"<<endl;
            fft_reset(); 
            return false;
        }
        /* ----- run cepstral analysis ----- */
        if (spro_process_audiofile(audioFileName.c_str(), tmpPrmFileName)) {
            cerr<<"SpkDetServer error -- cepstral analysis failed for file "<<tmpPrmFileName<<endl;
            fft_reset(); dct_reset();
            return false;
        }
        
        /* ----- add the resulting feature file to the feature server ----- */
        lstFeatureFile.addElement(tmpPrmFileBasename);
        delete _fs;
        _fs = new FeatureServer(*_config, lstFeatureFile);
                
        return true;
    }
    catch (Exception& e) {
        cout <<"SpkDetServer error -- Parameterization failed: "<< e.toString().c_str() << endl;
        return false;
    }
#else
    cerr<<"SpkDetServer error -- Parameterization requested, but the software was not compiled with SPro."<<endl;
    return false;
#endif
}


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
        cerr<<"not enough data read: read only : "<<ss<<" instead of "<<(*size)<<endl;
    }
    
    return ss;
}

/*! \fn bool G_Reset(const int &isockfd, String filename)
 *  \brief  reset all servers
 *
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      filename        config filename
 *
 *  \return true if servers reseted, otherwise false
 */
bool G_Reset(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        _config->load(filename);
        _fs->reset();
#if defined(SPRO)
        initSpro();
#endif //SPRO
        if (_config->existsParam("inputWorldFilename")) 
            _ms->loadMixtureGD(_config->getParam("inputWorldFilename"));
        else
            _ms->reset();
        _ss->reset();
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
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if no exception throwing, otherwise false
 */
bool G_Status(const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try{
        uint32_t featureCount, spkCount;
        uint8_t UBMLoaded=1;
        long idx;
        
        featureCount = htonl(_fs->getFeatureCount());
        
        spkCount = _ms->getMixtureCount();
        idx = _ms->getMixtureIndex("UBM");
        if(idx==-1) {
            UBMLoaded=0;
            spkCount = htonl(spkCount);
        }
        else
            spkCount = htonl(spkCount-1);

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
        for (unsigned long lcptr=0 ; lcptr<_ms->getMixtureCount() ; lcptr++) {
            String stmp(_ms->getMixture(lcptr).getId());
            if (stmp!="UBM")
                write(isockfd, stmp.c_str(), stmp.length()+1);
        }
        write(isockfd, &UBMLoaded, 1);

        cerr<<"audioServer size = 0";
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      opt         option name
 *  \param[in]      optValue        option value
 *
 *  \return true if no exception throwing, otherwise false
 */
bool G_SendOpt(const int &isockfd, String opt, String optValue) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        _config->setParam(opt, optValue);
#if defined(SPRO)
        initSpro();
#endif //SPRO
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
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if no exception throwing, otherwise false
 */
bool A_Reset(const int &isockfd){
    uint8_t cc = RSD_NO_ERROR;
    try{
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      filename        filename to save audio server
 *
 *  \return true if no exception throwing, otherwise false
 */
bool A_Save(const int &isockfd, String filename) {                      // Not yet implemented
    uint8_t cc = RSD_NO_ERROR;
    try {
        ofstream fout(filename.c_str(), ofstream::binary);
        
        fout.close();
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      filename        filename of the acoustic parameters to load
 *
 *  \return true if no exception throwing, otherwise false
 */
bool A_Load(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;
    try {
        if (parameterize_audio(filename)) {
            write(isockfd, &cc, 1);
            return true;
        }
        else {
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
    }
    catch (Exception& e) {
        cc = RSD_UNDEFINED_ERROR;
        write(isockfd, &cc, 1);
        cout <<"A_Load Exception"<< e.toString().c_str() << endl;
        return false;
    }
}

/*! \fn bool A_Send(const int &isockfd, uint32_t &size, unit8_t *data)
 *  \brief  Receive an audio signal from the client, parameterize it and add it to the feature server
 *
 *  \param[in]      isockfd     socket where send data
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
        else
            write(isockfd, &cc, 1);
        bytesread += size;
        fout.write((char*)data, size);
        read_command(isockfd, &loccommand, &locsize, &locdata);
        while (locsize!=0 && loccommand==A_SEND) {
            write(isockfd, &cc, 1);
            bytesread += size;
            fout.write((char*)locdata, locsize);
            free(locdata);
            locdata=NULL;
            read_command(isockfd, &loccommand, &locsize, &locdata);
        }
        fout.close();
        
        if (parameterize_audio(tmpAudioFileName)) {
            write(isockfd, &cc, 1);
            return true;
        }
        else {
            cerr<<"Failed to parameterize audio file: "<<tmpAudioFileName<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
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
 *  \param[in]      isockfd     socket where send data
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Reset(const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        delete _fs;
        _fs = new FeatureServer();
        lstFeatureFile.reset();
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      filename        name of the feature file to save the feature server
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Save(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        if (!(_config->existsParam("saveFeatureFileFormat"))) {
            cerr<<"saveFeatureFileFormat set to : \'SPRO4\'"<<endl;
            _config->setParam("saveFeatureFileFormat", "SPRO4");
        }
        if (!(_config->existsParam("saveFeatureFileExtension"))) {
            cerr<<"saveFeatureFileExtension set to : \'.prm\'"<<endl;
            _config->setParam("saveFeatureFileExtension", ".prm");
        }
        if (!(_config->existsParam("featureFlags"))) {
            cerr<<"featureFlags set to : \'"<<_fs->getFeatureFlags().getString()<<"\'"<<endl;
            _config->setParam("featureFlags", _fs->getFeatureFlags().getString());
        }
        if (!(_config->existsParam("sampleRate"))) {
            String s;
            String ss = s.valueOf(_fs->getSampleRate());
            cerr<<"sampleRate set to : \'"<<ss<<"\'"<<endl;
            _config->setParam("sampleRate", ss);
        }       
        FeatureFileWriter w(filename, *_config);  
        outputFeatureFile(*_config, *_fs, 0, _fs->getFeatureCount(), w) ;
        w.close();
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      filename        name of the feature file
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Load(const int &isockfd, String filename) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        if (!(_config->existsParam("loadFeatureFileFormat"))) {
            cerr<<"loadFeatureFileFormat set to : \'SPRO4\'"<<endl;
            _config->setParam("loadFeatureFileFormat", "SPRO4");
        }
        if (!(_config->existsParam("loadFeatureFileExtension"))) {
            cerr<<"loadFeatureFileExtension set to : \'\'"<<endl;
            _config->setParam("loadFeatureFileExtension", "");
        }
        lstFeatureFile.addElement(filename);
        delete _fs;
        _fs = new FeatureServer(*_config, lstFeatureFile);

        if (!(_config->existsParam("loadFeatureFileVectSize"))) {
            String s;
            String ss = s.valueOf(_fs->getVectSize());
            cerr<<"loadFeatureFileVectSize set to : \'"<<ss<<"\'"<<endl;
            _config->setParam("loadFeatureFileVectSize", ss);
        }
        if (!(_config->existsParam("featureFlags"))) {
            cerr<<"featureFlags set to : \'"<<_fs->getFeatureFlags().getString()<<"\'"<<endl;
            _config->setParam("featureFlags", _fs->getFeatureFlags().getString());
        }
        if (!(_config->existsParam("sampleRate"))) {
            String s;
            String ss = s.valueOf(_fs->getSampleRate());
            cerr<<"sampleRate set to : \'"<<ss<<"\'"<<endl;
            _config->setParam("sampleRate", ss);
        }
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      size        size of first paquet (read until a paquet with a size at 0 is receive)
 *  \param[in]      data        feature data
 *
 *  \return true if no exception throwing, otherwise false
 */
bool F_Send(const int &isockfd, uint32_t &size, uint8_t *data) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
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
        char file[40];
        
        bzero(file, 40);
        snprintf(file, 39, "./prm/%02d%02d%02d_%02d%02d%02d_tmp.prm", tt->tm_year%100, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
        fout.open(file, ofstream::binary);
        fout.write((char*)data, size);
        bytesread += size;
        read_command(isockfd, &loccommand, &locsize, &locdata);
        while (locsize!=0 && loccommand==F_SEND) {
            bytesread += size;
            fout.write((char*)locdata, locsize);
            free(locdata);
            locdata=NULL;
            read_command(isockfd, &loccommand, &locsize, &locdata);
        }
        cout<<bytesread<<endl;
        fout.close();
        
        FeatureServer lfs(*_config, file);
        _config->setParam("saveFeatureFileFormat", sloadFeatureFileFormat);             // temporary featureServer loaded, config re-initialized
        if (!(_config->existsParam("saveFeatureFileExtension"))) {
            cerr<<"saveFeatureFileExtension set to : \'\'"<<endl;
            _config->setParam("saveFeatureFileExtension", "");
        }
        bzero(file, 40);
        snprintf(file, 39, "./prm/%02d%02d%02d_%02d%02d%02d.prm", tt->tm_year%100, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
        FeatureFileWriter w(file, *_config);  
        outputFeatureFile(*_config, lfs, 0, lfs.getFeatureCount(), w) ;
        w.close();

        lstFeatureFile.addElement(file);
        delete _fs;
        _fs = new FeatureServer(*_config, lstFeatureFile);

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
 *  \brief  reset all USER mixture not the world
 *
 *  \param[in]      isockfd     socket where send data

*
 *  \return true if no exception throwing, otherwise false
 */
bool M_Reset(const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        long idx = _ms->getMixtureIndex("UBM");
        if(idx==-1) {
            cerr<<"no world model"<<endl;
            _ms->reset();
        }
        else {
            MixtureGD& m=_ms->getMixtureGD(idx);
            _ms->reset();
            _ms->setMixtureId(m,"UBM"); 
        }
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      fileName        filename to save the mixture
 *  \param[in]      uId         user_id of the mixture to save
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Save (const int &isockfd, String fileName, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        if (!(_config->existsParam("saveMixtureFileFormat"))) {
            cerr<<"saveMixtureFileFormat set to : \'RAW\'"<<endl;
            _config->setParam("saveMixtureFileFormat", "RAW");
        }
        if (!(_config->existsParam("saveMixtureFileExtension"))) {
            cerr<<"saveMixtureFileExtension set to : \'.gmm\'"<<endl;
            _config->setParam("saveMixtureFileExtension", ".gmm");
        }
        long idx = _ms->getMixtureIndex(uId);
        if(idx==-1) {
            cerr<<"Mixture not found : "<<uId<<endl;        
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        _ms->getMixture(idx).save(fileName,*_config);
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      fileName        filename to load the mixture
 *  \param[in]      uId         user_id of the mixture to load
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Load (const int &isockfd, String fileName, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        if (!(_config->existsParam("loadMixtureFileFormat"))) {
            cerr<<"loadMixtureFileFormat set to : \'RAW\'"<<endl;
            _config->setParam("loadMixtureFileFormat", "RAW");
        }
        MixtureGD &m=_ms->loadMixtureGD(fileName);
        _ms->setMixtureId(m,uId);
        if(_ms->getMixtureCount()<100) {
            cumulScore[_ms->getMixtureCount()-1].uId = uId;
            cumulScore[_ms->getMixtureCount()-1].score = 0.0;
            cumulScore[_ms->getMixtureCount()-1].nbFrame = 0;
        }
        else {
            cerr<<"Not enough memory for Cumul"<<endl;
            exit(-1);
        }
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      fileName        filename to load the mixture
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_WLoad (const int &isockfd, String fileName) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        if (!(_config->existsParam("loadMixtureFileFormat"))) {
            cerr<<"loadMixtureFileFormat set to : \'RAW\'"<<endl;
            _config->setParam("loadMixtureFileFormat", "RAW");
        }
        long idx = _ms->getMixtureIndex("UBM");             //check if a previous UBM exists and if true delete it
        if(idx!=-1) {
            _ms->deleteMixtures(idx, idx);
            _ms->deleteUnusedDistribs();
        }       
        MixtureGD &m=_ms->loadMixtureGD(fileName);          //load UBM
        _ms->setMixtureId(m, "UBM"); 
            if(_ms->getMixtureCount()<100) {
            cumulScore[_ms->getMixtureCount()-1].uId = "UBM";
            cumulScore[_ms->getMixtureCount()-1].score = 0.0;
            cumulScore[_ms->getMixtureCount()-1].nbFrame = 0;
        }
        else {
            cerr<<"Not enough memory for Cumul"<<endl;
            exit(-1);
        }

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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      uId         user_id of the model to delete
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Del (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        long idx = _ms->getMixtureIndex(uId);
        if(idx==-1) {
            cerr<<"Mixture not found : "<<uId<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        _ms->deleteMixtures(idx, idx);  
        _ms->deleteUnusedDistribs();
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      uId         user_id of the model to adapt
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Adapt (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        StatServer ss(*_config, *_ms);
        long idx = _ms->getMixtureIndex(uId);
        if(idx==-1) {
            cerr<<"Mixture not found : "<<uId<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD& m=_ms->getMixtureGD(idx);
        idx = _ms->getMixtureIndex("UBM");
        if(idx==-1) {
            cerr<<"Mixture not found : UBM"<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD& world = _ms->getMixtureGD(idx);
        SegServer fakeSegServer;
        fakeSegServer.createCluster(0);
        SegCluster& fakeSeg=fakeSegServer.getCluster(0);
        fakeSeg.add(fakeSegServer.createSeg(0,_fs->getFeatureCount(),0,"",""));
        adaptModel(*_config,ss,*_ms,*_fs,fakeSeg,world,m);
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      uId         user_id of the model to train
 *
 *  \return true if no exception throwing, otherwise false
 */
bool M_Train (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        StatServer ss(*_config, *_ms);
        long idx = _ms->getMixtureIndex("UBM");
        if(idx==-1) {
            cerr<<"Mixture not found : "<<uId<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD& world = _ms->getMixtureGD(idx);
        MixtureGD& m = _ms->duplicateMixture(world,DUPL_DISTRIB);
        _ms->setMixtureId(m,uId);
        SegServer fakeSegServer;
        fakeSegServer.createCluster(0);
        SegCluster& fakeSeg=fakeSegServer.getCluster(0);
        fakeSeg.add(fakeSegServer.createSeg(0,_fs->getFeatureCount(),0,"",""));
        adaptModel(*_config,ss,*_ms,*_fs,fakeSeg,world,m);
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      uId         user_id of the model to check
 *
 *  \return true if no exception throwing, otherwise false
 */
bool I_Det (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        //StatServer _ss(*_config, *_ms);
        uint8_t decision;
        float score;
        double threshold=0.0;       
        int idx = _ms->getMixtureIndex(uId);
        if(idx==-1) {
            cerr<<"Mixture not found : "<<uId<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD &client = _ms->getMixtureGD(idx);
        idx = _ms->getMixtureIndex("UBM");
        if(idx==-1) {
            cerr<<"World model not found : "<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD &world = _ms->getMixtureGD(idx);
        _ss->resetLLK(world);
        _ss->resetLLK(client);
        _fs->seekFeature(0); 
        Feature f;
        cerr<<_fs->getFeatureCount()<<endl;
        for (unsigned long idxFrame=0;idxFrame<_fs->getFeatureCount();idxFrame++) {
            _fs->readFeature(f);
            _ss->computeAndAccumulateLLK(world, f,DETERMINE_TOP_DISTRIBS);
            _ss->computeAndAccumulateLLK(client,f,USE_TOP_DISTRIBS);
        }
        cout<<endl<<endl;
        cerr<<"world="<<_ss->getMeanLLK(world)<<" client="<<_ss->getMeanLLK(client)<<endl;
        score = _ss->getMeanLLK(client)-_ss->getMeanLLK(world);
        cerr<<"Score = "<<score<<endl;
        write(isockfd, &cc, 1);
        write(isockfd, &score, sizeof(float));
        if (_config->existsParam("threshold")) 
            threshold=(_config->getParam("threshold")).toDouble();
        if (score>threshold) 
            decision=RSD_ACCEPT; 
        else 
            decision=RSD_REJECT;
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
 *  \param[in]      isockfd     socket where send data
 *
 *  \return 
 */
bool I_Id (const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        long nbSpk = _ms->getMixtureCount();
        long idxW, idx;
        uint8_t decision;
        Feature f;
    
        idxW = _ms->getMixtureIndex("UBM");
        if(idxW==-1) {
            cerr<<"World model not found : "<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD &world = _ms->getMixtureGD(idxW);  

        float *ftscore = new float[nbSpk-1];
        unsigned long *ulIdx = new unsigned long[nbSpk-1], lcptr;

        for (lcptr=0, idx=0 ; idx<nbSpk ; idx++) {
            if (idx!=idxW) {
                MixtureGD &m = _ms->getMixtureGD(idx);  
                _ss->resetLLK(m);
                _ss->resetLLK(world);
                _fs->seekFeature(0); 
                for (unsigned long idxFrame=0;idxFrame<_fs->getFeatureCount();idxFrame++) {
                    _fs->readFeature(f);
                    _ss->computeAndAccumulateLLK(world, f,DETERMINE_TOP_DISTRIBS);
                    _ss->computeAndAccumulateLLK(m,f,USE_TOP_DISTRIBS);
                    }
                cout<<endl<<endl;
                ulIdx[lcptr] = idx;
                ftscore[lcptr++] = _ss->getMeanLLK(m)-_ss->getMeanLLK(world);
                cout<<"model "<<_ss->getMeanLLK(m)<<" world "<<_ss->getMeanLLK(world)<<" score "<<ftscore[lcptr-1]<<endl;
            }
        }
        double threshold=0.0;
        if (_config->existsParam("threshold")) 
            threshold=(_config->getParam("threshold")).toDouble();
        
        write(isockfd, &cc, 1);
        int32_t tmp = htonl(nbSpk-1);
        write(isockfd, &tmp, 4);
        for (idx=0, lcptr=0 ; idx<nbSpk ; idx++) {
            if (idx!=idxW) {
                String stmp(_ms->getMixture(idx).getId());
                cerr<<stmp.c_str()<<" ";
                write(isockfd, stmp.c_str(), stmp.length()+1);
                cerr<<ftscore[lcptr]<<endl;
                write(isockfd, &ftscore[lcptr], sizeof(float));
                if (ftscore[lcptr++]>threshold) 
                    decision=RSD_ACCEPT; 
                else 
                    decision=RSD_REJECT;
                write(isockfd, &decision, 1);        
            }
        }
        delete[] ftscore;
        delete[] ulIdx;
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
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      uId         user_id of the model to check
 *
 *  \return true if no exception throwing, otherwise false
 */
bool I_DetCum (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        uint8_t decision;
        float score;
        double threshold=0.0;       
        int idx = _ms->getMixtureIndex(uId);
        if(idx==-1) {
            cerr<<"Mixture not found : "<<uId<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD &client = _ms->getMixtureGD(idx);
        idx = _ms->getMixtureIndex("UBM");
        if(idx==-1) {
            cerr<<"World model not found : "<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD &world = _ms->getMixtureGD(idx);
        _ss->resetLLK(world);
        _ss->resetLLK(client);
        _fs->seekFeature(0); 
        Feature f;
        cerr<<_fs->getFeatureCount()<<endl;
        for (unsigned long idxFrame=0;idxFrame<_fs->getFeatureCount();idxFrame++) {
            _fs->readFeature(f);
            _ss->computeAndAccumulateLLK(world, f,DETERMINE_TOP_DISTRIBS);
            _ss->computeAndAccumulateLLK(client,f,USE_TOP_DISTRIBS);
        }
        cout<<endl<<endl;
        cerr<<"world="<<_ss->getMeanLLK(world)<<" client="<<_ss->getMeanLLK(client)<<endl;
        score = _ss->getMeanLLK(client)-_ss->getMeanLLK(world);
        cerr<<"Score = "<<score<<endl;
        for (unsigned long idxSpk=0 ; idxSpk<_ms->getMixtureCount() ; idxSpk++ ) {
            if (uId==cumulScore[idxSpk].uId) {
                float ratio= _fs->getFeatureCount()/(_fs->getFeatureCount()+cumulScore[idxSpk].nbFrame);
                cumulScore[idxSpk].score = ratio*score + (1-ratio)*cumulScore[idxSpk].score;
                cumulScore[idxSpk].nbFrame+=_fs->getFeatureCount();
                score = cumulScore[idxSpk].score;
            }
        }
        write(isockfd, &cc, 1);
        write(isockfd, &score, sizeof(float));
        if (_config->existsParam("threshold")) 
            threshold=(_config->getParam("threshold")).toDouble();
        if (score>threshold) 
            decision=RSD_ACCEPT; 
        else 
            decision=RSD_REJECT;
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
 *  \param[in]      isockfd     socket where send data
 *
 *  \return 
 */
bool I_IdCum (const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        long nbSpk = _ms->getMixtureCount();
        long idxW, idx;
        uint8_t decision;
        Feature f;
    
        idxW = _ms->getMixtureIndex("UBM");
        if(idxW==-1) {
            cerr<<"World model not found : "<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        MixtureGD &world = _ms->getMixtureGD(idxW);  

        float *ftscore = new float[nbSpk-1];
        unsigned long *ulIdx = new unsigned long[nbSpk-1], lcptr;

        for (lcptr=0, idx=0 ; idx<nbSpk ; idx++) {
            if (idx!=idxW) {
                MixtureGD &m = _ms->getMixtureGD(idx);  
                _ss->resetLLK(m);
                _ss->resetLLK(world);
                _fs->seekFeature(0); 
                for (unsigned long idxFrame=0;idxFrame<_fs->getFeatureCount();idxFrame++) {
                    _fs->readFeature(f);
                    _ss->computeAndAccumulateLLK(world, f,DETERMINE_TOP_DISTRIBS);
                    _ss->computeAndAccumulateLLK(m,f,USE_TOP_DISTRIBS);
                }
                cout<<endl<<endl;
                ulIdx[lcptr] = idx;
                ftscore[lcptr++] = _ss->getMeanLLK(m)-_ss->getMeanLLK(world);
                cout<<"model "<<_ss->getMeanLLK(m)<<" world "<<_ss->getMeanLLK(world)<<" score "<<ftscore[lcptr-1]<<endl;

                for (long idxCum=0 ; idxCum<nbSpk ; idxCum++) {
                    if (_ms->getMixture(idxCum).getId()==cumulScore[idxCum].uId) {
                        float ratio= _fs->getFeatureCount()/(_fs->getFeatureCount()+cumulScore[idxCum].nbFrame);
                        cumulScore[idxCum].score = ratio*ftscore[idx] + (1-ratio)*cumulScore[idxCum].score;
                        cumulScore[idxCum].nbFrame+=_fs->getFeatureCount();
                        ftscore[idx] = cumulScore[idxCum].score;
                    }
                }
            }
        }

        double threshold=0.0;
        if (_config->existsParam("threshold")) 
            threshold=(_config->getParam("threshold")).toDouble();
        
        write(isockfd, &cc, 1);
        int32_t tmp = htonl(nbSpk-1);
        write(isockfd, &tmp, 4);
        for (idx=0, lcptr=0 ; idx<nbSpk ; idx++) {
            if (idx!=idxW) {
                String stmp(_ms->getMixture(idx).getId());
                cerr<<stmp.c_str()<<" ";
                write(isockfd, stmp.c_str(), stmp.length()+1);
                cerr<<ftscore[lcptr]<<endl;
                write(isockfd, &ftscore[lcptr], sizeof(float));
                if (ftscore[lcptr++]>threshold) 
                    decision=RSD_ACCEPT; 
                else 
                    decision=RSD_REJECT;
                write(isockfd, &decision, 1);        
            }
        }
        delete[] ftscore;
        delete[] ulIdx;
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
 *  \brief  Compute the score between features in memory and a given user and accumulate it
 *
 *  \param[in]      isockfd     socket where send data
 *  \param[in]      uId         user_id of the model to check
 *
 *  \return true if no exception throwing, otherwise false
 */
bool I_DetCumR (const int &isockfd, String uId) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        for (unsigned long idxSpk=0 ; idxSpk<_ms->getMixtureCount() ; idxSpk++ ) {
            if (uId==cumulScore[idxSpk].uId) {
                cumulScore[idxSpk].score = 0.0;
                cumulScore[idxSpk].nbFrame =0;
            }
        }
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
 *  \brief  reset accumulator for all speakers
 *
 *  \param[in]      isockfd     socket where send data
 *
 *  \return 
 */
bool I_IdCumR (const int &isockfd) {
    uint8_t cc = RSD_NO_ERROR;                                          // server answer by default (when all is OK)
    try {
        long idxW = _ms->getMixtureIndex("UBM");
        if(idxW==-1) {
            cerr<<"World model not found : "<<endl;
            cc = RSD_UNDEFINED_ERROR;
            write(isockfd, &cc, 1);
            return false;
        }
        for (unsigned long idx=0 ; idx<_ms->getMixtureCount() ; idx++ ) {
            if (idx!=(unsigned long)idxW) {
                cumulScore[idx].score = 0.0;
                cumulScore[idx].nbFrame =0;
            }
        }
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
    
    // Check the existence of directories where audio and features will be written
    if (access("./audio", R_OK|W_OK|W_OK) != 0) {
        if (errno == ENOENT) {
            if (mkdir("./audio",0777) != 0) {
                cerr<<"Directory ./audio does not exist and could not be created. Please create it and relaunch the server."<<endl;
                exit(EXIT_FAILURE);
            }
        } else {
            cerr<<"Directory ./audio cannot be accessed. Please make it accessible for reading and writing and relaunch the server."<<endl;
            exit(EXIT_FAILURE);
        }
    }
    if (access("./prm", R_OK|W_OK|W_OK) != 0) {
        if (errno == ENOENT) {
            if (mkdir("./prm",0777) != 0) {
                cerr<<"Directory ./prm does not exist and could not be created. Please create it and relaunch the server."<<endl;
                exit(EXIT_FAILURE);
            }
        } else {
            cerr<<"Directory ./prm cannot be accessed. Please make it accessible for reading and writing and relaunch the server."<<endl;
            exit(EXIT_FAILURE);
        }
    }
            
    
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
    
    if (bind(listenfd, (struct sockaddr *) &servaddr, sizeof(struct sockaddr_in)) < 0) {
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

    _fs=new FeatureServer(*_config);  
    _ms=new MixtureServer(*_config);
    _ss=new StatServer(*_config);
#if defined(SPRO)
    initSpro();
#endif //SPRO
    
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
                delete _fs ;  
                delete _ms; 
                delete _ss;
                exit(EXIT_SUCCESS);
            default :                                               // only with a client using a different version of the protocol
                cout<<"unrecognized command : "<<(int)command<<endl;
                close(connfd);
                close(listenfd);
                delete _fs ;  
                delete _ms; 
                delete _ss;
                exit(EXIT_FAILURE);
        }
        if (data) {
            free(data); 
            data=NULL;
        }
    }
}

