// SimpleSpkDetSystem.cpp
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
 * \file SimpleSpkDetSystem.cpp
 * \version 1.1
 *
 * \brief Simple (but complete) speaker detection system
 *
**/


#include <unistd.h>
#include <sys/stat.h>

#include "SimpleSpkDetSystem.h"
#include "GeneralTools.h"
#include "NormFeat.h"
#include "EnergyDetector.h"

#if defined(SPRO)
extern "C" {
#include "spro.h"
}
#endif 

using namespace alize;
using namespace std;


#if defined(SPRO)
void SimpleSpkDetSystem::initSpro() {
	
	// Default settings
	SPRO_format = SPRO_SIG_PCM16_FORMAT;   ///< signal file format
	SPRO_sampleRate = 8000.0;              ///< input signal sample rate
	SPRO_channel = 1;                      ///< channel to process
	SPRO_lswap = 0;                        ///< change input sample byte order
	SPRO_ibs = 10000000;                   ///< input buffer size
	SPRO_obs = 10000000;                   ///< output buffer size (in bytes)
	SPRO_emphco = 0.95;                    ///< pre-emphasis coefficient
	SPRO_fm_l = 20.0;                      ///< frame length in ms
	SPRO_fm_d = 10.0;                      ///< frame shift in ms
	SPRO_win = SPRO_HAMMING_WINDOW;        ///< weighting window
	SPRO_nfilters = 24;                    ///< number of filters
	SPRO_alpha = 0.0;                      ///< frequency deformation parameter
	SPRO_usemel = 0;                       ///< use Mel scale?
	SPRO_f_min = 0.0;                      ///< lower frequency bound
	SPRO_f_max = 0.0;                      ///< higher frequency bound
	SPRO_fftnpts = 512;                    ///< FFT length
	SPRO_numceps = 12;                     ///< number of cepstral coefficients
	SPRO_lifter = 0;                       ///< liftering value
	SPRO_flag = 0;                         ///< output stream description
	SPRO_winlen = 0;                       ///< length in frames of the CMS window
	SPRO_escale = 0.0;                     ///< energy scale factor
	SPRO_trace = 1;                        ///< trace level
	
    if (_config->existsParam("SPRO_format")) {
        if (_config->getParam("SPRO_format") == "SPRO_SIG_PCM16_FORMAT")
            SPRO_format = SPRO_SIG_PCM16_FORMAT;
        else if (_config->getParam("SPRO_format") == "SPRO_SIG_WAVE_FORMAT")
            SPRO_format = SPRO_SIG_WAVE_FORMAT;
        else if (_config->getParam("SPRO_format") == "SPRO_SIG_SPHERE_FORMAT")
#if defined(SPHERE)
            SPRO_format = SPRO_SIG_SPHERE_FORMAT;
#else
            cerr<<"SPRO_format set to Sphere in config file, but the system was not compiled with Sphere support."<<endl;
#endif //SPHERE
    }
    if (_config->existsParam("SPRO_sampleRate"))
        SPRO_sampleRate = _config->getParam("SPRO_sampleRate").toDouble();
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
        SPRO_usemel =  _config->getParam("SPRO_usemel").toBool();
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
int SimpleSpkDetSystem::spro_cepstral_analysis(sigstream_t *is, spfstream_t *os, unsigned long * frameCount) {
    unsigned short *idx;
    unsigned long l, d, j;
    float *w = NULL, *r = NULL;
    sample_t *buf;
    spsig_t *s;
    spf_t *e, *c;
    double energy;
    int status;
    
    *frameCount = 0;
    
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
        (*frameCount)++;
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
int SimpleSpkDetSystem::spro_process_audiofile(const char *ifn, char *ofn, unsigned long *frameCount) {
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
    if ((is = sig_stream_open(ifn, SPRO_format, SPRO_sampleRate, SPRO_ibs, SPRO_lswap)) == NULL) {
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
    if ((status = spro_cepstral_analysis(is, os, frameCount)) != 0) {
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


void SimpleSpkDetSystem::normalizeFeatures(String tmpPrmFileBasename) {
    // Setup config for energy normalization
    Config energyNormCfg(*_config);
    energyNormCfg.setParam("loadFeatureFileExtension",".init.prm");
    energyNormCfg.setParam("saveFeatureFileExtension",".en.prm");
    energyNormCfg.setParam("labelSelectedFrames","all");
    energyNormCfg.setParam("addDefaultLabel","true");
    energyNormCfg.setParam("defaultLabel","all");
    energyNormCfg.setParam("segmentalMode","false");
    energyNormCfg.setParam("vectSize","1");
    energyNormCfg.setParam("writeAllFeatures","false");
#if defined(SPRO)
    energyNormCfg.setParam("featureServerMask",String::valueOf(SPRO_numceps));
#endif
    energyNormCfg.setParam("inputFeatureFilename",tmpPrmFileBasename);
	// Normalize energy
	normFeat(energyNormCfg);
	// Register the resulting file for future deletion
	tmpFeatureFiles.push_back(energyNormCfg.getParam_featureFilesPath() + tmpPrmFileBasename + ".en.prm");

    // Setup config for energy-based speech detection
    Config energyDetectorCfg(*_config);
    energyDetectorCfg.setParam("loadFeatureFileExtension",".en.prm");
    energyDetectorCfg.setParam("featureServerMask","0");
	energyDetectorCfg.setParam("vectSize","1");
    energyDetectorCfg.setParam("minLLK","-200");
    energyDetectorCfg.setParam("maxLLK","1000");
    energyDetectorCfg.setParam("saveFeatureFileSPro3DataKind","FBCEPSTRA");
    energyDetectorCfg.setParam("addDefaultLabel","true");
    energyDetectorCfg.setParam("defaultLabel","all");
    energyDetectorCfg.setParam("labelSelectedFrames","all");
    energyDetectorCfg.setParam("labelOutputFrames","speech");
    energyDetectorCfg.setParam("saveLabelFileExtension",".lbl");
    energyDetectorCfg.setParam("nbTrainIt","8");
    energyDetectorCfg.setParam("varianceFlooring","0.0001");
    energyDetectorCfg.setParam("varianceCeiling","1.5");
    energyDetectorCfg.setParam("mixtureDistribCount","3");
    energyDetectorCfg.setParam("baggedFrameProbabilityInit","0.001");
    energyDetectorCfg.setParam("thresholdMode","meanStd");
    energyDetectorCfg.setParam("alpha","2");
    energyDetectorCfg.setParam("segmentalMode","file");
	String lblPath= energyDetectorCfg.getParam("labelFilesPath");
	SegServer segServer; // Create the segment server for dealing with selected/unselected segments
	// Run energy detector
	SegCluster &outputSeg = energyDetector(energyDetectorCfg, segServer, tmpPrmFileBasename);
	outputLabelFile(outputSeg, lblPath+tmpPrmFileBasename+".lbl", energyDetectorCfg);
	// Register the resulting file for future deletion
	tmpFeatureFiles.push_back(lblPath+tmpPrmFileBasename+".lbl");

    // Setup config for feature normalization
    Config featureNormCfg;
    featureNormCfg.setParam("featureFilesPath", _config->getParam("featureFilesPath"));
    featureNormCfg.setParam("loadFeatureFileFormat", _config->getParam("loadFeatureFileFormat"));
    featureNormCfg.setParam("saveFeatureFileFormat", _config->getParam("saveFeatureFileFormat"));
    featureNormCfg.setParam("labelFilesPath", _config->getParam("labelFilesPath"));
    featureNormCfg.setParam("featureServerMode", _config->getParam("featureServerMode"));
    featureNormCfg.setParam("featureServerMemAlloc", _config->getParam("featureServerMemAlloc"));
    featureNormCfg.setParam("bigEndian", _config->getParam("bigEndian"));
    featureNormCfg.setParam("frameLength", _config->getParam("frameLength"));
    featureNormCfg.setParam("sampleRate", _config->getParam("sampleRate"));
    featureNormCfg.setParam("segmentalMode", _config->getParam("segmentalMode"));
    featureNormCfg.setParam("loadFeatureFileBigEndian", _config->getParam("loadFeatureFileBigEndian"));
        
    featureNormCfg.setParam("featureServerBufferSize","ALL_FEATURES");
    featureNormCfg.setParam("loadFeatureFileExtension",".init.prm");
    featureNormCfg.setParam("saveFeatureFileExtension",".norm.prm");
	featureNormCfg.setParam("addDefaultLabel","false");
	featureNormCfg.setParam("defaultLabel","speech");
    featureNormCfg.setParam("labelSelectedFrames","speech");
    featureNormCfg.setParam("writeAllFeatures","true");
    featureNormCfg.setParam("inputFeatureFilename",tmpPrmFileBasename);
	// Normalize features
    normFeat(featureNormCfg);
	// Register the resulting file for future deletion
	tmpFeatureFiles.push_back(featureNormCfg.getParam_featureFilesPath() + tmpPrmFileBasename + ".norm.prm");
}


bool SimpleSpkDetSystem::parameterizeAudio(String audioFileName, const char *basename) {
#if defined(SPRO)       
    try {
		const char *featureDir;
        char tmpPrmFileBasename[40];
        char tmpPrmFileName[255];
        unsigned long frameCount;
		
		if (basename == NULL) {
			generateTmpBasename(tmpPrmFileBasename);
		} else {
			snprintf(tmpPrmFileBasename, 40, "%s", basename);
		}
		
		featureDir = _config->getParam_featureFilesPath().c_str();
        snprintf(tmpPrmFileName, 254, "%s/%s.init.prm", featureDir, tmpPrmFileBasename);

        /* ----- initialize necessary stuff ----- */
        if (fft_init(SPRO_fftnpts)) {
            cerr<<"SimpleSpkDetSystem error -- cannot initialize FFT with "<<SPRO_fftnpts<<" points"<<endl;
            return false;
        }
        if (dct_init(SPRO_nfilters, SPRO_numceps)) {
            cerr<<"SimpleSpkDetSystem error -- cannot initialize "<<SPRO_nfilters<<"x"<<SPRO_numceps<<" DCT kernel"<<endl;
            fft_reset(); 
            return false;
        }
        /* ----- run cepstral analysis ----- */
        if (spro_process_audiofile(audioFileName.c_str(), tmpPrmFileName, &frameCount)) {
            cerr<<"SimpleSpkDetSystem error -- cepstral analysis failed for file "<<tmpPrmFileName<<endl;
            fft_reset(); dct_reset();
            return false;
        }
		
		tmpFeatureFiles.push_back(tmpPrmFileName);

        /* ----- normalize the prm file ----- */
        normalizeFeatures(tmpPrmFileBasename);
		
        /* ----- add the resulting feature file to the feature server ----- */
        lstFeatureFile.addElement(tmpPrmFileBasename);
        delete _fs;
        _fs = new FeatureServer(*_config, lstFeatureFile);
        featureCounts.push_back(frameCount);
        return true;
    }
    catch (Exception& e) {
        cout <<"SimpleSpkDetSystem error -- Parameterization failed: "<< e.toString().c_str() << endl;
        return false;
    }
#else
    cerr<<"SimpleSpkDetSystem error -- Parameterization requested, but the software was not compiled with SPro."<<endl;
    return false;
#endif
}


long SimpleSpkDetSystem::featureCount() {
    return _fs->getFeatureCount();
}

long SimpleSpkDetSystem::speakerCount() {
    if (_ms->getMixtureIndex("UBM") == -1)
        return _ms->getMixtureCount();
    else
        return _ms->getMixtureCount() - 1; // Not counting the UBM
}

bool SimpleSpkDetSystem::isUBMLoaded() {
    return (_ms->getMixtureIndex("UBM") != -1);
}

vector<String> SimpleSpkDetSystem::speakerIDs() {
    vector<String> spkIDs;
    for (unsigned long i=0 ; i<_ms->getMixtureCount(); i++) {
        String aSpkID(_ms->getMixture(i).getId());
        if (aSpkID!="UBM")
            spkIDs.push_back(aSpkID);
    }
    return spkIDs;
}

/*! \fn void SimpleSpkDetSystem::setOption(const int &isockfd, String opt, String optValue)
 *  \brief  set an option of the configuration to optValue
 *
 *  \param[in]      opt         option name
 *  \param[in]      optValue        option value
 */
void SimpleSpkDetSystem::setOption(String opt, String optValue) {
    _config->setParam(opt, optValue);
#if defined(SPRO)
	initSpro();
#endif //SPRO
//	setupDirs();
	
	//TODO: see if it is worth keeping this method.
	//		The problem is that, in order to handle some new options correctly, we would have
	//		to reset the feature and/or mixture servers. What do we do with the data previously
	//		loaded in the various servers?
	//		Probably not worth the trouble.
}

/*! \fn void SimpleSpkDetSystem::decisionThreshold()
 *  \brief  Returns the theshold used in verifySpeaker() and identifySpeaker()
 *          Only speakers with scores higher than the decision threshold will be
 *          considered a match for the target speaker.
 */
double SimpleSpkDetSystem::decisionThreshold() {
	return _decisionThreshold;
}

/*! \fn void SimpleSpkDetSystem::setDecisionThreshold(double newValue)
 *  \brief  Sets the theshold used in verifySpeaker() and identifySpeaker()
 *          Only speakers with scores higher than the decision threshold will be
 *          considered a match for the target speaker.
 *
 *  \param[in]      newValue        the new decision threshold
 */
void SimpleSpkDetSystem::setDecisionThreshold(double newValue) {
	_decisionThreshold = newValue;
}




/*! \fn void SimpleSpkDetSystem::resetAudio()
 *  \brief  Remove all the audio data
 */
void SimpleSpkDetSystem::resetAudio() {
	for (int i=0; i<tmpAudioFiles.size(); i++)
		unlink(tmpAudioFiles[i].c_str());
	tmpAudioFiles.clear();
}

/*! \fn void SimpleSpkDetSystem::saveAudio(string filename)
 *  \brief  TODO :: save the audio server
 *
 *  \param[in]      filename        filename to save audio server
 */
void SimpleSpkDetSystem::saveAudio(String filename) {
	ofstream fout(filename.c_str(), ofstream::binary);
	//TODO: implement it
	fout.close();
}

/*! \fn void SimpleSpkDetSystem::addAudio(string filename)
 *  \brief  Load an audio file, parameterize it and add it to the feature server
 *
 *  \param[in]      filename        filename of the acoustic parameters to load
 */
void SimpleSpkDetSystem::addAudio(String filename) {
    if (!parameterizeAudio(filename))
        throw IOException("Cannot parameterize audio file", __FILE__, __LINE__,filename);
}

/*! \fn void SimpleSpkDetSystem::addAudio(uint32_t dataSize, void *data)
 *  \brief  Receive an audio signal (following the format specified in the configuration), parameterize it and add it to the feature server
 *
 *  \param[in]      dataSize    number of data bytes
 *  \param[in]      data        audio data
 */
void SimpleSpkDetSystem::addAudio(uint32_t dataSize, void *data) {
    if (dataSize == 0) {
        return;
    }
    
    ofstream fout;
    uint32_t bytesread=0, locsize;
    uint8_t *locdata;
    
	const char *audioDir;
	audioDir = _config->getParam_audioFilesPath().c_str();
	
	char basename[32];
	generateTmpBasename(basename);
	
    char tmpAudioFileName[255];
    snprintf(tmpAudioFileName, 254, "%s/%s.audio", audioDir, basename);
    fout.open(tmpAudioFileName, ofstream::binary);
    if (!fout.is_open())
        throw IOException("Failed to open audio file for writing", __FILE__, __LINE__,tmpAudioFileName);
	
    fout.write((char*)data, dataSize);
    fout.close();
	
	tmpAudioFiles.push_back(String(tmpAudioFileName));
    
    if (!parameterizeAudio(tmpAudioFileName, basename))
        throw IOException("Failed to parameterize audio file", __FILE__, __LINE__,tmpAudioFileName);
}

/*! \fn void SimpleSpkDetSystem::addAudio(uint32_t sampleCount, int16_t *samples)
 *  \brief  Receive an audio signal as 16-bit signed integer linear PCM, parameterize it and add it to the feature server
 *
 *  \param[in]      sampleCount   number of audio samples
 *  \param[in]      samples       audio data, as 16-bit signed integer linear PCM
 */
void SimpleSpkDetSystem::addAudio(uint32_t sampleCount, int16_t *samples) {
#if defined(SPRO)
	int currentFormat = SPRO_format;
	SPRO_format = SPRO_SIG_PCM16_FORMAT;
#endif
	addAudio(sampleCount*2, (void*)samples);
#if defined(SPRO)
	SPRO_format = currentFormat;
#endif
}



/*! \fn void SimpleSpkDetSystem::resetFeatures()
 *  \brief  Remove all the features from the feature server
 */
void SimpleSpkDetSystem::resetFeatures() {
	delete _fs;
	_fs = new FeatureServer();
	lstFeatureFile.reset();
	featureCounts = vector<unsigned long>();
	for (int i=0; i<tmpFeatureFiles.size(); i++)
		unlink(tmpFeatureFiles[i].c_str());
	tmpFeatureFiles.clear();
}

/*! \fn void SimpleSpkDetSystem::saveFeatures(String filename)
 *  \brief  Save the content of the feature server to a file
 *
 *  \param[in]      filename        name of the feature file to save the feature server
 */
void SimpleSpkDetSystem::saveFeatures(String filename) {
	if (!(_config->existsParam("saveFeatureFileFormat"))) {
		_config->setParam("saveFeatureFileFormat", "SPRO4");
	}
	if (!(_config->existsParam("saveFeatureFileExtension"))) {
		_config->setParam("saveFeatureFileExtension", ".prm");
	}
	if (!(_config->existsParam("featureFlags"))) {
		_config->setParam("featureFlags", _fs->getFeatureFlags().getString());
	}
	if (!(_config->existsParam("sampleRate"))) {
		String s;
		String ss = s.valueOf(_fs->getSampleRate());
		_config->setParam("sampleRate", ss);
	}
	FeatureFileWriter w(filename, *_config);
	outputFeatureFile(*_config, *_fs, 0, _fs->getFeatureCount(), w) ;
	w.close();
}

/*! \fn void SimpleSpkDetSystem::addFeatures(String filename)
 *  \brief  add features from a feature file
 *
 *  \param[in]      filename        name of the feature file
 */
void SimpleSpkDetSystem::addFeatures(String filename) {
	if (!(_config->existsParam("loadFeatureFileFormat"))) {
		_config->setParam("loadFeatureFileFormat", "SPRO4");
	}
	if (!(_config->existsParam("loadFeatureFileExtension"))) {
		_config->setParam("loadFeatureFileExtension", "");
	}
	lstFeatureFile.reset(); // Temporary: we work with only 1 feature file in this mode
	// TODO: handle multiple files
	lstFeatureFile.addElement(filename);
	delete _fs;
	_fs = new FeatureServer(*_config, lstFeatureFile);
	featureCounts = vector<unsigned long>(); // Temporary: same as above
	featureCounts.push_back(_fs->getFeatureCount()); // Temporary: same as above
	
	if (!(_config->existsParam("loadFeatureFileVectSize"))) {
		String s;
		String ss = s.valueOf(_fs->getVectSize());
		_config->setParam("loadFeatureFileVectSize", ss);
	}
	if (!(_config->existsParam("featureFlags"))) {
		_config->setParam("featureFlags", _fs->getFeatureFlags().getString());
	}
	if (!(_config->existsParam("sampleRate"))) {
		String s;
		String ss = s.valueOf(_fs->getSampleRate());
		_config->setParam("sampleRate", ss);
	}
}

/*! \fn void SimpleSpkDetSystem::addFeatures(uint32_t &dataSize, uint8_t *data)
 *  \brief  receive features from the client
 *
 *  \param[in]      dataSize    size of the feature data in bytes
 *  \param[in]      data        feature data
 *
 *  \return true if no exception throwing, otherwise false
 */
void SimpleSpkDetSystem::addFeatures(uint32_t dataSize, uint8_t *data) {
	ofstream fout;
	uint32_t locsize;
	uint8_t *locdata;
	
	const String sloadFeatureFileFormat(_config->getParam("loadFeatureFileFormat"));
	_config->setParam("loadFeatureFileFormat", "RAW");
	if (!(_config->existsParam("loadFeatureFileExtension"))) {
		_config->setParam("loadFeatureFileExtension", "");
	}
	
	char basename[32];
	generateTmpBasename(basename);
	
	const char *featureDir;
	featureDir = _config->getParam_featureFilesPath().c_str();
	
	char tmpFeatureFileName[255];
	snprintf(tmpFeatureFileName, 254, "%s/%s_tmp.prm", featureDir, basename);
	fout.open(tmpFeatureFileName, ofstream::binary);
	if (!fout.is_open())
		throw IOException("Failed to open temporary feature file for writing", __FILE__, __LINE__,tmpFeatureFileName);
	fout.write((char*)data, dataSize);
	fout.close();
	
	FeatureServer lfs(*_config, tmpFeatureFileName);
	_config->setParam("saveFeatureFileFormat", sloadFeatureFileFormat);             // temporary featureServer loaded, config re-initialized
	if (!(_config->existsParam("saveFeatureFileExtension"))) {
		_config->setParam("saveFeatureFileExtension", "");
	}
	FeatureFileWriter fileWriter(basename, *_config);
	outputFeatureFile(*_config, lfs, 0, lfs.getFeatureCount(), fileWriter) ;
	fileWriter.close();
	
	lstFeatureFile.reset(); // Temporary: we work with only 1 feature file in this mode
	// TODO: handle multiple files
	lstFeatureFile.addElement(tmpFeatureFileName);
	delete _fs;
	_fs = new FeatureServer(*_config, lstFeatureFile);
	featureCounts = vector<unsigned long>(); // Temporary: same as above
	featureCounts.push_back(_fs->getFeatureCount()); // Temporary: same as above
}


/*! \fn void SimpleSpkDetSystem::removeAllSpeakers()
 *  \brief  reset all USER mixtures, leaving only the world model
 *
 */
void SimpleSpkDetSystem::removeAllSpeakers() { //TODO debug this
	long idx = _ms->getMixtureIndex("UBM");
	if(idx==-1) {
		_ms->reset();
	}
	else {
		MixtureGD& m=_ms->getMixtureGD(idx);
		_ms->reset();
		_ms->setMixtureId(m,"UBM");
	}
}

/*! \fn void SimpleSpkDetSystem::saveSpeakerModel(String uId, String fileName)
 *  \brief  save the \a uId mixture to the \a filename file
 *
 *  \param[in]      uId         user_id of the mixture to save
 *  \param[in]      fileName        filename to save the mixture
 */
void SimpleSpkDetSystem::saveSpeakerModel(String uId, String fileName) {
	if (!(_config->existsParam("saveMixtureFileFormat"))) {
		_config->setParam("saveMixtureFileFormat", "RAW");
	}
	if (!(_config->existsParam("saveMixtureFileExtension"))) {
		_config->setParam("saveMixtureFileExtension", ".gmm");
	}
	long idx = _ms->getMixtureIndex(uId);
	if(idx==-1)
		throw Exception("Mixture not found", __FILE__, __LINE__);
	_ms->getMixture(idx).save(fileName,*_config);
}

/*! \fn void SimpleSpkDetSystem::loadSpeakerModel(String uId, String fileName)
 *  \brief  load the \a uId mixture from the \a filename file
 *
 *  \param[in]      uId         user_id of the mixture to load
 *  \param[in]      fileName        filename to load the mixture
 */
void SimpleSpkDetSystem::loadSpeakerModel(String uId, String fileName) {
	if (!(_config->existsParam("loadMixtureFileFormat"))) {
		_config->setParam("loadMixtureFileFormat", "RAW");
	}
	MixtureGD &m=_ms->loadMixtureGD(fileName);
	_ms->setMixtureId(m,uId);
	ScoreAcc tmp;
	tmp.uId = uId;
	tmp.score = 0.0;
	tmp.frameCount = 0;
	accumulatedScores.push_back(tmp);
}

/*! \fn void SimpleSpkDetSystem::loadBackgroundModel(String fileName)
 *  \brief  load the world model from the \a filename file
 *
 *  \param[in]      fileName        filename to load the mixture
 */
void SimpleSpkDetSystem::loadBackgroundModel(String fileName) {
	if (!(_config->existsParam("loadMixtureFileFormat"))) {
		_config->setParam("loadMixtureFileFormat", "RAW");
	}
	long idx = _ms->getMixtureIndex("UBM");             //check if a previous UBM exists and if true delete it
	if(idx!=-1) {
		_ms->deleteMixtures(idx, idx);
		_ms->deleteUnusedDistribs();
	}
	
	MixtureGD &m=_ms->loadMixtureGD(fileName);          //load UBM
	_ms->setMixtureId(m, "UBM");
	ScoreAcc tmp;
	tmp.uId = "UBM";
	tmp.score = 0.0;
	tmp.frameCount = 0;
	accumulatedScores.push_back(tmp);
}

/*! \fn void SimpleSpkDetSystem::removeSpeaker(String uId)
 *  \brief  delete the mixture with user_id to \a uId
 *
 *  \param[in]      uId         user_id of the model to delete
 */
void SimpleSpkDetSystem::removeSpeaker(String uId) {
	long idx = _ms->getMixtureIndex(uId);
	if(idx==-1)
		throw Exception("Mixture not found", __FILE__, __LINE__);
	_ms->deleteMixtures(idx, idx);
	_ms->deleteUnusedDistribs();
}

/*! \fn void SimpleSpkDetSystem::adaptSpeakerModel(String uId)
 *  \brief  adapt the \a uId mixture with the features currently in memory
 *
 *  \param[in]      uId         user_id of the model to adapt
 */
void SimpleSpkDetSystem::adaptSpeakerModel(String uId) {
	StatServer ss(*_config, *_ms);
	long idx = _ms->getMixtureIndex(uId);
	if(idx==-1)
		throw Exception("Mixture not found", __FILE__, __LINE__);
	MixtureGD& m=_ms->getMixtureGD(idx);
	idx = _ms->getMixtureIndex("UBM");
	if(idx==-1)
		throw Exception("UBM not found", __FILE__, __LINE__);
	MixtureGD& world = _ms->getMixtureGD(idx);
	
	// Create a temporary model by adapting the UBM, the same way as in M_Train
	MixtureGD& tmpModel = _ms->duplicateMixture(world,DUPL_DISTRIB);
	_ms->setMixtureId(tmpModel,"mtmp");
	SegServer fakeSegServer;
	fakeSegServer.createCluster(0);
	SegCluster& fakeSeg=fakeSegServer.getCluster(0);
	for (unsigned long i=0; i < lstFeatureFile.getElementCount(); i++) {
		fakeSeg.add(fakeSegServer.createSeg(0,featureCounts[i],0,"",lstFeatureFile.getElement(i,false)));
	}
	adaptModel(*_config,ss,*_ms,*_fs,fakeSeg,world,m);
	
	// Compute a linear interpolation between the initial speaker model and the temporary model created above
	unsigned long vectSize = world.getVectSize();
	unsigned long distribCount = world.getDistribCount();
	double alpha = 0.5;
	if (_config->existsParam("MAPAlpha")) {
		alpha = _config->getParam("MAPAlpha").toDouble();
	}
	for ( unsigned long indC=0; indC < distribCount; indC++) {
		DistribGD& t = m.getDistrib(indC);          // A priori data for a component (from the previous speaker model)
		DistribGD& c = tmpModel.getDistrib(indC);   // Statistics for the component estimated on the new data
		for (unsigned long coef=0;coef<vectSize;coef++) {
			double res=(alpha*t.getMean(coef)) +((1-alpha)*c.getMean(coef));
			t.setMean(res, coef);
		}
	}
	long idx1 = _ms->getMixtureIndex("mtmp");
	_ms->deleteMixtures(idx1, idx1);
	_ms->deleteUnusedDistribs();
}

/*! \fn void SimpleSpkDetSystem::createSpeakerModel(String uId)
 *  \brief  train a mixture with the features in memory and assigne it the user_id \a uId
 *
 *  \param[in]      uId         user_id of the model to train
 */
void SimpleSpkDetSystem::createSpeakerModel(String uId) {
	StatServer ss(*_config, *_ms);
	long idx = _ms->getMixtureIndex("UBM");
	if(idx==-1)
		throw Exception("UBM not found", __FILE__, __LINE__);
	MixtureGD& world = _ms->getMixtureGD(idx);
	MixtureGD& m = _ms->duplicateMixture(world,DUPL_DISTRIB);
	_ms->setMixtureId(m,uId);
	SegServer fakeSegServer;
	fakeSegServer.createCluster(0);
	SegCluster& fakeSeg=fakeSegServer.getCluster(0);
	for (unsigned long i=0; i < lstFeatureFile.getElementCount(); i++) {
		fakeSeg.add(fakeSegServer.createSeg(0,featureCounts[i],0,"",lstFeatureFile.getElement(i,false)));
	}
	adaptModel(*_config,ss,*_ms,*_fs,fakeSeg,world,m);
}


/*! \fn bool SimpleSpkDetSystem::verifySpeaker(String targetSpeakerId, float &resultingScore, bool withScoreAccumulation = false)
 *  \brief  Check the features in memory against a given user
 *
 *  \param[in]      targetSpeakerId   		user_id of the model to check
 *  \param[out]     resultingScore    		the score obtained by the features
 *  \param[in]      withScoreAccumulation   combine the score with previous scores for this speaker
 *
 *  \return true if the features match the given speaker, otherwise false
 */
bool SimpleSpkDetSystem::verifySpeaker(String targetSpeakerId, float &resultingScore, bool withScoreAccumulation) {
	int idx = _ms->getMixtureIndex(targetSpeakerId);
	if(idx==-1)
		throw Exception("Mixture not found: "+targetSpeakerId, __FILE__, __LINE__);
	MixtureGD &client = _ms->getMixtureGD(idx);
	idx = _ms->getMixtureIndex("UBM");
	if(idx==-1)
		throw Exception("UBM not found", __FILE__, __LINE__);
	MixtureGD &world = _ms->getMixtureGD(idx);
	_ss->resetLLK(world);
	_ss->resetLLK(client);
	_fs->seekFeature(0);
	Feature f;
	for (unsigned long idxFrame=0;idxFrame<_fs->getFeatureCount();idxFrame++) {
		_fs->readFeature(f);
		_ss->computeAndAccumulateLLK(world, f,DETERMINE_TOP_DISTRIBS);
		_ss->computeAndAccumulateLLK(client,f,USE_TOP_DISTRIBS);
	}
	resultingScore = _ss->getMeanLLK(client)-_ss->getMeanLLK(world);
	
	if (withScoreAccumulation) {
		for (unsigned long idxSpk=0 ; idxSpk<_ms->getMixtureCount() ; idxSpk++) {
			String uId = _ms->getMixtureGD(idxSpk).getId();
			if (uId==accumulatedScores[idxSpk].uId) {
				float ratio= _fs->getFeatureCount()/(_fs->getFeatureCount()+accumulatedScores[idxSpk].frameCount);
				accumulatedScores[idxSpk].score = ratio*resultingScore + (1-ratio)*accumulatedScores[idxSpk].score;
				accumulatedScores[idxSpk].frameCount +=_fs->getFeatureCount();
				resultingScore = accumulatedScores[idxSpk].score;
				break;
			}
		}
	}
	
	return (resultingScore > _decisionThreshold);
}


/*! \fn bool SimpleSpkDetSystem::identifySpeaker(String &foundSpeakerId, float &resultingScore, bool withScoreAccumulation = false)
 *  \brief  Check the features in memory against all the users and find the best match
 *
 *  \param[out]     foundSpeakerId          user_id of the best-matching speaker
 *  \param[out]     resultingScore          the score obtained by the features
 *  \param[in]      withScoreAccumulation   combine the scores with the previous scores obtained for the speakers
 *
 *  \return true if the features were found to match a speaker, otherwise false
 */
bool SimpleSpkDetSystem::identifySpeaker(String &foundSpeakerId, float &resultingScore, bool withScoreAccumulation) {
	long nbSpk = _ms->getMixtureCount();
	long idxW, idx;
	bool decision;
	Feature f;
	
	idxW = _ms->getMixtureIndex("UBM");
	if(idxW==-1)
		throw Exception("UBM not found", __FILE__, __LINE__);
	MixtureGD &world = _ms->getMixtureGD(idxW);
	
	float bestScore = -INFINITY;
	unsigned long bestScoreIndex = idxW;
	
	for (idx=0 ; idx<nbSpk ; idx++) {
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
			float score = _ss->getMeanLLK(m)-_ss->getMeanLLK(world);
			
			if (withScoreAccumulation) {
				for (unsigned long idxSpk=0 ; idxSpk<_ms->getMixtureCount() ; idxSpk++) {
					String uId = _ms->getMixtureGD(idxSpk).getId();
					if (uId==accumulatedScores[idxSpk].uId) {
						float ratio= _fs->getFeatureCount()/(_fs->getFeatureCount()+accumulatedScores[idxSpk].frameCount);
						accumulatedScores[idxSpk].score = ratio*score + (1-ratio)*accumulatedScores[idxSpk].score;
						accumulatedScores[idxSpk].frameCount +=_fs->getFeatureCount();
						score = accumulatedScores[idxSpk].score;
						break;
					}
				}
			}
			
			if (score > bestScore) {
				bestScore = score;
				bestScoreIndex = idx;
			}
		}
	}
	
	foundSpeakerId = _ms->getMixture(bestScoreIndex).getId();
	resultingScore = bestScore;
	
	return (bestScore > _decisionThreshold);
}


/*! \fn void SimpleSpkDetSystem::resetAccumulatedScore(String uId)
 *  \brief  Reset the score accumulator for the given speaker id
 *
 *  \param[in]      uId         user_id for which to reset the accumulator
 */
void SimpleSpkDetSystem::resetAccumulatedScore(String uId) {
	for (unsigned long idxSpk=0 ; idxSpk<_ms->getMixtureCount() ; idxSpk++ ) {
		if (uId==accumulatedScores[idxSpk].uId) {
			accumulatedScores[idxSpk].score = 0.0;
			accumulatedScores[idxSpk].frameCount = 0;
			break;
		}
	}
}


/*! \fn void SimpleSpkDetSystem::resetAllAccumulatedScores()
 *  \brief  Reset the score accumulator for all the speakers
 */
void SimpleSpkDetSystem::resetAllAccumulatedScores() {
	for (unsigned long idx=0 ; idx < accumulatedScores.size() ; idx++ ) {
		accumulatedScores[idx].score = 0.0;
		accumulatedScores[idx].frameCount = 0;
	}
}


/*! \fn void SimpleSpkDetSystem::setupDir(String parameterName, String defaultPath)
 *  \brief  Checks that the parameter exists in the config, otherwise sets it using the provided default path.
 *			All paths are considered relative to the working directory given to the constructor.
 *
 *  \param[in]      parameterName    Typically, one of "audioFilesPath", "featureFilesPath", etc.
 *  \param[in]      defaultPath      Default path to use if no value found in the configuration.
 */
void SimpleSpkDetSystem::setupDir(String parameterName, String defaultPath) {
	String path;
	if (_config->existsParam(parameterName))
		path = _workdirPath + _config->getParam(parameterName);
	else
		path = _workdirPath + defaultPath;
	if (!path.endsWith("/"))
		path += "/";
	_config->setParam(parameterName, path);
	
	if (access(path.c_str(), R_OK|W_OK) != 0) {
		if (errno == ENOENT) {
			if (mkdir(path.c_str(),0770) != 0) {
				throw FileNotFoundException("Directory does not exist and could not be created.", __FILE__, __LINE__, path);
			}
		} else {
			throw FileNotFoundException("Directory cannot be accessed.", __FILE__, __LINE__, path);
		}
	}
}


/*! \fn SimpleSpkDetSystem::SimpleSpkDetSystem(Config &config)
 *  \brief  Constructor
 *
 *  \param[in]      config      configuration
 */
SimpleSpkDetSystem::SimpleSpkDetSystem(Config &config, String workdirPath) {
	_config = new Config(config);
	_workdirPath = workdirPath;
	if (!_workdirPath.endsWith("/"))
		_workdirPath += "/";
	
	// Set default paths for model storage and temporary audio and feature files
	setupDir("audioFilesPath", "audio/");
	setupDir("featureFilesPath", "prm/");
	setupDir("mixtureFilesPath", "gmm/");
	setupDir("labelFilesPath", "lbl/");
	
	debug = _config->getParam_debug();
	verbose = (_config->existsParam("verbose")) ? _config->getParam("verbose").toBool() : false;
	if (_config->existsParam("verboseLevel")) {
		verboseLevel = _config->getParam("verboseLevel").toLong();
		verbose = (verboseLevel>0);
	} else {
		verboseLevel = verbose ? 1 : 0;
	}

	if (_config->existsParam("threshold"))
		_decisionThreshold = (_config->getParam("threshold")).toDouble();
	else
		_decisionThreshold = 0.0;

	_fs=new FeatureServer(*_config);
	_ms=new MixtureServer(*_config);
	_ss=new StatServer(*_config);
	if (_config->existsParam("inputWorldFilename"))
		loadBackgroundModel(_config->getParam("inputWorldFilename"));
#if defined(SPRO)
	initSpro();
#endif //SPRO
}

SimpleSpkDetSystem::~SimpleSpkDetSystem() {
	resetAudio();
	resetFeatures();
	delete _ss;
	delete _ms;
	delete _fs;
    delete _config;
}

String SimpleSpkDetSystem::getClassName() const { return "SimpleSpkDetSystem"; }


/*! \fn void SimpleSpkDetSystem::generateTmpBasename(char *buffer)
 *  \brief  Fills the provided character buffer with a "unique-enough" basename for a temporary file.
 *			The name is unique enough only in the context of the application, within an application-owned folder.
 *			I.e, not to be used for files in /tmp.
 *
 *  \param[out]     buffer    A pre-allocated buffer, long enough (at least 32 bytes) to store the generated name.
 */
void SimpleSpkDetSystem::generateTmpBasename(char *buffer) {
	time_t t; time(&t);
	struct tm *tt= gmtime(&t);
	sprintf(buffer, "%02d%02d%02d_%02d%02d%02d_%lx", tt->tm_year%100, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec, random());
}

