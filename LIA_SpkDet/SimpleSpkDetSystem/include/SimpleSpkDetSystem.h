/**
 * \file SimpleSpkDetSystem.h
 * \version 1.0
 *
 * \brief Simple (but complete) speaker detection system
 *
 **/

#ifndef __SimpleSpkDetSystem_H__
#define __SimpleSpkDetSystem_H__

#include "alize.h"
#if defined(SPRO)
extern "C" {
#include "spro.h"
}
#endif

namespace alize
{
	class SimpleSpkDetSystem:public Object
	{
	public:
        explicit SimpleSpkDetSystem(Config &config);
        ~SimpleSpkDetSystem();
		
		virtual String getClassName() const;
        
        long featureCount();
        long speakerCount();
        bool isUBMLoaded();
        vector<String> speakerIDs();
		void setOption(String opt, String optValue);
		
		void addAudio(uint32_t dataSize, uint8_t *data);
		void addAudio(String filename);
		void saveAudio(String filename);
		void resetAudio();
		
		void addFeatures(uint32_t dataSize, uint8_t *data);
		void addFeatures(String filename);
		void saveFeatures(String filename);
		void resetFeatures();
        
        void loadBackgroundModel(String fileName);
        void loadSpeakerModel(String uId, String fileName);
        void saveSpeakerModel(String uId, String fileName);
        void removeSpeaker(String uId);
        void removeAllSpeakers();
        void createSpeakerModel(String uId);
        void adaptSpeakerModel(String uId);
        
        bool verifySpeaker(String targetSpeakerId, float &resultingScore, bool withScoreAccumulation = false);
        bool identifySpeaker(String &foundSpeakerId, float &resultingScore, bool withScoreAccumulation = false);
        
        void resetAccumulatedScore(String uId);
        void resetAllAccumulatedScores();

	private:
		FeatureServer* _fs;                     ///< feature server
		MixtureServer* _ms;                     ///< mixture server
		StatServer* _ss;                        ///< stat server
		Config* _config;                        ///< configuration file
		XLine lstFeatureFile;                   ///< list of feature files loaded in the feature server
		vector<unsigned long> featureCounts;    ///< size of each feature file in the feature server
        struct ScoreAcc;
		vector<ScoreAcc> accumulatedScores;
		
#if defined(SPRO)
		int SPRO_format;
		float SPRO_Fs;
		int SPRO_channel;
		int SPRO_lswap;
		size_t SPRO_ibs;
		size_t SPRO_obs;
		float SPRO_emphco;
		float SPRO_fm_l;
		float SPRO_fm_d;
		int SPRO_win;
		unsigned short SPRO_nfilters;
		float SPRO_alpha;
		int SPRO_usemel;
		float SPRO_f_min;
		float SPRO_f_max;
		int SPRO_fftnpts;
		unsigned short SPRO_numceps;
		int SPRO_lifter;
		int SPRO_flag;
		unsigned long SPRO_winlen;
		float SPRO_escale;
		int SPRO_trace;
        
        void initSpro();
        int spro_cepstral_analysis(sigstream_t *is, spfstream_t *os, unsigned long * frameCount);
        int spro_process_audiofile(const char *ifn, char *ofn, unsigned long *frameCount);
#endif //SPRO
        
        bool parameterize_audio(String audioFileName);
		void setupTmpDirs();
		
	};
	
} // namespace alize

#endif
