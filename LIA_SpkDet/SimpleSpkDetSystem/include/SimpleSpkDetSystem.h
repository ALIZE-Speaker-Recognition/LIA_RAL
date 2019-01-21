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
#include <vector>


namespace alize
{
	class SimpleSpkDetSystem:public Object
	{
	public:
        explicit SimpleSpkDetSystem(Config &config, String workdirPath = ".");
        ~SimpleSpkDetSystem();
		
		virtual String getClassName() const;
        
        long featureCount();
        long speakerCount();
        bool isUBMLoaded();
        vector<String> speakerIDs();
		void setOption(String opt, String optValue);
	    double decisionThreshold();
	    void setDecisionThreshold(double newValue);
		
	    void addAudio(uint32_t sampleCount, int16_t *samples);  // for 16 bit linear PCM
		void addAudio(uint32_t dataSize, void *data);           // for other formats (following the format specified in the configuration file)
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
	    String _workdirPath;                    ///< working directory (for model storage + temp files)
	    double _decisionThreshold;              ///< for score comparison in verifySpeaker and identifySpeaker
		XLine lstFeatureFile;                   ///< list of feature files loaded in the feature server
		vector<unsigned long> featureCounts;    ///< size of each feature file in the feature server
		struct ScoreAcc {
			String uId;
			float score;
			unsigned long frameCount;
		};

		vector<ScoreAcc> accumulatedScores;
	    vector<String> tmpAudioFiles;
	    vector<String> tmpFeatureFiles;
		
#if defined(SPRO)
		int SPRO_format;
		float SPRO_sampleRate;
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
		spflag_t SPRO_flag;
		unsigned long SPRO_winlen;
		float SPRO_escale;
		int SPRO_trace;
        
        void initSpro();
        int spro_cepstral_analysis(sigstream_t *is, spfstream_t *os, unsigned long * frameCount);
        int spro_process_audiofile(const char *ifn, char *ofn, unsigned long *frameCount);
#endif //SPRO

        void normalizeFeatures(String tmpPrmFileBasename);
        bool parameterizeAudio(String audioFileName, const char *basename = NULL);
		void setupDir(String parameterName, String defaultPath);
	    void generateTmpBasename(char *buffer);
		
	};
	
} // namespace alize

#endif
