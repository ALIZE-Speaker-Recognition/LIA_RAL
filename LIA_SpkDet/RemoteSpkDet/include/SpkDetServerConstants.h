/**
 * \file SpkDetServerConstants.h
 * \author Christophe LEVY, Alexandre PRETI
 * \version 1.0
 * \date september, 13th 2007
 *
 * \brief Server for Alize use
 *
**/

#ifndef __SpkDet_ServerConstants_H__
#define __SpkDet_ServerConstants_H__

/* Codes for the commands */

#define G_QUIT          0      ///< close the socket and quit
#define G_LIST          1      ///< print the list of supported commands
#define G_RESET         2      ///< Restart LIA_SpkDet, reading the configuration file
#define G_STATUS        3      ///< Request status information
#define G_SENDOPT       4      ///< Set an option of the configuration to a value.

#define A_RESET         10     ///< not really implemented
#define A_SAVE          11     ///< not really implemented
#define A_LOAD          12     ///< not really implemented
#define A_SEND          13     ///< not really implemented

#define F_RESET         30     ///< Reset the feature server
#define F_SAVE          31     ///< Save the feature in a file (if the buffer is large enough)
#define F_LOAD          32     ///< Read a feature file
#define F_SEND          33     ///< Send feature data to the server

#define M_RESET         50     ///< Suppress all the user models (not the world)
#define M_SAVE          51     ///< Save a model in a file
#define M_LOAD          52     ///< Read a model from a file
#define M_WLOAD         53     ///< Read the world model from a file
#define M_DEL           54     ///< Remove from the memory the user model
#define M_ADAPT         55     ///< Adapt an existing model, using the features in the memory
#define M_TRAIN         56     ///< Create a model for an user, using the feature in memory

#define I_DET           70     ///< Compute the score between the features in memory and a given user
#define I_ID            71     ///< Compute the score between the features in memory and all the users in memory
#define I_DETCUM        72     ///< 
#define I_IDCUM         73     ///< 
#define I_DETCUMR       74     ///< 
#define I_IDCUMR        75     ///< 
#define I_IDCUMGETLIST  76     ///< 


/* Return values for an identification decision */

#define RSD_ACCEPT          1
#define RSD_REJECT          0


/* Error codes */

#define RSD_NO_ERROR        0
#define RSD_UNDEFINED_ERROR 1

#endif
