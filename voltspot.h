#ifndef __VOLTSPOT_H_
#define __VOLTSPOT_H_

#include "util.h"

/* global configuration parameters for HotSpot	*/
typedef struct global_config_t_st
{
	/* floorplan input file */
	char flp_file[STR_SIZE];
	/* input power trace file */
	char p_infile[STR_SIZE];
	/* output file for the voltage trace */
	char v_outfile[STR_SIZE];
	/* input configuration parameters from file	*/
	char config[STR_SIZE];
	/* output configuration parameters to file	*/
	char dump_config[STR_SIZE];
}global_config_t;

/* 
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void global_config_from_strs(global_config_t *config, str_pair *table, int size);

#endif
