#ifndef __FLP_H_
#define __FLP_H_

#include "util.h"

#define MAX_UNITS		8192

/* 
 * chip edge has true dead space, which is 
 * modeled by the following blocks
 */
#define RIM_LEFT		1
#define RIM_RIGHT		2
#define RIM_TOP			4
#define RIM_BOTTOM		8
#define RIM_PREFIX		"RIM"
#define RIM_LEFT_STR	RIM_PREFIX"_left"
#define RIM_RIGHT_STR	RIM_PREFIX"_right"
#define RIM_TOP_STR		RIM_PREFIX"_top"
#define RIM_BOTTOM_STR	RIM_PREFIX"_bottom"

/* prefix denoting dead block	*/
#define DEAD_PREFIX		"_"

/* flags denoting orientation	*/
/* rotated orientations	*/
#define	ROT_0		0x01	/* normal	*/
#define	ROT_90		0x02	/* 90 degrees anticlockwise */
#define	ROT_180		0x04	/* 180 degrees anticlockwise */
#define	ROT_270		0x08	/* 270 degrees anticlockwise */
/* flipped + rotated orientations	*/
#define	FLIP_0		0x10	/* flip about y axis of ROT_0	*/
#define	FLIP_90		0x20	/* flip about y axis of ROT_90	*/
#define	FLIP_180	0x40	/* flip about y axis of ROT_180	*/
#define	FLIP_270	0x80	/* flip about y axis of ROT_270	*/
#define ORIENTS_N	8		/* total no. of orientations	*/
#define ORIENT_NDEF	0xDF	/* undefined orientation	*/

/* type for holding the above flags	*/
//typedef unsigned char orient_t;

/* placed functional unit */
typedef struct PDN_unit_t_st
{
	char name[STR_SIZE];
	int    domain; // For multi-domain PDN
	double width;
	double height;
	double leftx;
	double bottomy;
}PDN_unit_t;

/* floorplan data structure	*/
typedef struct PDN_flp_t_st
{
	PDN_unit_t *units;
	int n_units;
  int n_domain;
  /* density of wires between units	*/
  //double **wire_density;
} PDN_flp_t;

/* skip floorplanning and read floorplan directly from file */
PDN_flp_t *PDN_read_flp(char *file);
/* 
 * print the floorplan in a FIG like format 
 * that can be read by tofig.pl to produce 
 * an xfig output 
 */
//void print_flp_fig (PDN_flp_t *flp);
///* debug print	*/
//void print_flp (PDN_flp_t *flp);
///* translate the floorplan to new origin (x,y)	*/
void PDN_flp_translate(PDN_flp_t *flp, double x, double y);
///* scale the floorplan by a factor 'factor'	*/
//void flp_scale(PDN_flp_t *flp, double factor);
/* 
 * change the orientation of the floorplan by
 * rotating and/or flipping. the target orientation
 * is specified in 'target'. 'width', 'height', 'xorig'
 * and 'yorig' are those of 'flp' respectively.
 */
//void flp_change_orient(PDN_flp_t *flp, double xorig, double yorig,
//					   double width, double height, orient_t target);

/* dump the floorplan onto a file	*/
//void dump_flp(PDN_flp_t *flp, char *file, int dump_connects);
/* memory uninitialization	*/
void PDN_free_flp(PDN_flp_t *flp);

/* placed floorplan access routines	*/

///* get unit index from its name	*/
int PDN_get_blk_index(PDN_flp_t *flp, char *name);
///* are the units horizontally adjacent?	*/
//int is_horiz_adj(PDN_flp_t *flp, int i, int j);
///* are the units vertically adjacent?	*/
//int is_vert_adj (PDN_flp_t *flp, int i, int j);
///* shared length between units	*/
//double get_shared_len(PDN_flp_t *flp, int i, int j);
///* total chip width	*/
double PDN_get_total_width(PDN_flp_t *flp);
///* total chip height */
double PDN_get_total_height(PDN_flp_t *flp);
///* x and y origins	*/
//double get_minx(PDN_flp_t *flp);
//double get_miny(PDN_flp_t *flp);
///* precondition: L2 should have been wrapped around	*/
//double get_core_width(PDN_flp_t *flp, char *l2_label);
///* precondition: L2 should have been wrapped around	*/
//double get_core_height(PDN_flp_t *flp, char *l2_label);
///* other queries	*/
//double get_manhattan_dist(PDN_flp_t *flp, int i, int j);
//double get_total_area(PDN_flp_t *flp);
//double get_core_area(PDN_flp_t *flp, char *l2_label);
//double get_core_occupied_area(PDN_flp_t *flp, char *l2_label);

#endif
