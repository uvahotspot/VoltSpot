#ifndef __PDN_SIM_H_
#define __PDN_SIM_H_

#include "flp.h"
#include "slu_ddefs.h"

#define LAYER_VDD		    1
#define LAYER_GND	    	0

// location configurations
#define CONFIG_CUS  		0
#define CONFIG_ALL  		1
// Pad
#define PADCONFIGS  		2
#define PAD_NUM_UNKNOWN	0
#define PAD_NUM_KNOWN  	1
// TSV
#define TSVCONFIGS  		2
// IVR
#define IVRLOC  		    1

#define STEADY          0
#define TRANSIENT       1

// for steady state
#define PDN_STEADY_EXTRA 2
#define PKG_GND          0
#define PKG_VDD          1

// TSV location definition
#define VDDTSV           3
#define GNDTSV           5
#define VTHTSV           7
#define GTHTSV           9

#define DEFAULT_MLAYERS	4	
/* metal layer configuration file constants */
#define MLCF_NPARAMS		6	/* no. of parameters per layer	*/
#define MLCF_SNO				0	/* serial number */
#define MLCF_PITCH      1	/* metal pitch */
#define MLCF_WIDTH      2	/* metal width */
#define MLCF_THICK      3	/* metal thickness */
#define MLCF_RHO        4	/* metal resistivity */
#define MLCF_DIREC      5	/* metal layer direction*/

#define MLCF_X          0	/* Horizontal wire */
#define MLCF_Y          1	/* Vertical wire */

/* 3D layer configuration file constants */
#define LCF_3D_NPARAMS	4	/* no. of parameters per layer	*/
#define LCF_3D_SNO			0	/* serial number */
#define LCF_3D_FLP      1	/* floorplan file */
#define LCF_3D_MLCF     2	/* metal layer spec file */
#define LCF_3D_TSVL     3	/* TSV loc file */

// for non-uniform de-cap placement, portion of decaps located in hot area
#define LOCAL_CAP_RATIO  0.9

// maximum number of cores allowed
#define MAX_CORE_NUM     32

// maximum grid dimension
#define MAX_DIM	    	   1000

/* block list: block to grid mapping data structure.
 * list of blocks mapped to a grid cell	
 */
typedef struct blist_t_st
{
	/* index of the mapped block	*/
	int idx;
	/* ratio of this block's area within the grid cell 
	 * to the total area of the grid cell
	 */
	double occupancy;
	/* next block mapped to the same cell	*/
	struct blist_t_st *next;
}blist_t;

typedef struct PDN_C4_t_st 
{
	/* pad grid size	*/
	int pad_grid_col;
	int pad_grid_row;

	/* pad count	*/
	int vdd_num;
	int gnd_num;

	/* pad R and L	*/
  double pad_r;
  double pad_l;

	/* pad current limit */
  double Ith;

	/* pad location	*/
	/* bit mark, see PGPAD, IOPAD, IORAN */
	int **vdd_loc;
	int **gnd_loc;

  /* Array stores pad index */
  int *vdd_padidx;
  int *gnd_padidx;

}PDN_C4_t;

typedef struct metal_geo_t_st 
{
	double pitch;
	double width;
	double thick;
	double rho;
	int    direc;
}metal_geo_t;

typedef struct metal_gridRL_t_st 
{
	double r;	/* onchip resistors	*/
	double l;	/* onchip inductance	*/
}metal_gridRL_t;

typedef struct PDN_metal_t_st 
{
	int n_metal;

  metal_geo_t    *geo;
  metal_gridRL_t *gridRL;

}PDN_metal_t;

typedef struct IVR_t_st 
{
	int num_IVR;

  double R_drop; // static model

  double freq;
  double Ron_top;
  double Ron_bottom;
  double c;

  int **loc;
}IVR_t;

typedef struct PDN_TSV_t_st 
{
	double w;
	double h;
	double t;

  double r;
  double l;
  //double c;

  int num_vdd;
  int num_gnd;

	char file[STR_SIZE];
	int **loc; // TSVs connecting ith and i+1th layer
}PDN_TSV_t;

typedef struct PDN_layer_t_st 
{
	/* configuration parameters	*/
	int no;				  /* serial number	*/
	/* floorplan */
	PDN_flp_t *flp;
	/* metal layer information	*/
	PDN_metal_t metal_layers;
	/* TSV info */
  PDN_TSV_t tsv;
	/* On-Chip capacitance */
	double **cap_c;
	/* block-grid map - 2-d array of block lists	*/
	blist_t ***b2gmap;
}PDN_layer_t;

/* PDN config */
typedef struct PDN_config_t_st
{
  int run_PDN;
	int PDN_grid_intv;

	int PDN_padconfig;
	int padloc_format;
	double PDN_padpitch;
	double PDN_padD;
	double PDN_padR;

  int PDN_pkgLC;
  int PDN_gridL;
  int vgradient_analyse;
	double vdd;
	double gnd;

	double proc_clock_freq;
	int PDN_step_percycle;
	int ptrace_sampling_intvl;
	int PDN_ptrace_warmup;	
  double PDN_decap_dense;
  double PDN_decap_ratio;
  int PDN_decap_unifm;
  double PDN_padL;
  double PDN_pkg_sL;
  double PDN_pkg_sR;
  double PDN_pkg_C;
  double PDN_pkg_pR;
  double PDN_pkg_pL;
  double PDN_pkg_scale;

	int v_stacking;
	int TSV_config;
	double TSV_R;
	double TSV_L;
	char layer_file_3D[STR_SIZE];
	char IVR_loc_file[STR_SIZE];
  double SC_freq;
  double SC_totcap;
  double SC_Rontop;
  double SC_Ronbtm;

	char mlayer_spec_file[STR_SIZE];
	char padloc_file_in[STR_SIZE];
	char padloc_file_out[STR_SIZE];
	char padcur_file[STR_SIZE];
	char tsvcur_file[STR_SIZE];
	char gridvol_file[STR_SIZE];
	char vio_file[STR_SIZE];
	char trans_vgradient_file[STR_SIZE];
	char senloc_file[STR_SIZE];
	char node_viotrace_file[STR_SIZE];

	double PDN_cur_dense;
	double PDN_noise_th;

	int reserve_io;
  int MC_pads;
  double IO_dense;

  int    PDN_multi_dom;
  int    animation;
  int    frame_intv;
  double legend_lwr;
  double legend_upr;
  double legend_curupr;
  int    PDN_sin_pattern;
  int    PDN_sin_totstep;
  double PDN_sin_freq;

}PDN_config_t;

/* PDN model's internal vector datatype	*/
typedef struct model_vector_t_st
{
	/* array of 3-d grid of nodes	*/
	double ***cuboid;

}model_vector_t;

/* PDN model	*/
typedef struct model_t_st
{
	/* configuration */
	PDN_config_t config;

	/* dimensions	*/
	double width;
	double height;
	/* grid resolution	*/
	int rows;
	int cols;

	/* is PDN 3D */
	int is_3D;

	/* layer information	*/
	PDN_layer_t *layers;
	int n_layers;

	/* C4 Pads information	*/
	PDN_C4_t *c4;

	/* IVR information	*/
	IVR_t *sc_converter;

	/* sum total of the functional blocks of all floorplans	*/
	int total_n_blocks;

	/* internal state - most recently computed 
	 * steady state voltages
	 */
	double *last_steady;

	/* for trans compute */
  int trans_matrix_dim;
  double *last_trans;
  double *last_power;

}model_t;

/* config setup */
PDN_config_t default_PDN_config(void);
void PDN_config_add_from_strs(PDN_config_t *config, str_pair *table, int size);

/* allocate/initiate/destruct */
model_vector_t *new_model_vector(model_t *model);
model_t *alloc_model(PDN_config_t *config, PDN_flp_t *flp_default);
PDN_C4_t *alloc_C4_PDN(model_t *model, int pad_grid_col, int pad_grid_row);
IVR_t *alloc_IVR(model_t *model);
double ***dcuboid_PDN(int nr, int nc, int nl);
void populate_R_model_PDN(model_t *model);
void populate_LC_model_PDN(model_t *model);
void populate_C4_PDN(model_t *model);
void populate_TSV_PDN(model_t *model);
void populate_IVR(model_t *model);
void alloc_layers_PDN(model_t *model, PDN_flp_t *flp_default);
PDN_metal_t alloc_metal_layers(int n_metal);
void populate_single_layer_PDN(model_t *model, PDN_flp_t *flp_default);
void populate_default_mlayers(model_t *model, int layer);
void parse_layer_file_PDN(model_t *model, FILE *fp);
void parse_tsv_loc(model_t *model, int layer);
void parse_IVR_loc(model_t *model, char *file);
void delete_model(model_t *model);
void set_current_limit(model_t *model);
void free_layers_PDN(model_t *model);
void free_metal_layers(PDN_metal_t l);
void free_TSV(PDN_TSV_t t);
void free_IVR(IVR_t *v);
void free_PDN_C4(PDN_C4_t *p);
void free_model_vector(model_vector_t *v);

/* block, grid and flp */
void PDN_xlate_vector_b2g(model_t *model, double *b, model_vector_t *g);
void set_bgmap_PDN(model_t *model, PDN_layer_t *layer);
void reset_b2gmap_PDN(model_t *model, PDN_layer_t *layer);
blist_t ***new_b2gmap_PDN(int rows, int cols);
blist_t *new_blist_PDN(int idx, double occupancy);
double blist_avg_PDN(blist_t *ptr, PDN_flp_t *flp, double *v);
void blist_append_PDN(blist_t *head, int idx, double occupancy);
void delete_b2gmap_PDN(blist_t ***b2gmap, int rows, int cols);

/* Multi branch */
void parse_metal_layer_file(model_t *model, char *file, int layer);

/* decap distribution */
int decap_under_blk(model_t *model, int idx, int layer);
double get_onchip_cap(model_t *model, int idx, int layer);

/* Steady state */
void steady_state_PDN(model_t *model, double *power);
void steady_state_PDN_regular(model_t *model, double *power);
void steady_state_PDN_vs(model_t *model, double *power);
SuperMatrix build_steady_grid_matrix(model_t *model);
SuperMatrix build_steady_rhs_vector(model_t *model, model_vector_t *power, double **rhs);
SuperMatrix build_steady_grid_matrix_vs(model_t *model, model_vector_t *power);
SuperMatrix build_steady_rhs_vector_vs(model_t *model, double **rhs);
void set_heuristic_vol_PDN_regular(model_t *model, double *temp, model_vector_t *power);
void set_heuristic_vol_PDN_vs(model_t *model, double *temp, model_vector_t *power);
void compute_pkg_vol(model_t *model, model_vector_t *temp, model_vector_t *power);

/* Transient state */
int PDN_trans_matrix_dim(model_t *model);
void PDN_init_trans_vector(model_t *model, double *g);
void trans_SLU_init(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r);
void trans_SLU_init_nogridL(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r);
void trans_matrix_build_3D(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int **perm_c, int **perm_r);
void compute_PDN_SLU(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r);
void Finalize_rhs(model_t *model, model_vector_t *power);

/* Matrix Operation */
int PadIdx_to_GridIdx(model_t *model, int idx, int layer);
int GridIdx_to_PadIdx(model_t *model, int idx, int layer);
char *GridIdx_to_UnitName(model_t *model, int idx, int layer);

/* Domain */
int compare_domain(model_t *model, int idx, int offset);
int does_branch_exist(model_t *model, int idx);
void set_flp_domain(PDN_flp_t *flp, PDN_config_t *config);
int get_core_id(char *blk_name);

/* SC converters */
void switch_SCconverters(model_t *model, int counter, int intvl);

#endif
