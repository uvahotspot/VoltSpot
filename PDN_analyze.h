#ifndef __PDN_ANALYZE_H_
#define __PDN_ANALYZE_H_

#include "PDN_sim.h"

// sensor location mark
#define SENSOR          1

/* PDN model's per node stats */
typedef struct node_stats_t_st
{
  int **counter_2D;
  double **max_2D;
  double **prev_2D;
  double **integral_2D;
  int *counter_1D;
  double *max_1D;
  double *integral_1D;

}node_stats_t;

typedef struct status_st 
{
  int trans_counter;
  int draw_counter;

  double maxdroop;
  double *layer_maxdroop;
  double *maxIR;
  double PDN_power;

  double *TSVcur;
  double *maxTSVcur;

  double sumcur;
  double curIc, avgIc, dIc;
  double curVc, avgVc;

  double vio_area_ratio;
  double pkgDrop;

  double prevIc;
  double prevVc;

	/* pad current	*/
	double **vdd_cur;
	double **gnd_cur;
  double maxcur_vdd, maxcur_gnd;
  double max_pad_dense;
  int    curIth; //# of pads that current exceed threshold
  int    curIth_vdd, curIth_gnd;

	/* sensor location */
  int **sensor_loc;
  double **sensor_noise;

	/* For averaging results per cycle */
  double *cycle_avg;

	/* voltage stats holder */
  node_stats_t gridstats;
  node_stats_t blkstats;
  node_stats_t corestats;
	/* grid voltage gradient record */
  node_stats_t vgradient;

}status_t;

/* Init */
status_t *alloc_status(model_t *model);
void populate_status(model_t *model, status_t *status);
void free_status(model_t *model, status_t *status);

/* Analyze data in PDN */
void PDN_steady_analyze(model_t *model, status_t *status);
void PDN_trans_analyze(model_t *model, status_t *status);
void print_step_singlenode(model_t *model, FILE *fp, int pi, int pj);
void print_status_singlenode(model_t *model, status_t *status, FILE *fp, int pi, int pj);
void print_trans_analyze(model_t *model, status_t *status, FILE *fp);
void print_trans_header(model_t *model, status_t *status, FILE *fp);
void print_steady_analyze(model_t *model, status_t *status);
void get_pad_current(model_t *model, status_t *status);
void get_TSV_current(model_t *model, status_t *status);
double get_steady_PDN_power(model_t *model);
void get_maxIR(model_t *model, double *maxIR);
void onchipV_analyze_2D(model_t *model, status_t *status);
void onchipV_analyze_3D(model_t *model, status_t *status);
void perstep_droop_3D(model_t *model, FILE *fp);
void trans_vgradient_analyze(model_t *model, status_t *status);

/* voltage sensing */
void parse_sensor_loc(model_t *model, status_t *status, char *file);

/* dump data */
void dump_files(model_t *model, status_t *status, int mode);
void dump_PDN_config(model_t *model, FILE *fp);
void dump_grid_drop(model_t *model, char *file);
void dump_violation(model_t *model, status_t *status, char *file);
void dump_trans_vgradient(model_t *model,  status_t *status, char *file);
void dump_tsv_cur(model_t *model, status_t *status, char *file);

/* Animation */
void draw_single_gif(model_t *model, status_t *status, int draw_counter, int mode);
void create_animation(int draw_counter, int animation);

#endif
