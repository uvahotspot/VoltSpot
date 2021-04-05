#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include <math.h>
#include <assert.h>

/* Lib for SuperLU */
#include "slu_ddefs.h"

#include "PDN_sim.h"
#include "flp.h"
#include "util.h"
#include "pad.h"

/* default PDN configuration parameters	*/
PDN_config_t default_PDN_config(void)
{
  PDN_config_t config;

  config.run_PDN = 0; 

  config.PDN_padconfig     = 1;
  config.padloc_format     = 1;
  config.reserve_io        = 0;
  config.PDN_grid_intv     = 5;
  config.PDN_pkgLC         = 1;
  config.PDN_gridL         = 1;
  config.vgradient_analyse = 0;
  config.PDN_padpitch      = 285e-6;
  config.PDN_padD          = 130e-6;
  config.PDN_padR          = 10e-3;
  config.vdd = 1;
  config.gnd = 0;

  config.proc_clock_freq = 3.7e9;
  config.PDN_step_percycle = 5;
  config.ptrace_sampling_intvl = 1;
  config.PDN_ptrace_warmup  = 0;
  config.PDN_decap_dense  = 12;
  config.PDN_decap_ratio  = 0.3;
  config.PDN_decap_unifm  = 0;
  config.PDN_padL   = 72e-12;
  config.PDN_pkg_sL = 120e-12;
  config.PDN_pkg_sR = 1e-3;
  config.PDN_pkg_C  = 26e-6;
  config.PDN_pkg_pR = 0.54e-3;
  config.PDN_pkg_pL = 5.61e-12;
  config.PDN_pkg_scale = 1;

  config.v_stacking      = 0;
  config.TSV_config      = 1;
  config.TSV_R           = 46.85e-3;
  config.TSV_L           = 34.262e-12;
  strcpy(config.layer_file_3D, NULLFILE);
  strcpy(config.IVR_loc_file, NULLFILE);
  config.SC_freq = 50e6;
  config.SC_totcap = 8e-9;
  config.SC_Rontop = 1.052;
  config.SC_Ronbtm = 1.17;

  config.PDN_multi_dom   = 0; 
  config.animation       = 0; 
  config.frame_intv      = 5; 
  config.legend_lwr      = 0.05; 
  config.legend_upr      = 0.05; 
  config.legend_curupr   = 0.02; 
  config.PDN_sin_pattern = 0; 
  config.PDN_sin_totstep = 10000; 
  config.PDN_sin_freq    = 100e6; 

  strcpy(config.mlayer_spec_file, NULLFILE);
  strcpy(config.padloc_file_in, NULLFILE);
  strcpy(config.padloc_file_out, NULLFILE);
  strcpy(config.vio_file, NULLFILE);
  strcpy(config.node_viotrace_file, NULLFILE);
  strcpy(config.trans_vgradient_file, NULLFILE);
  strcpy(config.senloc_file, NULLFILE);
  strcpy(config.padcur_file, NULLFILE);
  strcpy(config.tsvcur_file, NULLFILE);
  strcpy(config.gridvol_file, NULLFILE);

  config.MC_pads      = 80; 
  config.IO_dense     = 0.8; 
  config.PDN_cur_dense= 8.5e7; 

  config.PDN_noise_th = 5; 

  return config;
}

/* 
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void PDN_config_add_from_strs(PDN_config_t *config, str_pair *table, int size)
{
  int idx;
  if ((idx = get_str_index(table, size, "run_PDN")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->run_PDN) != 1)
      fatal("invalid format for configuration  parameter run_PDN\n");
  if ((idx = get_str_index(table, size, "PDN_padconfig")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_padconfig) != 1)
      fatal("invalid format for configuration  parameter PDN_padconfig\n");
  if ((idx = get_str_index(table, size, "padloc_format")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->padloc_format) != 1)
      fatal("invalid format for configuration  parameter padloc_format\n");
  if ((idx = get_str_index(table, size, "reserve_io")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->reserve_io) != 1)
      fatal("invalid format for configuration  parameter reserve_io\n");
  if ((idx = get_str_index(table, size, "PDN_grid_intv")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_grid_intv) != 1)
      fatal("invalid format for configuration  parameter PDN_grid_intv\n");
  if ((idx = get_str_index(table, size, "PDN_pkgLC")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_pkgLC) != 1)
      fatal("invalid format for configuration  parameter PDN_pkgLC\n");
  if ((idx = get_str_index(table, size, "PDN_gridL")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_gridL) != 1)
      fatal("invalid format for configuration  parameter PDN_gridL\n");
  if ((idx = get_str_index(table, size, "vgradient_analyse")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->vgradient_analyse) != 1)
      fatal("invalid format for configuration  parameter vgradient_analyse\n");
  if ((idx = get_str_index(table, size, "PDN_padpitch")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_padpitch) != 1)
      fatal("invalid format for configuration  parameter PDN_padpitch\n");
  if ((idx = get_str_index(table, size, "PDN_padD")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_padD) != 1)
      fatal("invalid format for configuration  parameter PDN_padD\n");
  if ((idx = get_str_index(table, size, "PDN_padR")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_padR) != 1)
      fatal("invalid format for configuration  parameter PDN_padR\n");
  if ((idx = get_str_index(table, size, "vdd")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->vdd) != 1)
      fatal("invalid format for configuration  parameter vdd\n");
  if ((idx = get_str_index(table, size, "gnd")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->gnd) != 1)
      fatal("invalid format for configuration  parameter gnd\n");
  if ((idx = get_str_index(table, size, "padloc_file_in")) >= 0)
    if(sscanf(table[idx].value, "%s", config->padloc_file_in) != 1)
      fatal("invalid format for configuration  parameter padloc_file_in\n");
  if ((idx = get_str_index(table, size, "padloc_file_out")) >= 0)
    if(sscanf(table[idx].value, "%s", config->padloc_file_out) != 1)
      fatal("invalid format for configuration  parameter padloc_file_out\n");
  if ((idx = get_str_index(table, size, "vio_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->vio_file) != 1)
      fatal("invalid format for configuration  parameter vio_file\n");
  if ((idx = get_str_index(table, size, "node_viotrace_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->node_viotrace_file) != 1)
      fatal("invalid format for configuration  parameter node_viotrace_file\n");
  if ((idx = get_str_index(table, size, "trans_vgradient_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->trans_vgradient_file) != 1)
      fatal("invalid format for configuration  parameter trans_vgradient_file\n");
  if ((idx = get_str_index(table, size, "senloc_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->senloc_file) != 1)
      fatal("invalid format for configuration  parameter senloc_file\n");
  if ((idx = get_str_index(table, size, "mlayer_spec_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->mlayer_spec_file) != 1)
      fatal("invalid format for configuration  parameter mlayer_spec_file\n");
  if ((idx = get_str_index(table, size, "layer_file_3D")) >= 0)
    if(sscanf(table[idx].value, "%s", config->layer_file_3D) != 1)
      fatal("invalid format for configuration  parameter layer_file_3D\n");
  if ((idx = get_str_index(table, size, "v_stacking")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->v_stacking) != 1)
      fatal("invalid format for configuration  parameter v_stacking\n");
  if ((idx = get_str_index(table, size, "TSV_config")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->TSV_config) != 1)
      fatal("invalid format for configuration  parameter TSV_config\n");
  if ((idx = get_str_index(table, size, "TSV_R")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->TSV_R) != 1)
      fatal("invalid format for configuration  parameter TSV_R\n");
  if ((idx = get_str_index(table, size, "TSV_L")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->TSV_L) != 1)
      fatal("invalid format for configuration  parameter TSV_L\n");
  if ((idx = get_str_index(table, size, "IVR_loc_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->IVR_loc_file) != 1)
      fatal("invalid format for configuration  parameter IVR_loc_file\n");
  if ((idx = get_str_index(table, size, "padcur_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->padcur_file) != 1)
      fatal("invalid format for configuration  parameter padcur_file\n");
  if ((idx = get_str_index(table, size, "tsvcur_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->tsvcur_file) != 1)
      fatal("invalid format for configuration  parameter tsvcur_file\n");
  if ((idx = get_str_index(table, size, "gridvol_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->gridvol_file) != 1)
      fatal("invalid format for configuration  parameter gridvol_file\n");
  if ((idx = get_str_index(table, size, "PDN_cur_dense")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_cur_dense) != 1)
      fatal("invalid format for configuration  parameter PDN_cur_dense\n");
  if ((idx = get_str_index(table, size, "PDN_multi_dom")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_multi_dom) != 1)
      fatal("invalid format for configuration  parameter PDN_multi_dom\n");
  if ((idx = get_str_index(table, size, "animation")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->animation) != 1)
      fatal("invalid format for configuration  parameter animation\n");
  if ((idx = get_str_index(table, size, "frame_intv")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->frame_intv) != 1)
      fatal("invalid format for configuration  parameter frame_intv\n");
  if ((idx = get_str_index(table, size, "legend_lwr")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->legend_lwr) != 1)
      fatal("invalid format for configuration  parameter legend_lwr\n");
  if ((idx = get_str_index(table, size, "legend_upr")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->legend_upr) != 1)
      fatal("invalid format for configuration  parameter legend_upr\n");
  if ((idx = get_str_index(table, size, "legend_curupr")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->legend_curupr) != 1)
      fatal("invalid format for configuration  parameter legend_curupr\n");
  if ((idx = get_str_index(table, size, "PDN_sin_pattern")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_sin_pattern) != 1)
      fatal("invalid format for configuration  parameter PDN_sin_pattern\n");
  if ((idx = get_str_index(table, size, "PDN_sin_totstep")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_sin_totstep) != 1)
      fatal("invalid format for configuration  parameter PDN_sin_totstep\n");
  if ((idx = get_str_index(table, size, "PDN_sin_freq")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_sin_freq) != 1)
      fatal("invalid format for configuration  parameter PDN_sin_freq\n");
  if ((idx = get_str_index(table, size, "PDN_noise_th")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_noise_th) != 1)
      fatal("invalid format for configuration  parameter PDN_noise_th\n");
  if ((idx = get_str_index(table, size, "MC_pads")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->MC_pads) != 1)
      fatal("invalid format for configuration  parameter MC_pads\n");
  if ((idx = get_str_index(table, size, "IO_dense")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->IO_dense) != 1)
      fatal("invalid format for configuration  parameter IO_dense\n");
  if ((idx = get_str_index(table, size, "PDN_step_percycle")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_step_percycle) != 1)
      fatal("invalid format for configuration  parameter PDN_step_percycle\n");
  if ((idx = get_str_index(table, size, "proc_clock_freq")) >= 0)
    if(sscanf(table[idx].value, "%le", &config->proc_clock_freq) != 1)
      fatal("invalid format for configuration  parameter proc_clock_freq\n");
  if ((idx = get_str_index(table, size, "ptrace_sampling_intvl")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->ptrace_sampling_intvl) != 1)
      fatal("invalid format for configuration  parameter ptrace_sampling_intvl\n");
  if ((idx = get_str_index(table, size, "PDN_ptrace_warmup")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_ptrace_warmup) != 1)
      fatal("invalid format for configuration  parameter PDN_ptrace_warmup\n");
  if ((idx = get_str_index(table, size, "PDN_decap_dense")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_decap_dense) != 1)
      fatal("invalid format for configuration  parameter PDN_decap_dense\n");
  if ((idx = get_str_index(table, size, "PDN_decap_ratio")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_decap_ratio) != 1)
      fatal("invalid format for configuration  parameter PDN_decap_ratio\n");
  if ((idx = get_str_index(table, size, "PDN_decap_unifm")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->PDN_decap_unifm) != 1)
      fatal("invalid format for configuration  parameter PDN_decap_unifm\n");
  if ((idx = get_str_index(table, size, "PDN_padL")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_padL) != 1)
      fatal("invalid format for configuration  parameter PDN_padL\n");
  if ((idx = get_str_index(table, size, "PDN_pkg_sL")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_pkg_sL) != 1)
      fatal("invalid format for configuration  parameter PDN_pkg_sL\n");
  if ((idx = get_str_index(table, size, "PDN_pkg_sR")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_pkg_sR) != 1)
      fatal("invalid format for configuration  parameter PDN_pkg_sR\n");
  if ((idx = get_str_index(table, size, "PDN_pkg_C")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_pkg_C) != 1)
      fatal("invalid format for configuration  parameter PDN_pkg_C\n");
  if ((idx = get_str_index(table, size, "PDN_pkg_pR")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_pkg_pR) != 1)
      fatal("invalid format for configuration  parameter PDN_pkg_pR\n");
  if ((idx = get_str_index(table, size, "PDN_pkg_pL")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_pkg_pL) != 1)
      fatal("invalid format for configuration  parameter PDN_pkg_pL\n");
  if ((idx = get_str_index(table, size, "PDN_pkg_scale")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->PDN_pkg_scale) != 1)
      fatal("invalid format for configuration  parameter PDN_pkg_scale\n");
  if ((idx = get_str_index(table, size, "SC_freq")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->SC_freq) != 1)
      fatal("invalid format for configuration  parameter SC_freq\n");
  if ((idx = get_str_index(table, size, "SC_totcap")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->SC_totcap) != 1)
      fatal("invalid format for configuration  parameter SC_totcap\n");
  if ((idx = get_str_index(table, size, "SC_Rontop")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->SC_Rontop) != 1)
      fatal("invalid format for configuration  parameter SC_Rontop\n");
  if ((idx = get_str_index(table, size, "SC_Ronbtm")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->SC_Ronbtm) != 1)
      fatal("invalid format for configuration  parameter SC_Ronbtm\n");

  if (config->PDN_grid_intv < 1)
    fatal("PDN grid intv option should not be zero or negative!!\n");
  if (config->run_PDN < 0 || config->run_PDN > 2 )
    fatal("run_PDN should be 0, 1 or 2\n");
  if (config->PDN_padD <= 0)
    fatal("PDN pad diameter should be greater than 0\n");
  if (config->PDN_padR <= 0)
    fatal("PDN pad resistivity should be greater than 0\n");
  if (config->PDN_pkgLC != 0 && config->PDN_pkgLC != 1)
    fatal("PDN_pkgLC shoud be either 0 or 1\n");
  if (config->PDN_gridL != 0 && config->PDN_gridL != 1)
    fatal("PDN_gridL shoud be either 0 or 1\n");
  if (config->vgradient_analyse != 0 && config->vgradient_analyse != 1)
    fatal("vgradient_analyse shoud be either 0 or 1\n");
  if (config->PDN_padconfig < 0 || config->PDN_padconfig >= PADCONFIGS)
    fatal("Invalid PDN_padconfig input!\n");
  if (config->TSV_config < 0 || config->TSV_config >= TSVCONFIGS)
    fatal("Invalid TSV_config input!\n");
  if (config->TSV_R < 0)
    fatal("TSV_R should be larger than 0!\n");
  if (config->TSV_L < 0)
    fatal("TSV_L should be larger than 0!\n");
  if (config->v_stacking < 0 || config->v_stacking > 1)
    fatal("Invalid v_stacking input!\n");
  if (config->padloc_format < 0 || config->padloc_format > 3)
    fatal("Invalid padloc_format input!\n");
  if (config->reserve_io < 0 || config->reserve_io > 3)
    fatal("reserve_io should be 1/0!\n");
  if (config->PDN_cur_dense <= 0 )
    fatal("PDN_cur_dense should be larger than 0\n");
  if (config->PDN_multi_dom < -1)
    fatal("PDN_multi_dom should be greater than -1\n");
  if (config->PDN_decap_unifm < 0)
    fatal("PDN_decap_unifm should be greater than or equal to 0\n");
  if (config->animation < 0 || config->animation > 3)
    fatal("animation should be either 0 or 1\n");
  if (config->frame_intv < 0)
    fatal("frame_intv should be greater than 0\n");
  if (config->PDN_sin_pattern < 0 || config->PDN_sin_pattern > 1)
    fatal("PDN_sin_pattern should be either 0 or 1\n");
  if (config->PDN_sin_totstep < 0)
    fatal("PDN_sin_totstep should be larger than 0 \n");
  if (config->PDN_sin_freq < 0)
    fatal("PDN_sin_freq should be greater than 0\n");
  if (config->vdd < config->gnd)
    fatal("vdd should be larger than gnd\n");
  if (config->PDN_noise_th < 0)
    fatal("PDN_noise_th should be greater than 0!\n");
  if (config->MC_pads < 0)
    fatal("MC_pads should be equal or larger than 0!\n");
  if (config->IO_dense <= 0 || config->IO_dense > 1)
    fatal("IO_dense should be between 0 and 1 (!=0)!\n");
  if (config->ptrace_sampling_intvl <= 0)
    fatal("PDN ptrace interval should be greater than 0\n");
  if (config->PDN_step_percycle<= 0)
    fatal("PDN transient timestep should be greater than 0\n");
  if (config->proc_clock_freq<= 0)
    fatal("PDN clock freq should be greater than 0\n");
  if (config->PDN_ptrace_warmup < 0)
    fatal("PDN warm up time should be larger than 0\n");
  if (config->PDN_decap_dense <= 0)
    fatal("PDN decap density should be greater than 0\n");
  if (config->PDN_decap_ratio <= 0 || config->PDN_decap_ratio > 1)
    fatal("PDN decap ratio should within (0,1] \n");
  if (config->PDN_padL <= 0)
    fatal("PDN pad L should be greater than 0\n");
  if (config->PDN_pkg_sL <= 0)
    fatal("PDN package L should be greater than 0\n");
  if (config->PDN_pkg_sR < 0)
    fatal("PDN package R should be greater than 0\n");
  if (config->PDN_pkg_C <= 0)
    fatal("PDN package C should be greater than 0\n");
  if (config->PDN_pkg_pR <= 0)
    fatal("PDN package parallel R should be greater than 0\n");
  if (config->PDN_pkg_pL <= 0)
    fatal("PDN package parallel L should be greater than 0\n");
  if (config->PDN_pkg_scale <= 0)
    fatal("PDN package scale should be greater than 0\n");
  if (config->SC_freq <= 0)
    fatal("SC frequency should be greater than 0\n");
  if (config->SC_totcap <= 0)
    fatal("SC total cap should be greater than 0\n");
  if (config->SC_Rontop <= 0)
    fatal("SC Ron_top should be greater than 0\n");
  if (config->SC_Ronbtm <= 0)
    fatal("SC Ron_btm should be greater than 0\n");
}

/* constructor also initilize parameters	*/ 
model_t *alloc_model(PDN_config_t *config, PDN_flp_t *flp_default)
{
  int i;
  int pad_grid_col, pad_grid_row;

  model_t *model;

  model = (model_t *) calloc (1, sizeof(model_t));
  if (!model)
    fatal("memory allocation error\n");

  model->config = *config;

  model->width = PDN_get_total_width(flp_default);
  model->height = PDN_get_total_height(flp_default);

  pad_grid_col = floor(model->width / (model->config.PDN_padpitch));
  pad_grid_row = floor(model->height / (model->config.PDN_padpitch));
  if(!pad_grid_row | !pad_grid_col)
    fatal("Pad pithch is lager than chip scale, please fix!\n");

  model->rows = model->config.PDN_grid_intv * (pad_grid_row - 1) + 1;
  model->cols = model->config.PDN_grid_intv * (pad_grid_col - 1) + 1;

  // Sanity check on problem size
  // Force exit if grid size is too large
  if((model->rows < 0) || (model->cols < 0))
    fatal("On-chip grid size calculation overflow! Check your floorplan unit (Should all in meter, not mm or um).\n");
  if((pad_grid_col == 1) && (pad_grid_row == 1))
    fatal("This design only has one C4 pad! Please double-check your floorplan and/or C4 pad pitch\n");
  if((model->rows > MAX_DIM) || (model->cols > MAX_DIM))
    fatal("Problem size exceeding limit!\nCheck your floorplan unit (Should all in meter, not mm or um)\nIf you do want to specify a problem that large, please override this assertion in function alloc_model\n");

  /* set C4 pads */
  model->c4 = alloc_C4_PDN(model, pad_grid_col, pad_grid_row);

  if(strcmp(model->config.layer_file_3D, NULLFILE))
    model->is_3D = TRUE;
  else
    model->is_3D = FALSE;

  /* get layer information	*/
  alloc_layers_PDN(model, flp_default);

  if(model->config.v_stacking)
    model->sc_converter = alloc_IVR(model);

  /* count the total no. of blocks */
  model->total_n_blocks = 0;
  for(i=0; i<model->n_layers; i++){
      model->total_n_blocks += model->layers[i].flp->n_units;
  }

  /* allocate internal state	*/
  /* Two extra nodes for package, each layer has two grid */
  model->last_steady = dvector(2*model->rows*model->cols*model->n_layers + PDN_STEADY_EXTRA);
  model->last_power = dvector(model->rows*model->cols*model->n_layers);

  return model;
}

PDN_metal_t alloc_metal_layers(int n_metal)
{
  PDN_metal_t metal_layers;

  metal_layers.n_metal = n_metal;
  metal_layers.geo     = (metal_geo_t *) calloc (n_metal, sizeof(metal_geo_t));
  metal_layers.gridRL  = (metal_gridRL_t *) calloc (n_metal, sizeof(metal_gridRL_t));

  return metal_layers;
}

void set_current_limit(model_t *model)
{
  double padD, padR;
  double dense;

  padD = model->config.PDN_padD;
  padR = padD / 2;
  dense= model->config.PDN_cur_dense;

  model->c4->Ith = dense * PI * padR * padR;
}

void alloc_layers_PDN(model_t *model, PDN_flp_t *flp_default)
{
  int i;
	FILE *fp = NULL;
	char str[STR_SIZE];
  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;

  if (model->is_3D){
		if (!strcasecmp(model->config.layer_file_3D, "stdin"))
			fp = stdin;
		else
			fp = fopen (model->config.layer_file_3D, "r");
		if (!fp) {
			sprintf(str, "error opening lcf file %s\n", model->config.layer_file_3D);
			fatal(str);
		}
  }

  /* compute the no. of layers	*/
  if (model->is_3D){
      model->n_layers = count_significant_lines(fp);
      if (model->n_layers % LCF_3D_NPARAMS)
        fatal("wrong no. of lines in layer file\n");
      model->n_layers /= LCF_3D_NPARAMS;
  }
  else
    model->n_layers = 1;

  /* allocate initial memory */
  model->layers = (PDN_layer_t *) calloc (model->n_layers, sizeof(PDN_layer_t));
  if (!model->layers)
    fatal("memory allocation error\n");

	/* populate layers */
  if (model->is_3D){
		parse_layer_file_PDN(model, fp);
		warning("3D layer configuration file specified. overriding default floorplan with those in lcf file, but requires default flp and flps in lcf file have same width and height...\n");
  }
  else
    populate_single_layer_PDN(model, flp_default);

  /* allocate cap and set bgmap*/
  for(i=0; i<model->n_layers; i++){
      model->layers[i].cap_c = dmatrix(nr, nc);
      zero_dvector(model->layers[i].cap_c[0], nr*nc);
      set_bgmap_PDN(model, &model->layers[i]);
  }

  if (model->is_3D && fp != stdin)
    fclose(fp);
}

IVR_t *alloc_IVR(model_t *model)
{
  IVR_t *ivr;

  ivr = (IVR_t *) calloc (1, sizeof(IVR_t));
  ivr->loc = imatrix(model->rows, model->cols);
  zero_ivector(ivr->loc[0], model->rows*model->cols);

  return ivr;
}

void populate_IVR(model_t *model)
{
  int i, j;
  int counter = 0;
  int nr = model->rows;
  int nc = model->cols;

  IVR_t *v = model->sc_converter;

  if(strcmp(model->config.IVR_loc_file, NULLFILE))
    parse_IVR_loc(model, model->config.IVR_loc_file);

  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++) {
        if(IVRLOC == v->loc[i][j])
          counter++;
    }
  v->num_IVR = counter;
  v->freq = 2*2*model->config.SC_freq;
}

void populate_single_layer_PDN(model_t *model, PDN_flp_t *flp_default)
{
  int num_layer;
  FILE *fp = NULL;
  char str[STR_SIZE];

  if(strcmp(model->config.mlayer_spec_file, NULLFILE)){
      fp = fopen(model->config.mlayer_spec_file, "r");
      if(!fp) {
          sprintf(str, "error opening mlayer_spec file %s\n", model->config.mlayer_spec_file);
          fatal(str);
      }
      num_layer = count_significant_lines(fp);
      if (num_layer % MLCF_NPARAMS)
        fatal("wrong no. of lines in metal layer file\n");
      num_layer /= MLCF_NPARAMS;
      fclose(fp);
  }
  else
    num_layer = DEFAULT_MLAYERS;

  model->layers[0].no = 0;
  model->layers[0].flp = flp_default;
  model->layers[0].metal_layers = alloc_metal_layers(num_layer);
  model->layers[0].b2gmap = new_b2gmap_PDN(model->rows, model->cols);
  if(strcmp(model->config.mlayer_spec_file, NULLFILE))
    parse_metal_layer_file(model, model->config.mlayer_spec_file, 0);
  else
    populate_default_mlayers(model, 0);
}

void populate_default_mlayers(model_t *model, int layer)
{
  model->layers[0].metal_layers.geo[0].pitch = 30e-6;
  model->layers[0].metal_layers.geo[0].width = 10e-6;
  model->layers[0].metal_layers.geo[0].thick = 3.5e-6;
  model->layers[0].metal_layers.geo[0].rho   = 1.68e-8;
  model->layers[0].metal_layers.geo[0].direc = MLCF_X;

  model->layers[0].metal_layers.geo[1].pitch = 30e-6;
  model->layers[0].metal_layers.geo[1].width = 10e-6;
  model->layers[0].metal_layers.geo[1].thick = 3.5e-6;
  model->layers[0].metal_layers.geo[1].rho   = 1.68e-8;
  model->layers[0].metal_layers.geo[1].direc = MLCF_Y;

  model->layers[0].metal_layers.geo[2].pitch = 810e-9;
  model->layers[0].metal_layers.geo[2].width = 400e-9;
  model->layers[0].metal_layers.geo[2].thick = 720e-9;
  model->layers[0].metal_layers.geo[2].rho   = 1.68e-8;
  model->layers[0].metal_layers.geo[2].direc = MLCF_X;

  model->layers[0].metal_layers.geo[3].pitch = 810e-9;
  model->layers[0].metal_layers.geo[3].width = 400e-9;
  model->layers[0].metal_layers.geo[3].thick = 720e-9;
  model->layers[0].metal_layers.geo[3].rho   = 1.68e-8;
  model->layers[0].metal_layers.geo[3].direc = MLCF_Y;
}

void populate_R_model_PDN(model_t *model)
{
  int i, l;
  int metal_grid_col, metal_grid_row;
  double pitch, width, thick, rho;
  int direc;

  /* short cuts */
  int nl = model->n_layers;
  double cw = model->width;
  double ch = model->height;
  double nr = model->rows;
  double nc = model->cols;
  int vs = model->config.v_stacking;

  /* layer specific resistances	
   * R = rho * lenth / Cross-sectional area
   */
  for(l=0; l<nl; l++)
    for(i=0; i<model->layers[l].metal_layers.n_metal; i++){
        pitch = model->layers[l].metal_layers.geo[i].pitch;
        width = model->layers[l].metal_layers.geo[i].width;
        thick = model->layers[l].metal_layers.geo[i].thick;
        rho   = model->layers[l].metal_layers.geo[i].rho;
        direc = model->layers[l].metal_layers.geo[i].direc;

        metal_grid_col = floor(cw/(2*pitch));
        metal_grid_row = floor(ch/(2*pitch));

        if(MLCF_X == direc){
            model->layers[l].metal_layers.gridRL[i].r = 
              (1/(nc-1))*(nr/metal_grid_row) * (rho*cw/(width*thick));
        }
        else{
            model->layers[l].metal_layers.gridRL[i].r = 
              (1/(nr-1))*(nc/metal_grid_col) * (rho*ch/(width*thick));
        }
    }

  // TSV resistance
  if(model->is_3D){
      for(l=0; l<nl; l++){
          model->layers[l].tsv.r = model->config.TSV_R;
      }
  }

  // C4 resistance
  model->c4->pad_r = model->config.PDN_padR; 

  // IVR resistance
  if(vs){
      model->sc_converter->R_drop = 0.46;
      model->sc_converter->Ron_top = 2*model->config.SC_Rontop;
      model->sc_converter->Ron_bottom = 2*model->config.SC_Ronbtm;
  }

  // scale pkg parameter for mtdmn "block" model
  model->config.PDN_pkg_sR *= model->config.PDN_pkg_scale; 
}

void populate_LC_model_PDN(model_t *model)
{
  int i, j, l, d;
  double w, t, s, p;
  double mgc, mgr;
  int decap_count_hot = 0;
  int decap_count_cold = 0;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int vs = model->config.v_stacking;
  double cw = model->width;
  double ch = model->height;
  double tot_decap, tot_ldcap, decap_percell_hot, decap_percell_cold, ldcap_percell;

  /* package LC's	*/
  // scale pkg parameter for mtdmn "block" model
  model->config.PDN_pkg_sL *= model->config.PDN_pkg_scale; 
  model->config.PDN_pkg_C  /= model->config.PDN_pkg_scale; 
  model->config.PDN_pkg_pR *= model->config.PDN_pkg_scale; 
  model->config.PDN_pkg_pL *= model->config.PDN_pkg_scale; 

  /* Total load cap is estimated from equation W = alpha*f*C*V2 */
  tot_ldcap = 16.4e-8 / model->config.PDN_pkg_scale;
  //tot_ldcap = 5e-8 / model->config.PDN_pkg_scale;
  ldcap_percell = (2*tot_ldcap) / (nr*nc); // all cell has load cap, x2 for two grid

  /* Total decap is calculated based on user input */
  tot_decap = model->config.PDN_decap_dense * \
              model->config.PDN_decap_ratio * \
              model->height * model->width * \
              (1.0e-9 / 1.0e-6);

  for(l=0; l<nl; l++)
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++){
          if(decap_under_blk(model, (i*nc+j), l))
            decap_count_hot++;
          else
            decap_count_cold++;
      }

  if (model->config.PDN_decap_unifm > 0){
          decap_percell_hot = (2*tot_decap*LOCAL_CAP_RATIO ) / (decap_count_hot);
          decap_percell_cold = (2*tot_decap*(1-LOCAL_CAP_RATIO)) / (decap_count_cold);
  }

  /* layer specific capacitances	*/
  for(l=0; l<nl; l++){
      // Place all decap around designated area
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++) {
            if(0 == model->config.PDN_decap_unifm){
                if(model->is_3D)
                  model->layers[l].cap_c[i][j] = (tot_decap+tot_ldcap)/(nr*nc);
                else
                  model->layers[l].cap_c[i][j] = (2*(tot_decap+tot_ldcap))/(nr*nc);
            }
            else{
                if(decap_under_blk(model, (i*nc+j), l)){
                    model->layers[l].cap_c[i][j] = decap_percell_hot + ldcap_percell;
                }
                else{
                    model->layers[l].cap_c[i][j] = decap_percell_cold + ldcap_percell;
                }
            }
        }
  }

  // TSV LC
  if(model->is_3D){
      for(l=0; l<nl; l++){
          model->layers[l].tsv.l = model->config.TSV_L;
      }
  }

  /* IVR C */
  if(vs){
      model->sc_converter->c = model->config.SC_totcap/2;
  }

  /* pad L */
  model->c4->pad_l = model->config.PDN_padL;

  /* on-chip L */
  for(l=0; l<nl; l++)
    for(i=0; i<model->layers[l].metal_layers.n_metal; i++) {
        p  = model->layers[l].metal_layers.geo[i].pitch;
        w  = model->layers[l].metal_layers.geo[i].width;
        t  = model->layers[l].metal_layers.geo[i].thick;
        s  = model->layers[l].metal_layers.geo[i].pitch - 
          model->layers[l].metal_layers.geo[i].width;
        d  = model->layers[l].metal_layers.geo[i].direc;

        mgc = floor(cw/(2*p));
        mgr = floor(ch/(2*p));

        if(MLCF_X == d){
            model->layers[l].metal_layers.gridRL[i].l = 
              (cw/(nc-1))*(nr/mgr)*(1/PI)*VAPER*(1.5+log(2/PI)+log((w+s)/(w+t)));
        }
        else{
            model->layers[l].metal_layers.gridRL[i].l = 
              (ch/(nr-1))*(nc/mgc)*(1/PI)*VAPER*(1.5+log(2/PI)+log((w+s)/(w+t)));
        }
    }
}

int decap_under_blk(model_t *model, int idx, int layer)
{
  int i,j;
  int nc = model->cols;
  char *name;
  blist_t *ptr;

  if(0 == model->config.PDN_decap_unifm)
    return 1;

  j = idx % nc;
  i = (int) floor(idx/nc);

  ptr = model->layers[layer].b2gmap[i][j];
  name = model->layers[layer].flp->units[ptr->idx].name;

  // A quick and dirty way to do identify places to put decap
  if((name[0] == 'A' &&
      name[1] == 'L' &&
      name[2] == 'U')){
      return 1;
  }
  else if((name[0] == 'R' &&
           name[1] == 'O' &&
           name[2] == 'B')){
      return 1;
  }
  else if((name[0] == 'I' &&
           name[1] == 'n' &&
           name[2] == 't' &&
           name[3] == 'I' &&
           name[4] == 'W') &&
          model->config.PDN_decap_unifm > 1){
      return 1;
  }
  else if((name[0] == 'I' &&
           name[1] == 'n' &&
           name[2] == 't' &&
           name[3] == 'R' &&
           name[4] == 'F') &&
          model->config.PDN_decap_unifm > 1){
      return 1;
  }
  else if((name[0] == 'F' &&
           name[1] == 'P' &&
           name[2] == 'U') &&
          model->config.PDN_decap_unifm > 2){
      return 1;
  }
  else if((name[0] == 'C' &&
           name[1] == 'p' &&
           name[2] == 'l' &&
           name[3] == 'A' &&
           name[4] == 'L' &&
           name[5] == 'U') &&
          model->config.PDN_decap_unifm > 2){
      return 1;
  }
  else
    return 0;
}

double get_onchip_cap(model_t *model, int idx, int layer)
{
  int i,j;
  int nc = model->cols;

  j = idx % nc;
  i = (int) floor(idx/nc);

  return model->layers[layer].cap_c[i][j];
}

void steady_state_PDN(model_t *model, double *power)
{
  if(model->config.v_stacking)
    steady_state_PDN_vs(model, power);
  else
    steady_state_PDN_regular(model, power);
}

void steady_state_PDN_regular(model_t *model, double *power)
{
  model_vector_t *p;
  SuperMatrix A, L, U, B;
  double   *rhs;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      info;
  superlu_options_t options;
  SuperLUStat_t stat;

  int          i;
  DNformat     *Astore;
  double       *dp;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;

  p = new_model_vector(model);

  /* map the block power numbers to the grid	*/
  PDN_xlate_vector_b2g(model, power, p);

  /* Calculate pkg vdd/gnd  */
  set_heuristic_vol_PDN_regular(model, model->last_steady, p);

  A = build_steady_grid_matrix(model);
  B = build_steady_rhs_vector(model, p, &rhs);

  if ( !(perm_r = intMalloc(2*nl*nr*nc)) ) fatal("Malloc fails for perm_r[].\n");
  if ( !(perm_c = intMalloc(2*nl*nr*nc)) ) fatal("Malloc fails for perm_c[].\n");

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = 0.01;
  options.SymmetricMode = YES;
  options.Equil = YES;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve the linear system. */
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

  Astore = (DNformat *) B.Store;
  dp = (double *) Astore->nzval;
  //copy results back to last_steady
  for(i=0; i<2*nl*nr*nc; ++i){
      model->last_steady[i] = dp[i];
  }

  free_model_vector(p);
  SUPERLU_FREE (rhs);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);
}

void steady_state_PDN_vs(model_t *model, double *power)
{
  model_vector_t *p;
  SuperMatrix A, L, U, B;
  double   *rhs;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      info;
  superlu_options_t options;
  SuperLUStat_t stat;

  int          i;
  DNformat     *Astore;
  double       *dp;

  double       max_error, cur_error;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;

  p = new_model_vector(model);

  /* map the block power numbers to the grid	*/
  PDN_xlate_vector_b2g(model, power, p);

  /* Calculate pkg vdd/gnd  */
  set_heuristic_vol_PDN_vs(model, model->last_steady, p);

  A = build_steady_grid_matrix_vs(model, p);
  B = build_steady_rhs_vector_vs(model, &rhs);

  if ( !(perm_r = intMalloc(2*nl*nr*nc)) ) fatal("Malloc fails for perm_r[].\n");
  if ( !(perm_c = intMalloc(2*nl*nr*nc)) ) fatal("Malloc fails for perm_c[].\n");

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = 0.01;
  options.SymmetricMode = YES;
  options.Equil = YES;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve the linear system. */
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

  Astore = (DNformat *) B.Store;
  dp = (double *) Astore->nzval;
  //copy results back to last_steady
  for(i=0; i<2*nl*nr*nc; ++i){
      model->last_steady[i] = dp[i];
  }

  StatFree(&stat);
  SUPERLU_FREE (rhs);
  free_model_vector(p);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
}

/* translate power/voltage between block and grid vectors	*/
void PDN_xlate_vector_b2g(model_t *model, double *b, model_vector_t *g)
{
  int i, j, l, base = 0;
  double area;

  /* area of a single grid cell	*/
  area = (model->width * model->height) / (model->cols * model->rows);

  for(l=0; l<model->n_layers; l++){
      for(i=0; i<model->rows; i++)
        for(j=0; j<model->cols; j++) {
            /* for each grid cell, the power density / voltage are 
             * the average of the power densities / voltage of the 
             * blocks in it weighted by their occupancies
             */
            /* convert power density to power	*/ 
            g->cuboid[l][i][j] = blist_avg_PDN(model->layers[l].b2gmap[i][j], 
                                               model->layers[l].flp, &b[base]) * area;
        }
      base+= model->layers[l].flp->n_units;
  }
}

/* constructor	*/
model_vector_t *new_model_vector(model_t *model)
{
  model_vector_t *v;

  v = (model_vector_t *) calloc (1, sizeof(model_vector_t));
  if (!v)
    fatal("memory allocation error\n");

  v->cuboid = dcuboid_PDN(model->rows, model->cols, model->n_layers);
  return v;
}

/* create a linked list node and append it at the end	*/
void blist_append_PDN(blist_t *head, int idx, double occupancy)
{
  blist_t *tail = NULL;

  if(!head)
    fatal("blist_append_PDN called with empty list\n");

  /* traverse till the end	*/
  for(; head; head = head->next)
    tail = head;

  /* append */
  tail->next =  new_blist_PDN(idx, occupancy);
}

/* destructor	*/
void delete_b2gmap_PDN(blist_t ***b2gmap, int rows, int cols)
{
  int i, j;
  blist_t *ptr, *temp;

  /* free the linked list	*/
  for(i=0; i < rows; i++)
    for(j=0; j < cols; j++) {
        ptr = b2gmap[i][j];
        while(ptr) {
            temp = ptr->next;
            free(ptr);
            ptr = temp;
        }
    }

  /* free the array space	*/
  free(b2gmap[0]);
  free(b2gmap);
}

/* setup the block and grid mapping data structures	*/
void set_bgmap_PDN(model_t *model, PDN_layer_t *layer)
{
  /* i1, i2, j1 and j2 are indices of the boundary grid cells	*/
  int i, j, u, i1, i2, j1, j2;

  /* shortcuts for cell width(cw) and cell height(ch)	*/
  double cw = model->width / model->cols;
  double ch = model->height / model->rows;

  /* initialize	*/
  reset_b2gmap_PDN(model, layer);

  /* for each functional unit	*/
  for(u=0; u<layer->flp->n_units; u++) {
      /* shortcuts for unit boundaries	*/
      double lu = layer->flp->units[u].leftx;
      double ru = lu + layer->flp->units[u].width;
      double bu = layer->flp->units[u].bottomy;
      double tu = bu + layer->flp->units[u].height;

      /* top index (lesser row) = rows - ceil (topy / cell height)	*/
      i1 = model->rows - tolerant_ceil(tu/ch);
      /* bottom index (greater row) = rows - floor (bottomy / cell height)	*/
      i2 = model->rows - tolerant_floor(bu/ch);
      /* left index = floor (leftx / cell width)	*/
      j1 = tolerant_floor(lu/cw);
      /* right index = ceil (rightx / cell width)	*/
      j2 = tolerant_ceil(ru/cw);
      /* sanity check	*/
      if((i1 < 0) || (j1 < 0))
        fatal("negative grid cell start index!\n");
      if((i2 > model->rows) || (j2 > model->cols))
        fatal("grid cell end index out of bounds!\n");
      if((i1 >= i2) || (j1 >= j2))
        fatal("invalid floorplan spec or grid resolution\n");

      /* setup b2gmap	*/
      /* for each grid cell in this unit	*/
      for(i=i1; i<i2; i++)
        for(j=j1; j<j2; j++)
          /* grid cells fully overlapped by this unit	*/
          if ((i > i1) && (i < i2-1) && (j > j1) && (j < j2-1)) {
              /* first unit in the list	*/
              if (!layer->b2gmap[i][j])
                layer->b2gmap[i][j] = new_blist_PDN(u, 1.0);
              else {
                  /* this should not occur since the grid cell is 
                   * fully covered and hence, no other unit should 
                   * be sharing it
                   */
                  blist_append_PDN(layer->b2gmap[i][j], u, 0.0); // We ignore the second block
                  warning("overlap of functional blocks? ignoring the second one\n");
              }
              /* boundary grid cells partially overlapped by this unit	*/
          } else {
              /* shortcuts for cell boundaries	*/
              double lc = j * cw, rc = (j+1) * cw;
              double tc = model->height - i * ch;
              double bc = model->height - (i+1) * ch;

              /* shortcuts for overlap width and height	*/
              double oh = (MIN(tu, tc) - MAX(bu, bc));
              double ow = (MIN(ru, rc) - MAX(lu, lc));
              double occupancy;

              /* overlap tolerance	*/
              if (eq(oh/ch, 0))
                oh = 0;
              else if (eq(oh/ch, 1))
                oh = ch;

              if (eq(ow/cw, 0))
                ow = 0;
              else if (eq(ow/cw, 1))
                ow = cw;

              occupancy = (oh * ow) / (ch * cw);
              if (oh < 0 || ow < 0)
                fatal("negative overlap!\n");

              /* first unit in the list	*/
              if (!layer->b2gmap[i][j])
                layer->b2gmap[i][j] = new_blist_PDN(u, occupancy);
              else
                /* append at the end	*/
                blist_append_PDN(layer->b2gmap[i][j], u, occupancy);
          }
  }

  /* 
   * sanity check	
   test_b2gmap(model, layer);
   */
}

/* re-initialize */
void reset_b2gmap_PDN(model_t *model, PDN_layer_t *layer)
{
  int i, j;
  blist_t *ptr, *temp;

  /* free the linked list	*/
  for(i=0; i<model->rows; i++)
    for(j=0; j<model->cols; j++) {
        ptr = layer->b2gmap[i][j];
        while(ptr) {
            temp = ptr->next;
            free(ptr);
            ptr = temp;
        }
        layer->b2gmap[i][j] = NULL;
    }
}

/* a simple way to set initial voltage*/
void set_heuristic_vol_PDN_regular(model_t *model, double *temp, model_vector_t *power)
{
  int i, j, l;
  int idx;
  /* current sum of pads */
  double pad_sum = 0;

  double r_pkg_s = model->config.PDN_pkg_sR;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;

  for(l=0; l<nl; l++)
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++) {
          // VDD layer idx
          idx = l*nr*nc + i*nc + j ;
          temp[idx] = vdd;
          // GND layer idx
          idx = nl*nr*nc + l*nr*nc + i*nc + j ;
          temp[idx] = gnd;

          pad_sum += power->cuboid[l][i][j];
      }
  pad_sum /= (vdd-gnd);

  temp[2*nl*nc*nr+PKG_VDD] = vdd - pad_sum * r_pkg_s; 
  temp[2*nl*nc*nr+PKG_GND] = gnd + pad_sum * r_pkg_s;
}

void set_heuristic_vol_PDN_vs(model_t *model, double *temp, model_vector_t *power)
{
  int i, j, l;
  int idx;

  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;

  for(l=0; l<nl; l++)
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++) {
          // GND layer idx
          idx = 2*l*nr*nc + i*nc + j ;
          temp[idx] = l*vdd;
          // Vdd layer idx
          idx = (2*l+1)*nr*nc + i*nc + j ;
          temp[idx] = (l+1)*vdd;
      }
}

blist_t ***new_b2gmap_PDN(int rows, int cols)
{
  int i;
  blist_t ***b2gmap;

  b2gmap = (blist_t ***) calloc (rows, sizeof(blist_t **));
  b2gmap[0] = (blist_t **) calloc (rows * cols, sizeof(blist_t *));
  if (!b2gmap || !b2gmap[0])
    fatal("memory allocation error\n");

  for(i=1; i<rows; i++) 
    b2gmap[i] = b2gmap[0] + cols * i;

  return b2gmap;	
}

/* constructors	*/
blist_t *new_blist_PDN(int idx, double occupancy)
{
  blist_t *ptr = (blist_t *) calloc (1, sizeof(blist_t));
  if (!ptr)
    fatal("memory allocation error\n");
  ptr->idx = idx;
  ptr->occupancy = occupancy;
  ptr->next = NULL;
  return ptr;
}

/* compute the power/voltage average weighted by occupancies	*/
double blist_avg_PDN(blist_t *ptr, PDN_flp_t *flp, double *v)
{
  double  val = 0.0;

  for(; ptr; ptr = ptr->next) {
      val += ptr->occupancy * v[ptr->idx] / (flp->units[ptr->idx].width * 
                                             flp->units[ptr->idx].height);
  }		

  return val;		   
}

double ***dcuboid_PDN(int nr, int nc, int nl)
{
  int i, j;
  double ***m;

  /* 1-d array of pointers to the rows of the 2-d array below	*/
  m = (double ***) calloc (nl, sizeof(double **));
  assert(m != NULL);
  /* 2-d array of pointers denoting (layer, row)	*/
  m[0] = (double **) calloc (nl * nr, sizeof(double *));
  assert(m[0] != NULL);
  /* the actual 3-d data array	*/
  m[0][0] = (double *) calloc (nl * nr * nc, sizeof(double));
  assert(m[0][0] != NULL);

  /* remaining pointers of the 1-d pointer array	*/
  for(i=1; i<nl; i++)
    m[i] = m[0] + nr*i;

  /* remaining pointers of the 2-d pointer array	*/
  for(i=0; i<nl; i++)
    for(j=0; j<nr; j++)
      /* to reach the jth row in the ith layer, 
       * one has to cross i layers i.e., i*(nr*nc)
       * values first and then j rows i.e., j*nc 
       * values next
       */
      m[i][j] = m[0][0] + (nr*nc)*i + nc*j;

  return m;
}

PDN_C4_t *alloc_C4_PDN(model_t *model, int pad_grid_col, int pad_grid_row)
{
  PDN_C4_t *c4;

  c4 = (PDN_C4_t *) calloc (1, sizeof(PDN_C4_t));
  if (!c4)
    fatal("memory allocation error\n");

  c4->gnd_num = 0;
  c4->vdd_num = 0;

  c4->vdd_loc = imatrix(model->rows, model->cols);
  c4->gnd_loc = imatrix(model->rows, model->cols);
  zero_ivector(c4->vdd_loc[0], model->rows*model->cols);
  zero_ivector(c4->gnd_loc[0], model->rows*model->cols);

  c4->pad_grid_row = pad_grid_row;
  c4->pad_grid_col = pad_grid_col;

  c4->vdd_padidx = ivector(pad_grid_row * pad_grid_col);
  c4->gnd_padidx = ivector(pad_grid_row * pad_grid_col);

  return c4;
}

void populate_C4_PDN(model_t *model)
{
  int i, j, nr, nc;
  int pl;
  int rpg, cpg;//pad grid
  int r_cordt, c_cordt;//cordt for pad
  int itv_row, itv_col;
  int vdd_idx, gnd_idx;

  /* shortcuts */
  nr = model->rows;
  nc = model->cols;
  pl = model->config.PDN_padconfig;
  rpg = model->c4->pad_grid_row;
  cpg = model->c4->pad_grid_col;
  itv_row = (int)(nr - 1)/(rpg - 1);
  itv_col = (int)(nc - 1)/(cpg - 1);

  /* initalize */
  for(i=0; i<nr; i++){
      for(j=0; j<nc; j++){
          model->c4->gnd_loc[i][j] = UNDEF;
          model->c4->vdd_loc[i][j] = UNDEF;
      }
  }

  /* pad config */
  if (pl == CONFIG_CUS){
      parse_pad_loc(model, model->config.padloc_file_in);
  }
  else if (pl == CONFIG_ALL){
      // All seats are filled with pads
      for(i=0; i<rpg; i++) {
        for(j=0; j<cpg; j++){
            r_cordt = i * itv_row;
            c_cordt = j * itv_col;
            if ( (i%2) == (j%2) ){
                model->c4->vdd_loc[r_cordt][c_cordt] |= PGPAD;
            }
            else {
                model->c4->gnd_loc[r_cordt][c_cordt] |= PGPAD;
            }
        }
      }
  }
  else{
      fatal("Unexcepted pad configuration setup\n");
  }

  if(0 == model->config.reserve_io){
      mark_all(model, IORAN);
  }
  else if(1 == model->config.reserve_io){
      mark_pads(model, IOPAD);
      remove_pg_for_io(model);
  }
  else if(2 == model->config.reserve_io){
      mark_pads(model, IORAN);
  }

  set_current_limit(model);


  // Assign pad ID, map pad ID to grid index
  for(i=0; i<model->c4->pad_grid_row*model->c4->pad_grid_col; i++){
      model->c4->vdd_padidx[i] = -1;
      model->c4->gnd_padidx[i] = -1;
  }

  vdd_idx = 0;
  gnd_idx = 0;
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++){
        if (PGPAD & model->c4->vdd_loc[i][j]){
            model->c4->vdd_padidx[vdd_idx] = i*nc + j;
            vdd_idx++;
        }
        if (PGPAD & model->c4->gnd_loc[i][j]){
            model->c4->gnd_padidx[gnd_idx] = i*nc + j;
            gnd_idx++;
        }
    }

  // power and ground pad number
  model->c4->vdd_num = vdd_idx;
  model->c4->gnd_num = gnd_idx;
}

void populate_TSV_PDN(model_t *model)
{
  int i, j, l;

  /* shortcuts */
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int tl = model->config.TSV_config;

  /* initalize */
  for(l=0; l<nl; l++)
    for(i=0; i<nr; i++){
        for(j=0; j<nc; j++){
            model->layers[l].tsv.loc[i][j] = UNDEF;
        }
    }

  /* TSV config */
  if (tl == CONFIG_CUS){
      for(l=0; l<nl-1; l++)
        parse_tsv_loc(model, l);
  }
  else if (tl == CONFIG_ALL){
      for(l=0; l<nl-1; l++)
        for(i=0; i<nr; i++)
          for(j=0; j<nc; j++){
              if( (i%2) == (j%2) ){
                  model->layers[l].tsv.loc[i][j] = VDDTSV;
              }
              else{
                  model->layers[l].tsv.loc[i][j] = GNDTSV;
              }
          }
  }
  else{
      fatal("Unexcepted pad configuration setup\n");
  }

  //sanity check rules about TSV
  //we force all layers have same TSV distribution for now
  for(l=1; l<nl-1; l++)
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++){
          if(model->layers[0].tsv.loc[i][j] !=
             model->layers[l].tsv.loc[i][j])
            fatal("Different TSV distribution between layers detected!\n");
      }

  //count TSV num
  for(l=0; l<nl-1; l++){
      model->layers[l].tsv.num_vdd = 0;
      model->layers[l].tsv.num_gnd = 0;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            if(model->layers[l].tsv.loc[i][j] == VDDTSV){
                model->layers[l].tsv.num_vdd += 1;
            }
            else if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                model->layers[l].tsv.num_gnd += 1;
            }
        }
  }
}

void delete_model(model_t *model)
{
  free_layers_PDN(model);
  free_dvector(model->last_steady);
  free_dvector(model->last_trans);
  free_dvector(model->last_power);
  free_PDN_C4(model->c4);
  if(model->config.v_stacking)
    free_IVR(model->sc_converter);
  free(model);
}

void free_layers_PDN(model_t *model)
{
  int i;
  for(i=0; i<model->n_layers; i++){
      free_metal_layers(model->layers[i].metal_layers);
      if(model->is_3D){
          free_TSV(model->layers[i].tsv);
          PDN_free_flp(model->layers[i].flp);
      }
      free_dmatrix(model->layers[i].cap_c);
      delete_b2gmap_PDN(model->layers[i].b2gmap, model->rows, model->cols);
  }
  free(model->layers);
}

void free_TSV(PDN_TSV_t t)
{
  free_imatrix(t.loc);
}

void free_IVR(IVR_t *v)
{
  free_imatrix(v->loc);
  free(v);
}

void free_metal_layers(PDN_metal_t l)
{
  free(l.geo);
  free(l.gridRL);
}

/* free the C4 pads */
void free_PDN_C4(PDN_C4_t *p)
{
  free_imatrix(p->gnd_loc);
  free_imatrix(p->vdd_loc);
  free_ivector(p->vdd_padidx);
  free_ivector(p->gnd_padidx);
  free(p);
}

/* destructor	*/
void free_model_vector(model_vector_t *v)
{
  free_dcuboid(v->cuboid);
  free(v);
}

void compute_PDN_SLU(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r)
{
  SuperMatrix B;
  SuperLUStat_t stat;
  int info;

  model_vector_t *p;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int m = model->trans_matrix_dim;

  // L * Yn, stored in last_trans
  // with package, L = (E + J * delta_t / 2)
  SparseMatrix_mul_Vector(A, model->last_trans);

  p = new_model_vector(model);

  /* map the block power/temp numbers to the grid	*/
  PDN_xlate_vector_b2g(model, power, p);

  Finalize_rhs(model, p);

  dCreate_Dense_Matrix(&B, m, 1, model->last_trans, m, SLU_DN, SLU_D, SLU_GE);

  StatInit(&stat);
  dgstrs(NOTRANS, L, U, perm_c, perm_r, &B, &stat, &info);

  // save previous power trace for next computation
  copy_dvector(model->last_power, p->cuboid[0][0], nr*nc*nl);

  //We free B after set B's nzval ptr to NULL
  //Because last_trans is in B
  DNformat    *Xstore;
  Xstore = (DNformat *) B.Store;
  Xstore->nzval = NULL;
  Destroy_SuperMatrix_Store(&B);

  free_model_vector(p);
  StatFree(&stat);
}

int PDN_trans_matrix_dim(model_t *model)
{
  int m, l;
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int np = model->c4->vdd_num + model->c4->gnd_num;
  int nml = model->layers[0].metal_layers.n_metal/2;
  int nbr = 2*nr*nc - nr - nc;
  int nvtsv = model->layers[0].tsv.num_vdd;
  int ngtsv = model->layers[0].tsv.num_gnd;
  int nIVR;

  if (!model->is_3D){
      m = 2*nr*nc + np;
      if(model->config.PDN_pkgLC)
        m += 2;
      if(model->config.PDN_gridL)
        m += 2*nbr*nml;
  }
  else{
      m = 2*nr*nc*nl + np;

      if(model->config.v_stacking){
          nIVR = model->sc_converter->num_IVR;
          m += ngtsv*(nl-1);
          m += nIVR*(nl-1)*4;
      }
      else{
          m += (nvtsv+ngtsv)*(nl-1);
      }

      if(model->config.PDN_pkgLC)
        m += 2;
  }
  return m;
}

int PadIdx_to_GridIdx(model_t *model, int idx, int layer)
{
  if(LAYER_VDD == layer)
    return model->c4->vdd_padidx[idx];
  else if(LAYER_GND == layer)
    return model->c4->gnd_padidx[idx];
  else
    fatal("Wrong arguemnt layer received by PadIdx_to_GridIdx!\n");

  return 0;
}

int GridIdx_to_PadIdx(model_t *model, int idx, int layer)
{
  int size = model->c4->pad_grid_row*model->c4->pad_grid_col;
  int i;

  if(LAYER_VDD == layer){
      for(i=0; i<size; i++){
          if(idx == model->c4->vdd_padidx[i])
            return i;
      }
  }
  else if(LAYER_GND == layer){
      for(i=0; i<size; i++){
          if(idx == model->c4->gnd_padidx[i])
            return i;
      }
  }
  else
    fatal("Wrong arguemnt layer received by GridIdx_to_PadIdx!\n");

  return -1;
}

char *GridIdx_to_UnitName(model_t *model, int idx, int layer)
{
  int i,j;
  int nc = model->cols;
  blist_t *ptr;

  j = idx % nc;
  i = (int) floor(idx/nc);

  ptr = model->layers[layer].b2gmap[i][j];

  return model->layers[layer].flp->units[ptr->idx].name;
}

int compare_domain(model_t *model, int idx, int offset)
{
  int ia,ja,ib,jb;
  int nc = model->cols;
  blist_t *ptra, *ptrb;
  int mtdmn  = model->config.PDN_multi_dom;

  if(!mtdmn)
    return 1;

  ja = idx % nc;
  ia = (int) floor(idx/nc);
  jb = (idx+offset) % nc;
  ib = (int) floor((idx+offset)/nc);

  ptra = model->layers[0].b2gmap[ia][ja];
  ptrb = model->layers[0].b2gmap[ib][jb];

  if(model->layers[0].flp->units[ptra->idx].domain == 
     model->layers[0].flp->units[ptrb->idx].domain)
    return 1;
  else
    return 0;
}

// when multidomain, return whether a grid branch exist
int does_branch_exist(model_t *model, int idx)
{
  int ia,ja,ib,jb;
  int nr = model->rows;
  int nc = model->cols;
  blist_t *ptra, *ptrb;
  int mtdmn  = model->config.PDN_multi_dom;

  if(!mtdmn)
    return 1;

  if(idx < ((nc-1)*nr)){
      ja = idx % (nc-1);
      ia = (int) floor(idx/(nc-1));
      jb = ja + 1;
      ib = ia;
  }
  else{
      ja = (idx-(nc-1)*nr) % nc;
      ia = (int) floor((idx-(nc-1)*nr)/nc);
      jb = ja;
      ib = ia + 1;
  }

  ptra = model->layers[0].b2gmap[ia][ja];
  ptrb = model->layers[0].b2gmap[ib][jb];

  if(model->layers[0].flp->units[ptra->idx].domain == 
     model->layers[0].flp->units[ptrb->idx].domain){
      return 1;
  }
  else{
      return 0;
  }
}

void set_flp_domain(PDN_flp_t *flp, PDN_config_t *config)
{
  int i;
  int core_id;
  int num_dom = 0;
  char *name;

  for(i=0; i<flp->n_units; i++){
      name = flp->units[i].name;

      core_id = get_core_id(name);

      if(core_id < 0){
          flp->units[i].domain = 0;
      }
      else{
          if(config->PDN_multi_dom > 0){
              // assumes that core id starts with 1
              flp->units[i].domain = floor((core_id-1) / config->PDN_multi_dom) + 1;
          }
          else if(-1 == config->PDN_multi_dom){
              // non-uniform domain size for 16 core
              if(core_id > 15)
                flp->units[i].domain = 5;
              else if(core_id > 14)
                flp->units[i].domain = 4;
              else if(core_id > 12)
                flp->units[i].domain = 3;
              else if(core_id > 8)
                flp->units[i].domain = 2;
              else
                flp->units[i].domain = 1;
          }
          else
            fatal("Unsupported multi_dom config");
      }

      if(flp->units[i].domain > num_dom)
        num_dom = flp->units[i].domain; 
  }

  /* # of conceptual domains, domains with same id 
   * are not guaranteed to connect together */
  flp->n_domain = num_dom + 1;
}

int get_core_id(char *blk_name)
{
  int core_id;
  int length = strlen(blk_name);

  // MC does not belong to any core
  if(blk_name[0] == 'M' &&
     blk_name[1] == 'C'){
      if(blk_name[length-2] <= '9' &&
         blk_name[length-2] >= '0')
        core_id = -1 * (10 * (blk_name[length-2]-'0') + (blk_name[length-1]-'0'));
      else
        core_id = -1 * (blk_name[length-1]-'0');
  }
  else{
      if(blk_name[length-2] <= '9' &&
         blk_name[length-2] >= '0')
        core_id = 10*(blk_name[length-2]-'0') + (blk_name[length-1]-'0');
      else
        core_id = (blk_name[length-1]-'0');
  }

  if(0 == core_id)
    fatal("Core ID should start with 1\n");

  if(core_id > MAX_CORE_NUM)
    fatal("Core id exceeds maximum value\n");

  return core_id;
}

/* parse the layer file open for reading	*/
void parse_layer_file_PDN(model_t *model, FILE *fp)
{
  char line[LINE_SIZE], *ptr;
  FILE *mfp = NULL;
  char str[STR_SIZE];
  int count, i = 0;
  int field = LCF_3D_SNO, ival;
  int num_layer;

  fseek(fp, 0, SEEK_SET);
  count = 0;
  while (!feof(fp) && count < (model->n_layers * LCF_3D_NPARAMS)) {
      fgets(line, LINE_SIZE, fp);
      if (feof(fp))
        break;

      /* ignore comments and empty lines	*/
      ptr = strtok(line, " \r\t\n");
      if (!ptr || ptr[0] == '#')
        continue;

      switch (field) 
        {
        case LCF_3D_SNO:
          if (sscanf(ptr, "%d", &ival) != 1)
            fatal("invalid layer number\n");
          if(ival >= model->n_layers || ival < 0)
            fatal("layer number must be >= 0 and < no. of layers\n");
          if (model->layers[ival].no != 0)
            fatal("layer numbers must be unique\n");
          i = ival;
          model->layers[i].no = ival;
          break;
        case LCF_3D_FLP:
          model->layers[i].flp = PDN_read_flp(ptr);
          if(!eq(model->width, PDN_get_total_width(model->layers[i].flp)) || 
             !eq(model->height, PDN_get_total_height(model->layers[i].flp)))
            fatal("width and height differ across layers\n");
          break;
        case LCF_3D_MLCF:
          mfp = fopen(ptr, "r");
          if(!mfp) {
              sprintf(str, "error opening mlcf file %s in lcf file\n", ptr);
              fatal(str);
          }
          num_layer = count_significant_lines(mfp);
          if (num_layer % MLCF_NPARAMS)
            fatal("wrong no. of lines in metal layer file\n");
          num_layer /= MLCF_NPARAMS;
          model->layers[i].metal_layers = alloc_metal_layers(num_layer);
          parse_metal_layer_file(model, ptr, i);
          fclose(mfp);
          break;
        case LCF_3D_TSVL:
          strcpy(model->layers[i].tsv.file, ptr);
          model->layers[i].tsv.loc = imatrix(model->rows, model->cols);
          zero_ivector(model->layers[i].tsv.loc[0], model->rows*model->cols);
          break;
        default:
          fatal("invalid field id\n");
          break;
        }
      field = (field + 1) % LCF_3D_NPARAMS;
      count++;
  }

  /* allocate the block-grid maps */
  for(i=0; i < model->n_layers; i++) {
      model->layers[i].b2gmap = new_b2gmap_PDN(model->rows, model->cols);
  }
}

/* parse the layer file open for reading	*/
void parse_metal_layer_file(model_t *model, char *file, int layer)
{
  FILE *fp;
  char str[LINE_SIZE], line[LINE_SIZE], *ptr;
  int count, i = 0, j, field = MLCF_SNO, ival;
  double dval;
  int n = model->layers[layer].metal_layers.n_metal;
  int *set_flag = ivector(n);
  zero_ivector(set_flag, n);

  fp = fopen (file, "r");
  if (!fp) {
      sprintf(str, "error opening mlcf file %s\n", file);
      fatal(str);
  }

  fseek(fp, 0, SEEK_SET);
  count = 0;
  while (!feof(fp) && count < (n * MLCF_NPARAMS)) {
      fgets(line, LINE_SIZE, fp);
      if (feof(fp))
        break;

      /* ignore comments and empty lines	*/
      ptr = strtok(line, " \r\t\n");
      if (!ptr || ptr[0] == '#')
        continue;

      switch (field) {
        case MLCF_SNO:
          if (sscanf(ptr, "%d", &ival) != 1)
            fatal("invalid layer number\n");
          if(ival >= n || ival < 0)
            fatal("layer number must be >= 0 and < no. of metal layers\n");
          if (set_flag[ival] == TRUE)
            fatal("layer numbers must be unique\n");
          i = ival;
          set_flag[i] = TRUE;
          break;
        case MLCF_PITCH:
          if (sscanf(ptr, "%lf", &dval) != 1)
            fatal("invalid pitch\n");
          model->layers[layer].metal_layers.geo[i].pitch = dval;
          break;
        case MLCF_WIDTH:
          if (sscanf(ptr, "%lf", &dval) != 1)
            fatal("invalid width\n");
          model->layers[layer].metal_layers.geo[i].width = dval;
          break;
        case MLCF_THICK:
          if (sscanf(ptr, "%lf", &dval) != 1)
            fatal("invalid thick\n");
          model->layers[layer].metal_layers.geo[i].thick = dval;
          break;
        case MLCF_RHO:
          if (sscanf(ptr, "%lf", &dval) != 1)
            fatal("invalid rho\n");
          model->layers[layer].metal_layers.geo[i].rho = dval;
          break;
        case MLCF_DIREC:
          if (sscanf(ptr, "%d", &ival) != 1)
            fatal("invalid direction\n");
          if ((MLCF_X != ival) && (MLCF_Y != ival))
            fatal("direction input should either 0 or 1\n");
          model->layers[layer].metal_layers.geo[i].direc = ival;
          break;
        default:
          fatal("invalid field id\n");
          break;
      }
      field = (field + 1) % MLCF_NPARAMS;
      count++;
  }
  fclose(fp);

  for(j=0; j<n; j++){
      if(set_flag[j] != TRUE){
          sprintf(str, "Metal layer %d was left uninitialized\n", j);
          fatal(str);
      }
  }
  free_ivector(set_flag);

  //sanity check
  //for simplicity, we force that n_metal must be even
  //and adjacent layers cannot be in the same direction
  //for now, we force that layer 0 must be x-direc
  if(n%2)
    fatal("Please specify even number of metal layers\n");

  for(i=0; i<n; i++){
      if(i%2 != model->layers[layer].metal_layers.geo[i].direc)
        fatal("Please make sure that layer 0 is x, 1 is y, etc. \n");
  }
}

void parse_tsv_loc(model_t *model, int layer)
{
  char str[LINE_SIZE], copy[LINE_SIZE]; 
  char s[STR_SIZE];
  char name[STR_SIZE];
  int grid_x, grid_y;
  char *ptr;
  FILE *fp;

  /* short cuts */  
  int nr = model->rows;
  int nc = model->cols;
  int **tloc = model->layers[layer].tsv.loc;

  if (!strcasecmp(model->layers[layer].tsv.file, "stdin"))
    fp = stdin;
  else
    fp = fopen (model->layers[layer].tsv.file, "r");

  if (!fp) {
      sprintf(s, "error opening tsvloc file %s\n", model->layers[layer].tsv.file);
      fatal(s);
  }

  fseek(fp, 0, SEEK_SET);
  while(!feof(fp)) {
      fgets(str, LINE_SIZE, fp);
      if (feof(fp))
        break;
      strcpy(copy, str);

      /* ignore comments and empty lines */
      ptr = strtok(str, " \r\t\n");
      if (!ptr || ptr[0] == '#')
        continue;

      if (sscanf(copy, "%s%d%d", name, &grid_x, &grid_y) == 3) {
          if ((grid_x >= model->cols) || (grid_y >= model->rows) ||
              (grid_x < 0) || (grid_y < 0)){
              printf("x = %d, y = %d\n", grid_x, grid_y);
              printf("grid_size = %d*%d\n", model->cols, model->rows);
              fatal("TSV location does not fit in current grid!\n");
          }

          if ((VDDTSV == tloc[nr - grid_y - 1][grid_x]) ||
              (GNDTSV == tloc[nr - grid_y - 1][grid_x]))
            warning("Duplication exists in tsv location file .tsvloc\n");

          if (!strcmp(name, "V")){
              tloc[nr - grid_y - 1][grid_x] = VDDTSV;
          }
          else if (!strcmp(name, "G")){
              tloc[nr - grid_y - 1][grid_x] = GNDTSV;
          }
          else
            fatal("TSV location is neither V(vdd) nor G(ground)\n");
      }
      else
        fatal("invalid tsv location file format\n");
  }

  if(fp != stdin)
    fclose(fp);
}

void parse_IVR_loc(model_t *model, char *file)
{
  char str[LINE_SIZE], copy[LINE_SIZE]; 
  char s[STR_SIZE];
  int grid_x, grid_y;
  char *ptr;
  FILE *fp;

  /* short cuts */  
  int nr = model->rows;
  int nc = model->cols;

  if (!strcasecmp(file, "stdin"))
    fp = stdin;
  else
    fp = fopen (file, "r");

  if (!fp) {
      sprintf(s, "error opening IVR file %s\n", file);
      fatal(s);
  }

  fseek(fp, 0, SEEK_SET);
  while(!feof(fp)) {
      fgets(str, LINE_SIZE, fp);
      if (feof(fp))
        break;
      strcpy(copy, str);

      /* ignore comments and empty lines */
      ptr = strtok(str, " \r\t\n");
      if (!ptr || ptr[0] == '#')
        continue;

      // parse index as metal grid
      if (sscanf(copy, "%d%d", &grid_x, &grid_y) == 2) {
          if ((grid_x >= nc) || (grid_y >= nr) || (grid_x < 0) || (grid_y < 0)){
              printf("x = %d, y = %d\n", grid_x, grid_y);
              printf("grid_size = %d*%d\n", nc, nr);
              fatal("IVR location does not fit in current grid!\n");
          }
          if (IVRLOC == model->sc_converter->loc[nr - grid_y - 1][grid_x])
            warning("Duplication exists in tsv location file\n");
          model->sc_converter->loc[nr - grid_y - 1][grid_x] = IVRLOC;
      }
      else
        fatal("invalid pad location file format\n");
  }

  if(fp != stdin)
    fclose(fp);
}

void switch_SCconverters(model_t *model, int counter, int intvl)
{
  int i, j, l;
  int ivr_counter = 0;
  int num_intvl;
  int grididx;
  int head_idx, foot_idx;
  int t1_idx, t2_idx; int b1_idx, b2_idx;
  double Vt1_p, Vt2_p; double Vb1_p, Vb2_p;
  double Vhead, Vfoot;

  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int ngtsv = model->layers[0].tsv.num_gnd;
  int ntsv = ngtsv*(nl-1);
  int **ivr_loc = model->sc_converter->loc;
  double *g = model->last_trans;

  if(!(counter % intvl)){
      num_intvl = counter / intvl;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            if(IVRLOC == ivr_loc[i][j]){
                grididx = i*nc + j;
                for(l=0; l<nl-1; l++){
                    t1_idx = nvp+ngp+nl*2*nr*nc+ntsv+4*ivr_counter+0;
                    t2_idx = nvp+ngp+nl*2*nr*nc+ntsv+4*ivr_counter+1;
                    b1_idx = nvp+ngp+nl*2*nr*nc+ntsv+4*ivr_counter+2;
                    b2_idx = nvp+ngp+nl*2*nr*nc+ntsv+4*ivr_counter+3;
                    head_idx = nvp+ngp+(l+1)*2*nr*nc+grididx;
                    foot_idx = nvp+ngp+l*2*nr*nc+nr*nc+grididx;

                    Vt1_p = g[t1_idx]; Vt2_p = g[t2_idx];
                    Vb1_p = g[b1_idx]; Vb2_p = g[b2_idx];
                    Vhead = g[head_idx]; Vfoot = g[foot_idx];

                    if(num_intvl % 2){
                        g[t1_idx] = Vhead + Vfoot - Vb1_p;
                        g[b1_idx] = Vhead + Vfoot - Vt1_p;
                    }
                    else{
                        g[t2_idx] = Vhead + Vfoot - Vb2_p;
                        g[b2_idx] = Vhead + Vfoot - Vt2_p;
                    }
                    ivr_counter++;
                }
            }
        }
  }
}
