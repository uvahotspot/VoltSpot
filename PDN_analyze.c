#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>

#include "PDN_sim.h"
#include "PDN_analyze.h"
#include "util.h"
#include "pad.h"

static int stat_counter =0;

status_t *alloc_status(model_t *model)
{
  int l, i;
  status_t *status;

  status = (status_t *) calloc(1, sizeof(status_t));

  status->trans_counter = 0;
  status->draw_counter = 0;

  status->maxcur_vdd = 0;
  status->maxcur_gnd = 0;
  status->max_pad_dense = 0;
  status->curIth = 0;
  status->curIth_vdd = 0;
  status->curIth_gnd = 0;
  status->maxdroop = 0;
  status->curIc = 0;
  status->avgIc = 0;
  status->dIc = 0;
  status->curVc = 0;
  status->avgVc = 0;
  status->vio_area_ratio = 0;
  status->pkgDrop = 0;
  status->PDN_power = 0;
  status->sumcur = 0;

  status->prevIc = 0;
  status->prevVc = 1;

  /* for averaging resutls*/
  int m = model->trans_matrix_dim;
  status->cycle_avg  = dvector(m);
  zero_dvector(status->cycle_avg, m);

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nbr = 2*nr*nc - nr - nc;

  status->maxIR = dvector(nl);
  status->layer_maxdroop = dvector(nl);

	/* tsv current */
  int ntsv = 0;
  for(l=0; l<nl-1; l++){
      ntsv += model->layers[l].tsv.num_gnd;
      ntsv += model->layers[l].tsv.num_vdd;
  }
  status->TSVcur = dvector(ntsv);
  status->maxTSVcur = dvector(nl-1);

	/* sensor location */
  status->sensor_loc = imatrix(nr, nc);
  zero_ivector(status->sensor_loc[0], nr*nc);
  status->sensor_noise = dmatrix(nr, nc);
  zero_dvector(status->sensor_noise[0], nr*nc);

	/* pad current */
  status->vdd_cur = dmatrix(nr, nc);
  status->gnd_cur = dmatrix(nr, nc);
  zero_dvector(status->vdd_cur[0], nr*nc);
  zero_dvector(status->gnd_cur[0], nr*nc);

	/* tsv current */

  /* allocate voltage record */
  status->gridstats.counter_2D  = imatrix(nr, nc);
  status->gridstats.max_2D      = dmatrix(nr, nc);
  status->gridstats.integral_2D = dmatrix(nr, nc);
  status->blkstats.counter_1D   = ivector(model->total_n_blocks);
  status->blkstats.max_1D       = dvector(model->total_n_blocks);
  status->corestats.max_1D      = dvector(MAX_CORE_NUM);
  zero_ivector(status->gridstats.counter_2D[0], nr*nc);
  negL_dvector(status->gridstats.max_2D[0], nr*nc);
  zero_dvector(status->gridstats.integral_2D[0], nr*nc);
  zero_ivector(status->blkstats.counter_1D, model->total_n_blocks);
  negL_dvector(status->blkstats.max_1D, model->total_n_blocks);
  negL_dvector(status->corestats.max_1D, MAX_CORE_NUM);

  /* allocate current record */
  //for trans WP
  if(model->config.vgradient_analyse){
      status->vgradient.max_1D = dvector(9*nbr);
      zero_dvector(status->vgradient.max_1D, 6*nbr);
      for(i = 0; i<3*nbr; i++)
        status->vgradient.max_1D[i] = -1000.0;
      for(i = 3*nbr; i<6*nbr; i++)
        status->vgradient.max_1D[i] = 1000.0;
  }
  //if(model->config.vgradient_analyse){
  //    status->vgradient.max_1D = dvector(2*nbr);
  //    zero_dvector(status->vgradient.max_1D, 2*nbr);
  //}

  return status;
}

void free_status(model_t *model, status_t *status)
{
  free_dvector(status->cycle_avg);
  free_imatrix(status->sensor_loc);
  free_dmatrix(status->sensor_noise);
  free_dvector(status->maxIR);
  free_dvector(status->layer_maxdroop);
  free_dvector(status->TSVcur);
  free_dvector(status->maxTSVcur);
  free_imatrix(status->gridstats.counter_2D);
  free_dmatrix(status->gridstats.max_2D);
  free_dmatrix(status->gridstats.integral_2D);
  free_ivector(status->blkstats.counter_1D);
  free_dvector(status->blkstats.max_1D);
  free_dvector(status->corestats.max_1D);
  free_dmatrix(status->vdd_cur);
  free_dmatrix(status->gnd_cur);
  if(model->config.vgradient_analyse)
    free_dvector(status->vgradient.max_1D);
  free(status);
}

void populate_status(model_t *model, status_t *status)
{
  if (strcmp(model->config.senloc_file, NULLFILE)){
      parse_sensor_loc(model, status, model->config.senloc_file);
  }
}

void parse_sensor_loc(model_t *model, status_t *status, char *file)
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
      sprintf(s, "error opening sensor_loc file %s\n", file);
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

      if (sscanf(copy, "%d%d", &grid_x, &grid_y) == 2) {
          if ((grid_x >= model->cols) || (grid_y >= model->rows) ||
              (grid_x < 0) || (grid_y < 0)){
              printf("x = %d, y = %d\n", grid_x, grid_y);
              printf("grid_size = %d*%d\n", nc,nr);
              fatal("Sensor location does not fit in current grid!\n");
          }
          if (SENSOR == status->sensor_loc[nr - grid_y - 1][grid_x])
            warning("Duplication exists in sensor location file\n");
          status->sensor_loc[nr - grid_y - 1][grid_x] = SENSOR;
      }
      else
        fatal("invalid pad location file format\n");
  }

  if(fp != stdout)
    fclose(fp);	
}

void PDN_steady_analyze(model_t *model, status_t *status)
{
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;

  //pad current
  get_pad_current(model, status);

  //TSV current
  if(model->is_3D)
    get_TSV_current(model, status);

  /* IR drop*/
  get_maxIR(model, status->maxIR);

  //calculate static PDN power loss
  status->PDN_power = get_steady_PDN_power(model);
}

void print_steady_analyze(model_t *model, status_t *status)
{
  int l;
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;

  int p_count = 1;
  printf("%d: chip dimension (width/height, meter):  %lf/%lf\n", 
         p_count, model->width, model->height);p_count++;
  printf("%d: supply voltage vdd (Volt):             %.2lf\n", 
         p_count, model->config.vdd);p_count++;
  printf("%d: virtual grid size (num_cols/num_rows): %d / %d\n",
         p_count, nc, nr);p_count++;
  printf("%d: pad grid size (num_cols/num_rows):     %d / %d\n",
         p_count, model->c4->pad_grid_col, model->c4->pad_grid_row);p_count++;
  printf("%d: number of pads (vdd/gnd):              %d / %d\n",
         p_count, model->c4->vdd_num, model->c4->gnd_num);p_count++;
  printf("%d: sum of current (A):                    %lf\n",
         p_count, status->sumcur);p_count++;
  printf("%d: Max pad current (A):                   %.4lf / %.4lf\n",
         p_count, status->maxcur_vdd, status->maxcur_gnd);p_count++;
  printf("%d: max on-chip IR drop (%%Vdd):            ", 
         p_count);p_count++;
  for(l=0; l<nl; l++)
    printf("%.3lf / ", status->maxIR[l]);
  printf("\n");
  if(model->is_3D){
      printf("%d: TSV Count:                             ", 
             p_count);p_count++;
      for(l=0; l<nl-1; l++)
        printf("%d/", model->layers[l].tsv.num_gnd + model->layers[l].tsv.num_vdd);
      printf("\n");

      printf("%d: max TSV current (A):                   ", 
             p_count);p_count++;
      for(l=0; l<nl-1; l++)
        printf("%.4lf/", status->maxTSVcur[l]);
      printf("\n");
  }
  printf("%d: PDN static power loss (W):             %lf\n", 
         p_count, status->PDN_power);p_count++;
}

void PDN_trans_analyze(model_t *model, status_t *status)
{
  int i;
  int m = model->trans_matrix_dim;
  double l_pkg_p = model->config.PDN_pkg_pL;
  double r_pkg_p = model->config.PDN_pkg_pR;
  double delta_t = 1/model->config.proc_clock_freq;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  for(i=0; i<m; i++)
    status->cycle_avg[i] += model->last_trans[i];

  if(!(status->trans_counter % model->config.PDN_step_percycle)){
      for(i=0; i<m; i++){
          status->cycle_avg[i] /= model->config.PDN_step_percycle;
      }

      if(model->config.PDN_pkgLC){
          /* Package node droop */
          status->curIc   = status->cycle_avg[m-1];
          status->curVc   = status->cycle_avg[m-2];
          status->avgIc   = (status->curIc + status->prevIc) / 2;
          status->avgVc   = (status->curVc + status->prevVc) / 2;
          status->dIc     = status->curIc - status->prevIc;

          status->pkgDrop = status->avgIc*r_pkg_p + status->avgVc + l_pkg_p*status->dIc/delta_t;
          status->pkgDrop = 100 * ((vdd-gnd) - status->pkgDrop) / (vdd-gnd);

          status->prevIc = status->curIc;
          status->prevVc = status->curVc;
      }
      /* Onchip node droop */
      if(!model->is_3D)
        onchipV_analyze_2D(model, status);
      else
        onchipV_analyze_3D(model, status);

      zero_dvector(status->cycle_avg, m);
  }

  if((model->config.animation) &&
     !(status->trans_counter % model->config.frame_intv)){
      draw_single_gif(model, status, status->draw_counter, TRANSIENT);
      status->draw_counter++;
  }
}

void print_trans_header(model_t *model, status_t *status, FILE *fp)
{
  int i, j, l;
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int **sloc = status->sensor_loc;

  if(model->config.PDN_pkgLC){
      //fprintf(fp, "pkgCur\t");
      fprintf(fp, "pkgDrop\t");
  }
  fprintf(fp, "max_onchip_drop\t");
  //fprintf(fp, "vio_area_ratio\t");

  //for(i=0; i<MAX_CORE_NUM; i++)
  //  if(-LARGENUM != status->corestats.max_1D[i])
  //    fprintf(fp, "C_%d\t", i);

  if(model->is_3D){
      for(l=0; l<nl; l++)
        fprintf(fp, "Layer_%d\t", l);
  }

  if(strcmp(model->config.senloc_file, NULLFILE)){
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++)
          if(SENSOR == sloc[i][j])
            fprintf(fp, "S_%d_%d\t", j, nr-i-1);
  }

  fprintf(fp, "\n");
}

void print_trans_analyze(model_t *model, status_t *status, FILE *fp)
{
  int i, j, l;
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int **sloc = status->sensor_loc;

  if(model->config.PDN_pkgLC){
      //fprintf(fp, "%.4lf\t", status->curIc);
      fprintf(fp, "%.4lf\t", status->pkgDrop);
  }
  fprintf(fp, "%.4lf\t", status->maxdroop);
  //fprintf(fp, "%.4lf\t", status->vio_area_ratio);

  //for(i=0; i<MAX_CORE_NUM; i++){
  //    if(-LARGENUM != status->corestats.max_1D[i]){
  //        fprintf(fp, "%lf\t", status->corestats.max_1D[i]);
  //    }
  //}

  if(model->is_3D){
      for(l=0; l<nl; l++)
        fprintf(fp, "%.4lf\t", status->layer_maxdroop[l]);
  }

  if (strcmp(model->config.senloc_file, NULLFILE)){
      for(i=0; i < nr; i++){
          for(j=0; j < nc; j++){
              if (SENSOR == sloc[i][j]){
                  fprintf(fp, "%.4lf\t", status->sensor_noise[i][j]);
              }
          }
      }
  }
  fprintf(fp, "\n");
}

void print_step_singlenode(model_t *model, FILE *fp, int pi, int pj)
{
  double vdroop, vnoise;

  int nr = model->rows;
  int nc = model->cols;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nbr = 2*nr*nc - nr - nc;
  int nml = model->layers[0].metal_layers.n_metal/2;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  int cur_idx = pi*nc + pj;
  if(model->config.PDN_gridL)
    vdroop = model->last_trans[nvp+cur_idx] - model->last_trans[nvp+nr*nc+nbr*nml+ngp+cur_idx];
  else
    vdroop = model->last_trans[nvp+cur_idx] - model->last_trans[nvp+nr*nc+ngp+cur_idx];

  vnoise = 100*(vdd-gnd-vdroop)/(vdd-gnd);

  fprintf(fp, "%lf\t%lf\n", model->last_trans[nvp+cur_idx], vnoise);
}

void print_status_singlenode(model_t *model, status_t *status, FILE *fp, int pi, int pj)
{
  double vdroop, vnoise;

  int nr = model->rows;
  int nc = model->cols;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nbr = 2*nr*nc - nr - nc;
  int nml = model->layers[0].metal_layers.n_metal/2;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  if((pi>=nr) || (pj>=nc)){
      fatal("out of bound! func. print_status_singlenode\n");
  }

  int cur_idx = pi*nc + pj;
  if(model->config.PDN_gridL)
    vdroop = status->cycle_avg[nvp+cur_idx] - status->cycle_avg[nvp+nr*nc+nbr*nml+ngp+cur_idx];
  else
    vdroop = status->cycle_avg[nvp+cur_idx] - status->cycle_avg[nvp+nr*nc+ngp+cur_idx];

  vnoise = 100*(vdd-gnd-vdroop)/(vdd-gnd);

  fprintf(fp, "%lf\t%lf\n", status->cycle_avg[nvp+cur_idx], vnoise);

}

// get current for pad and top level metal
void get_pad_current(model_t *model, status_t *status)
{
  int i, j, k;
  double r_tsv_pad;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int vs = model->config.v_stacking;
  double **vcur = status->vdd_cur;
  double **gcur = status->gnd_cur;
  int **vloc = model->c4->vdd_loc;
  int **gloc = model->c4->gnd_loc;
  double *v = model->last_steady;
  double r_pad = model->c4->pad_r;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;
  double pkg_vdd = model->last_steady[2*nl*nr*nc+PKG_VDD];
  double pkg_gnd = model->last_steady[2*nl*nr*nc+PKG_GND];
  double padD = model->config.PDN_padD/2;
  double pad_area = PI * padD * padD;

  /* pad current */
  if(!vs){
      for(i = 0; i < nr; i++)
        for(j = 0; j < nc; j++){
            if (vloc[i][j] & PGPAD)
              vcur[i][j] = (pkg_vdd - v[i*nc + j]) / r_pad;
            else
              vcur[i][j] = 0;

            if (gloc[i][j] & PGPAD)
              gcur[i][j] = (v[nl*nr*nc + i*nc + j] - pkg_gnd) / r_pad;
            else
              gcur[i][j] = 0;
        }
  }
  else{
      r_tsv_pad = r_pad;
      for(k=0; k<nl-1; k++)
        r_tsv_pad += model->layers[k].tsv.r;

      for(i = 0; i < nr; i++)
        for(j = 0; j < nc; j++){
            if (vloc[i][j] & PGPAD)
              vcur[i][j] = (nl*vdd - v[(2*nl-1)*nr*nc + i*nc + j]) / r_tsv_pad;
            else
              vcur[i][j] = 0;

            if (gloc[i][j] & PGPAD)
              gcur[i][j] = (v[i*nc + j] - gnd) / r_pad;
            else
              gcur[i][j] = 0;
        }
  }

  /* Vdd pad current */
  status->sumcur = sum_dmatrix(vcur, nr, nc);
  status->maxcur_vdd = max_dmatrix_skip0(vcur, nr, nc);
  status->curIth_vdd = above_threshold_dmatrix(vcur, nr, nc, model->c4->Ith);
  /* Gnd pad current */
  status->maxcur_gnd = max_dmatrix_skip0(gcur, nr, nc);
  status->curIth_gnd = above_threshold_dmatrix(gcur, nr, nc, model->c4->Ith);

  if(status->maxcur_vdd > status->maxcur_gnd)
    status->max_pad_dense = status->maxcur_vdd;
  else
    status->max_pad_dense = status->maxcur_gnd;

  status->max_pad_dense /= pad_area;
}

void get_TSV_current(model_t *model, status_t *status)
{
  int i, j, l;
  int a_grid_idx, b_grid_idx;
  double cur_val = 0;
  double layer_max, layer_sum;
  double Rv;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int vs = model->config.v_stacking;
  double *v = model->last_steady;

  int counter = 0;
  for(l=0; l<nl-1; l++){
      layer_max = 0;
      //layer_sum = 0;
      Rv = model->layers[l].tsv.r;

      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            if(!vs){
                if((VDDTSV == model->layers[l].tsv.loc[i][j])){
                    a_grid_idx = l*nr*nc + i*nc + j;
                    b_grid_idx = a_grid_idx + nr*nc;
                    cur_val = (v[a_grid_idx] - v[b_grid_idx]) / Rv;
                    status->TSVcur[counter] = cur_val;
                    counter++;
                }
                else if((GNDTSV == model->layers[l].tsv.loc[i][j])){
                    a_grid_idx = nl*nr*nc + l*nr*nc + i*nc + j;
                    b_grid_idx = a_grid_idx + nr*nc;
                    cur_val = (v[b_grid_idx] - v[a_grid_idx]) / Rv;
                    status->TSVcur[counter] = cur_val;
                    counter++;
                    //layer_sum += cur_val;
                }
            }
            else{
                a_grid_idx = l*2*nr*nc + nr*nc + i*nc + j;
                b_grid_idx = a_grid_idx + nr*nc;
                if((GNDTSV == model->layers[l].tsv.loc[i][j])){
                    cur_val = (v[b_grid_idx] - v[a_grid_idx]) / Rv;
                    status->TSVcur[counter] = cur_val;
                    counter++;
                    //layer_sum += cur_val;
                }
            }
            if(cur_val > layer_max)
              layer_max = cur_val;
        }
      status->maxTSVcur[l] = layer_max;
      //printf("%.2lf\n", layer_sum);
  }
}

double get_steady_PDN_power(model_t *model)
{
  int i, j, l, k, grid_idx, direc;
  int a_grid_idx, b_grid_idx;
  double r_tsv_pad;
  double pad_cur, branch_cur, tsv_cur;
  double cur_sum, cur_sq_sum, vddpad_cur_sq_sum ;
  double pad_power, pkg_power, onchip_power, tsv_power;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int vs = model->config.v_stacking;
  int **vloc = model->c4->vdd_loc;
  int **gloc = model->c4->gnd_loc;
  double *v = model->last_steady;
  double r_pad = model->c4->pad_r;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;
  double pkg_vdd = model->last_steady[2*nl*nr*nc+PKG_VDD];
  double pkg_gnd = model->last_steady[2*nl*nr*nc+PKG_GND];
  double pkg_r = model->config.PDN_pkg_sR;
  double rv; double rx; double ry;

  if(!vs){
      cur_sum = 0;
      cur_sq_sum = 0;
      /* pad power */
      l = 0;
      for(i = 0; i < nr; i++)
        for(j = 0; j < nc; j++){
            grid_idx = l*nr*nc + i*nc + j;
            if (vloc[i][j] & PGPAD){
                pad_cur = (pkg_vdd - v[grid_idx]) / r_pad;
                cur_sq_sum += pad_cur*pad_cur;
                cur_sum += pad_cur;
            }
            if (gloc[i][j] & PGPAD){
                pad_cur = (v[nl*nr*nc + grid_idx] - pkg_gnd) / r_pad;
                cur_sq_sum += pad_cur*pad_cur;
            }
        }
      pad_power = cur_sq_sum * r_pad;

      /* on-chip power */
      onchip_power = 0;
      for(l=0; l<nl; l++){
          rx = 0; ry = 0;
          //calculate grid R
          for(i=0; i<model->layers[l].metal_layers.n_metal; i++){
              direc = model->layers[l].metal_layers.geo[i].direc;
              if(MLCF_X == direc)
                rx += 1/model->layers[l].metal_layers.gridRL[i].r;
              else
                ry += 1/model->layers[l].metal_layers.gridRL[i].r;
          }
          rx = 1/rx; ry = 1/ry;

          cur_sq_sum = 0;
          for(i=0; i<nr; i++)
            for(j=0; j<nc-1; j++){
                grid_idx = l*nr*nc + i*nc + j;
                branch_cur = (v[grid_idx] - v[grid_idx+1])/rx;
                cur_sq_sum += branch_cur*branch_cur;
                branch_cur = (v[nl*nr*nc+grid_idx] - v[nl*nr*nc+grid_idx+1])/rx;
                cur_sq_sum += branch_cur*branch_cur;
            }
          onchip_power += cur_sq_sum * rx;

          cur_sq_sum = 0;
          for(i=0; i<nr-1; i++)
            for(j=0; j<nc; j++){
                grid_idx = l*nr*nc + i*nc + j;
                branch_cur = (v[grid_idx] - v[grid_idx+nc])/ry;
                cur_sq_sum += branch_cur*branch_cur;
                branch_cur = (v[nl*nr*nc+grid_idx] - v[nl*nr*nc+grid_idx+nc])/ry;
                cur_sq_sum += branch_cur*branch_cur;
            }
          onchip_power += cur_sq_sum * ry;
      }

      cur_sq_sum = 0;
      tsv_power = 0;
      /* TSV power */
      for(l=0; l<nl-1; l++){
          rv = model->layers[l].tsv.r;

          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                if((VDDTSV == model->layers[l].tsv.loc[i][j])){
                    a_grid_idx = l*nr*nc + i*nc + j;
                    b_grid_idx = a_grid_idx + nr*nc;
                    tsv_cur = (v[a_grid_idx] - v[b_grid_idx]) / rv;
                    cur_sq_sum += tsv_cur*tsv_cur;
                }
                else if((GNDTSV == model->layers[l].tsv.loc[i][j])){
                    a_grid_idx = nl*nr*nc + l*nr*nc + i*nc + j;
                    b_grid_idx = a_grid_idx + nr*nc;
                    tsv_cur = (v[b_grid_idx] - v[a_grid_idx]) / rv;
                    cur_sq_sum += tsv_cur*tsv_cur;
                }
            }
          tsv_power += cur_sq_sum * rv; 
      }

      pkg_power = 2 * cur_sum * cur_sum * pkg_r;
  }
  else{
      /* pad power */
      cur_sq_sum = 0;
      vddpad_cur_sq_sum = 0;
      r_tsv_pad = r_pad;
      for(k=0; k<nl-1; k++)
        r_tsv_pad += model->layers[k].tsv.r;

      for(i = 0; i < nr; i++)
        for(j = 0; j < nc; j++){
            if (vloc[i][j] & PGPAD){
                pad_cur = (nl*vdd - v[(2*nl-1)*nr*nc + i*nc + j]) / r_tsv_pad;
                cur_sq_sum += pad_cur*pad_cur;
                vddpad_cur_sq_sum += pad_cur*pad_cur;
            }
            if (gloc[i][j] & PGPAD){
                pad_cur = (v[i*nc + j] - gnd) / r_pad;
                cur_sq_sum += pad_cur*pad_cur;
            }
        }
      pad_power = cur_sq_sum * r_pad;

      for(l=0; l<nl; l++){
          rx = 0; ry = 0;
          //calculate grid R
          for(i=0; i<model->layers[l].metal_layers.n_metal; i++){
              direc = model->layers[l].metal_layers.geo[i].direc;
              if(MLCF_X == direc)
                rx += 1/model->layers[l].metal_layers.gridRL[i].r;
              else
                ry += 1/model->layers[l].metal_layers.gridRL[i].r;
          }
          rx = 1/rx; ry = 1/ry;

          cur_sq_sum = 0;
          /* on-chip power */
          for(i=0; i<nr; i++)
            for(j=0; j<nc-1; j++){
                grid_idx = (2*l)*nr*nc + i*nc + j;
                branch_cur = (v[grid_idx] - v[grid_idx+1])/rx;
                cur_sq_sum += branch_cur*branch_cur;
                grid_idx = (2*l+1)*nr*nc + i*nc + j;
                branch_cur = (v[grid_idx] - v[grid_idx+1])/rx;
                cur_sq_sum += branch_cur*branch_cur;
            }
          onchip_power += cur_sq_sum * rx;

          cur_sq_sum = 0;
          for(i=0; i<nr-1; i++)
            for(j=0; j<nc; j++){
                grid_idx = (2*l)*nr*nc + i*nc + j;
                branch_cur = fabs((v[grid_idx] - v[grid_idx+nc])/ry);
                cur_sq_sum += branch_cur*branch_cur;
                grid_idx = (2*l+1)*nr*nc + i*nc + j;
                branch_cur = fabs((v[grid_idx] - v[grid_idx+nc])/ry);
                cur_sq_sum += branch_cur*branch_cur;
            }
          onchip_power += cur_sq_sum * ry;
      }

      cur_sq_sum = 0;
      tsv_power = 0;
      /* TSV power */
      for(l=0; l<nl-1; l++){
          rv = model->layers[l].tsv.r;

          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                if((GNDTSV == model->layers[l].tsv.loc[i][j])){
                    a_grid_idx = l*2*nr*nc + nr*nc + i*nc + j;
                    b_grid_idx = a_grid_idx + nr*nc;
                    tsv_cur = (v[b_grid_idx] - v[a_grid_idx]) / rv;
                    cur_sq_sum += tsv_cur*tsv_cur;
                }
            }
          tsv_power += cur_sq_sum * rv; 
      }
      tsv_power += vddpad_cur_sq_sum * (r_tsv_pad - r_pad); 
  }

  return pad_power + pkg_power + onchip_power + tsv_power;
}

void get_maxIR(model_t *model, double *maxIR)
{
  int i, j, l;
  int v_grid_idx, g_grid_idx;
  int cur_id;
  double temp;
  double min;
  blist_t *ptr;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int vs = model->config.v_stacking;
  double *v = model->last_steady;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  if(!vs){
      for(l = 0; l < nl; l++){
          min = LARGENUM;
          for(i = 0; i < nr; i++)
            for(j = 0; j < nc; j++){
                v_grid_idx = l*nr*nc + i*nc + j;
                g_grid_idx = nl*nr*nc + v_grid_idx;
                temp = v[v_grid_idx] - v[g_grid_idx];
                if (temp < min){
                    min = temp;
                }
            }
          maxIR[l] = 100 * (1-min/(vdd-gnd));
      }
  }
  else{
      for(l = 0; l < nl; l++){
          min = LARGENUM;
          for(i = 0; i < nr; i++)
            for(j = 0; j < nc; j++){
                g_grid_idx = l*2*nr*nc + i*nc + j;
                v_grid_idx = g_grid_idx + nr*nc;
                temp = v[v_grid_idx] - v[g_grid_idx];
                if (temp < min){
                    min = temp;
                }
            }
          maxIR[l] = 100 * (1-min/(vdd-gnd));
      }
  }

  return ;
}

/* Analyze onchip voltage droop, return adjusted max percent */
void onchipV_analyze_2D(model_t *model, status_t *status)
{
  int i, j;
  double temp, vnoise;
  double max = -LARGENUM;
  int blk_idx, core_id;
  int cur;
  int print_nodetrace = 0;
  blist_t *ptr;
  FILE *fp_t, *fp_v;
  char str[STR_SIZE];

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nml = model->layers[0].metal_layers.n_metal/2;
  int nbr = 2*nr*nc - nr - nc;
  int **sloc = status->sensor_loc;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  if (strcmp(model->config.node_viotrace_file, NULLFILE)){
      print_nodetrace = 1;
      fp_t = fopen (model->config.node_viotrace_file, "a");
      if (!fp_t) {
          sprintf (str,"error: %s could not be opened for writing\n", model->config.node_viotrace_file);
          fatal(str);
      }
  }

  status->vio_area_ratio = 0;
  negL_dvector(status->corestats.max_1D, MAX_CORE_NUM);
  for(i = 0; i < nr; i++)
    for(j = 0; j < nc; j++){
        cur = i*nc + j;
        ptr = model->layers[0].b2gmap[i][j];
        blk_idx = ptr->idx;
        core_id = get_core_id(model->layers[0].flp->units[blk_idx].name);

        if(!model->is_3D){
            if(model->config.PDN_gridL)
              temp = status->cycle_avg[nvp+cur] - status->cycle_avg[nvp+nr*nc+nbr*nml+ngp+cur];
            else
              temp = status->cycle_avg[nvp+cur] - status->cycle_avg[nvp+nr*nc+ngp+cur];
        }
        else{
            temp = status->cycle_avg[nvp+ngp+cur] - status->cycle_avg[nvp+ngp+nr*nc+cur];
        }
        vnoise = 100*(vdd-gnd-temp)/(vdd-gnd);

        if (vnoise > max){ max = vnoise; }

        status->gridstats.integral_2D[i][j] += vnoise;
        /* per cycle noise amplitude analysis */
        if(vnoise > status->gridstats.max_2D[i][j]){
            status->gridstats.max_2D[i][j] = vnoise;
        }
        if(vnoise > status->blkstats.max_1D[blk_idx]){
            status->blkstats.max_1D[blk_idx] = vnoise;
        }
        // per-core max droop of this analysis run
        if((core_id >= 0) && (vnoise > status->corestats.max_1D[core_id]))
          status->corestats.max_1D[core_id] = vnoise;

        /* per cycle violation analysis */
        if(vnoise > (model->config.PDN_noise_th)){
            /* grid violation counter */
            status->gridstats.counter_2D[i][j]++;
            /* per block violation counter */
            if(status->blkstats.counter_1D[blk_idx] >= 0){
                // at most increment once in a cycle
                status->blkstats.counter_1D[blk_idx]++;
                status->blkstats.counter_1D[blk_idx] *= -1;
            }
            /* get vio grid # */
            status->vio_area_ratio += 1;
        }

        /* Update sensor record */
        if (SENSOR == sloc[i][j]){
            status->sensor_noise[i][j] = vnoise;
        }

        /* per cycle grid violation trace */
        if(print_nodetrace){
            if(vnoise > (model->config.PDN_noise_th)){ fprintf(fp_t, "1"); }
            else{ fprintf(fp_t, "0"); }
        }
    }

  /* convert vio area to ratio */
  status->vio_area_ratio /= nr*nc;

  abs_ivector(status->blkstats.counter_1D, model->total_n_blocks);

  if(print_nodetrace){
      fprintf(fp_t, "\n");
      fclose(fp_t);	
  }

  //print grid voltage trace
  if (strcmp(model->config.gridvol_file, NULLFILE)){
      fp_v = fopen (model->config.gridvol_file, "a");
      if (!fp_v) {
          sprintf (str,"error: %s could not be opened for writing\n", model->config.gridvol_file);
          fatal(str);
      }
      for(i = 0; i < nr; i++)
        for(j = 0; j < nc; j++){
            cur = i*nc + j;
            fprintf(fp_v, "%.7f ", status->cycle_avg[nvp+cur]);
        }
      fprintf(fp_v, "\n");
      fclose(fp_v);	
  }

  status->maxdroop = max;

  //transient grid voltage gradient analysis
  if(model->config.vgradient_analyse){
      trans_vgradient_analyze(model, status);
  }
}

void onchipV_analyze_3D(model_t *model, status_t *status)
{
  int i, j, l;
  int vdd_idx, gnd_idx;
  double rail_v, vnoise;
  double max, l_max;
  int cur;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  max = -LARGENUM;
  for(l=0; l<nl; l++){
      l_max = -LARGENUM;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            cur = i*nc + j;
            vdd_idx = nvp + ngp + l*2*nr*nc + cur; 
            gnd_idx = nvp + ngp + l*2*nr*nc + nr*nc + cur; 

            rail_v = status->cycle_avg[vdd_idx] - status->cycle_avg[gnd_idx];
            vnoise = 100*(vdd-gnd-rail_v)/(vdd-gnd);

            if (vnoise > l_max){ l_max = vnoise; }
            if (vnoise > max){ max = vnoise; }
        }
      status->layer_maxdroop[l] = l_max;
  }

  status->maxdroop = max;
}

void perstep_droop_3D(model_t *model, FILE *fp)
{
  int i, j, l;
  int vdd_idx, gnd_idx;
  double rail_v, vnoise;
  double max, l_max;
  int cur;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  max = -LARGENUM;
  for(l=0; l<nl; l++){
      l_max = -LARGENUM;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            cur = i*nc + j;
            vdd_idx = nvp + ngp + l*2*nr*nc + cur; 
            gnd_idx = nvp + ngp + l*2*nr*nc + nr*nc + cur; 

            rail_v = model->last_trans[vdd_idx] - model->last_trans[gnd_idx];
            vnoise = 100*(vdd-gnd-rail_v)/(vdd-gnd);

            if (vnoise > l_max){ l_max = vnoise; }
            if (vnoise > max){ max = vnoise; }
        }
      fprintf(fp, "%.4lf\t", l_max);
  }
  fprintf(fp, "%.4lf\n", max);
}

//for trans WP
void trans_vgradient_analyze(model_t *model, status_t *status)
{
  int i,j;
  double v_grad;
  double *rst_v, *rst_g;
  int brc_idx, nod_idx;

  /* shortcuts */
  int nr = model->rows;
  int nc = model->cols;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nml = model->layers[0].metal_layers.n_metal/2;
  int nbr = 2*nr*nc - nr - nc;

  // max
  //////////////////////////////////////////////////
  //VDD
  rst_v = &status->cycle_avg[nvp];
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_v[nod_idx + 1]);
          if(v_grad > status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_v[nod_idx + nc]);
          if(v_grad > status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }

  //GND
  if(model->config.PDN_gridL)
    rst_g = &status->cycle_avg[nvp+ngp+nr*nc+nml*nbr];
  else
    rst_g = &status->cycle_avg[nvp+ngp+nr*nc];

  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr + i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_g[nod_idx] - rst_g[nod_idx + 1]);
          if(v_grad > status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr + (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_g[nod_idx] - rst_g[nod_idx + nc]);
          if(v_grad > status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }

  // VDD-GND
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr*2 + i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_g[nod_idx]) - (rst_v[nod_idx + 1]-rst_g[nod_idx + 1]);
          if(v_grad > status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr*2 + (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_g[nod_idx]) - (rst_v[nod_idx + nc] - rst_g[nod_idx + nc]);
          if(v_grad > status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }


  // min
  //////////////////////////////////////////////////
  //VDD
  rst_v = &status->cycle_avg[nvp];
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr*3 + i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_v[nod_idx + 1]);
          if(v_grad < status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr*3 + (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_v[nod_idx + nc]);
          if(v_grad < status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }

  //GND
  if(model->config.PDN_gridL)
    rst_g = &status->cycle_avg[nvp+ngp+nr*nc+nml*nbr];
  else
    rst_g = &status->cycle_avg[nvp+ngp+nr*nc];

  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr*4 + i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_g[nod_idx] - rst_g[nod_idx + 1]);
          if(v_grad < status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr*4 + (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_g[nod_idx] - rst_g[nod_idx + nc]);
          if(v_grad < status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }

  // VDD-GND
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr*5 + i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_g[nod_idx]) - (rst_v[nod_idx + 1]-rst_g[nod_idx + 1]);
          if(v_grad < status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr*5 + (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_g[nod_idx]) - (rst_v[nod_idx + nc] - rst_g[nod_idx + nc]);
          if(v_grad < status->vgradient.max_1D[brc_idx])
            status->vgradient.max_1D[brc_idx] = v_grad;
      }
  }

  // Root mean square
  //////////////////////////////////////////////////
  //VDD
  rst_v = &status->cycle_avg[nvp];
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr*6 + i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_v[nod_idx + 1]);
          if(stat_counter == 0)
            status->vgradient.max_1D[brc_idx] = v_grad*v_grad;
          else
            status->vgradient.max_1D[brc_idx] += v_grad*v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr*6 + (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_v[nod_idx] - rst_v[nod_idx + nc]);
          if(stat_counter == 0)
            status->vgradient.max_1D[brc_idx] = v_grad*v_grad;
          else
            status->vgradient.max_1D[brc_idx] += v_grad*v_grad;
      }
  }

  //GND
  if(model->config.PDN_gridL)
    rst_g = &status->cycle_avg[nvp+ngp+nr*nc+nml*nbr];
  else
    rst_g = &status->cycle_avg[nvp+ngp+nr*nc];

  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr*7 + i*(nc-1) + j;
          nod_idx = i*nc + j;
          v_grad = (rst_g[nod_idx] - rst_g[nod_idx + 1]);
          if(stat_counter == 0)
            status->vgradient.max_1D[brc_idx] = v_grad*v_grad;
          else
            status->vgradient.max_1D[brc_idx] += v_grad*v_grad;
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr*7 + (nc-1)*nr + i*nc + j;
          nod_idx = i*nc + j;
          v_grad = (rst_g[nod_idx] - rst_g[nod_idx + nc]);
          if(stat_counter == 0)
            status->vgradient.max_1D[brc_idx] = v_grad*v_grad;
          else
            status->vgradient.max_1D[brc_idx] += v_grad*v_grad;
      }
  }

  stat_counter++;
}

//void trans_vgradient_analyze(model_t *model, status_t *status)
//{
//  int i,j;
//  double v_grad;
//  double *rst;
//  int brc_idx, nod_idx;
//
//  /* shortcuts */
//  int nr = model->rows;
//  int nc = model->cols;
//  int nvp = model->c4->vdd_num;
//  int ngp = model->c4->gnd_num;
//  int nml = model->layers[0].metal_layers.n_metal/2;
//  int nbr = 2*nr*nc - nr - nc;
//
//  //VDD
//  rst = &status->cycle_avg[nvp];
//  for(i=0; i<nr; i++){
//      for(j=0; j<nc-1; j++){
//          brc_idx = i*(nc-1) + j;
//          nod_idx = i*nc + j;
//          v_grad = (rst[nod_idx] - rst[nod_idx + 1]);
//          if(fabs(v_grad) > fabs(status->vgradient.max_1D[brc_idx]))
//            status->vgradient.max_1D[brc_idx] = v_grad;
//      }
//  }
//  for(i=0; i<nr-1; i++){
//      for(j=0; j<nc; j++){
//          nod_idx = i*nc + j;
//          brc_idx = (nc-1)*nr + i*nc + j;
//          v_grad = (rst[nod_idx] - rst[nod_idx + nc]);
//          if(fabs(v_grad) > fabs(status->vgradient.max_1D[brc_idx]))
//            status->vgradient.max_1D[brc_idx] = v_grad;
//      }
//  }
//
//  //GND
//  if(model->config.PDN_gridL)
//    rst = &status->cycle_avg[nvp+ngp+nr*nc+nml*nbr];
//  else
//    rst = &status->cycle_avg[nvp+ngp+nr*nc];
//
//  for(i=0; i<nr; i++){
//      for(j=0; j<nc-1; j++){
//          brc_idx = nbr + i*(nc-1) + j;
//          nod_idx = i*nc + j;
//          v_grad = (rst[nod_idx] - rst[nod_idx + 1]);
//          if(fabs(v_grad) > fabs(status->vgradient.max_1D[brc_idx]))
//            status->vgradient.max_1D[brc_idx] = v_grad;
//      }
//  }
//  for(i=0; i<nr-1; i++){
//      for(j=0; j<nc; j++){
//          nod_idx = i*nc + j;
//          brc_idx = nbr + (nc-1)*nr + i*nc + j;
//          v_grad = (rst[nod_idx] - rst[nod_idx + nc]);
//          if(fabs(v_grad) > fabs(status->vgradient.max_1D[brc_idx]))
//            status->vgradient.max_1D[brc_idx] = v_grad;
//      }
//  }
//}

void dump_files(model_t *model, status_t *status, int mode)
{
  if (strcmp(model->config.padloc_file_out, NULLFILE)){
      dump_anypadloc(model, model->config.padloc_file_out, PGPAD);
  }

  if(STEADY == mode){
      if (strcmp(model->config.gridvol_file, NULLFILE)){
          dump_grid_drop(model, model->config.gridvol_file);
      }
      if (strcmp(model->config.padcur_file, NULLFILE)){
          dump_cur_with_cordt(model, status, model->config.padcur_file);
      }
      if ((model->is_3D) && (strcmp(model->config.tsvcur_file, NULLFILE))){
          dump_tsv_cur(model, status, model->config.tsvcur_file);
      }
  }

  if(TRANSIENT == mode){
      if (strcmp(model->config.vio_file, NULLFILE)){
          dump_violation(model, status, model->config.vio_file);
      }
      if(model->config.vgradient_analyse)
        if (strcmp(model->config.trans_vgradient_file, NULLFILE)){
            dump_trans_vgradient(model, status, model->config.trans_vgradient_file);
        }
  }
}

void dump_PDN_config(model_t *model, FILE *fp)
{
  int i;

  fprintf(fp, "PDN_grid_intv=%d\n", model->config.PDN_grid_intv);
  for(i=0; i<model->layers[0].metal_layers.n_metal; i++){
      fprintf(fp, "Metal Layer %d: p%e\tw%e\tt%e\n", 
              i, model->layers[0].metal_layers.geo[i].pitch,
              model->layers[0].metal_layers.geo[i].width, 
              model->layers[0].metal_layers.geo[i].thick);
      fprintf(fp, "R:%lf\tL:%e\n", 
              model->layers[0].metal_layers.gridRL[i].r, 
              model->layers[0].metal_layers.gridRL[i].l);
  }

  fprintf(fp, "PDN_padpitch=%e\n", model->config.PDN_padpitch);
  fprintf(fp, "PDN_padD=%e\n", model->config.PDN_padD);
  fprintf(fp, "PDN_padR=%e\n", model->config.PDN_padR);
  fprintf(fp, "PDN_cur_dense=%e\n", model->config.PDN_cur_dense);

  fprintf(fp, "PDN_padconfig=%d\n", model->config.PDN_padconfig);
  fprintf(fp, "padloc_format=%d\n", model->config.padloc_format);

  fprintf(fp, "PDN_decap_dense=%e\n", model->config.PDN_decap_dense);
  fprintf(fp, "PDN_decap_ratio=%e\n", model->config.PDN_decap_ratio);
  fprintf(fp, "PDN_decap_unifm=%d\n", model->config.PDN_decap_unifm);
  fprintf(fp, "PDN_padL=%e\n", model->config.PDN_padL);
  fprintf(fp, "PDN_pkg_sL=%e\n", model->config.PDN_pkg_sL);
  fprintf(fp, "PDN_pkg_sR=%e\n", model->config.PDN_pkg_sR);
  fprintf(fp, "PDN_pkg_C=%e\n", model->config.PDN_pkg_C);
  fprintf(fp, "PDN_pkg_pR=%e\n", model->config.PDN_pkg_pR);
  fprintf(fp, "PDN_pkg_pL=%e\n", model->config.PDN_pkg_pL);
  fprintf(fp, "PDN_pkg_scale=%e\n", model->config.PDN_pkg_scale);

  fprintf(fp, "vdd=%e\n", model->config.vdd);
  fprintf(fp, "gnd=%e\n", model->config.gnd);

  fprintf(fp, "proc_clock_freq=%e\n", model->config.proc_clock_freq);
  fprintf(fp, "PDN_step_percycle=%d\n", model->config.PDN_step_percycle);
  fprintf(fp, "ptrace_sampling_intvl=%d\n", model->config.ptrace_sampling_intvl);
  fprintf(fp, "PDN_ptrace_warmup=%d\n", model->config.PDN_ptrace_warmup);
  fprintf(fp, "PDN_multi_dom=%d\n", model->config.PDN_multi_dom);
  //fprintf(fp, "animation=%d\n", model->config.animation);
  //fprintf(fp, "frame_intv=%d\n", model->config.frame_intv);
  //fprintf(fp, "legend_lwr=%lf\n", model->config.legend_lwr);
  //fprintf(fp, "legend_upr=%lf\n", model->config.legend_upr);
  //fprintf(fp, "legend_curupr=%lf\n", model->config.legend_curupr);

  //fprintf(fp, "reserve_io=%d\n", model->config.reserve_io);
  //fprintf(fp, "MC_pads=%d\n", model->config.MC_pads);
  //fprintf(fp, "IO_dense=%lf\n", model->config.IO_dense);

  fprintf(fp, "v_stacking=%d\n", model->config.v_stacking);
  fprintf(fp, "TSV_config=%d\n", model->config.TSV_config);
  fprintf(fp, "TSV_R=%lf\n", model->config.TSV_R);
  fprintf(fp, "TSV_L=%lf\n", model->config.TSV_L);
  fprintf(fp, "SC_freq=%lf\n", model->config.SC_freq);
  fprintf(fp, "SC_totcap=%lf\n", model->config.SC_totcap);
  fprintf(fp, "SC_Rontop=%lf\n", model->config.SC_Rontop);
  fprintf(fp, "SC_Ronbtm=%lf\n", model->config.SC_Ronbtm);

  fprintf(fp, "run_PDN=%d\n", model->config.run_PDN);
  fprintf(fp, "PDN_pkgLC=%d\n", model->config.PDN_pkgLC);
  fprintf(fp, "PDN_gridL=%d\n", model->config.PDN_gridL);
  fprintf(fp, "PDN_sin_pattern=%d\n", model->config.PDN_sin_pattern);
  fprintf(fp, "PDN_sin_totstep=%d\n", model->config.PDN_sin_totstep);
  fprintf(fp, "PDN_sin_freq=%e\n", model->config.PDN_sin_freq);
  //fprintf(fp, "vgradient_analyse=%d\n", model->config.vgradient_analyse);

  fprintf(fp, "mlayer_spec_file=%s\n", model->config.mlayer_spec_file);
  fprintf(fp, "layer_file_3D=%s\n", model->config.layer_file_3D);
  if(0 == model->config.PDN_padconfig)
    fprintf(fp, "padloc_file_in=%s\n", model->config.padloc_file_in);
  fprintf(fp, "padloc_file_out=%s\n", model->config.padloc_file_out);
  fprintf(fp, "vio_file=%s\n", model->config.vio_file);
  fprintf(fp, "node_viotrace_file=%s\n", model->config.node_viotrace_file);
  fprintf(fp, "trans_vgradient_file=%s\n", model->config.trans_vgradient_file);
  fprintf(fp, "senloc_file=%s\n", model->config.senloc_file);
  fprintf(fp, "padcur_file=%s\n", model->config.padcur_file);
  fprintf(fp, "tsvcur_file=%s\n", model->config.tsvcur_file);
  fprintf(fp, "gridvol_file=%s\n", model->config.gridvol_file);

  fprintf(fp, "END_OF_CONFIGS\n");
}

/* dump all grid voltage to a data file */
void dump_grid_drop(model_t *model, char *file)
{
  int i, j, l;
  int v_grid_idx, g_grid_idx;
  char str[STR_SIZE];
  FILE *fp;
  double max_cur_dense;
  double drop_ratio;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int vs = model->config.v_stacking;
  double *v = model->last_steady;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  if (!strcasecmp(file, "stdout"))
    fp = stdout;
  else if (!strcasecmp(file, "stderr"))
    fp = stderr;
  else 	
    fp = fopen (file, "w");
  if (!fp) {
      sprintf (str,"error: %s could not be opened for writing\n", file);
      fatal(str);
  }

  for(l=0; l < nl; l++){
      fprintf(fp, "#Layer:%d\n", l);
      for(i=0; i < nr; i++){
          for(j=0; j < nc; j++){
              if(!vs){
                  v_grid_idx = l*nr*nc + i*nc + j;
                  g_grid_idx = nl*nr*nc + v_grid_idx;
              }
              else{
                  g_grid_idx = l*2*nr*nc + i*nc + j;
                  v_grid_idx = g_grid_idx + nr*nc;
              }
              drop_ratio = 100*(1-(v[v_grid_idx] - v[g_grid_idx])/(vdd-gnd));
              fprintf(fp, "%d\t%d\t%lf\n", j, nr-i-1, drop_ratio);
          }
      }
      fprintf(fp, "\n");
  }

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

void dump_violation(model_t *model, status_t *status, char *file)
{
  int i, j;
  char str[STR_SIZE];
  FILE *fp;
  char *name;

  /* shortcuts */
  int nr = model->rows;
  int nc = model->cols;
  int nb = model->layers[0].flp->n_units;

  if (!strcasecmp(file, "stdout"))
    fp = stdout;
  else if (!strcasecmp(file, "stderr"))
    fp = stderr;
  else 	
    fp = fopen (file, "w");
  if (!fp) {
      sprintf (str,"error: %s could not be opened for writing\n", file);
      fatal(str);
  }

  for(i=0; i < nr; i++)
    for(j=0; j < nc; j++){
        fprintf(fp, "%d\t%d\t%d\t%lf\t%lf", j, nr-i-1, 
                status->gridstats.counter_2D[i][j], 
                status->gridstats.integral_2D[i][j],
                status->gridstats.max_2D[i][j]);
        fprintf(fp, "\n");
    }

  fprintf(fp, "-_-b\n");
  for(i=0; i < nb; i++){
      name = model->layers[0].flp->units[i].name;
      fprintf(fp, "%s\t%d\t%lf\n", name, status->blkstats.counter_1D[i],
              status->blkstats.max_1D[i]);
  }

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

// for transient walking pad
void dump_trans_vgradient(model_t *model,  status_t *status, char *file)
{
  int i, j;
  char str[STR_SIZE];
  FILE *fp;
  FILE *fp_info;
  char *name;
  int brc_idx;

  /* shortcuts */
  int nr = model->rows;
  int nc = model->cols;
  int nbr = 2*nr*nc - nr - nc;

  /* Output important information */
  fp_info = fopen ("trans.info", "w");
  if (!fp_info) 
    {
      printf("Error in opening trans.info\n");
      /* sprintf(str,"error: %s could not be opened for writing\n", "trans.info"); */
      /* fatal(str); */
    }
  fprintf(fp_info, "NC %d\n", nc);
  fprintf(fp_info, "NR %d\n", nr);
  fprintf(fp_info, "GRID_INTV %d\n", model->config.PDN_grid_intv);
  /* fprintf(fp_info, "RX %d\n", nc); */
  /* fprintf(fp_info, "RY %d\n", nc); */
  fclose(fp_info);	

  if (!strcasecmp(file, "stdout"))
    fp = stdout;
  else if (!strcasecmp(file, "stderr"))
    fp = stderr;
  else 	
    fp = fopen (file, "w");
  if (!fp) {
      sprintf (str,"error: %s could not be opened for writing\n", file);
      fatal(str);
  }

  if(stat_counter == 0) stat_counter = -1;

  //VDD
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = i*(nc-1) + j;
          fprintf(fp, "VV\t%d_%d\t%d_%d\t%e\t%e\t%e\n", j, nr-i, j+1, nr-i, status->vgradient.max_1D[brc_idx], status->vgradient.max_1D[brc_idx+nbr*3], sqrt(status->vgradient.max_1D[brc_idx+nbr*6]/(double)(stat_counter)));
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = (nc-1)*nr + i*nc + j;
          fprintf(fp, "VV\t%d_%d\t%d_%d\t%e\t%e\t%e\n", j, nr-i, j, nr-(i+1), status->vgradient.max_1D[brc_idx], status->vgradient.max_1D[brc_idx+nbr*3], sqrt(status->vgradient.max_1D[brc_idx+nbr*6]/(double)(stat_counter)));
      }
  }

  //GND
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr + i*(nc-1) + j;
          fprintf(fp, "GG\t%d_%d\t%d_%d\t%e\t%e\t%e\n", j, nr-i, j+1, nr-i, status->vgradient.max_1D[brc_idx], status->vgradient.max_1D[brc_idx+nbr*3], sqrt(status->vgradient.max_1D[brc_idx+nbr*6]/(double)(stat_counter)));
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr + (nc-1)*nr + i*nc + j;
          fprintf(fp, "GG\t%d_%d\t%d_%d\t%e\t%e\t%e\n", j, nr-i, j, nr-(i+1), status->vgradient.max_1D[brc_idx], status->vgradient.max_1D[brc_idx+nbr*3], sqrt(status->vgradient.max_1D[brc_idx+nbr*6]/(double)(stat_counter)));
      }
  }

  //VDD-GND
  for(i=0; i<nr; i++){
      for(j=0; j<nc-1; j++){
          brc_idx = nbr*2 + i*(nc-1) + j;
          fprintf(fp, "VG\t%d_%d\t%d_%d\t%e\t%e\n", j, nr-i, j+1, nr-i, status->vgradient.max_1D[brc_idx], status->vgradient.max_1D[brc_idx+nbr*3]);
      }
  }
  for(i=0; i<nr-1; i++){
      for(j=0; j<nc; j++){
          brc_idx = nbr*2 + (nc-1)*nr + i*nc + j;
          fprintf(fp, "VG\t%d_%d\t%d_%d\t%e\t%e\n", j, nr-i, j, nr-(i+1), status->vgradient.max_1D[brc_idx], status->vgradient.max_1D[brc_idx+nbr*3]);
      }
  }

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

//void dump_trans_vgradient(model_t *model,  status_t *status, char *file)
//{
//  int i, j;
//  char str[STR_SIZE];
//  FILE *fp;
//  int brc_idx;
//
//  /* shortcuts */
//  int nr = model->rows;
//  int nc = model->cols;
//  int nbr = 2*nr*nc - nr - nc;
//
//  if (!strcasecmp(file, "stdout"))
//    fp = stdout;
//  else if (!strcasecmp(file, "stderr"))
//    fp = stderr;
//  else 	
//    fp = fopen (file, "w");
//  if (!fp) {
//      sprintf (str,"error: %s could not be opened for writing\n", file);
//      fatal(str);
//  }
//
//  //VDD
//  for(i=0; i<nr; i++){
//      for(j=0; j<nc-1; j++){
//          brc_idx = i*(nc-1) + j;
//          fprintf(fp, "V\t%d_%d\t%d_%d\t%e\n", i, j, i, j+1, status->vgradient.max_1D[brc_idx]);
//      }
//  }
//  for(i=0; i<nr-1; i++){
//      for(j=0; j<nc; j++){
//          brc_idx = (nc-1)*nr + i*nc + j;
//          fprintf(fp, "V\t%d_%d\t%d_%d\t%e\n", i, j, i+1, j, status->vgradient.max_1D[brc_idx]);
//      }
//  }
//
//  //GND
//  for(i=0; i<nr; i++){
//      for(j=0; j<nc-1; j++){
//          brc_idx = nbr + i*(nc-1) + j;
//          fprintf(fp, "G\t%d_%d\t%d_%d\t%e\n", i, j, i, j+1, status->vgradient.max_1D[brc_idx]);
//      }
//  }
//  for(i=0; i<nr-1; i++){
//      for(j=0; j<nc; j++){
//          brc_idx = nbr + (nc-1)*nr + i*nc + j;
//          fprintf(fp, "G\t%d_%d\t%d_%d\t%e\n", i, j, i+1, j, status->vgradient.max_1D[brc_idx]);
//      }
//  }
//
//  if(fp != stdout && fp != stderr)
//    fclose(fp);	
//}

void draw_single_gif(model_t *model, status_t *status, int draw_counter, int mode)
{
  int i, j, l, dum;
  int cur_idx;
  char str[STR_SIZE];
  char posfix[STR_SIZE];
  char file1[STR_SIZE], file2[STR_SIZE];
  FILE *fp1, *fp2;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nvp = model->c4->vdd_num;
  double vdd = model->config.vdd;

  double upr_ratio = model->config.legend_upr;
  double lwr_ratio = model->config.legend_lwr;
  double upr_cur   = model->config.legend_curupr;

  fp1 = fopen ("temp1.dat", "w");
  if (!fp1) {
      sprintf(str,"error: %s could not be opened for writing\n", "temp1.dat");
      fatal(str);
  }
  if(3 == model->config.animation){
      fp2 = fopen ("temp2.dat", "w");
      if (!fp2) {
          sprintf(str,"error: %s could not be opened for writing\n", "temp2.dat");
          fatal(str);
      }
  }

  if(TRANSIENT == mode){
      if(3 == model->config.animation){
          fatal("Do not support current plot for transient yet!\n");
      }

      if(!model->config.PDN_pkgLC)
        fatal("Do not support plotting without pkg yet!\n");

      for(i = 0; i < nr; i++)
        for(j = 0; j < nc; j++){
            fprintf(fp1, "%d\t%d\t%lf\n",j, nr-i-1, 
                    model->last_trans[nvp+i*nc+j]);
        }
  }
  else if(STEADY == mode){
      if(3 == model->config.animation){
          fatal("No longer support steady-state on-chip current plot");
      }
      else{
          l = 0;
          for(i = 0; i < nr; i++)
            for(j = 0; j < nc; j++){
                cur_idx = l*nr*nc + i*nc + j;
                fprintf(fp1, "%d\t%d\t%lf\n",j, nr-i-1, model->last_steady[cur_idx]);
            }
      }
  }
  else
    fatal("Function Call error in single_draw_gif");

  fclose(fp1);	
  if(3 == model->config.animation)
    fclose(fp2);	

  fp1 = fopen ("temp1.gpi", "w");
  if (!fp1) {
      sprintf(str,"error: %s could not be opened for writing\n", "temp1.gpi");
      fatal(str);
  }
  if(3 == model->config.animation){
      fp2 = fopen ("temp2.gpi", "w");
      if (!fp2) {
          sprintf(str,"error: %s could not be opened for writing\n", "temp2.gpi");
          fatal(str);
      }
  }

  int_to_str(draw_counter, 8, posfix);
  sprintf(file1,"frame_%s.1.gif", posfix);
  sprintf(file2,"frame_%s.2.gif", posfix);

  if(1 == model->config.animation){
      fprintf(fp1, "set style rectangle back fc lt -3 fillstyle  solid 1.00 border -1\n");
      fprintf(fp1, "unset key\n");
      fprintf(fp1, "set view map\n");
      fprintf(fp1, "set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n");
      fprintf(fp1, "set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n");
      fprintf(fp1, "set xrange [ 0 : %d ];\n", nc);
      fprintf(fp1, "set yrange [ 0 : %d ];\n", nr);
      fprintf(fp1, "set ylabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n");
      fprintf(fp1, "set cblabel \"Voltage\" \n");
      fprintf(fp1, "set cblabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n");
      fprintf(fp1, "set cbrange [ %lf : %lf ] noreverse nowriteback\n", vdd*(1-lwr_ratio), vdd*(1+upr_ratio));
      fprintf(fp1, "set palette rgbformulae 21, 22, 23 \n");
      fprintf(fp1, "set term gif\n");
      fprintf(fp1, "set output \"%s\"\n", file1);
      fprintf(fp1, "plot \"temp1.dat\" using 1:2:3 with image\n");
  }
  else if(2 == model->config.animation){
      fprintf(fp1, "set view 60, 30, 0.85\n");
      fprintf(fp1, "set contour base\n");
      fprintf(fp1, "set dgrid3d %d %d\n", nc, nr);
      fprintf(fp1, "set hidden3d\n");
      fprintf(fp1, "set cntrparam levels auto 10\n");
      fprintf(fp1, "unset key\n");
      fprintf(fp1, "unset xtics\n");
      fprintf(fp1, "unset ytics\n");
      fprintf(fp1, "unset ztics\n");
      fprintf(fp1, "set xrange [ 0 : %d ];\n", nc);
      fprintf(fp1, "set yrange [ 0 : %d ];\n", nr);
      fprintf(fp1, "set zrange [ %lf : %lf ];\n", vdd*(1-lwr_ratio), vdd*(1+upr_ratio));
      fprintf(fp1, "set term gif\n");
      fprintf(fp1, "set output \"%s\"\n", file1);
      fprintf(fp1, "splot \"temp1.dat\" using 1:2:3 with lines\n");
  }
  else if(3 == model->config.animation){
      fprintf(fp1, "set style rectangle back fc lt -3 fillstyle  solid 1.00 border -1\n");
      fprintf(fp1, "unset key\n");
      fprintf(fp1, "set view map\n");
      fprintf(fp1, "set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n");
      fprintf(fp1, "set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n");
      fprintf(fp1, "set xrange [ 0 : %d ];\n", nc);
      fprintf(fp1, "set yrange [ 0 : %d ];\n", nr);
      fprintf(fp1, "set ylabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n");
      fprintf(fp1, "set cblabel \"Current\" \n");
      fprintf(fp1, "set cblabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n");
      fprintf(fp1, "set cbrange [ 0 : %lf ] noreverse nowriteback\n", upr_cur);
      fprintf(fp1, "set palette rgbformulae 21, 22, 23 \n");
      fprintf(fp1, "set term gif\n");
      fprintf(fp1, "set output \"%s\"\n", file1);
      fprintf(fp1, "plot \"temp1.dat\" using 1:2:3 with image\n");

      fprintf(fp2, "set style rectangle back fc lt -3 fillstyle  solid 1.00 border -1\n");
      fprintf(fp2, "unset key\n");
      fprintf(fp2, "set view map\n");
      fprintf(fp2, "set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n");
      fprintf(fp2, "set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n");
      fprintf(fp2, "set xrange [ 0 : %d ];\n", nc);
      fprintf(fp2, "set yrange [ 0 : %d ];\n", nr);
      fprintf(fp2, "set ylabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n");
      fprintf(fp2, "set cblabel \"Current\" \n");
      fprintf(fp2, "set cblabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n");
      fprintf(fp2, "set cbrange [ 0 : %lf ] noreverse nowriteback\n", upr_cur);
      fprintf(fp2, "set palette rgbformulae 21, 22, 23 \n");
      fprintf(fp2, "set term gif\n");
      fprintf(fp2, "set output \"%s\"\n", file2);
      fprintf(fp2, "plot \"temp2.dat\" using 1:2:3 with image\n");
  }

  fclose(fp1);	
  if(3 == model->config.animation)
    fclose(fp2);	

  dum = system("gnuplot temp1.gpi");
  dum = remove("temp1.dat");
  dum = remove("temp1.gpi");
  if(3 == model->config.animation){
      dum = system("gnuplot temp2.gpi");
      dum = remove("temp2.dat");
      dum = remove("temp2.gpi");
  }
}

void create_animation(int draw_counter, int animation)
{
  int i, dum;
  int delay;
  char file1[STR_SIZE], file2[STR_SIZE];
  char cmd[STR_SIZE], posfix[STR_SIZE];

  delay = 10;

  sprintf(cmd,"gifsicle --delay=%d --loop frame_*.1.gif > anim.1.gif", delay);
  dum = system(cmd);
  for(i=0; i<draw_counter; i++){
      int_to_str(i, 8, posfix);
      sprintf(file1,"frame_%s.1.gif", posfix);
      dum = remove(file1);
  }

  if(3 == animation){
      sprintf(cmd,"gifsicle --delay=%d --loop frame_*.2.gif > anim.2.gif", delay);
      dum = system(cmd);
      for(i=0; i<draw_counter; i++){
          int_to_str(i, 8, posfix);
          sprintf(file2,"frame_%s.2.gif", posfix);
          dum = remove(file2);
      }
  }
}

void dump_tsv_cur(model_t *model, status_t *status, char *file)
{
  int l, i, j;
  char str[STR_SIZE];
  FILE *fp;

  /* shortcuts */
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int vs = model->config.v_stacking;
  int **vloc = model->c4->vdd_loc;
  int ntsv = 0;
  for(l=0; l<nl-1; l++){
      ntsv += model->layers[l].tsv.num_gnd;
      ntsv += model->layers[l].tsv.num_vdd;
  }

  if (!strcasecmp(file, "stdout"))
    fp = stdout;
  else if (!strcasecmp(file, "stderr"))
    fp = stderr;
  else 	
    fp = fopen (file, "w");
  if (!fp) {
      sprintf (str,"error: %s could not be opened for writing\n", file);
      fatal(str);
  }

  for(i=0; i<ntsv; i++)
    fprintf(fp, "%e\n", status->TSVcur[i]);

  if(vs){
      fprintf(fp, "\n");
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++)
          for(l=0; l<nl-1; l++)
            if (vloc[i][j] & PGPAD){
                fprintf(fp, "%e\n", status->vdd_cur[i][j]);
            }
  }

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}
