/* 
 * This is a trace-level thermal simulator. It reads power values 
 * from an input trace file and outputs the corresponding instantaneous 
 * voltage values to an output trace file. It also outputs the steady 
 * state voltage/current values to designated files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "PDN_sim.h"
#include "PDN_analyze.h"
#include "flp.h"
#include "util.h"
#include "voltspot.h"

/* Lib for SuperLU */
#include "slu_ddefs.h"

void usage(int argc, char **argv)
{
  fprintf(stdout, "Usage: %s -f <file> -p <file> [-c <file>] [options]\n", argv[0]);
  fprintf(stdout, "Options:(may be specified in any order, within \"[]\" means optional)\n");
  fprintf(stdout, "   -f <file>\tfloorplan input file (e.g. example.flp) - overridden by the\n");
  fprintf(stdout, "            \tlayer configuration file (e.g. 3D.lcf) when the\n");
  fprintf(stdout, "            \tlatter is specified\n");
  fprintf(stdout, "   -p <file>\tpower trace input file (e.g. example.ptrace)\n");
  fprintf(stdout, "  [-c <file>]\tinput configuration parameters from file (e.g. pdn.config)\n");
  fprintf(stdout, "  [-v <file>]\ttransient PDN output file - will skip transient simulation\n");
  fprintf(stdout, "            \tif not provided\n");
  fprintf(stdout, "  [options]\tzero or more options of the form \"-<name> <value>\",\n");
  fprintf(stdout, "           \toverride the options from config file.\n");
}

/* 
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void global_config_from_strs(global_config_t *config, str_pair *table, int size)
{
  int idx;
  if((idx = get_str_index(table, size, "f")) >= 0) {
      if(sscanf(table[idx].value, "%s", config->flp_file) != 1)
        fatal("invalid format for configuration  parameter flp_file\n");
  } else {
      fatal("required parameter flp_file missing. check usage\n");
  }
  if((idx = get_str_index(table, size, "p")) >= 0) {
      if(sscanf(table[idx].value, "%s", config->p_infile) != 1)
        fatal("invalid format for configuration  parameter p_infile\n");
  } else {
      fatal("required parameter p_infile missing. check usage\n");
  }
  if((idx = get_str_index(table, size, "c")) >= 0) {
      if(sscanf(table[idx].value, "%s", config->config) != 1)
        fatal("invalid format for configuration  parameter config\n");
  } else {
      strcpy(config->config, NULLFILE);
  }
  if((idx = get_str_index(table, size, "d")) >= 0) {
      if(sscanf(table[idx].value, "%s", config->dump_config) != 1)
        fatal("invalid format for configuration  parameter dump_config\n");
  } else {
      strcpy(config->dump_config, NULLFILE);
  }
  if((idx = get_str_index(table, size, "v")) >= 0) {
      if(sscanf(table[idx].value, "%s", config->v_outfile) != 1)
        fatal("invalid format for configuration  parameter v_outfile\n");
  } else {
      strcpy(config->v_outfile, NULLFILE);
  }
}

/* 
 * read a single line of trace file containing names
 * of functional blocks
 */
int read_names(FILE *fp, char **names)
{
  char line[LINE_SIZE], temp[LINE_SIZE], *src;
  int i;

  /* skip empty lines	*/
  do {
      /* read the entire line	*/
      fgets(line, LINE_SIZE, fp);
      if(feof(fp))
        fatal("not enough names in trace file\n");
      strcpy(temp, line);
      src = strtok(temp, " \r\t\n");
  } while (!src);

  /* new line not read yet	*/	
  if(line[strlen(line)-1] != '\n')
    fatal("line too long\n");

  /* chop the names from the line read	*/
  for(i=0,src=line; *src && i < MAX_UNITS; i++) {
      if(!sscanf(src, "%s", names[i]))
        fatal("invalid format of names\n");
      src += strlen(names[i]);
      while (isspace((int)*src))
        src++;
  }
  if(*src && i == MAX_UNITS)
    fatal("no. of units exceeded limit\n");

  return i;
}

/* read a single line of power trace numbers	*/
int read_vals(FILE *fp, double *vals)
{
  char line[LINE_SIZE], temp[LINE_SIZE], *src;
  int i;

  /* skip empty lines	*/
  do {
      /* read the entire line	*/
      fgets(line, LINE_SIZE, fp);
      if(feof(fp))
        return 0;
      strcpy(temp, line);
      src = strtok(temp, " \r\t\n");
  } while (!src);

  /* new line not read yet	*/	
  if(line[strlen(line)-1] != '\n')
    fatal("line too long\n");

  /* chop the power values from the line read	*/
  for(i=0,src=line; *src && i < MAX_UNITS; i++) {
      if(!sscanf(src, "%s", temp) || !sscanf(src, "%lf", &vals[i]))
        fatal("invalid format of values\n");
      src += strlen(temp);
      while (isspace((int)*src))
        src++;
  }
  if(*src && i == MAX_UNITS)
    fatal("no. of entries exceeded limit\n");

  return i;
}

int compare_names(char **n1, char **n2, int size)
{
  int i;
  int flag = 0;

  for(i=0; i < size-1; i++)
    if(strcmp(n1[i], n2[i])){
        flag = 1;
        break;
    }

  return flag;
}

char **alloc_names(int nr, int nc)
{
  int i;
  char **m;

  m = (char **) calloc (nr, sizeof(char *));
  assert(m != NULL);
  m[0] = (char *) calloc (nr * nc, sizeof(char));
  assert(m[0] != NULL);

  for (i = 1; i < nr; i++)
    m[i] =  m[0] + nc * i;

  return m;
}

void free_names(char **m)
{
  free(m[0]);
  free(m);
}

int main(int argc, char **argv)
{
  int i, j, m, idx, base = 0, count = 0, n = 0;
  int num, size, lines = 0;
  int do_transient = TRUE;
  int intv, trans_intvl, switch_intvl;
  int warmup_steps;

  char **names; double *vals;
  /* trace file pointers	*/
  FILE *pin; 
  FILE *vout = NULL;
  /* floorplan	*/
  PDN_flp_t *flp;

  /* PDN model */
  model_t *model;
  status_t *status;

  /* instantaneous power values	*/
  double *power;
  double *cur_power, *delta_power;
  double total_power = 0.0;

  /* steady state power values	*/
  double *overall_power;
  /* PDN model configuration parameters	*/
  PDN_config_t config;
  /* global configuration parameters	*/
  global_config_t global_config;
  /* table to hold options and configuration */
  str_pair table[MAX_ENTRIES];

  /* SLU Matrix for Transient Solving */
  SuperMatrix A, L, U;
  int         *perm_r; /* row permutations from partial pivoting */
  int         *perm_c; /* column permutation vector */

  if(!(argc >= 6 && argc % 2)) {
      usage(argc, argv);
      return 1;
  }

  size = parse_cmdline(table, MAX_ENTRIES, argc, argv);
  global_config_from_strs(&global_config, table, size);

  /* no PDN transient simulation, only steady state	*/
  if(!strcmp(global_config.v_outfile, NULLFILE))
    do_transient = FALSE;

  /* read configuration file	*/
  if(strcmp(global_config.config, NULLFILE))
    size += read_str_pairs(&table[size], MAX_ENTRIES, global_config.config);

  /* 
   * earlier entries override later ones. so, command line options 
   * have priority over config file 
   */
  size = str_pairs_remove_duplicates(table, size);

  /* get defaults PDN config*/
  config = default_PDN_config();
  /* modify according to command line / config file	*/
  PDN_config_add_from_strs(&config, table, size);

  /* initialization: the flp_file global configuration 
   * parameter is overridden by the layer configuration 
   * file in the layer_file_3D when the latter is specified.
   */
  flp = PDN_read_flp(global_config.flp_file);

  /* Manually set domain for PDN */
  if(config.PDN_multi_dom)
    set_flp_domain(flp, &config);

  /* Simulation mode checking */
  int is_3D = FALSE;
  if(strcmp(config.layer_file_3D, NULLFILE))
    is_3D = TRUE;
  if((!is_3D) && config.v_stacking)
    fatal("Please turn on 3D if you want to do voltage stacking.\n");
  if(do_transient){
      if(1 == config.run_PDN){
          fatal("We do not allow running transient simulation along with\
                steady state simulation.\n");
      }
      else if(2 == config.run_PDN){
          fatal("Running transient simulation and per-line steady state simulation together \
                does not make much sense. Modify voltspot.c to override this interruption.\n");
      }
  }
  if(config.PDN_multi_dom && config.v_stacking)
    fatal("Multi_dom and v_stacking do not work together yet.\n");


  /* allocate and initialize the R model for PDN	*/
  model = alloc_model(&config, flp);

  populate_C4_PDN(model);
  if(is_3D)
    populate_TSV_PDN(model);
  if(config.v_stacking)
    populate_IVR(model);
  populate_R_model_PDN(model);
  if(do_transient){
      populate_LC_model_PDN(model);

      m = PDN_trans_matrix_dim(model);
      model->trans_matrix_dim = m;
      model->last_trans = dvector(m);
      if(!is_3D){
          if( !(perm_r = (int *) calloc(m, sizeof(int))) ) fatal("Malloc fails for perm_r[].\n");
          if( !(perm_c = (int *) calloc(m, sizeof(int))) ) fatal("Malloc fails for perm_c[].\n");
      }
  }

  status = alloc_status(model);
  populate_status(model, status);

  /* n is the sum total of the number of functional blocks
   * of all the floorplans */
  n = model->total_n_blocks;

  /* allocate the power arrays	*/
  power = dvector(n);

  /* For trans PDN step input elimination */
  cur_power = dvector(n);
  delta_power = dvector(n);
  overall_power = dvector(n);

  if(!(pin = fopen(global_config.p_infile, "r")))
    fatal("unable to open power trace input file\n");
  if(do_transient && !(vout = fopen(global_config.v_outfile, "w")))
    fatal("unable to open PDN trace file for output\n");

  /* names of functional units	*/
  names = alloc_names(MAX_UNITS, STR_SIZE);

  if(read_names(pin, names) != n)
    fatal("no. of units in floorplan and trace file differ\n");


  vals = dvector(MAX_UNITS);
  /* Initialize Matrix, power trace*/
  if(do_transient){
      status->trans_counter = 0;
      status->draw_counter = 0;

      /* Zero power before the first line */
      zero_dvector(cur_power, n);
      /* Create J matrix, LU decomposition */
      if(!is_3D){
          if(model->config.PDN_gridL){
              trans_SLU_init(model, cur_power, &A, &L, &U, perm_c, perm_r);
          }
          else{
              trans_SLU_init_nogridL(model, cur_power, &A, &L, &U, perm_c, perm_r);
          }
      }

      /* Init vector  */
      PDN_init_trans_vector(model, model->last_trans);

      if(model->config.PDN_sin_pattern)
        trans_intvl = model->config.PDN_sin_totstep;
      else
        trans_intvl = model->config.ptrace_sampling_intvl * model->config.PDN_step_percycle;

      warmup_steps = model->config.PDN_ptrace_warmup * trans_intvl;

      if(model->config.v_stacking){
          double Fsw_ratio = model->config.proc_clock_freq*model->config.PDN_step_percycle/model->sc_converter->freq;
          if(fabs(floor(Fsw_ratio+0.5)-Fsw_ratio) < DELTA)
            switch_intvl = (int) Fsw_ratio;
          else
            fatal("Please use an integer ratio between SC converters switching T and solver timestep!\n");
      }
      else{
          switch_intvl = 0;
      }

      dump_PDN_config(model, vout);
      print_trans_header(model, status, vout);
  }

  /* read the instantaneous power trace	*/
  while ((num=read_vals(pin, vals)) != 0) {
      if(num != n){
          fatal("invalid trace file format\n");
      }

      /* permute the power numbers according to the floorplan order	*/
      for(i=0, base=0, count=0; i<model->n_layers; i++) {
          for(j=0; j<model->layers[i].flp->n_units; j++) {
              idx = PDN_get_blk_index(model->layers[i].flp, names[count+j]);
              power[base+idx] = vals[count+j];
          }
          count += model->layers[i].flp->n_units;
          base += model->layers[i].flp->n_units;	
      }

      /* compute PDN voltage */
      if(do_transient){
          if(is_3D){
              if(model->config.PDN_sin_pattern){
                  for(intv=0; intv<trans_intvl; intv++){
                      status->trans_counter++;
                      sin_dvector(cur_power, power, 
                                  intv/(model->config.proc_clock_freq*model->config.PDN_step_percycle), 
                                  model->config.PDN_sin_freq, n);
                      trans_matrix_build_3D(model, cur_power, &A, &L, &U, &perm_c, &perm_r);
                      compute_PDN_SLU(model, cur_power, &A, &L, &U, perm_c, perm_r);

                      if(model->config.v_stacking)
                        switch_SCconverters(model, status->trans_counter, switch_intvl);

                      perstep_droop_3D(model, vout);

                      SUPERLU_FREE (perm_r);
                      SUPERLU_FREE (perm_c);
                      Destroy_CompCol_Matrix(&A);
                      Destroy_SuperNode_Matrix(&L);
                      Destroy_CompCol_Matrix(&U);
                  }
                  lines++;
                  break;
              }
              else{
                  trans_matrix_build_3D(model, power, &A, &L, &U, &perm_c, &perm_r);
                  for(intv=0; intv<trans_intvl; intv++){
                      status->trans_counter++;
                      compute_PDN_SLU(model, power, &A, &L, &U, perm_c, perm_r);

                      if(model->config.v_stacking)
                        switch_SCconverters(model, status->trans_counter, switch_intvl);

                      if(status->trans_counter > warmup_steps){ // ignore warm up
                          PDN_trans_analyze(model, status);
                          if(!(status->trans_counter % model->config.PDN_step_percycle))
                            print_trans_analyze(model, status, vout);
                      }
                  }
                  SUPERLU_FREE (perm_r);
                  SUPERLU_FREE (perm_c);
                  Destroy_CompCol_Matrix(&A);
                  Destroy_SuperNode_Matrix(&L);
                  Destroy_CompCol_Matrix(&U);
              }
          }
          else{
              if(model->config.PDN_sin_pattern){
                  for(intv=0; intv<trans_intvl; intv++){
                      sin_dvector(cur_power, power, 
                                  intv/(model->config.proc_clock_freq*model->config.PDN_step_percycle), 
                                  model->config.PDN_sin_freq, n);
                      compute_PDN_SLU(model, cur_power, &A, &L, &U, perm_c, perm_r);
                      print_step_singlenode(model, vout, 0, 0);
                  }
                  lines++;
                  break;
              }
              else{
                  sub_dvector(delta_power, power, cur_power, n);
                  mul_val_dvector(delta_power, ((double)1/trans_intvl), n);
                  for(intv=0; intv<trans_intvl; intv++){
                      status->trans_counter++;
                      add_dvector(cur_power, cur_power, delta_power, n);

                      compute_PDN_SLU(model, cur_power, &A, &L, &U, perm_c, perm_r);

                      if(status->trans_counter > warmup_steps){ // ignore warm up
                          PDN_trans_analyze(model, status);
                          if(!(status->trans_counter % model->config.PDN_step_percycle))
                            print_trans_analyze(model, status, vout);
                      }
                  }
              }
          }
          printf("Done with ptrace line %d\n", lines+1);
      }

      /* per-line steady state */
      if(2 == model->config.run_PDN){
          steady_state_PDN(model, power);
          PDN_steady_analyze(model, status);
          if(model->config.animation){
              draw_single_gif(model, status, status->draw_counter, STEADY);
              status->draw_counter++;
          }
          printf("Done with ptrace line %d\n", lines+1);
      }

      /* for computing average	*/
      for(i=0, base=0; i<model->n_layers; i++) {
          for(j=0; j<model->layers[i].flp->n_units; j++)
            overall_power[base+j] += power[base+j];
          base += model->layers[i].flp->n_units;	
      }
      lines++;
  }

  if(!lines)
    fatal("no power numbers in trace file\n");

  int tot_cycles;
  if(do_transient){
      /* Mark end in vtrace */
      fprintf(vout, "#END\n");

      // calculate per cycle average gridnode noise
      tot_cycles = (lines - model->config.PDN_ptrace_warmup) * model->config.ptrace_sampling_intvl;
      mul_val_dvector(status->gridstats.integral_2D[0], ((double)1/tot_cycles), model->rows*model->cols);

      dump_files(model, status, TRANSIENT);
  }

  if(2 == model->config.run_PDN){
      dump_files(model, status, STEADY);
  }

  if(model->config.animation)
    create_animation(status->draw_counter, model->config.animation);

  /* for computing average	*/
  for(i=0, base=0; i<model->n_layers; i++) {
      for(j=0; j<model->layers[i].flp->n_units; j++) {
          overall_power[base+j] /= lines;
          total_power += overall_power[base+j];
      }
      base += model->layers[i].flp->n_units;	
  }

  if(1 == model->config.run_PDN){
      /* steady state PDN */
      steady_state_PDN(model, overall_power);
      /* analyze and print PDN status */
      PDN_steady_analyze(model, status);
      print_steady_analyze(model, status);
      dump_files(model, status, STEADY);
  }

  /* cleanup	*/
  fclose(pin);
  PDN_free_flp(flp);
  free_dvector(power);
  free_dvector(cur_power);
  free_dvector(delta_power);
  free_dvector(overall_power);
  free_names(names);
  free_dvector(vals);
  if(do_transient && (!is_3D)){
      fclose(vout);
      SUPERLU_FREE (perm_r);
      SUPERLU_FREE (perm_c);
      Destroy_CompCol_Matrix(&A);
      Destroy_SuperNode_Matrix(&L);
      Destroy_CompCol_Matrix(&U);
  }
  free_status(model, status);
  delete_model(model);

  return 0;
}
