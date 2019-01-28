#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include <math.h>
#include <assert.h>

#include "pad.h"

/* parse pad location information	
 * initialize padloc matrix in PDN struct
 */
void parse_pad_loc(model_t *model, char *file)
{
  char str[LINE_SIZE], copy[LINE_SIZE]; 
  char s[STR_SIZE];
  char name[STR_SIZE];
  int grid_x, grid_y;
  int pad_x, pad_y;
  char *ptr;
  FILE *fp;

  /* short cuts */  
  int nr = model->rows;
  int nc = model->cols;
  int rpg = model->c4->pad_grid_row;
  int cpg = model->c4->pad_grid_col;
  int itv_row = (int)(nr - 1)/(rpg - 1);
  int itv_col = (int)(nc - 1)/(cpg - 1);

  if (!strcasecmp(file, "stdin"))
    fp = stdin;
  else
    fp = fopen (file, "r");

  if (!fp) {
      sprintf(s, "error opening padloc file %s\n", file);
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

      if((0 == model->config.padloc_format) || 
         (1 == model->config.padloc_format)){
          // parse index as metal grid
          if (sscanf(copy, "%s%d%d", name, &grid_x, &grid_y) == 3) {
              if ((grid_x >= model->cols) || (grid_y >= model->rows) ||
                  (grid_x < 0) || (grid_y < 0)){
                  printf("x = %d, y = %d\n", grid_x, grid_y);
                  printf("grid_size = %d*%d\n", model->cols, model->rows);
                  fatal("Pad location does not fit in current grid!\n");
              }
              if (!strcmp(name, "V")){
                  if (PGPAD & model->c4->vdd_loc[nr - grid_y - 1][grid_x])
                    warning("Duplication exists in pad location file .padloc\n");
                  model->c4->vdd_loc[nr - grid_y - 1][grid_x] |= PGPAD;
              }
              else if (!strcmp(name, "G")){
                  if (PGPAD & model->c4->gnd_loc[nr - grid_y - 1][grid_x])
                    warning("Duplication exists in pad location file .padloc\n");
                  model->c4->gnd_loc[nr - grid_y - 1][grid_x] |= PGPAD;
              }
              else
                fatal("Pad location is neither V(vdd) nor G(ground)\n");
          }
          else
            fatal("invalid pad location file format\n");
      }
      else{
          // parse index as pad grid
          if (sscanf(copy, "%s%d%d", name, &pad_x, &pad_y) == 3) {
              grid_x = pad_x * itv_col;
              grid_y = (rpg - pad_y - 1) * itv_row;

              if ((pad_x >= cpg) || (pad_y >= rpg) ||
                  (pad_x < 0) || (pad_y < 0)){
                  printf("pad_x = %d, pad_y = %d\n", pad_x, pad_y);
                  printf("grid_size = %d*%d\n", cpg, rpg);
                  fatal("Pad location does not fit in current grid!\n");
              }
              if (!strcmp(name, "V")){
                  if (PGPAD & model->c4->vdd_loc[grid_y][grid_x])
                    warning("Duplication exists in pad location file .padloc\n");
                  model->c4->vdd_loc[grid_y][grid_x] |= PGPAD;
              }
              else if (!strcmp(name, "G")){
                  if (PGPAD & model->c4->gnd_loc[grid_y][grid_x])
                    warning("Duplication exists in pad location file .padloc\n");
                  model->c4->gnd_loc[grid_y][grid_x] |= PGPAD;
              }
              else
                fatal("Pad location is neither V(vdd) nor G(ground)\n");
          }
          else
            fatal("invalid pad location file format\n");
      }
  }

  if(fp != stdin)
    fclose(fp);
}

/* check whether a grid node is allowed to place IO */
int node_allow_io(model_t *model, int r, int c)
{
  int allow = 0;
  char *name;
  blist_t *ptr;

  ptr = model->layers[0].b2gmap[r][c];
  name = model->layers[0].flp->units[ptr->idx].name;

  if(strlen(name) < 3){
      warning("Some block's name is shorter than 3 char, not going put IO under it\n");
      allow = 0;
  }

  if(name[0] == 'M' &&
     name[1] == 'C'){
      allow = 1;
  }
  else if(name[0] == 'L' &&
          name[1] == '2'){
      allow = 1;
  }
  else if(name[0] == 'N' &&
          name[1] == 'o' &&
          name[2] == 'C'){
      allow = 1;
  }
  else
    allow = 0;

  return allow;
}

/* check whether a grid node is pad candidate */
int is_grid_pad(model_t *model, int r, int c, int layer)
{
  int rpg = model->c4->pad_grid_row;
  int cpg = model->c4->pad_grid_col;
  int nr = model->rows;
  int nc = model->cols;
  int itv_r, itv_c;

  itv_r = (int)((nr - 1)/(rpg - 1));
  itv_c = (int)((nc - 1)/(cpg - 1));

  if(r % itv_r)
    return 0;
  if(c % itv_c)
    return 0;

  if(LAYER_VDD == layer){
      if(((r/itv_r)%2) == ((c/itv_c)%2))
        return 1;
  }
  else if(LAYER_GND == layer){
      if(((r/itv_r)%2) != ((c/itv_c)%2))
        return 1;
  }
  else
    fatal("Wrong arguemnt layer received by is_grid_pad!\n");

  return 0;
}

/* Mark I/O pads as IOPAD or mark range as IORAN
 * Firt mark pads under MCs, if need to mark more,
 * expand along y axis.
 * */
void mark_pads(model_t *model, int pad_type)
{
  int i, m, n;
  PDN_flp_t *flp = model->layers[0].flp;
  int flp_count = flp->n_units;
  double w = model->width;
  double h = model->height;
  double lx, rx, by, ty;
  int lx_cordt, rx_cordt, by_cordt, ty_cordt;
  int **vp = model->c4->vdd_loc;
  int **gp = model->c4->gnd_loc;
  int rpg = model->c4->pad_grid_row;
  int cpg = model->c4->pad_grid_col;
  int nr = model->rows;
  int nc = model->cols;
  int top, btm, lft, rgt;
  int mc_count = 0;
  int padleft_v = 0;
  int padleft_g = 0;

  for(i = 0; i < flp_count; i++){
      if((flp->units[i].name[0] == 'M') && (flp->units[i].name[1] == 'C')){
          if(IOPAD & pad_type){
              padleft_v = model->config.MC_pads / (2 * 0.5);
              padleft_g = model->config.MC_pads / (2 * 0.5);
          }
          else if(IORAN & pad_type){
              padleft_v = model->config.MC_pads / (2 * model->config.IO_dense);
              padleft_g = model->config.MC_pads / (2 * model->config.IO_dense);
          }

          lx = flp->units[i].leftx;
          rx = lx + flp->units[i].width;
          by = flp->units[i].bottomy;
          ty = by + flp->units[i].height;
          lx_cordt = ceil(lx * (nc-1) / w);
          rx_cordt = floor(rx * (nc-1) / w);
          ty_cordt = nr - ceil(by * (nr-1) / h) - 1;
          by_cordt = nr - floor(ty * (nr-1) / h) - 1;

          //under MC
          for(m=by_cordt; m<=ty_cordt; m++)
            for(n=lx_cordt; n<=rx_cordt; n++){
                if ((padleft_v > 0) && 
                    (1 == is_grid_pad(model, m, n, LAYER_VDD))){
                    if (!(vp[m][n] & pad_type)){
                        vp[m][n] |= pad_type;
                        padleft_v--;
                    }
                }

                if ((padleft_g > 0) && 
                    (1 == is_grid_pad(model, m, n, LAYER_GND))){
                    if (!(gp[m][n] & pad_type)){
                        gp[m][n] |= pad_type;
                        padleft_g--;
                    }
                }
            }

          top = ty_cordt;
          btm = by_cordt;
          lft = lx_cordt;
          rgt = rx_cordt;
          while(padleft_v > 0){
              // top
              for(n=lft; n<=rgt; n++){
                  if(padleft_v == 0)
                    break;
                  if(top > nr-1)
                    break;
                  if (1 == is_grid_pad(model, top, n, LAYER_VDD)){
                      if (!(vp[top][n] & pad_type) && node_allow_io(model, top, n)){
                          vp[top][n] |= pad_type;
                          padleft_v--;
                      }
                  }
              }

              // bottom
              for(n=lft; n<=rgt; n++){
                  if(padleft_v == 0)
                    break;
                  if(btm < 0)
                    break;
                  if (1 == is_grid_pad(model, btm, n, LAYER_VDD)){
                      if (!(vp[btm][n] & pad_type) && node_allow_io(model, btm, n)){
                          vp[btm][n] |= pad_type;
                          padleft_v--;
                      }
                  }
              }

              // left
              for(m=btm; m<=top; m++){
                  if(padleft_v == 0)
                    break;
                  if(lft < 0)
                    break;
                  if (1 == is_grid_pad(model, m, lft, LAYER_VDD)){
                      if (!(vp[m][lft] & pad_type) && node_allow_io(model, m, lft)){
                          vp[m][lft] |= pad_type;
                          padleft_v--;
                      }
                  }
              }

              // right
              for(m=btm; m<=top; m++){
                  if(padleft_v == 0)
                    break;
                  if(rgt > nc-1)
                    break;
                  if (1 == is_grid_pad(model, m, rgt, LAYER_VDD)){
                      if (!(vp[m][rgt] & pad_type) && node_allow_io(model, m, rgt)){
                          vp[m][rgt] |= pad_type;
                          padleft_v--;
                      }
                  }
              }

              if(top < nr-1)
                top++;
              else if(btm > 0)
                btm--;
              else if(lft > 0)
                lft--;
              else if(rgt < nc-1)
                rgt++;
              else{
                  warning("All VDD pads are marked!\n");
                  break;
              }
          }

          top = ty_cordt;
          btm = by_cordt;
          lft = lx_cordt;
          rgt = rx_cordt;

          while(padleft_g > 0){
              // top
              for(n=lft; n<=rgt; n++){
                  if(padleft_g == 0)
                    break;
                  if(top > nr-1)
                    break;
                  if (1 == is_grid_pad(model, top, n, LAYER_GND)){
                      if (!(gp[top][n] & pad_type) && node_allow_io(model, top, n)){
                          gp[top][n] |= pad_type;
                          padleft_g--;
                      }
                  }
              }

              // bottom
              for(n=lft; n<=rgt; n++){
                  if(padleft_g == 0)
                    break;
                  if(btm < 0)
                    break;
                  if (1 == is_grid_pad(model, btm, n, LAYER_GND)){
                      if (!(gp[btm][n] & pad_type) && node_allow_io(model, btm, n)){
                          gp[btm][n] |= pad_type;
                          padleft_g--;
                      }
                  }
              }

              // left
              for(m=btm; m<=top; m++){
                  if(padleft_g == 0)
                    break;
                  if(lft < 0)
                    break;
                  if (1 == is_grid_pad(model, m, lft, LAYER_GND)){
                      if (!(gp[m][lft] & pad_type) && node_allow_io(model, m, lft)){
                          gp[m][lft] |= pad_type;
                          padleft_g--;
                      }
                  }
              }

              // right
              for(m=btm; m<=top; m++){
                  if(padleft_g == 0)
                    break;
                  if(rgt > nc-1)
                    break;
                  if (1 == is_grid_pad(model, m, rgt, LAYER_GND)){
                      if (!(gp[m][rgt] & pad_type) && node_allow_io(model, m, rgt)){
                          gp[m][rgt] |= pad_type;
                          padleft_g--;
                      }
                  }
              }

              if(top < nr-1)
                top++;
              else if(btm > 0)
                btm--;
              else if(lft > 0)
                lft--;
              else if(rgt < nc-1)
                rgt++;
              else{
                  warning("All GND pads are marked!\n");
                  break;
              }

          }
          if(padleft_v)
            printf("There are %d Vdd pads not marked yet for MC%d\n", padleft_v, mc_count);
          if(padleft_g)
            printf("There are %d Gnd pads not marked yet for MC%d\n", padleft_g, mc_count);

          mc_count++;
      }
  }
}

void mark_all(model_t *model, int pad_type)
{
  int i, j, nr, nc;
  int rpg, cpg;//pad grid
  int r_cordt, c_cordt;//cordt for pad
  int itv_row, itv_col;

  /* shortcuts */
  nr = model->rows;
  nc = model->cols;
  rpg = model->c4->pad_grid_row;
  cpg = model->c4->pad_grid_col;
  itv_row = (int)(nr - 1)/(rpg - 1);
  itv_col = (int)(nc - 1)/(cpg - 1);

  for(i=0; i<rpg; i++)
    for(j=0; j<cpg; j++){
        r_cordt = i * itv_row;
        c_cordt = j * itv_col;
        if ( (i%2) == (j%2) ){
            model->c4->vdd_loc[r_cordt][c_cordt] |= pad_type;
        }
        else {
            model->c4->gnd_loc[r_cordt][c_cordt] |= pad_type;
        }
    }
}

//radom remove pads for IO
void remove_pg_for_io(model_t *model)
{
  int i, j;
  int vcounter, gcounter;

  vcounter = 0;
  gcounter = 0;

  /* shortcuts */
  int nr = model->rows;
  int nc = model->cols;
  int **vloc = model->c4->vdd_loc;
  int **gloc = model->c4->gnd_loc;

  for(i=0; i<nr; i++){
      for(j=0; j<nc; j++){
          if ((vloc[i][j] & PGPAD) && (vloc[i][j] & IOPAD) &&
              (!check_pad_neighbour(model, i, j, LAYER_VDD)) ){
              vloc[i][j] &= (~PGPAD);
              vcounter++;
          }
          if ((gloc[i][j] & PGPAD) && (gloc[i][j] & IOPAD) &&
              (!check_pad_neighbour(model, i, j, LAYER_GND)) ){
              gloc[i][j] &= (~PGPAD);
              gcounter++;
          }
      }
  }
}

/* Find the neighbor pads of a given pad node.
 * If all of them are connected to power pad, return 0,
 * if not, return 1
 */
int check_pad_neighbour(model_t *model, int r, int c, int grid_type)
{
  int nr, nc;
  int **loc;
  int rpg_r, rpg_c;
  int itv_r, itv_c;
  int rpg, cpg;//pad grid
  int cordt_r,cordt_c;
  int flag = 0;

  /* Shortcuts */
  if(LAYER_VDD == grid_type){
      loc = model->c4->vdd_loc;
  }
  else if(LAYER_GND == grid_type){
      loc = model->c4->gnd_loc;
  }
  else
    fatal("Func check_pad_neighbour received wrong argument grid_type!\n");

  nr = model->rows;
  nc = model->cols;
  rpg = model->c4->pad_grid_row;
  cpg = model->c4->pad_grid_col;

  itv_r = (int)((nr - 1)/(rpg - 1));
  itv_c = (int)((nc - 1)/(cpg - 1));

  if((r%itv_r) || (c%itv_c))
    fatal("Function check_pad_neighbour: got non-pad grid input!\n");

  rpg_r = r / itv_r;
  rpg_c = c / itv_c;

  // North
  if(rpg_r < rpg-2){
      cordt_r = itv_r * (rpg_r + 2);
      cordt_c = itv_c * rpg_c;
      if (!(PGPAD & loc[cordt_r][cordt_c]))
        flag = 1;
  }
  // South
  if(rpg_r > 1){
      cordt_r = itv_r * (rpg_r - 2);
      cordt_c = itv_c * rpg_c;
      if (!(PGPAD & loc[cordt_r][cordt_c]))
        flag = 1;
  }
  // West
  if(rpg_c > 1){
      cordt_r = itv_r * rpg_r;
      cordt_c = itv_c * (rpg_c - 2);
      if (!(PGPAD & loc[cordt_r][cordt_c]))
        flag = 1;
  }
  // East
  if(rpg_c < cpg-2){
      cordt_r = itv_r * rpg_r;
      cordt_c = itv_c * (rpg_c + 2);
      if (!(PGPAD & loc[cordt_r][cordt_c]))
        flag = 1;
  }
  return flag;

  /* An alternative way to do interleave
   * Tend to end up with uneven Vdd/Gnd number
  // North East
  if((rpg_r < rpg-1) && (rpg_c < cpg-1)){
  cordt_r = itv_r * (rpg_r + 1);
  cordt_c = itv_c * (rpg_c + 1);
  if (!(PGPAD & loc[cordt_r][cordt_c]))
  flag = 1;
  }
  // North West
  if((rpg_r < rpg-1) && (rpg_c > 0)){
  cordt_r = itv_r * (rpg_r + 1);
  cordt_c = itv_c * (rpg_c - 1);
  if (!(PGPAD & loc[cordt_r][cordt_c]))
  flag = 1;
  }
  // South West
  if((rpg_r > 0) && (rpg_c > 0)){
  cordt_r = itv_r * (rpg_r - 1);
  cordt_c = itv_c * (rpg_c - 1);
  if (!(PGPAD & loc[cordt_r][cordt_c]))
  flag = 1;
  }
  // South East
  if((rpg_r > 0) && (rpg_c < cpg-1)){
  cordt_r = itv_r * (rpg_r - 1);
  cordt_c = itv_c * (rpg_c + 1);
  if (!(PGPAD & loc[cordt_r][cordt_c]))
  flag = 1;
  }
  */
}

void dump_anypadloc(model_t *model, char *file, int type)
{
  int i, j;
  int pad_r, pad_c;
  char str[STR_SIZE];
  FILE *fp;

  /* shortcuts */
  int nr = model->rows;
  int nc = model->cols;
  int **vloc = model->c4->vdd_loc;
  int **gloc = model->c4->gnd_loc;
  int rpg = model->c4->pad_grid_row;
  int cpg = model->c4->pad_grid_col;
  int itv_row = (int)(nr - 1)/(rpg - 1);
  int itv_col = (int)(nc - 1)/(cpg - 1);

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

  if((0 == model->config.padloc_format) || 
     (2 == model->config.padloc_format)){
      // dump metal grid index
      for(i=0; i < nr; i++){
          for(j=0; j < nc; j++){
              if (vloc[i][j] & type){
                  fprintf(fp, "V\t%d\t%d\n", j, nr-i-1);
              }
              if (gloc[i][j] & type){
                  fprintf(fp, "G\t%d\t%d\n", j, nr-i-1);
              }
          }
      }
  }
  else{
      // dump pad grid index
      for(i=0; i < nr; i++){
          for(j=0; j < nc; j++){
              if (vloc[i][j] & type){
                  pad_r = i/itv_row;
                  pad_c = j/itv_col;
                  fprintf(fp, "V\t%d\t%d\n", pad_c, rpg-pad_r-1);
              }
              if (gloc[i][j] & type){
                  pad_r = i/itv_row;
                  pad_c = j/itv_col;
                  fprintf(fp, "G\t%d\t%d\n", pad_c, rpg-pad_r-1);
              }
          }
      }
  }

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

void dump_cur_with_cordt(model_t *model, status_t *status, char *file)
{
  int i, j;
  int pad_r, pad_c;
  char str[STR_SIZE];
  FILE *fp;

  /* shortcuts */
  int nr = model->rows;
  int nc = model->cols;
  int rpg = model->c4->pad_grid_row;
  int cpg = model->c4->pad_grid_col;
  int itv_row = (int)(nr - 1)/(rpg - 1);
  int itv_col = (int)(nc - 1)/(cpg - 1);
  int **vloc = model->c4->vdd_loc;
  int **gloc = model->c4->gnd_loc;
  double **vcur = status->vdd_cur;
  double **gcur = status->gnd_cur;
  double padD = model->config.PDN_padD/2;
  double pad_area = PI * padD * padD;

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

  if((0 == model->config.padloc_format) || 
     (2 == model->config.padloc_format)){
      // dump metal grid index
      for(i=0; i < nr; i++){
          for(j=0; j < nc; j++){
              if (vloc[i][j] & PGPAD){
                  fprintf(fp, "V\t%d\t%d\t%e\n", j, nr-i-1, vcur[i][j]/pad_area);
              }
              if (gloc[i][j] & PGPAD){
                  fprintf(fp, "G\t%d\t%d\t%e\n", j, nr-i-1, gcur[i][j]/pad_area);
              }
          }
      }
  }
  else{
      // dump pad grid index
      for(i=0; i < nr; i++){
          for(j=0; j < nc; j++){
              if (vloc[i][j] & PGPAD){
                  pad_r = i/itv_row;
                  pad_c = j/itv_col;
                  fprintf(fp, "V\t%d\t%d\t%e\n", pad_c, rpg-pad_r-1, vcur[i][j]/pad_area);
              }
              if (gloc[i][j] & PGPAD){
                  pad_r = i/itv_row;
                  pad_c = j/itv_col;
                  fprintf(fp, "G\t%d\t%d\t%e\n", pad_c, rpg-pad_r-1, gcur[i][j]/pad_area);
              }
          }
      }
  }

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

// count pads that at least have type pad_type
int count_pads(int **loc, int nr, int nc, int pad_type)
{
  int sum = 0;
  int i, j;

  for(i = 0; i < nr; i++)
    for(j = 0; j < nc; j++){
        if ((pad_type & loc[i][j]) == pad_type)
          sum++;
    }

  return sum;
}

