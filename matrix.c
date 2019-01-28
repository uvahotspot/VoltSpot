#include <stdio.h>
#include <stdlib.h>

/* Lib for SuperLU */
#include "slu_ddefs.h"

#include "PDN_sim.h"
#include "util.h"
#include "pad.h"

SuperMatrix build_steady_grid_matrix(model_t *model)
{
  SuperMatrix A;
  double   *a;
  int      *asub, *xa;
  double   *cooV;
  int      *cooX, *cooY;

  int      i, j, l, m, n, nnz;
  int      direc;
  double   dia_val;
  int      curidx, grididx;
  int      xoffset, yoffset;
  int      has_pad;
  int      *padloc;
  int      num_tsv;
  double   Rx, Ry, Rv;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  double Rp = model->c4->pad_r;
  int is_3D = model->is_3D;
  int mtdmn  = model->config.PDN_multi_dom;

  //count the total number of TSVs
  num_tsv = 0;
  if(is_3D){
      for(l=0; l<nl-1; l++){
          num_tsv += model->layers[l].tsv.num_vdd;
          num_tsv += model->layers[l].tsv.num_gnd;
      }
  }

  /* Initialize matrix A. */
  /* Layout: L0v L1v ... Ln-1v L0g L1g ... Ln-1g  */
  m = n = 2 * nl * nc * nr;
  /* Num of non-zeros
   * Five diagonal, c*r, c*(r-1), (c-1)*r, c*(r-1), (c-1)*r
   */
  nnz = 5*nr*nc- 2*(nr+nc);
  nnz *= nl;
  nnz *= 2;  // both Vdd and Gnd
  nnz += 2*num_tsv;

  if ( !(cooV = doubleMalloc(nnz)) ) fatal("Malloc fails for cooV[].\n");
  if ( !(cooX = intMalloc(nnz)) ) fatal("Malloc fails for cooX[].\n");
  if ( !(cooY = intMalloc(nnz)) ) fatal("Malloc fails for cooY[].\n");

  if ( !(a = doubleMalloc(nnz)) ) fatal("Malloc fails for a[].\n");
  if ( !(asub = intMalloc(nnz)) ) fatal("Malloc fails for asub[].\n");
  if ( !(xa = intMalloc(n+1)) ) fatal("Malloc fails for xa[].\n");

  curidx = 0;
  /* VDD grid diagonal block */
  padloc = model->c4->vdd_loc[0];
  for(l=0; l<nl; l++){
      //calculate grid R
      Rx = 0; Ry = 0;
      for(i=0; i<model->layers[l].metal_layers.n_metal; i++){
          direc = model->layers[l].metal_layers.geo[i].direc;
          if(MLCF_X == direc)
            Rx += 1/model->layers[l].metal_layers.gridRL[i].r;
          else
            Ry += 1/model->layers[l].metal_layers.gridRL[i].r;
      }
      Rx = 1/Rx; Ry = 1/Ry;
      Rv = model->layers[l].tsv.r; //assume uniform TSV R

      xoffset = l*nr*nc;
      yoffset = l*nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            if((0 == l) && (PGPAD & padloc[grididx]))
              has_pad = 1;
            else
              has_pad = 0;

            dia_val = 0;
            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;}
                else{cooV[curidx] = 0; curidx++;}
            }
            if(is_3D){
                if((l<nl-1) && (VDDTSV == model->layers[l].tsv.loc[i][j])){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1/Rv; curidx++; dia_val += 1/Rv;
                }
                if((l>0) && (VDDTSV == model->layers[l-1].tsv.loc[i][j])){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1/Rv; curidx++; dia_val += 1/Rv;
                }
            }

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = has_pad ? (dia_val+1/Rp) : dia_val; curidx++;
        }
  }

  /* GND grid diagonal block */
  padloc = model->c4->gnd_loc[0];
  for(l=0; l<nl; l++){
      //calculate grid R
      Rx = 0; Ry = 0;
      for(i=0; i<model->layers[l].metal_layers.n_metal; i++){
          direc = model->layers[l].metal_layers.geo[i].direc;
          if(MLCF_X == direc)
            Rx += 1/model->layers[l].metal_layers.gridRL[i].r;
          else
            Ry += 1/model->layers[l].metal_layers.gridRL[i].r;
      }
      Rx = 1/Rx; Ry = 1/Ry;
      Rv = model->layers[l].tsv.r; //assume uniform TSV R

      xoffset = nl*nr*nc + l*nr*nc;
      yoffset = nl*nr*nc + l*nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            dia_val = 0;
            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;
            }

            if(is_3D){
                if((l<nl-1) && (GNDTSV == model->layers[l].tsv.loc[i][j])){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1/Rv; curidx++; dia_val += 1/Rv;
                }
                if((l>0) && (GNDTSV == model->layers[l-1].tsv.loc[i][j])){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1/Rv; curidx++; dia_val += 1/Rv;
                }
            }

            // bumps
            if((0 == l) && (PGPAD & padloc[grididx]))
              dia_val += 1/Rp; 

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = dia_val; curidx++;
        }
  }

  if(curidx != nnz)
    fatal("Steady-state Matrix build error: less elements than nnz\n");

  coo2csc(n, nnz, cooX, cooY, cooV, asub, xa, a);

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

  free(cooV);
  free(cooX);
  free(cooY);

  return A;
}

SuperMatrix build_steady_rhs_vector(model_t *model, model_vector_t *power, double **rhs)
{
  SuperMatrix B;
  int      idx, nrhs;
  int      i, j, l;
  int      **loc;
  double   pkg_vdd, pkg_gnd;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int m  = 2 * nl * nr * nc;
  double Rp = model->c4->pad_r;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  pkg_vdd = model->last_steady[2*nl*nr*nc + PKG_VDD];
  pkg_gnd = model->last_steady[2*nl*nr*nc + PKG_GND];

  nrhs = 1;
  if ( !(*rhs = doubleMalloc(m*nrhs)) ) fatal("Malloc fails for rhs[].\n");

  //VDD layer
  for(l=0; l<nl; l++)
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++) {
          idx = l*nr*nc + i*nc + j ;
          loc = model->c4->vdd_loc;
          if((0 == l) && (PGPAD & loc[i][j])) // Bottom layer connecting to C4 pads
            (*rhs)[idx] = pkg_vdd/Rp - power->cuboid[l][i][j]/(vdd-gnd);
          else// I = P/v
            (*rhs)[idx] = -power->cuboid[l][i][j]/(vdd-gnd);
      }

  //GND layer
  for(l=0; l<nl; l++)
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++) {
          idx = nl*nr*nc + l*nr*nc + i*nc + j ;
          loc = model->c4->gnd_loc;
          if((0 == l) && (PGPAD & loc[i][j])) // Bottom layer connecting to C4 pads
            (*rhs)[idx] = pkg_gnd/Rp + power->cuboid[l][i][j]/(vdd-gnd);
          else// I = P/v
            (*rhs)[idx] = power->cuboid[l][i][j]/(vdd-gnd);
      }

  dCreate_Dense_Matrix(&B, m, nrhs, *rhs, m, SLU_DN, SLU_D, SLU_GE);

  return B;
}

SuperMatrix build_steady_grid_matrix_vs(model_t *model, model_vector_t *power)
{
  SuperMatrix A;
  double   *a;
  int      *asub, *xa;
  double   *cooV;
  int      *cooX, *cooY;

  int      i, j, k, l, m, n, nnz;
  int      cur_layer;
  int      direc;
  double   dia_val;
  int      curidx, grididx;
  int      xoffset, yoffset;
  int      num_tsv;
  double   Rx, Ry, Rp, Rv, Rivr, Rload;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int *padloc_v = model->c4->vdd_loc[0];
  int *padloc_g = model->c4->gnd_loc[0];
  int *ivrloc = model->sc_converter->loc[0];
  int num_ivr = model->sc_converter->num_IVR;
  double vdd = model->config.vdd - model->config.gnd;

  // Count the total number of TSVs
  // Only using gnd TSVs to connect layers.
  // Vdd TSVs are used for top layer supply, but only the ones connected to bumps are populated
  num_tsv = 0;
  for(l=0; l<nl-1; l++)
    num_tsv += model->layers[l].tsv.num_gnd;

  /* Initialize matrix A. */
  /* Layout: L0g L0v L1g L1v ... Ln-1g Ln-1v */
  m = n = 2 * nl * nc * nr;
  /* Num of non-zeros
   * Five diagonal, c*r, c*(r-1), (c-1)*r, c*(r-1), (c-1)*r
   */
  nnz = 5*nr*nc- 2*(nr+nc);
  nnz *= nl;
  nnz *= 2;  // both Vdd and Gnd
  nnz += 2*num_tsv;
  nnz += 2*nl*nr*nc; // load resistors
  nnz += 2*(nl-2)*num_ivr; // for IVRs

  if ( !(cooV = doubleMalloc(nnz)) ) fatal("Malloc fails for cooV[].\n");
  if ( !(cooX = intMalloc(nnz)) ) fatal("Malloc fails for cooX[].\n");
  if ( !(cooY = intMalloc(nnz)) ) fatal("Malloc fails for cooY[].\n");

  if ( !(a = doubleMalloc(nnz)) ) fatal("Malloc fails for a[].\n");
  if ( !(asub = intMalloc(nnz)) ) fatal("Malloc fails for asub[].\n");
  if ( !(xa = intMalloc(n+1)) ) fatal("Malloc fails for xa[].\n");

  curidx = 0;
  for(l=0; l<2*nl; l++){
      //calculate grid R
      cur_layer = l >> 1;
      Rx = 0; Ry = 0;
      for(i=0; i<model->layers[cur_layer].metal_layers.n_metal; i++){
          direc = model->layers[cur_layer].metal_layers.geo[i].direc;
          if(MLCF_X == direc)
            Rx += 1/model->layers[cur_layer].metal_layers.gridRL[i].r;
          else
            Ry += 1/model->layers[cur_layer].metal_layers.gridRL[i].r;
      }
      Rx = 1/Rx; Ry = 1/Ry;
      Rv = model->layers[cur_layer].tsv.r; //assume uniform TSV R
      Rp = model->c4->pad_r;
      Rivr = model->sc_converter->R_drop;
      if((2*nl-1) == l)
        for(k=0; k<nl-1; k++)
          Rp += model->layers[k].tsv.r;

      xoffset = l*nr*nc;
      yoffset = l*nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            if (power->cuboid[cur_layer][i][j] > 0)
              Rload = vdd*vdd/power->cuboid[cur_layer][i][j];
            else
              Rload = LARGEINT;

            dia_val = 0;
            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1/Ry; curidx++; dia_val += 1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1/Rx; curidx++; dia_val += 1/Rx;
            }
            if(l%2){//Vdd net
                // Load resistor
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                cooV[curidx] = -1/Rload; curidx++; dia_val += 1/Rload;
                // TSV
                if((cur_layer<nl-1) && (GNDTSV == model->layers[cur_layer].tsv.loc[i][j])){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1/Rv; curidx++; dia_val += 1/Rv;
                }
            }
            else{//Gnd net
                // Load resistor
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                cooV[curidx] = -1/Rload; curidx++; dia_val += 1/Rload;
                // TSV
                if((cur_layer>0) && (GNDTSV == model->layers[cur_layer-1].tsv.loc[i][j])){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1/Rv; curidx++; dia_val += 1/Rv;
                }
            }

            //add pads
            if(((0 == l) && (PGPAD & padloc_g[grididx])) || 
               ((2*nl-1 == l) && (PGPAD & padloc_v[grididx])))
              dia_val += 1/Rp;
            //add IVRs
            // head/foot room model
            if((l%2) && (cur_layer<nl-1) && (IVRLOC == ivrloc[grididx])){
                dia_val += 1/Rivr;
                if(cur_layer>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-2*nr*nc+yoffset;
                    cooV[curidx] = -1/(2*Rivr); curidx++;
                }
                if(cur_layer<nl-2){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+2*nr*nc+yoffset;
                    cooV[curidx] = -1/(2*Rivr); curidx++;
                }
            }

            // ideal voltage (l*Vdd) model
            //if((0 == l%2) && (l>0) && (IVRLOC == ivrloc[grididx]))
            //  dia_val += 1/Rivr;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = dia_val; curidx++;
        }
  }

  if(curidx != nnz)
    fatal("Steady-state Matrix build error: less elements than nnz\n");

  coo2csc(n, nnz, cooX, cooY, cooV, asub, xa, a);

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

  free(cooV);
  free(cooX);
  free(cooY);

  return A;
}

SuperMatrix build_steady_rhs_vector_vs(model_t *model, double **rhs)
{
  SuperMatrix B;
  int      idx, nrhs;
  int      i, j, l, k;
  double   Rp, Rivr;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int m  = 2 * nl * nr * nc;
  int **padloc_v = model->c4->vdd_loc;
  int **padloc_g = model->c4->gnd_loc;
  int **ivrloc = model->sc_converter->loc;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;
  double *v = model->last_steady;

  nrhs = 1;
  if ( !(*rhs = doubleMalloc(m*nrhs)) ) fatal("Malloc fails for rhs[].\n");

  for(l=0; l<2*nl; l++){
      Rp = model->c4->pad_r;
      Rivr = model->sc_converter->R_drop;
      if((2*nl-1) == l)
        for(k=0; k<nl-1; k++)
          Rp += model->layers[k].tsv.r;

      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++) {
            idx = l*nr*nc + i*nc + j ;
            if((0 == l) && (PGPAD & padloc_g[i][j]))
              (*rhs)[idx] = gnd/Rp;
            else if((2*nl-1 == l) && (PGPAD & padloc_v[i][j]))
              (*rhs)[idx] = nl*vdd/Rp; // Stacking voltage multiplication
            else if((1==l) && (IVRLOC == ivrloc[i][j])){
                (*rhs)[idx] = 0;
            }
            else if((2*nl-3==l) && (IVRLOC == ivrloc[i][j])){
                (*rhs)[idx] += nl*vdd/(2*Rivr);
            }
            //else if((0 == l%2) && (l>0) && (IVRLOC == ivrloc[i][j]))
            //    (*rhs)[idx] = (l/2)*vdd/Rivr;
            else
              (*rhs)[idx] = 0;
        }
  }

  dCreate_Dense_Matrix(&B, m, nrhs, *rhs, m, SLU_DN, SLU_D, SLU_GE);

  return B;
}

void trans_SLU_init(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r)
{
  SuperMatrix B, Z;
  double   *a, *b, *rhs;
  int      *asub, *xa, *bsub, *xb;
  double   *cooV;
  int      *cooX, *cooY;

  int      info, i, j, l, m, n, nnz;
  int      tmp_i, tmp_j;
  int      xoffset, yoffset;
  int      grididx, curidx;
  superlu_options_t options;
  SuperLUStat_t stat;
  model_vector_t *p;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nml = model->layers[0].metal_layers.n_metal/2;
  int mtdmn  = model->config.PDN_multi_dom;
  double Rp  = model->c4->pad_r;
  double Lp  = model->c4->pad_l;
  double Rsp = model->config.PDN_pkg_sR;
  double Lsp = model->config.PDN_pkg_sL;
  double Rpp = model->config.PDN_pkg_pR;
  double Lpp = model->config.PDN_pkg_pL;
  double Cpp = model->config.PDN_pkg_C;
  double delta_t = 1/(model->config.proc_clock_freq*model->config.PDN_step_percycle);
  double Cg;
  metal_gridRL_t *g = model->layers[0].metal_layers.gridRL;

  /* Intermediate values */
  double iA = 1 / (nvp*Lsp/Lp + 1);
  double iB = 1 / (ngp*Lsp/Lp + 1);
  double iM = 1 / ((iA+iB)*Lsp + Lpp);
  double iN = iM * ((iA+iB)*Rsp + Rpp);
  double iP = Lsp / Lp;
  double iQ = Rsp / Lp;
  double iS = Rsp - iP*Rp;
  double iT = Rp / Lp;

  int nbr = 2*nr*nc - nr - nc;

  p = new_model_vector(model);

  /* map the block power/temp numbers to the grid	*/
  PDN_xlate_vector_b2g(model, power, p);
  copy_dvector(model->last_power, p->cuboid[0][0], nr*nc*nl);

  /* Initialize matrix A. */
  m = n = model->trans_matrix_dim;

  /* Num of non-zeros
   * Matrix for onchipI/onchipI:              5*(2*nr*nc-nr-nc)*mlayers*2
   * two diagnal block of zeros for grid V:   nr * nc + nr * nc
   * two dense block for pads:                nvp * nvp, ngp * ngp
   * pad<->grid side blocks:                  2 * nvp, 2 * ngp
   * row for Vc (extra 0 in diaganal):        1 + 1
   * row for Ic:                              2*nvp + 2*ngp + 2
   * col for Ic:                              nvp + ngp
   */
  nnz = 2*5*nbr*nml; 
  nnz += 2*nr*nc;
  nnz += nvp*nvp + ngp*ngp;
  nnz += 2*(nvp+ngp);
  if(model->config.PDN_pkgLC){
      nnz += 2*nvp + 2*ngp + 2;
      nnz += nvp + ngp + 2;
  }

  if ( !(cooV = doubleMalloc(nnz)) ) fatal("Malloc fails for cooV[].\n");
  if ( !(cooX = intMalloc(nnz)) ) fatal("Malloc fails for cooX[].\n");
  if ( !(cooY = intMalloc(nnz)) ) fatal("Malloc fails for cooY[].\n");

  // Create A
  if ( !(a = doubleMalloc(nnz)) ) fatal("Malloc fails for a[].\n");
  if ( !(asub = intMalloc(nnz)) ) fatal("Malloc fails for asub[].\n");
  if ( !(xa = intMalloc(n+1)) ) fatal("Malloc fails for xa[].\n");

  curidx = 0;

  /* Vdd pad dense block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++)
    for(j=0; j<nvp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(i == j)
          cooV[curidx] = -iQ - iT;
        else
          cooV[curidx] = -iQ;
        curidx++;
    }

  /* Gnd pad dense block */
  xoffset = nvp + nr*nc + nbr*nml;
  yoffset = nvp + nr*nc + nbr*nml;
  for(i=0; i<ngp; i++)
    for(j=0; j<ngp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(i == j)
          cooV[curidx] = -iQ - iT;
        else
          cooV[curidx] = -iQ;
        curidx++;
    }

  /* Vdd grid voltage diagonal block */
  xoffset = nvp;
  yoffset = nvp;
  for(i=0; i<(nr*nc); i++){
      cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = 0;
      curidx++;
  }

  /* Gnd grid voltage diagonal block */
  xoffset = nvp + nr*nc + nbr*nml + ngp;
  yoffset = nvp + nr*nc + nbr*nml + ngp;
  for(i=0; i<(nr*nc); i++){
      cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = 0;
      curidx++;
  }

  /* Vdd grid current diagonal block */
  for(l=0; l<nml; l++){
      xoffset = nvp + nr*nc + nbr*l;
      yoffset = nvp + nr*nc + nbr*l;
      for(i=0; i<nbr; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
          if((!mtdmn) || does_branch_exist(model, i)){
              if(i < ((nc-1)*nr))
                cooV[curidx] = -g[2*l].r/g[2*l].l;
              else{
                  cooV[curidx] = -g[2*l+1].r/g[2*l+1].l;
              }
          }
          else{
              cooV[curidx] = 0;
          }
          curidx++;
      }
  }

  /* Gnd grid current diagonal block */
  for(l=0; l<nml; l++){
      xoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
      yoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
      for(i=0; i<nbr; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
          if(i < ((nc-1)*nr))
            cooV[curidx] = -g[2*l].r/g[2*l].l;
          else{
              cooV[curidx] = -g[2*l+1].r/g[2*l+1].l;
          }
          curidx++;
      }
  }

  /* Vdd grid - current block */
  for(l=0; l<nml; l++){
      xoffset = nvp;
      yoffset = nvp + nr*nc + nbr*l;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;
            Cg = get_onchip_cap(model, grididx, 0);

            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+((nc-1)*nr)+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -nc)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, nc)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, -1)){
                    cooV[curidx] = 1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                if((!mtdmn) || compare_domain(model, grididx, 1)){
                    cooV[curidx] = -1/Cg; curidx++;}
                else{cooV[curidx] = 0; curidx++;}
            }
        }
  }

  /* Gnd grid - current block */
  for(l=0; l<nml; l++){
      xoffset = nvp + nr*nc + nbr*nml + ngp;
      yoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;
            Cg = get_onchip_cap(model, grididx, 0);

            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                cooV[curidx] = 1/Cg; curidx++;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                cooV[curidx] = 1/Cg; curidx++;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = -1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+((nc-1)*nr)+yoffset;
                cooV[curidx] = -1/Cg; curidx++;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                cooV[curidx] = 1/Cg; curidx++;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                cooV[curidx] = 1/Cg; curidx++;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                cooV[curidx] = 1/Cg; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                cooV[curidx] = -1/Cg; curidx++;
            }
        }
  }

  /* Vdd current - voltage block */
  for(l=0; l<nml; l++){
      xoffset = nvp + nr*nc + nbr*l;
      yoffset = nvp;
      for(i=0; i<nbr; i++){
          if(i<((nc-1)*nr)){
              tmp_j = i % (nc-1);
              tmp_i = (int) floor(i/(nc-1));
              grididx = tmp_i*nc + tmp_j;
              cooX[curidx] = i+xoffset; cooY[curidx] = grididx+yoffset;
              if((!mtdmn) || does_branch_exist(model, i)){
                  cooV[curidx] = 1/g[2*l].l; curidx++;}
              else{cooV[curidx] = 0; curidx++;}

              cooX[curidx] = i+xoffset; cooY[curidx] = grididx+1+yoffset;
              if((!mtdmn) || does_branch_exist(model, i)){
                  cooV[curidx] = -1/g[2*l].l; curidx++;}
              else{cooV[curidx] = 0; curidx++;}
          }
          else{
              cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+yoffset;
              if((!mtdmn) || does_branch_exist(model, i)){
                  cooV[curidx] = 1/g[2*l+1].l; curidx++;}
              else{cooV[curidx] = 0; curidx++;}

              cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+nc+yoffset;
              if((!mtdmn) || does_branch_exist(model, i)){
                  cooV[curidx] = -1/g[2*l+1].l; curidx++;}
              else{cooV[curidx] = 0; curidx++;}
          }
      }
  }

  /* Gnd current - voltage block */
  for(l=0; l<nml; l++){
      xoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
      yoffset = nvp + nr*nc + nbr*nml + ngp;
      for(i=0; i<nbr; i++){
          if(i<((nc-1)*nr)){
              tmp_j = i % (nc-1);
              tmp_i = (int) floor(i/(nc-1));
              grididx = tmp_i*nc + tmp_j;
              cooX[curidx] = i+xoffset; cooY[curidx] = grididx+yoffset;
              cooV[curidx] = 1/g[2*l].l; curidx++;

              cooX[curidx] = i+xoffset; cooY[curidx] = grididx+1+yoffset;
              cooV[curidx] = -1/g[2*l].l; curidx++;
          }
          else{
              cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+yoffset;
              cooV[curidx] = 1/g[2*l+1].l; curidx++;

              cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+nc+yoffset;
              cooV[curidx] = -1/g[2*l+1].l; curidx++;
          }
      }
  }

  /* VDD pad two sparse block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
      if(-1 == grididx)
        fatal("Cannot find grid index for VDD pad\n");
      Cg = get_onchip_cap(model, grididx, 0);

      cooX[curidx] = grididx+nvp+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = 1/Cg; curidx++;

      cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+yoffset;
      cooV[curidx] = -1/Lp; curidx++;
  }

  /* GND pad two sparse block */
  xoffset = nvp + nr*nc + nbr*nml;
  yoffset = nvp + nr*nc + nbr*nml;
  for(i=0; i<ngp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
      if(-1 == grididx)
        fatal("Cannot find grid index for GND pad\n");
      Cg = get_onchip_cap(model, grididx, 0);

      cooX[curidx] = grididx+ngp+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = -1/Cg; curidx++;

      cooX[curidx] = i+xoffset; cooY[curidx] = grididx+ngp+yoffset;
      cooV[curidx] = 1/Lp; curidx++;
  }

  if(model->config.PDN_pkgLC){
      /* Bottom row, VDD part */
      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = -iM*iA*iS; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+yoffset;
          cooV[curidx] = iM*iA*iP; curidx++;
      }
      /* Bottom row, GND part */
      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = nvp + nr*nc + nbr*nml;
      for(i=0; i<ngp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = -iM*iB*iS; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+ngp+yoffset;
          cooV[curidx] = -iM*iB*iP; curidx++;
      }
      /* Right col, VDD part */
      xoffset = 0;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      for(i=0; i<nvp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iQ; curidx++;
      }
      /* Right col, GND part */
      xoffset = nvp + nr*nc + nbr*nml;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      for(i=0; i<ngp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iQ; curidx++;
      }

      /* Individual points */
      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = -iN; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = -iM; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 1/Cpp; curidx++;
  }

  if(curidx < nnz)
    fatal("COO Matrix has less elements than nnz!\n");
  else if(curidx > nnz)
    fatal("COO Matrix has more elements than nnz!\n");

  coo2csc(n, nnz, cooX, cooY, cooV, asub, xa, a);

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

  if(model->config.PDN_pkgLC){
      // Create B
      if ( !(b = doubleMalloc(nnz)) ) fatal("Malloc fails for b[].\n");
      if ( !(bsub = intMalloc(nnz)) ) fatal("Malloc fails for bsub[].\n");
      if ( !(xb = intMalloc(n+1)) ) fatal("Malloc fails for xb[].\n");

      curidx = 0;
      /* Vdd pad dense block */
      xoffset = 0;
      yoffset = 0;
      for(i=0; i<nvp; i++)
        for(j=0; j<nvp; j++){
            cooX[curidx] = i + xoffset;
            cooY[curidx] = j + yoffset;
            cooV[curidx] = -iP;
            curidx++;
        }

      /* Gnd pad dense block */
      xoffset = nvp + nr*nc + nbr*nml;
      yoffset = nvp + nr*nc + nbr*nml;
      for(i=0; i<ngp; i++)
        for(j=0; j<ngp; j++){
            cooX[curidx] = i + xoffset;
            cooY[curidx] = j + yoffset;
            cooV[curidx] = -iP;
            curidx++;
        }

      /* Vdd grid voltage diagonal block */
      xoffset = nvp;
      yoffset = nvp;
      for(i=0; i<(nr*nc); i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0;
          curidx++;
      }

      /* Gnd grid voltage diagonal block */
      xoffset = nvp + nr*nc + nbr*nml + ngp;
      yoffset = nvp + nr*nc + nbr*nml + ngp;
      for(i=0; i<(nr*nc); i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0;
          curidx++;
      }

      /* Vdd grid current diagonal block */
      for(l=0; l<nml; l++){
          xoffset = nvp + nr*nc + nbr*l;
          yoffset = nvp + nr*nc + nbr*l;
          for(i=0; i<nbr; i++){
              cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
              cooV[curidx] = 0;
              curidx++;
          }
      }

      /* Gnd grid current diagonal block */
      for(l=0; l<nml; l++){
          xoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
          yoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
          for(i=0; i<nbr; i++){
              cooX[curidx] = i+xoffset; cooY[curidx] = i+yoffset;
              cooV[curidx] = 0;
              curidx++;
          }
      }

      /* Vdd grid - current block */
      for(l=0; l<nml; l++){
          xoffset = nvp;
          yoffset = nvp + nr*nc + nbr*l;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;

                if(0 == grididx){//top left corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((nc-1) == grididx ){//top right corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((nr*nc-nc) == grididx){//bottom left corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((nr*nc-1) == grididx){//bottom right corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((grididx > 0) && (grididx < nc-1)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+((nc-1)*nr)+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if(0 == (grididx%nc)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if(0 == ((grididx+1)%nc)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else{
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
            }
      }

      /* Gnd grid - current block */
      for(l=0; l<nml; l++){
          xoffset = nvp + nr*nc + nbr*nml + ngp;
          yoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;

                if(0 == grididx){//top left corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((nc-1) == grididx ){//top right corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((nr*nc-nc) == grididx){//bottom left corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((nr*nc-1) == grididx){//bottom right corner
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((grididx > 0) && (grididx < nc-1)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+((nc-1)*nr)+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if(0 == (grididx%nc)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else if(0 == ((grididx+1)%nc)){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
                else{
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+(i-1)*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = (nc-1)*nr+i*nc+j+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j-1+yoffset;
                    cooV[curidx] = 0; curidx++;

                    cooX[curidx] = grididx+xoffset; cooY[curidx] = i*(nc-1)+j+yoffset;
                    cooV[curidx] = 0; curidx++;
                }
            }
      }

      /* Vdd current - voltage block */
      for(l=0; l<nml; l++){
          xoffset = nvp + nr*nc + nbr*l;
          yoffset = nvp;
          for(i=0; i<nbr; i++){
              if(i<((nc-1)*nr)){
                  tmp_j = i % (nc-1);
                  tmp_i = (int) floor(i/(nc-1));
                  grididx = tmp_i*nc + tmp_j;
                  cooX[curidx] = i+xoffset; cooY[curidx] = grididx+yoffset;
                  cooV[curidx] = 0; curidx++;

                  cooX[curidx] = i+xoffset; cooY[curidx] = grididx+1+yoffset;
                  cooV[curidx] = 0; curidx++;
              }
              else{
                  cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+yoffset;
                  cooV[curidx] = 0; curidx++;

                  cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+nc+yoffset;
                  cooV[curidx] = 0; curidx++;
              }
          }
      }

      /* Gnd current - voltage block */
      for(l=0; l<nml; l++){
          xoffset = nvp + nr*nc + nbr*nml + ngp + nr*nc + nbr*l;
          yoffset = nvp + nr*nc + nbr*nml + ngp;
          for(i=0; i<nbr; i++){
              if(i<((nc-1)*nr)){
                  tmp_j = i % (nc-1);
                  tmp_i = (int) floor(i/(nc-1));
                  grididx = tmp_i*nc + tmp_j;
                  cooX[curidx] = i+xoffset; cooY[curidx] = grididx+yoffset;
                  cooV[curidx] = 0; curidx++;

                  cooX[curidx] = i+xoffset; cooY[curidx] = grididx+1+yoffset;
                  cooV[curidx] = 0; curidx++;
              }
              else{
                  cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+yoffset;
                  cooV[curidx] = 0; curidx++;

                  cooX[curidx] = i+xoffset; cooY[curidx] = i-nr*(nc-1)+nc+yoffset;
                  cooV[curidx] = 0; curidx++;
              }
          }
      }

      /* VDD pad two sparse block */
      xoffset = 0;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = grididx+nvp+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+yoffset;
          cooV[curidx] = 0; curidx++;
      }

      /* GND pad two sparse block */
      xoffset = nvp + nr*nc + nbr*nml;
      yoffset = nvp + nr*nc + nbr*nml;
      for(i=0; i<ngp; i++){
          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for GND pad\n");

          cooX[curidx] = grididx+ngp+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+ngp+yoffset;
          cooV[curidx] = 0; curidx++;
      }

      /* Bottom row, VDD part */
      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+yoffset;
          cooV[curidx] = 0; curidx++;
      }
      /* Bottom row, GND part */
      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = nvp + nr*nc + nbr*nml;
      for(i=0; i<ngp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+ngp+yoffset;
          cooV[curidx] = 0; curidx++;
      }
      /* Right col, VDD part */
      xoffset = 0;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      for(i=0; i<nvp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iP; curidx++;
      }
      /* Right col, GND part */
      xoffset = nvp + nr*nc + nbr*nml;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      for(i=0; i<ngp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iP; curidx++;
      }

      /* Individual points */
      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml;
      yoffset = nvp + 2*nr*nc + ngp + 2*nbr*nml + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      if(curidx < nnz)
        fatal("COO Matrix has less elements than nnz!\n");
      else if(curidx > nnz)
        fatal("COO Matrix has more elements than nnz!\n");

      coo2csc(n, nnz, cooX, cooY, cooV, bsub, xb, b);

      dCreate_CompCol_Matrix(&B, m, n, nnz, b, bsub, xb, SLU_NC, SLU_D, SLU_GE);
  }

  // E - B - delta_t*A/2 for pkgLC
  // E - delta_t*A/2 for no pkgLC
  SparseMatrix_mul_SingleNum(A, -delta_t/2);
  if(model->config.PDN_pkgLC)
    SparseMatrix_add(A, &B, -1);
  SparseMatrix_add_Iden(A, 1);

  if ( !(rhs = doubleMalloc(m)) ) fatal("Malloc fails for rhs[].\n");
  for(i=0; i<m; ++i) rhs[i] = 1;
  dCreate_Dense_Matrix(&Z, m, 1, rhs, m, SLU_DN, SLU_D, SLU_GE);

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = 0.01;
  options.SymmetricMode = NO;
  options.Equil = YES;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve the linear system. */
  dgssv(&options, A, perm_c, perm_r, L, U, &Z, &stat, &info);

  // Create RHS
  // E - B + (delta_t/2)*J for pkgLC
  // E + (delta_t/2)*J for no pkgLC
  SparseMatrix_mul_SingleNum(A, -1);
  if(model->config.PDN_pkgLC)
    SparseMatrix_add(A, &B, -2);
  SparseMatrix_add_Iden(A, 2);

  /* De-allocate storage */
  SUPERLU_FREE (rhs);
  Destroy_SuperMatrix_Store(&Z);
  if(model->config.PDN_pkgLC)
    Destroy_CompCol_Matrix(&B);
  StatFree(&stat);
  free_model_vector(p);

  free(cooV);
  free(cooX);
  free(cooY);

}

void trans_SLU_init_nogridL(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r)
{
  SuperMatrix B, Z;
  double   *a, *b, *rhs;
  int      *asub, *xa, *bsub, *xb;
  double   *cooV;
  int      *cooX, *cooY;

  double   dia_val;
  int      info, i, j, m, n, nnz;
  int      direc;
  int      xoffset, yoffset;
  int      grididx, curidx;
  superlu_options_t options;
  SuperLUStat_t stat;
  model_vector_t *p;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int mtdmn  = model->config.PDN_multi_dom;
  double Rp  = model->c4->pad_r;
  double Lp  = model->c4->pad_l;
  double Rx  = 0; double Ry  = 0;
  double Rsp = model->config.PDN_pkg_sR;
  double Lsp = model->config.PDN_pkg_sL;
  double Rpp = model->config.PDN_pkg_pR;
  double Lpp = model->config.PDN_pkg_pL;
  double Cpp = model->config.PDN_pkg_C;
  double delta_t = 1/(model->config.proc_clock_freq*model->config.PDN_step_percycle);
  double Cg;

  //calculate grid R
  for(i=0; i<model->layers[0].metal_layers.n_metal; i++){
      direc = model->layers[0].metal_layers.geo[i].direc;
      if(MLCF_X == direc)
        Rx += 1/model->layers[0].metal_layers.gridRL[i].r;
      else
        Ry += 1/model->layers[0].metal_layers.gridRL[i].r;
  }
  Rx = 1/Rx;
  Ry = 1/Ry;

  /* Intermediate values */
  double iA = 1 / (nvp*Lsp/Lp + 1);
  double iB = 1 / (ngp*Lsp/Lp + 1);
  double iM = 1 / ((iA+iB)*Lsp + Lpp);
  double iN = iM * ((iA+iB)*Rsp + Rpp);
  double iP = Lsp / Lp;
  double iQ = Rsp / Lp;
  double iS = Rsp - iP*Rp;
  double iT = Rp / Lp;

  p = new_model_vector(model);

  /* map the block power/temp numbers to the grid	*/
  PDN_xlate_vector_b2g(model, power, p);
  copy_dvector(model->last_power, p->cuboid[0][0], nr*nc*nl);

  /* Initialize matrix A. */
  m = n = model->trans_matrix_dim;

  /* Num of non-zeros
   * Five diagonal of both vdd and gnd plane: c*r, c*(r-1), (c-1)*r, c*(r-1), (c-1)*r
   * two dense block for pads:                nvp * nvp, ngp * ngp
   * pad<->grid side blocks:                  2 * nvp, 2 * ngp
   * row for Vc (extra 0 in diaganal):        1 + 1
   * row for Ic:                              2*nvp + 2*ngp + 2
   * col for Ic:                              nvp + ngp
   */
  nnz = 2*(5*nr*nc - 2*(nr+nc)); 
  nnz += nvp*nvp + ngp*ngp;
  nnz += 2*(nvp+ngp);
  if(model->config.PDN_pkgLC){
      nnz += 2*nvp + 2*ngp + 2;
      nnz += nvp + ngp + 2;
  }

  if ( !(cooV = doubleMalloc(nnz)) ) fatal("Malloc fails for cooV[].\n");
  if ( !(cooX = intMalloc(nnz)) ) fatal("Malloc fails for cooX[].\n");
  if ( !(cooY = intMalloc(nnz)) ) fatal("Malloc fails for cooY[].\n");

  // Create A
  if ( !(a = doubleMalloc(nnz)) ) fatal("Malloc fails for a[].\n");
  if ( !(asub = intMalloc(nnz)) ) fatal("Malloc fails for asub[].\n");
  if ( !(xa = intMalloc(n+1)) ) fatal("Malloc fails for xa[].\n");

  curidx = 0;

  /* Vdd pad dense block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++)
    for(j=0; j<nvp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(i == j)
          cooV[curidx] = -iQ - iT;
        else
          cooV[curidx] = -iQ;
        curidx++;
    }

  /* Gnd pad dense block */
  xoffset = nvp + nr*nc;
  yoffset = nvp + nr*nc;
  for(i=0; i<ngp; i++)
    for(j=0; j<ngp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(i == j)
          cooV[curidx] = -iQ - iT;
        else
          cooV[curidx] = -iQ;
        curidx++;
    }

  /* Vdd grid diagonal block */
  xoffset = nvp;
  yoffset = nvp;
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++){
        grididx = i*nc + j;
        Cg = get_onchip_cap(model, grididx, 0);

        dia_val = 0;
        if(0 == grididx){//top left corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, 1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((nc-1) == grididx ){//top right corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((nr*nc-nc) == grididx){//bottom left corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, 1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((nr*nc-1) == grididx){//bottom right corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((grididx > 0) && (grididx < nc-1)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, 1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, 1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if(0 == (grididx%nc)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, 1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if(0 == ((grididx+1)%nc)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else{
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, nc)){
                cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, -1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            if((!mtdmn) || compare_domain(model, grididx, 1)){
                cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;}
            else{cooV[curidx] = 0; curidx++;}

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
    }

  /* Gnd grid diagonal block */
  xoffset = nvp + nr*nc + ngp;
  yoffset = nvp + nr*nc + ngp;
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++){
        grididx = i*nc + j;
        Cg = get_onchip_cap(model, grididx, 0);

        dia_val = 0;

        if(0 == grididx){//top left corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((nc-1) == grididx ){//top right corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((nr*nc-nc) == grididx){//bottom left corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((nr*nc-1) == grididx){//bottom right corner
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((grididx > 0) && (grididx < nc-1)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if(0 == (grididx%nc)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else if(0 == ((grididx+1)%nc)){
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
        else{
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
            cooV[curidx] = 1/(Cg*Ry); curidx++; dia_val+=1/Ry;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
            cooV[curidx] = 1/(Cg*Rx); curidx++; dia_val+=1/Rx;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cg; curidx++;
        }
    }

  /* VDD pad two sparse block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
      if(-1 == grididx)
        fatal("Cannot find grid index for VDD pad\n");
      Cg = get_onchip_cap(model, grididx, 0);

      cooX[curidx] = grididx+nvp+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = 1/Cg; curidx++;

      cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+yoffset;
      cooV[curidx] = -1/Lp; curidx++;
  }

  /* GND pad two sparse block */
  xoffset = nvp + nr*nc;
  yoffset = nvp + nr*nc;
  for(i=0; i<ngp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
      if(-1 == grididx)
        fatal("Cannot find grid index for GND pad\n");
      Cg = get_onchip_cap(model, grididx, 0);

      cooX[curidx] = grididx+ngp+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = -1/Cg; curidx++;

      cooX[curidx] = i+xoffset; cooY[curidx] = grididx+ngp+yoffset;
      cooV[curidx] = 1/Lp; curidx++;
  }

  if(model->config.PDN_pkgLC){
      /* Bottom row, VDD part */
      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = -iM*iA*iS; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+yoffset;
          cooV[curidx] = iM*iA*iP; curidx++;
      }
      /* Bottom row, GND part */
      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = nvp + nr*nc;
      for(i=0; i<ngp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = -iM*iB*iS; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+ngp+yoffset;
          cooV[curidx] = -iM*iB*iP; curidx++;
      }
      /* Right col, VDD part */
      xoffset = 0;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      for(i=0; i<nvp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iQ; curidx++;
      }
      /* Right col, GND part */
      xoffset = nvp + nr*nc;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      for(i=0; i<ngp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iQ; curidx++;
      }

      /* Individual points */
      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = -iN; curidx++;

      xoffset = nvp + 2*nr*nc + ngp;
      yoffset = nvp + 2*nr*nc + ngp;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = nvp + 2*nr*nc + ngp;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = -iM; curidx++;

      xoffset = nvp + 2*nr*nc + ngp;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 1/Cpp; curidx++;
  }

  if(curidx != nnz)
    fatal("COO Matrix has less elements than nnz!\n");

  coo2csc(n, nnz, cooX, cooY, cooV, asub, xa, a);

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

  if(model->config.PDN_pkgLC){
      // Create B
      if ( !(b = doubleMalloc(nnz)) ) fatal("Malloc fails for b[].\n");
      if ( !(bsub = intMalloc(nnz)) ) fatal("Malloc fails for bsub[].\n");
      if ( !(xb = intMalloc(n+1)) ) fatal("Malloc fails for xb[].\n");

      curidx = 0;
      /* Vdd pad dense block */
      xoffset = 0;
      yoffset = 0;
      for(i=0; i<nvp; i++)
        for(j=0; j<nvp; j++){
            cooX[curidx] = i + xoffset;
            cooY[curidx] = j + yoffset;
            cooV[curidx] = -iP;
            curidx++;
        }

      /* Gnd pad dense block */
      xoffset = nvp + nr*nc;
      yoffset = nvp + nr*nc;
      for(i=0; i<ngp; i++)
        for(j=0; j<ngp; j++){
            cooX[curidx] = i + xoffset;
            cooY[curidx] = j + yoffset;
            cooV[curidx] = -iP;
            curidx++;
        }

      /* Vdd grid diagonal block */
      xoffset = nvp;
      yoffset = nvp;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
        }

      /* Gnd grid diagonal block */
      xoffset = nvp + nr*nc + ngp;
      yoffset = nvp + nr*nc + ngp;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
        }

      /* VDD pad two sparse block */
      xoffset = 0;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = grididx+nvp+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+yoffset;
          cooV[curidx] = 0; curidx++;
      }

      /* GND pad two sparse block */
      xoffset = nvp + nr*nc;
      yoffset = nvp + nr*nc;
      for(i=0; i<ngp; i++){
          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for GND pad\n");

          cooX[curidx] = grididx+ngp+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+ngp+yoffset;
          cooV[curidx] = 0; curidx++;
      }

      /* Bottom row, VDD part */
      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+yoffset;
          cooV[curidx] = 0; curidx++;
      }
      /* Bottom row, GND part */
      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = nvp + nr*nc;
      for(i=0; i<ngp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+ngp+yoffset;
          cooV[curidx] = 0; curidx++;
      }
      /* Right col, VDD part */
      xoffset = 0;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      for(i=0; i<nvp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iP; curidx++;
      }
      /* Right col, GND part */
      xoffset = nvp + nr*nc;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      for(i=0; i<ngp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iP; curidx++;
      }

      /* Individual points */
      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp;
      yoffset = nvp + 2*nr*nc + ngp;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp + 1;
      yoffset = nvp + 2*nr*nc + ngp;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = nvp + 2*nr*nc + ngp;
      yoffset = nvp + 2*nr*nc + ngp + 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      if(curidx != nnz)
        fatal("COO Matrix has less elements than nnz!\n");

      coo2csc(n, nnz, cooX, cooY, cooV, bsub, xb, b);

      dCreate_CompCol_Matrix(&B, m, n, nnz, b, bsub, xb, SLU_NC, SLU_D, SLU_GE);
  }

  // E - B - delta_t*A/2 for pkgLC
  // E - delta_t*A/2 for no pkgLC
  SparseMatrix_mul_SingleNum(A, -delta_t/2);
  if(model->config.PDN_pkgLC)
    SparseMatrix_add(A, &B, -1);
  SparseMatrix_add_Iden(A, 1);

  if ( !(rhs = doubleMalloc(m)) ) fatal("Malloc fails for rhs[].\n");
  for(i=0; i<m; ++i) rhs[i] = 1;
  dCreate_Dense_Matrix(&Z, m, 1, rhs, m, SLU_DN, SLU_D, SLU_GE);

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = 0.01;
  options.SymmetricMode = NO;
  options.Equil = YES;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve the linear system. */
  dgssv(&options, A, perm_c, perm_r, L, U, &Z, &stat, &info);

  // Create RHS
  // E + (delta_t/2)*J
  SparseMatrix_mul_SingleNum(A, -1);
  if(model->config.PDN_pkgLC)
    SparseMatrix_add(A, &B, -2);
  SparseMatrix_add_Iden(A, 2);

  /* De-allocate storage */
  SUPERLU_FREE (rhs);
  Destroy_SuperMatrix_Store(&Z);
  if(model->config.PDN_pkgLC)
    Destroy_CompCol_Matrix(&B);
  StatFree(&stat);
  free_model_vector(p);

  free(cooV);
  free(cooX);
  free(cooY);
}

void trans_matrix_build_3D(model_t *model, double *power, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, int **perm_c, int **perm_r)
{
  SuperMatrix B, Z;
  double   *a, *b, *rhs;
  int      *asub, *xa, *bsub, *xb;
  double   *cooV;
  int      *cooX, *cooY;

  int      i, j, l, k, m, n, nnz;
  int      tmp_i, tmp_j;
  int      info, direc;
  int      xoffset, yoffset;
  int      grididx, curidx;
  int      ntsv, nIVR;
  int      **ivr_loc;
  double   dia_val;
  superlu_options_t options;
  SuperLUStat_t stat;
  model_vector_t *p;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nvtsv = model->layers[0].tsv.num_vdd;
  int ngtsv = model->layers[0].tsv.num_gnd;
  int tsv_count, ivr_idx;
  int vs = model->config.v_stacking;
  double Rpv  = model->c4->pad_r;
  double Lpv  = model->c4->pad_l;
  double Rpg  = model->c4->pad_r;
  double Lpg  = model->c4->pad_l;
  double Rsp = model->config.PDN_pkg_sR;
  double Lsp = model->config.PDN_pkg_sL;
  double Rpp = model->config.PDN_pkg_pR;
  double Lpp = model->config.PDN_pkg_pL;
  double Cpp = model->config.PDN_pkg_C;
  double delta_t = 1/(model->config.proc_clock_freq*model->config.PDN_step_percycle);
  double vdd = model->config.vdd;
  double Rx, Ry, Rload;
  double Cg, Cm;
  double Rtsv, Ltsv;
  double Ront, Ronb, Civr;

  /* Add TSV RL into pad RL */
  // This can be ignored if we can use larger TSVs to connect top-layer pads
  //if(vs){
  //    double Rtsv_sum = 0;
  //    double Ltsv_sum = 0;
  //    for(l=0; l<nl-1; l++){
  //        Rtsv_sum += model->layers[l].tsv.r;
  //        Ltsv_sum += model->layers[l].tsv.l;
  //    }
  //    Rpv += Rtsv_sum;
  //    Lpv += Ltsv_sum;
  //}

  /* Intermediate values */
  double iA = 1 / (nvp*Lsp/Lpv + 1);
  double iB = 1 / (ngp*Lsp/Lpg + 1);
  double iM = 1 / ((iA+iB)*Lsp + Lpp);
  double iN = iM * ((iA+iB)*Rsp + Rpp);
  double iPv = Lsp / Lpv;
  double iQv = Rsp / Lpv;
  double iSv = Rsp - iPv*Rpv;
  double iTv = Rpv / Lpv;
  double iPg = Lsp / Lpg;
  double iQg = Rsp / Lpg;
  double iSg = Rsp - iPg*Rpg;
  double iTg = Rpg / Lpg;

  /* map the block power numbers to the grid	*/
  p = new_model_vector(model);
  PDN_xlate_vector_b2g(model, power, p);

  if(vs){
      ntsv = ngtsv*(nl-1);
      nIVR = model->sc_converter->num_IVR;
      Ront = model->sc_converter->Ron_top;
      Ronb = model->sc_converter->Ron_bottom;
      Civr = model->sc_converter->c;
      ivr_loc = model->sc_converter->loc;
  }
  else{
      ntsv = (nvtsv+ngtsv)*(nl-1);
      nIVR = 0;
  }

  /* matrix dim */
  m = n = model->trans_matrix_dim;

  if( !(*perm_r = (int *) calloc(m, sizeof(int))) ) fatal("Malloc fails for perm_r[].\n");
  if( !(*perm_c = (int *) calloc(m, sizeof(int))) ) fatal("Malloc fails for perm_c[].\n");

  /* Num of non-zeros
   * Five diagonal of both vdd and gnd plane: c*r, c*(r-1), (c-1)*r, c*(r-1), (c-1)*r
   * single diagonal for Vdd-Gnd decap:       2*nr*nc
   * pad<->grid side blocks:                  2*nvp, 2*ngp
   * two dense block for pads:                nvp*nvp, ngp*ngp
   * two side blocks for IVR:                 2*8*(nl-1)*nIVR
   * diaganal for IVR:                        4*(nl-1)*nIVR
   * row for Vc (extra 0 in diaganal):        1+1
   * row for Ic:                              2*nvp + 2*ngp + 2
   * col for Ic:                              nvp + ngp
   */
  nnz = 2*(5*nr*nc - 2*(nr+nc)); 
  nnz += 2*nr*nc;
  nnz *= nl;
  nnz += 2*(nvp+ngp);
  nnz += nvp*nvp + ngp*ngp;
  nnz += 5*ntsv;
  nnz += (2*8*(nl-1)*nIVR + 4*(nl-1)*nIVR);
  if(model->config.PDN_pkgLC){
      nnz += 2*nvp + 2*ngp + 2;
      nnz += nvp + ngp + 2;
  }

  if ( !(cooV = doubleMalloc(nnz)) ) fatal("Malloc fails for cooV[].\n");
  if ( !(cooX = intMalloc(nnz)) ) fatal("Malloc fails for cooX[].\n");
  if ( !(cooY = intMalloc(nnz)) ) fatal("Malloc fails for cooY[].\n");

  // Create A
  if ( !(a = doubleMalloc(nnz)) ) fatal("Malloc fails for a[].\n");
  if ( !(asub = intMalloc(nnz)) ) fatal("Malloc fails for asub[].\n");
  if ( !(xa = intMalloc(n+1)) ) fatal("Malloc fails for xa[].\n");

  curidx = 0;

  /* Vdd pad dense block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++)
    for(j=0; j<nvp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(i == j)
          cooV[curidx] = -iQv - iTv;
        else
          cooV[curidx] = -iQv;
        curidx++;
    }

  /* Gnd pad dense block */
  xoffset = nvp;
  yoffset = nvp;
  for(i=0; i<ngp; i++)
    for(j=0; j<ngp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(i == j)
          cooV[curidx] = -iQg - iTg;
        else
          cooV[curidx] = -iQg;
        curidx++;
    }

  /* VDD pad two sparse block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
      if(-1 == grididx)
        fatal("Cannot find grid index for VDD pad\n");

      if(vs){
          tmp_i = floor(grididx/nc);
          tmp_j = grididx % nc;
          Cg = get_onchip_cap(model, grididx, nl-1);
          if((IVRLOC == ivr_loc[tmp_i][tmp_j]) && (nl > 1))
            Cm = Cg + 2*Civr;
          else
            Cm = Cg;

          cooX[curidx] = grididx+nvp+ngp+(nl-1)*2*nr*nc+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 1/Cm; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+ngp+(nl-1)*2*nr*nc+yoffset;
          cooV[curidx] = -1/Lpv; curidx++;
      }
      else{
          Cg = get_onchip_cap(model, grididx, 0);
          cooX[curidx] = grididx+nvp+ngp+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 1/Cg; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+ngp+yoffset;
          cooV[curidx] = -1/Lpv; curidx++;
      }
  }

  /* GND pad two sparse block */
  xoffset = nvp;
  yoffset = nvp;
  for(i=0; i<ngp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
      if(-1 == grididx)
        fatal("Cannot find grid index for GND pad\n");

      tmp_i = floor(grididx/nc);
      tmp_j = grididx % nc;
      Cg = get_onchip_cap(model, grididx, 0);
      if(vs && (IVRLOC == ivr_loc[tmp_i][tmp_j]))
        Cm = Cg + 2*Civr;
      else
        Cm = Cg;

      cooX[curidx] = grididx+ngp+nr*nc+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = -1/Cm; curidx++;

      cooX[curidx] = i+xoffset; cooY[curidx] = grididx+ngp+nr*nc+yoffset;
      cooV[curidx] = 1/Lpg; curidx++;
  }

  /* Onchip nodes */
  for(l=0; l<nl; l++){
      //calculate grid R
      Rx = 0; Ry = 0;
      for(i=0; i<model->layers[l].metal_layers.n_metal; i++){
          direc = model->layers[l].metal_layers.geo[i].direc;
          if(MLCF_X == direc)
            Rx += 1/model->layers[l].metal_layers.gridRL[i].r;
          else
            Ry += 1/model->layers[l].metal_layers.gridRL[i].r;
      }
      Rx = 1/Rx;
      Ry = 1/Ry;

      /* Vdd grid diagonal block */
      xoffset = nvp + ngp + l*2*nr*nc;
      yoffset = nvp + ngp + l*2*nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;
            Cg = get_onchip_cap(model, grididx, l);
            if(vs && (IVRLOC == ivr_loc[i][j]) && (l>0))
              Cm = Cg + 2*Civr;
            else
              Cm = Cg;

            if (p->cuboid[l][i][j] > 0)
              Rload = vdd*vdd/p->cuboid[l][i][j];
            else
              Rload = LARGEINT;

            dia_val = 0;
            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;
            }

            /* Vdd-Gnd grid diagonal for load resisotrs */
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
            cooV[curidx] = 1/(Cm*Rload); curidx++; dia_val += 1/Rload;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cm; curidx++;
        }

      /* Gnd grid diagonal block */
      xoffset = nvp + ngp + l*2*nr*nc + nr*nc;
      yoffset = nvp + ngp + l*2*nr*nc + nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;
            Cg = get_onchip_cap(model, grididx, l);
            if(vs && (IVRLOC == ivr_loc[i][j]) && (l<nl-1))
              Cm = Cg + 2*Civr;
            else
              Cm = Cg;

            if (p->cuboid[l][i][j] > 0)
              Rload = vdd*vdd/p->cuboid[l][i][j];
            else
              Rload = LARGEINT;

            dia_val = 0;
            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 1/(Cm*Ry); curidx++; dia_val+=1/Ry;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 1/(Cm*Rx); curidx++; dia_val+=1/Rx;
            }

            /* Gnd-Vdd grid diagonal for load resisotrs */
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
            cooV[curidx] = 1/(Cm*Rload); curidx++; dia_val += 1/Rload;

            if(vs && (IVRLOC == ivr_loc[i][j]) && (l>0))
              dia_val += 2/Ront + 2/Ronb;

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = -dia_val/Cm; curidx++;
        }
  }

  if(vs){
      /* Onchip nodes -> TSVs block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          xoffset = nvp + ngp + l*2*nr*nc;
          yoffset = nvp + ngp + nl*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    Cg = get_onchip_cap(model, grididx, l);
                    if((IVRLOC == ivr_loc[i][j]) && (l>0))
                      Cm = Cg + 2*Civr;
                    else
                      Cm = Cg;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 1/Cm; curidx++;

                    Cg = get_onchip_cap(model, grididx, l+1);
                    if((IVRLOC == ivr_loc[i][j]) && (l<nl-2))
                      Cm = Cg + 2*Civr;
                    else
                      Cm = Cg;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = -1/Cm; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSVs -> Onchip nodes block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          Ltsv = model->layers[l].tsv.l;

          xoffset = nvp + ngp + nl*2*nr*nc;
          yoffset = nvp + ngp + l*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+yoffset;
                    cooV[curidx] = -1/Ltsv; curidx++;
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 1/Ltsv; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSV diagonal block */
      xoffset = nvp + ngp + nl*2*nr*nc;
      yoffset = nvp + ngp + nl*2*nr*nc;
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          Rtsv = model->layers[l].tsv.r;
          Ltsv = model->layers[l].tsv.l;

          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = -Rtsv/Ltsv; curidx++;
                    tsv_count++;
                }
            }
      }
  }
  else{
      /* Onchip nodes -> TSVs block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          xoffset = nvp + ngp + l*2*nr*nc;
          yoffset = nvp + ngp + nl*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == VDDTSV){
                    Cg = get_onchip_cap(model, grididx, l);
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = -1/Cg; curidx++;
                    Cg = get_onchip_cap(model, grididx, l+1);
                    cooX[curidx] = grididx+2*nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 1/Cg; curidx++;
                    tsv_count++;
                }
            }
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    Cg = get_onchip_cap(model, grididx, l);
                    cooX[curidx] = grididx+nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = -1/Cg; curidx++;
                    Cg = get_onchip_cap(model, grididx, l+1);
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 1/Cg; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSVs -> Onchip nodes block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          Ltsv = model->layers[l].tsv.l;

          xoffset = nvp + ngp + nl*2*nr*nc;
          yoffset = nvp + ngp + l*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == VDDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+yoffset;
                    cooV[curidx] = 1/Ltsv; curidx++;
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+2*nr*nc+yoffset;
                    cooV[curidx] = -1/Ltsv; curidx++;
                    tsv_count++;
                }
            }
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = 1/Ltsv; curidx++;
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = -1/Ltsv; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSV diagonal block */
      xoffset = nvp + ngp + nl*2*nr*nc;
      yoffset = nvp + ngp + nl*2*nr*nc;
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          Rtsv = model->layers[l].tsv.r;
          Ltsv = model->layers[l].tsv.l;

          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                if(model->layers[l].tsv.loc[i][j] == VDDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = -Rtsv/Ltsv; curidx++;
                    tsv_count++;
                }
            }
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = -Rtsv/Ltsv; curidx++;
                    tsv_count++;
                }
            }
      }
  }

  /* IVR block */
  if(nIVR > 0){
      /* Onchip nodes -> IVR block */
      ivr_idx = 0;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            if(IVRLOC == ivr_loc[i][j]){
                for(l=0; l<nl-1; l++){
                    grididx = i*nc + j;
                    xoffset = nvp + ngp + l*2*nr*nc;
                    yoffset = nvp + ngp + nl*2*nr*nc + ntsv;

                    Cg = get_onchip_cap(model, grididx, l+1);
                    if(l<nl-2)
                      Cm = Cg + 2*Civr;
                    else
                      Cm = Cg;

                    cooX[curidx] = grididx+2*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 1/(Cm*Ront); curidx++;
                    ivr_idx++;
                    cooX[curidx] = grididx+2*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 1/(Cm*Ront); curidx++;
                    ivr_idx++;

                    cooX[curidx] = grididx+nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 1/(Cm*Ronb); curidx++;
                    ivr_idx++;
                    cooX[curidx] = grididx+nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 1/(Cm*Ronb); curidx++;
                    ivr_idx++;
                }
            }
        }

      /* IVR block -> Onchip nodes */
      ivr_idx = 0;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            if(IVRLOC == ivr_loc[i][j]){
                for(l=0; l<nl-1; l++){
                    grididx = i*nc + j;
                    xoffset = nvp + ngp + nl*2*nr*nc + ntsv;
                    yoffset = nvp + ngp + l*2*nr*nc;

                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+2*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 1/(Civr*Ront); curidx++;
                    ivr_idx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+2*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 1/(Civr*Ront); curidx++;
                    ivr_idx++;

                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 1/(Civr*Ronb); curidx++;
                    ivr_idx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 1/(Civr*Ronb); curidx++;
                    ivr_idx++;
                }
            }
        }

      /* Diagonal block */
      xoffset = nvp + ngp + nl*2*nr*nc + ntsv;
      yoffset = nvp + ngp + nl*2*nr*nc + ntsv;
      ivr_idx = 0;
      for(k=0; k<nIVR; k++)
        for(l=0; l<nl-1; l++){
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = -1/(Civr*Ront); curidx++; ivr_idx++;
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = -1/(Civr*Ront); curidx++; ivr_idx++;
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = -1/(Civr*Ronb); curidx++; ivr_idx++;
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = -1/(Civr*Ronb); curidx++; ivr_idx++;
        }
  }

  if(model->config.PDN_pkgLC){
      /* Bottom row, VDD part */
      xoffset = m - 1;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = -iM*iA*iSv; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          if(vs){
              cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+ngp+(nl-1)*2*nr*nc+yoffset;
          }
          else{
              cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+ngp+yoffset;
          }
          cooV[curidx] = iM*iA*iPv; curidx++;
      }
      /* Bottom row, GND part */
      xoffset = m - 1;
      yoffset = nvp;
      for(i=0; i<ngp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = -iM*iB*iSg; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+ngp+nr*nc+yoffset;
          cooV[curidx] = -iM*iB*iPg; curidx++;
      }
      /* Right col, VDD part */
      xoffset = 0;
      yoffset = m - 1;
      for(i=0; i<nvp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iQv; curidx++;
      }
      /* Right col, GND part */
      xoffset = nvp;
      yoffset = m - 1;
      for(i=0; i<ngp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iQg; curidx++;
      }

      /* Individual points */
      xoffset = m - 1;
      yoffset = m - 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = -iN; curidx++;

      xoffset = m - 2;
      yoffset = m - 2;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = m - 1;
      yoffset = m - 2;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = -iM; curidx++;

      xoffset = m - 2;
      yoffset = m - 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 1/Cpp; curidx++;
  }

  if(curidx != nnz){
      printf("%d\t%d\n", curidx, nnz);
      fatal("COO Matrix A has different element count than nnz!\n");
  }

  coo2csc(n, nnz, cooX, cooY, cooV, asub, xa, a);

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

  // Create B
  if ( !(b = doubleMalloc(nnz)) ) fatal("Malloc fails for b[].\n");
  if ( !(bsub = intMalloc(nnz)) ) fatal("Malloc fails for bsub[].\n");
  if ( !(xb = intMalloc(n+1)) ) fatal("Malloc fails for xb[].\n");

  curidx = 0;
  /* Vdd pad dense block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++)
    for(j=0; j<nvp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(model->config.PDN_pkgLC)
          cooV[curidx] = -iPv;
        else
          cooV[curidx] = 0;
        curidx++;
    }

  /* Gnd pad dense block */
  xoffset = nvp;
  yoffset = nvp;
  for(i=0; i<ngp; i++)
    for(j=0; j<ngp; j++){
        cooX[curidx] = i + xoffset;
        cooY[curidx] = j + yoffset;
        if(model->config.PDN_pkgLC)
          cooV[curidx] = -iPg;
        else
          cooV[curidx] = 0;
        curidx++;
    }

  /* VDD pad two sparse block */
  xoffset = 0;
  yoffset = 0;
  for(i=0; i<nvp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
      if(-1 == grididx)
        fatal("Cannot find grid index for VDD pad\n");

      if(vs){
          cooX[curidx] = grididx+nvp+ngp+(nl-1)*2*nr*nc+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+ngp+(nl-1)*2*nr*nc+yoffset;
          cooV[curidx] = 0; curidx++;
      }
      else{
          cooX[curidx] = grididx+nvp+ngp+xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          cooX[curidx] = i+xoffset; cooY[curidx] = grididx+nvp+ngp+yoffset;
          cooV[curidx] = 0; curidx++;
      }
  }

  /* GND pad two sparse block */
  xoffset = nvp;
  yoffset = nvp;
  for(i=0; i<ngp; i++){
      grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
      if(-1 == grididx)
        fatal("Cannot find grid index for GND pad\n");

      cooX[curidx] = grididx+ngp+nr*nc+xoffset; cooY[curidx] = i+yoffset;
      cooV[curidx] = 0; curidx++;

      cooX[curidx] = i+xoffset; cooY[curidx] = grididx+ngp+nr*nc+yoffset;
      cooV[curidx] = 0; curidx++;
  }

  /* Onchip nodes */
  for(l=0; l<nl; l++){
      /* Vdd grid diagonal block */
      xoffset = nvp + ngp + l*2*nr*nc;
      yoffset = nvp + ngp + l*2*nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = 0; curidx++;
        }

      /* Vdd-Gnd grid diagonal */
      xoffset = nvp + ngp + l*2*nr*nc;
      yoffset = nvp + ngp + l*2*nr*nc + nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = 0;
            curidx++;
        }

      /* Gnd grid diagonal block */
      xoffset = nvp + ngp + l*2*nr*nc + nr*nc;
      yoffset = nvp + ngp + l*2*nr*nc + nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = 0; curidx++;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = 0; curidx++;
            }

            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = 0; curidx++;
        }

      /* Gnd-Vdd grid  diagonal */
      xoffset = nvp + ngp + l*2*nr*nc + nr*nc;
      yoffset = nvp + ngp + l*2*nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;
            cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
            cooV[curidx] = 0;
            curidx++;
        }
  }

  if(vs){
      /* Onchip nodes -> TSVs block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          xoffset = nvp + ngp + l*2*nr*nc;
          yoffset = nvp + ngp + nl*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 0; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSVs -> Onchip nodes block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          xoffset = nvp + ngp + nl*2*nr*nc;
          yoffset = nvp + ngp + l*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSV diagonal block */
      xoffset = nvp + ngp + nl*2*nr*nc;
      yoffset = nvp + ngp + nl*2*nr*nc;
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 0; curidx++;
                    tsv_count++;
                }
            }
      }
  }
  else{
      /* Onchip nodes -> TSVs block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          xoffset = nvp + ngp + l*2*nr*nc;
          yoffset = nvp + ngp + nl*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == VDDTSV){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = grididx+2*nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 0; curidx++;
                    tsv_count++;
                }
            }
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = grididx+nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = tsv_count+yoffset;
                    cooV[curidx] = 0; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSVs -> Onchip nodes block */
      tsv_count = 0;
      for(l=0; l<nl-1; l++){
          xoffset = nvp + ngp + nl*2*nr*nc;
          yoffset = nvp + ngp + l*2*nr*nc;
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == VDDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+2*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    tsv_count++;
                }
            }
          for(i=0; i<nr; i++)
            for(j=0; j<nc; j++){
                grididx = i*nc + j;
                if(model->layers[l].tsv.loc[i][j] == GNDTSV){
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    cooX[curidx] = tsv_count+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    tsv_count++;
                }
            }
      }

      /* TSV diagonal block */
      xoffset = nvp + ngp + nl*2*nr*nc;
      yoffset = nvp + ngp + nl*2*nr*nc;
      for(k=0; k<ntsv; k++){
          cooX[curidx] = k+xoffset; cooY[curidx] = k+yoffset;
          cooV[curidx] = 0; curidx++;
      }
  }

  /* IVR block */
  if(nIVR > 0){
      /* Onchip nodes -> IVR block */
      ivr_idx = 0;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            if(IVRLOC == ivr_loc[i][j]){
                for(l=0; l<nl-1; l++){
                    grididx = i*nc + j;
                    xoffset = nvp + ngp + l*2*nr*nc;
                    yoffset = nvp + ngp + nl*2*nr*nc + ntsv;

                    Cg = get_onchip_cap(model, grididx, l+1);
                    Cm = Cg + 2*Civr;

                    cooX[curidx] = grididx+2*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = Civr/Cm; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;
                    cooX[curidx] = grididx+2*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = Civr/Cm; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;

                    Cg = get_onchip_cap(model, grididx, l);
                    Cm = Cg + 2*Civr;
                    cooX[curidx] = grididx+nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = Civr/Cm; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;
                    cooX[curidx] = grididx+nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = Civr/Cm; curidx++;
                    cooX[curidx] = grididx+3*nr*nc+xoffset; cooY[curidx] = ivr_idx+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;
                }
            }
        }

      /* IVR block -> Onchip nodes */
      ivr_idx = 0;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            if(IVRLOC == ivr_loc[i][j]){
                for(l=0; l<nl-1; l++){
                    grididx = i*nc + j;
                    xoffset = nvp + ngp + nl*2*nr*nc + ntsv;
                    yoffset = nvp + ngp + l*2*nr*nc;

                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+2*nr*nc+yoffset;
                    cooV[curidx] = 1; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+2*nr*nc+yoffset;
                    cooV[curidx] = 1; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;

                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = 1; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = 1; curidx++;
                    cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = grididx+3*nr*nc+yoffset;
                    cooV[curidx] = 0; curidx++;
                    ivr_idx++;
                }
            }
        }

      /* Diagonal block */
      xoffset = nvp + ngp + nl*2*nr*nc + ntsv;
      yoffset = nvp + ngp + nl*2*nr*nc + ntsv;
      ivr_idx = 0;
      for(k=0; k<nIVR; k++)
        for(l=0; l<nl-1; l++){
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = 0; curidx++; ivr_idx++;
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = 0; curidx++; ivr_idx++;
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = 0; curidx++; ivr_idx++;
            cooX[curidx] = ivr_idx+xoffset; cooY[curidx] = ivr_idx+yoffset;
            cooV[curidx] = 0; curidx++; ivr_idx++;
        }
  }

  if(model->config.PDN_pkgLC){
      /* Bottom row, VDD part */
      xoffset = m - 1;
      yoffset = 0;
      for(i=0; i<nvp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_VDD);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          if(vs){
              cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+ngp+(nl-1)*2*nr*nc+yoffset;
          }
          else{
              cooX[curidx] = xoffset; cooY[curidx] = grididx+nvp+ngp+yoffset;
          }
          cooV[curidx] = 0; curidx++;
      }
      /* Bottom row, GND part */
      xoffset = m - 1;
      yoffset = nvp;
      for(i=0; i<ngp; i++){
          cooX[curidx] = xoffset; cooY[curidx] = i+yoffset;
          cooV[curidx] = 0; curidx++;

          grididx = PadIdx_to_GridIdx(model, i, LAYER_GND);
          if(-1 == grididx)
            fatal("Cannot find grid index for VDD pad\n");

          cooX[curidx] = xoffset; cooY[curidx] = grididx+ngp+nr*nc+yoffset;
          cooV[curidx] = 0; curidx++;
      }
      /* Right col, VDD part */
      xoffset = 0;
      yoffset = m - 1;
      for(i=0; i<nvp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iPv; curidx++;
      }
      /* Right col, GND part */
      xoffset = nvp;
      yoffset = m - 1;
      for(i=0; i<ngp; i++){
          cooX[curidx] = i+xoffset; cooY[curidx] = yoffset;
          cooV[curidx] = -iPg; curidx++;
      }

      /* Individual points */
      xoffset = m - 1;
      yoffset = m - 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = m - 2;
      yoffset = m - 2;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = m - 1;
      yoffset = m - 2;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;

      xoffset = m - 2;
      yoffset = m - 1;
      cooX[curidx] = xoffset; cooY[curidx] = yoffset;
      cooV[curidx] = 0; curidx++;
  }

  if(curidx != nnz){
      printf("%d\t%d\n", curidx, nnz);
      fatal("COO Matrix B has different element count than nnz!\n");
  }

  coo2csc(n, nnz, cooX, cooY, cooV, bsub, xb, b);

  dCreate_CompCol_Matrix(&B, m, n, nnz, b, bsub, xb, SLU_NC, SLU_D, SLU_GE);

  // E - B - delta_t*A/2
  SparseMatrix_mul_SingleNum(A, -delta_t/2);
  SparseMatrix_add(A, &B, -1);
  SparseMatrix_add_Iden(A, 1);

  if ( !(rhs = doubleMalloc(m)) ) fatal("Malloc fails for rhs[].\n");
  for(i=0; i<m; ++i) rhs[i] = 1;
  dCreate_Dense_Matrix(&Z, m, 1, rhs, m, SLU_DN, SLU_D, SLU_GE);

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = 0.01;
  options.SymmetricMode = NO;
  options.Equil = YES;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve the linear system. */
  dgssv(&options, A, *perm_c, *perm_r, L, U, &Z, &stat, &info);

  // Create RHS
  // E + (delta_t/2)*J
  SparseMatrix_mul_SingleNum(A, -1);
  SparseMatrix_add(A, &B, -2);
  SparseMatrix_add_Iden(A, 2);

  /* De-allocate storage */
  SUPERLU_FREE (rhs);
  Destroy_SuperMatrix_Store(&Z);
  Destroy_CompCol_Matrix(&B);
  StatFree(&stat);
  free_model_vector(p);

  free(cooV);
  free(cooX);
  free(cooY);
}

void PDN_init_trans_vector(model_t *model, double *g)
{
  int i, l;
  int cur, ntsv;
  int nIVR;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nvtsv = model->layers[0].tsv.num_vdd;
  int ngtsv = model->layers[0].tsv.num_gnd;
  int nml = model->layers[0].metal_layers.n_metal/2;
  int nbr = 2*nr*nc - nr - nc;
  int vs = model->config.v_stacking;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;

  cur = 0;
  if(!model->is_3D){
      for(i=0; i<nvp; i++){
          g[cur] = 0; cur++;
      }
      for(i=0; i<(nr*nc); i++){
          g[cur] = vdd; cur++;
      }
      if(model->config.PDN_gridL){
          for(i=0; i<nbr*nml; i++){
              g[cur] = 0; cur++;
          }
      }
      for(i=0; i<ngp; i++){
          g[cur] = 0; cur++;
      }
      for(i=0; i<(nr*nc); i++){
          g[cur] = gnd; cur++;
      }
      if(model->config.PDN_gridL){
          for(i=0; i<nbr*nml; i++){
              g[cur] = 0; cur++;
          }
      }
      if(model->config.PDN_pkgLC){
          g[cur] = vdd-gnd; cur++;
          g[cur] = 0; cur++;
      }
  }
  else{
      if(vs)
        ntsv = ngtsv*(nl-1);
      else
        ntsv = (nvtsv+ngtsv)*(nl-1);

      for(i=0; i<nvp; i++){
          g[cur] = 0; cur++;
      }
      for(i=0; i<ngp; i++){
          g[cur] = 0; cur++;
      }
      for(l=0; l<nl; l++){
          if(vs){
              for(i=0; i<(nr*nc); i++){
                  g[cur] = (l+1)*vdd; cur++;
              }
              for(i=0; i<(nr*nc); i++){
                  g[cur] = l*vdd; cur++;
              }
          }
          else{
              for(i=0; i<(nr*nc); i++){
                  g[cur] = vdd; cur++;
              }
              for(i=0; i<(nr*nc); i++){
                  g[cur] = gnd; cur++;
              }
          }
      }
      for(i=0; i<ntsv; i++){
          g[cur] = 0;
          cur++;
      }
      if(vs){
          nIVR = model->sc_converter->num_IVR;
          for(i=0; i<nIVR; i++)
            for(l=0; l<nl-1; l++){
                g[cur] = (l+1)*vdd; cur++;
                g[cur] = (l+1)*vdd; cur++;
                g[cur] = (l+1)*vdd; cur++;
                g[cur] = (l+1)*vdd; cur++;
            }
      }
      if(model->config.PDN_pkgLC){
          if(vs)
            g[cur] = nl*(vdd-gnd);
          else
            g[cur] = vdd-gnd;
          cur++;
          g[cur] = 0;
          cur++;
      }
  }
}

void Finalize_rhs(model_t *model, model_vector_t *power)
{
  int i, l;
  int tmp_i, tmp_j;
  int cur;
  int ntsv, nIVR;
  int **ivr_loc;

  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int nvp = model->c4->vdd_num;
  int ngp = model->c4->gnd_num;
  int nvtsv = model->layers[0].tsv.num_vdd;
  int ngtsv = model->layers[0].tsv.num_gnd;
  int nbr = 2*nr*nc - nr - nc;
  int nml = model->layers[0].metal_layers.n_metal/2;
  int vs = model->config.v_stacking;

  double Cg, Cm, Civr;
  double Lpv = model->c4->pad_l;
  double Lpg = model->c4->pad_l;
  double Lsp = model->config.PDN_pkg_sL;
  double Lpp = model->config.PDN_pkg_pL;
  double vdd = model->config.vdd;
  double gnd = model->config.gnd;
  double delta_t = 1/(model->config.proc_clock_freq*model->config.PDN_step_percycle);

  /* Add TSV RL into pad RL */
  //if(model->is_3D && vs){
  //    double Ltsv_sum = 0;
  //    for(l=0; l<nl-1; l++){
  //        Ltsv_sum += model->layers[l].tsv.l;
  //    }
  //    Lpv += Ltsv_sum;
  //}

  double iA = 1 / (nvp*Lsp/Lpv + 1);
  double iB = 1 / (ngp*Lsp/Lpg + 1);
  double iM = 1 / ((iA+iB)*Lsp + Lpp);

  // Perform last_trans += delta_t(fn + f(n+1))/2
  cur = 0;
  if(!model->is_3D){
      for(i=0; i<nvp; i++){
          model->last_trans[cur] += delta_t * vdd / Lpv;
          cur++;
      }
      for(i=0; i<(nr*nc); i++){
          Cg = get_onchip_cap(model, i, 0);
          model->last_trans[cur] -= (delta_t/2) * 
            (power->cuboid[0][0][i] + model->last_power[i]) / (Cg * (vdd-gnd));
          cur++;
      }
      if(model->config.PDN_gridL){
          for(i=0; i<nbr*nml; i++){
              model->last_trans[cur] += 0;
              cur++;
          }
      }
      for(i=0; i<ngp; i++){
          model->last_trans[cur] -= delta_t * gnd / Lpg;
          cur++;
      }
      for(i=0; i<(nr*nc); i++){
          Cg = get_onchip_cap(model, i, 0);
          model->last_trans[cur] += (delta_t/2) * 
            (power->cuboid[0][0][i] + model->last_power[i]) / (Cg * (vdd-gnd));
          cur++;
      }
      if(model->config.PDN_gridL){
          for(i=0; i<nbr*nml; i++){
              model->last_trans[cur] += 0;
              cur++;
          }
      }
      if(model->config.PDN_pkgLC){
          model->last_trans[cur] += 0;
          cur++;
          model->last_trans[cur] += delta_t * (vdd*iM*iA - gnd*iM*iB);
          cur++;
      }
  }
  else{
      if(vs){
          ntsv = ngtsv*(nl-1);
          nIVR = model->sc_converter->num_IVR;
          Civr = model->sc_converter->c;
          ivr_loc = model->sc_converter->loc;
      }
      else{
          ntsv = (nvtsv+ngtsv)*(nl-1);
          nIVR = 0;
      }

      if(vs)
        for(i=0; i<nvp; i++){
            model->last_trans[cur] += delta_t * nl*vdd / Lpv;
            cur++;
        }
      else
        for(i=0; i<nvp; i++){
            model->last_trans[cur] += delta_t * vdd / Lpv;
            cur++;
        }

      for(i=0; i<ngp; i++){
          model->last_trans[cur] -= delta_t * gnd / Lpg;
          cur++;
      }
      for(l=0; l<nl; l++){
          for(i=0; i<(nr*nc); i++){
              model->last_trans[cur] += 0;
              cur++;
          }
          for(i=0; i<(nr*nc); i++){
              model->last_trans[cur] += 0;
              cur++;
          }
      }
      for(i=0; i<ntsv; i++){
          model->last_trans[cur] += 0;
          cur++;
      }
      for(i=0; i<nIVR*(nl-1)*4; i++){
          model->last_trans[cur] += 0;
          cur++;
      }
      if(model->config.PDN_pkgLC){
          model->last_trans[cur] += 0;
          cur++;
          if(vs)
            model->last_trans[cur] += delta_t * (nl*vdd*iM*iA - gnd*iM*iB);
          else
            model->last_trans[cur] += delta_t * (vdd*iM*iA - gnd*iM*iB);
          cur++;
      }
  }
}
