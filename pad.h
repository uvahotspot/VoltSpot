#ifndef __PAD_H_
#define __PAD_H_

#include "PDN_sim.h"
#include "PDN_analyze.h"

// pad location definition
#define UNDEF            0
#define PGPAD            1
#define IOPAD            (1<<1)
#define IORAN            (1<<2)

/* pad manipulation	*/
void mark_pads(model_t *model, int pad_type);
void mark_all(model_t *model, int pad_type);
void remove_pg_for_io(model_t *model);
int check_pad_neighbour(model_t *model, int r, int c, int grid_type);
int is_grid_pad(model_t *model, int r, int c, int layer);
int node_allow_io(model_t *model, int r, int c);
void parse_pad_loc(model_t *model, char *file);
int count_pads(int **loc, int nr, int nc, int pad_type);

void dump_anypadloc(model_t *model, char *file, int type);
void dump_cur_with_cordt(model_t *model, status_t *status, char *file);

#endif
