#include <stdio.h>
#include <string.h>
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif
#include <stdlib.h>
#include <math.h>

#include "flp.h"
#include "util.h"

/* translate the floorplan to new origin (x,y)	*/
void PDN_flp_translate(PDN_flp_t *flp, double x, double y)
{
	int i;
	double minx = flp->units[0].leftx;
	double miny = flp->units[0].bottomy;

	for (i=1; i < flp->n_units; i++) {
		if (minx > flp->units[i].leftx)
			minx = flp->units[i].leftx;
		if (miny > flp->units[i].bottomy)
			miny = flp->units[i].bottomy;
	}
	for (i=0; i < flp->n_units; i++) {
		flp->units[i].leftx += (x - minx);
		flp->units[i].bottomy += (y - miny);
	}
}

/* scale the floorplan by a factor 'factor'	*/
//void flp_scale(PDN_flp_t *flp, double factor)
//{
//	int i;
//	double minx = flp->units[0].leftx;
//	double miny = flp->units[0].bottomy;
//
//	for (i=1; i < flp->n_units; i++) {
//		if (minx > flp->units[i].leftx)
//			minx = flp->units[i].leftx;
//		if (miny > flp->units[i].bottomy)
//			miny = flp->units[i].bottomy;
//	}
//	for(i=0; i < flp->n_units; i++) {
//		flp->units[i].leftx = (flp->units[i].leftx - minx) * factor + minx;
//		flp->units[i].bottomy = (flp->units[i].bottomy - miny) * factor + miny;
//		flp->units[i].width *= factor;
//		flp->units[i].height *= factor;
//	}
//}

/* 
 * change the orientation of the floorplan by
 * rotating and/or flipping. the target orientation
 * is specified in 'target'. 'width', 'height', 'xorig'
 * and 'yorig' are those of 'flp' respectively.
 */
//void flp_change_orient(PDN_flp_t *flp, double xorig, double yorig,
//					   double width, double height, orient_t target)
//{
//	int i;
//
//	for(i=0; i < flp->n_units; i++) {
//		double leftx, bottomy, rightx, topy;
//		/* all co-ordinate calculations are 
//		 * done assuming (0,0) as the center. 
//		 * so, shift accordingly
//		 */
//		leftx = flp->units[i].leftx  - (xorig + width / 2.0);
//		bottomy = flp->units[i].bottomy - (yorig + height / 2.0);
//		rightx = leftx + flp->units[i].width;
//		topy = bottomy + flp->units[i].height;
//		/* when changing orientation, leftx and 
//		 * bottomy of a rectangle could change
//		 * to one of the other three corners. 
//		 * also, signs of the co-ordinates
//		 * change according to the rotation
//		 * or reflection. Further x & y are
//		 * swapped for rotations that are
//		 * odd multiples of 90 degrees
//		 */
//		switch(target) {
//			case ROT_0:
//					flp->units[i].leftx = leftx;
//					flp->units[i].bottomy = bottomy;
//					break;
//			case ROT_90:
//					flp->units[i].leftx = -topy;
//					flp->units[i].bottomy = leftx;
//					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
//					break;
//			case ROT_180:
//					flp->units[i].leftx = -rightx;
//					flp->units[i].bottomy = -topy;
//					break;
//			case ROT_270:
//					flp->units[i].leftx = bottomy;
//					flp->units[i].bottomy = -rightx;
//					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
//					break;
//			case FLIP_0:
//					flp->units[i].leftx = -rightx;
//					flp->units[i].bottomy = bottomy;
//					break;
//			case FLIP_90:
//					flp->units[i].leftx = bottomy;
//					flp->units[i].bottomy = leftx;
//					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
//					break;
//			case FLIP_180:
//					flp->units[i].leftx = leftx;
//					flp->units[i].bottomy = -topy;
//					break;
//			case FLIP_270:
//					flp->units[i].leftx = -topy;
//					flp->units[i].bottomy = -rightx;
//					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
//					break;
//			default:
//					fatal("unknown orientation\n");
//					break;
//		}
//		/* translate back to original origin	*/
//		flp->units[i].leftx += (xorig + width / 2.0);
//		flp->units[i].bottomy += (yorig + height / 2.0);
//	}
//}

/* 
 * find the number of units from the 
 * floorplan file
 */
int PDN_flp_count_units(FILE *fp)
{
    char str1[LINE_SIZE], str2[LINE_SIZE];
	char name[STR_SIZE];
	double leftx, bottomy, width, height;
	char *ptr;
    int count = 0;

	fseek(fp, 0, SEEK_SET);
	while(!feof(fp)) {
		fgets(str1, LINE_SIZE, fp);

		strcpy(str2, str1);
		
		/* ignore comments and empty lines	*/
		ptr = strtok(str1, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;

		/* functional block placement information	*/
		if (sscanf(str2, "%s%lf%lf%lf%lf", name, &leftx, &bottomy,
		  		   &width, &height) == 5)
			count++;

		if (feof(fp))
			break;

	}
	return count;
}

PDN_flp_t *PDN_flp_alloc_init_mem(int count)
{
	int i;
	PDN_flp_t *flp;
	flp = (PDN_flp_t *) calloc (1, sizeof(PDN_flp_t));
	if(!flp)
		fatal("memory allocation error\n");
	flp->units = (PDN_unit_t *) calloc(count, sizeof(PDN_unit_t));
	//flp->wire_density = (double **) calloc(count, sizeof(double *));
	//if (!flp->units || !flp->wire_density)
	if (!flp->units)
		fatal("memory allocation error\n");
	flp->n_units = count;

	//for (i=0; i < count; i++) {
	//  flp->wire_density[i] = (double *) calloc(count, sizeof(double));
	//  if (!flp->wire_density[i])
	//  	fatal("memory allocation error\n");
	//}
	return flp;
}

/* populate block information	*/
void PDN_flp_populate_blks(PDN_flp_t *flp, FILE *fp)
{
	int i=0;
	char str[LINE_SIZE], copy[LINE_SIZE]; 
	char name1[STR_SIZE], name2[STR_SIZE];
	double width, height, leftx, bottomy;
	//double wire_density;
	char *ptr;

	fseek(fp, 0, SEEK_SET);
	while(!feof(fp)) {		/* second pass	*/
		fgets(str, LINE_SIZE, fp);
		if (feof(fp))
			break;
		strcpy(copy, str);

		/* ignore comments and empty lines	*/
		ptr = strtok(str, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;

		if (sscanf(copy, "%s%lf%lf%lf%lf", name1, &width, &height, 
				   &leftx, &bottomy) == 5) {
			strcpy(flp->units[i].name, name1);
			flp->units[i].width = width;
			flp->units[i].height = height;
			flp->units[i].leftx = leftx;
			flp->units[i].bottomy = bottomy;
			i++;
			/* skip connectivity info	*/
    } else// if (sscanf(copy, "%s%s%lf", name1, name2, &wire_density) != 3) 
      fatal("invalid floorplan file format\n");
	}
	if (i != flp->n_units)
	  fatal("mismatch of number of units\n");
}

/* populate connectivity info	*/
//void flp_populate_connects(PDN_flp_t *flp, FILE *fp)
//{
//	char str1[LINE_SIZE], str2[LINE_SIZE]; 
//	char name1[STR_SIZE], name2[STR_SIZE];
//	/* dummy fields	*/
//	double f1, f2, f3, f4, f5, f6;
//	double wire_density;
//	char *ptr;
//	int x, y, temp;
//
//	/* initialize wire_density	*/
//	for(x=0; x < flp->n_units; x++)
//		for(y=0; y < flp->n_units; y++)
//			flp->wire_density[x][y] = 0.0;
//
//	fseek(fp, 0, SEEK_SET);
//	while(!feof(fp)) {
//		fgets(str1, LINE_SIZE, fp);
//		if (feof(fp))
//			break;
//		strcpy(str2, str1);
//
//		/* ignore comments and empty lines	*/
//		ptr = strtok(str1, " \r\t\n");
//		if (!ptr || ptr[0] == '#')
//			continue;
//
//		/* lines with unit positions	*/
//		if (sscanf(str2, "%s%lf%lf%lf%lf%lf%lf", name1, &f1, &f2, &f3, &f4, &f5, &f6) == 7 ||
//		  	/* flp_desc like lines. ignore them	*/
//		  	sscanf(str2, "%s%lf%lf%lf%d", name1, &f1, &f2, &f3, &temp) == 5)
//			continue;
//
//		 /* lines with connectivity info	*/
//		else if (sscanf(str2, "%s%s%lf", name1, name2, &wire_density) == 3) {
//			x = PDN_get_blk_index(flp, name1);
//			y = PDN_get_blk_index(flp, name2);
//
//			if (x == y)
//				fatal("block connected to itself?\n");
//
//			if (!flp->wire_density[x][y] && !flp->wire_density[y][x])
//				flp->wire_density[x][y] = flp->wire_density[y][x] = wire_density;
//			else if((flp->wire_density[x][y] != flp->wire_density[y][x]) ||
//			        (flp->wire_density[x][y] != wire_density)) {
//				sprintf(str2, "wrong connectivity information for blocks %s and %s\n", 
//				        name1, name2);
//				fatal(str2);
//			}
//		} else 
//		  	fatal("invalid floorplan file format\n");
//	} /* end while	*/
//}

PDN_flp_t *PDN_read_flp(char *file)
{
	char str[STR_SIZE];
	FILE *fp;
	PDN_flp_t *flp;
	int count, i, j;

	if (!strcasecmp(file, "stdin"))
		fp = stdin;
	else
		fp = fopen (file, "r");

	if (!fp) {
		sprintf(str, "error opening file %s\n", file);
		fatal(str);
	}

	/* 1st pass - find n_units	*/
	count = PDN_flp_count_units(fp);

	if(!count)
		fatal("no units specified in the floorplan file\n");

	/* allocate initial memory */
	flp = PDN_flp_alloc_init_mem(count);

	/* 2nd pass - populate block info	*/
  PDN_flp_populate_blks(flp, fp);

  //for (i=0; i < flp->n_units; i++)
  //  for (j=0; j < flp->n_units; j++)
  //    flp->wire_density[i][j] = 1.0;

  if(fp != stdin)
    fclose(fp);	

  /* make sure the origin is (0,0)	*/
  PDN_flp_translate(flp, 0, 0);	
  return flp;
}

//void dump_flp(PDN_flp_t *flp, char *file, int dump_connects)
//{
//	char str[STR_SIZE];
//	int i, j;
//	FILE *fp;
//
//	if (!strcasecmp(file, "stdout"))
//		fp = stdout;
//	else if (!strcasecmp(file, "stderr"))
//		fp = stderr;
//	else 	
//		fp = fopen (file, "w");
//
//	if (!fp) {
//		sprintf(str, "error opening file %s\n", file);
//		fatal(str);
//	}
//	/* functional unit placement info	*/
//	for(i=0; i < flp->n_units; i++)  {
//		fprintf(fp, "%s\t%.11f\t%.11f\t%.11f\t%.11f\n",
//				flp->units[i].name, flp->units[i].width, flp->units[i].height,
//				flp->units[i].leftx, flp->units[i].bottomy);
//	}
//
//	if (dump_connects) {
//		fprintf(fp, "\n");
//		/* connectivity information	*/
//		for(i=1; i < flp->n_units; i++)
//			for(j=0; j < i; j++)
//				if (flp->wire_density[i][j])
//					fprintf(fp, "%s\t%s\t%.3f\n", flp->units[i].name,
//							flp->units[j].name, flp->wire_density[i][j]);
//	}
//	
//	if(fp != stdout && fp != stderr)
//		fclose(fp);
//}

void PDN_free_flp(PDN_flp_t *flp)
{
  int i;
  //for (i=0; i < flp->n_units + compacted; i++) {
  //    free(flp->wire_density[i]);
  //}
  free(flp->units);
  //free(flp->wire_density);
  free(flp);
}

//void print_flp_fig (PDN_flp_t *flp)
//{
//	int i;
//	double leftx, bottomy, rightx, topy;
//
//	fprintf(stdout, "FIG starts\n");
//	for (i=0; i< flp->n_units; i++) {
//		leftx = flp->units[i].leftx;
//		bottomy = flp->units[i].bottomy;
//		rightx = flp->units[i].leftx + flp->units[i].width;
//		topy = flp->units[i].bottomy + flp->units[i].height;
//		fprintf(stdout, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n", 
//			    leftx, bottomy, leftx, topy, rightx, topy, rightx, bottomy, 
//				leftx, bottomy);
//		fprintf(stdout, "%s\n", flp->units[i].name);
//	}
//	fprintf(stdout, "FIG ends\n");
//}

/* debug print	*/
//void print_flp (PDN_flp_t *flp)
//{
//	int i, j;
//
//	fprintf(stdout, "printing floorplan information for %d blocks\n", flp->n_units);
//	fprintf(stdout, "name\tarea\twidth\theight\tleftx\tbottomy\trightx\ttopy\n");
//	for (i=0; i< flp->n_units; i++) {
//		double area, width, height, leftx, bottomy, rightx, topy;
//		char *name;
//		name = flp->units[i].name;
//		width = flp->units[i].width;
//		height = flp->units[i].height;
//		area = width * height;
//		leftx = flp->units[i].leftx;
//		bottomy = flp->units[i].bottomy;
//		rightx = flp->units[i].leftx + flp->units[i].width;
//		topy = flp->units[i].bottomy + flp->units[i].height;
//		fprintf(stdout, "%s\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", 
//			    name, area, width, height, leftx, bottomy, rightx, topy);
//	}
//	fprintf(stdout, "printing connections:\n");
//	for (i=0; i< flp->n_units; i++)
//		for (j=i+1; j < flp->n_units; j++)
//			if (flp->wire_density[i][j])
//				fprintf(stdout, "%s\t%s\t%lg\n", flp->units[i].name, 
//						flp->units[j].name, flp->wire_density[i][j]);
//}

int PDN_get_blk_index(PDN_flp_t *flp, char *name)
{
  int i;
  char msg[STR_SIZE];

  if (!flp)
    fatal("null pointer in PDN_get_blk_index\n");

  for (i = 0; i < flp->n_units; i++) {
      if (!strcasecmp(name, flp->units[i].name)) {
          return i;
      }
  }

  sprintf(msg, "block %s not found\n", name);
  fatal(msg);
  return -1;
}

//int is_horiz_adj(PDN_flp_t *flp, int i, int j)
//{
//	double x1, x2, x3, x4;
//	double y1, y2, y3, y4;
//
//	if (i == j) 
//		return FALSE;
//
//	x1 = flp->units[i].leftx;
//	x2 = x1 + flp->units[i].width;
//	x3 = flp->units[j].leftx;
//	x4 = x3 + flp->units[j].width;
//
//	y1 = flp->units[i].bottomy;
//	y2 = y1 + flp->units[i].height;
//	y3 = flp->units[j].bottomy;
//	y4 = y3 + flp->units[j].height;
//
//	/* diagonally adjacent => not adjacent */
//	if (eq(x2,x3) && eq(y2,y3))
//		return FALSE;
//	if (eq(x1,x4) && eq(y1,y4))
//		return FALSE;
//	if (eq(x2,x3) && eq(y1,y4))
//		return FALSE;
//	if (eq(x1,x4) && eq(y2,y3))
//		return FALSE;
//
//	if (eq(x1,x4) || eq(x2,x3))
//		if ((y3 >= y1 && y3 <= y2) || (y4 >= y1 && y4 <= y2) ||
//		    (y1 >= y3 && y1 <= y4) || (y2 >= y3 && y2 <= y4))
//			return TRUE;
//
//	return FALSE;
//}

//int is_vert_adj (PDN_flp_t *flp, int i, int j)
//{
//	double x1, x2, x3, x4;
//	double y1, y2, y3, y4;
//
//	if (i == j)
//		return FALSE;
//
//	x1 = flp->units[i].leftx;
//	x2 = x1 + flp->units[i].width;
//	x3 = flp->units[j].leftx;
//	x4 = x3 + flp->units[j].width;
//
//	y1 = flp->units[i].bottomy;
//	y2 = y1 + flp->units[i].height;
//	y3 = flp->units[j].bottomy;
//	y4 = y3 + flp->units[j].height;
//
//	/* diagonally adjacent => not adjacent */
//	if (eq(x2,x3) && eq(y2,y3))
//		return FALSE;
//	if (eq(x1,x4) && eq(y1,y4))
//		return FALSE;
//	if (eq(x2,x3) && eq(y1,y4))
//		return FALSE;
//	if (eq(x1,x4) && eq(y2,y3))
//		return FALSE;
//
//	if (eq(y1,y4) || eq(y2,y3))
//		if ((x3 >= x1 && x3 <= x2) || (x4 >= x1 && x4 <= x2) ||
//		    (x1 >= x3 && x1 <= x4) || (x2 >= x3 && x2 <= x4))
//			return TRUE;
//
//	return FALSE;
//}

//double get_shared_len(PDN_flp_t *flp, int i, int j)
//{
//	double p11, p12, p21, p22;
//	p11 = p12 = p21 = p22 = 0.0;
//
//	if (i==j) 
//		return FALSE;
//
//	if (is_horiz_adj(flp, i, j)) {
//		p11 = flp->units[i].bottomy;
//		p12 = p11 + flp->units[i].height;
//		p21 = flp->units[j].bottomy;
//		p22 = p21 + flp->units[j].height;
//	}
//
//	if (is_vert_adj(flp, i, j)) {
//		p11 = flp->units[i].leftx;
//		p12 = p11 + flp->units[i].width;
//		p21 = flp->units[j].leftx;
//		p22 = p21 + flp->units[j].width;
//	}
//
//	return (MIN(p12, p22) - MAX(p11, p21));
//}

double PDN_get_total_width(PDN_flp_t *flp)
{	
  int i;
  double min_x = flp->units[0].leftx;
  double max_x = flp->units[0].leftx + flp->units[0].width;

  for (i=1; i < flp->n_units; i++) {
      if (flp->units[i].leftx < min_x)
        min_x = flp->units[i].leftx;
      if (flp->units[i].leftx + flp->units[i].width > max_x)
        max_x = flp->units[i].leftx + flp->units[i].width;
  }

  return (max_x - min_x);
}

double PDN_get_total_height(PDN_flp_t *flp)
{	
  int i;
  double min_y = flp->units[0].bottomy;
  double max_y = flp->units[0].bottomy + flp->units[0].height;

  for (i=1; i < flp->n_units; i++) {
      if (flp->units[i].bottomy < min_y)
        min_y = flp->units[i].bottomy;
      if (flp->units[i].bottomy + flp->units[i].height > max_y)
        max_y = flp->units[i].bottomy + flp->units[i].height;
  }

  return (max_y - min_y);
}

//double get_minx(PDN_flp_t *flp)
//{
//	int i;
//	double min_x = flp->units[0].leftx;
//	
//	for (i=1; i < flp->n_units; i++)
//		if (flp->units[i].leftx < min_x)
//			min_x = flp->units[i].leftx;
//
//	return min_x;
//}

//double get_miny(PDN_flp_t *flp)
//{
//	int i;
//	double min_y = flp->units[0].bottomy;
//	
//	for (i=1; i < flp->n_units; i++)
//		if (flp->units[i].bottomy < min_y)
//			min_y = flp->units[i].bottomy;
//
//	return min_y;
//}

/* precondition: L2 should have been wrapped around	*/
//double get_core_width(PDN_flp_t *flp, char *l2_label)
//{
//	int i;
//	double min_x = LARGENUM;
//	double max_x = -LARGENUM;
//	
//	for (i=0; i < flp->n_units; i++) {
//		/* core is that part of the chip excluding the l2 and rim	*/
//		if (strstr(flp->units[i].name, l2_label) != flp->units[i].name &&
//			strstr(flp->units[i].name, RIM_PREFIX) != flp->units[i].name) {
//			if (flp->units[i].leftx < min_x)
//				min_x = flp->units[i].leftx;
//			if (flp->units[i].leftx + flp->units[i].width > max_x)
//				max_x = flp->units[i].leftx + flp->units[i].width;
//		}		
//	}
//
//	return (max_x - min_x);
//}

/* precondition: L2 should have been wrapped around	*/
//double get_core_height(PDN_flp_t *flp, char *l2_label)
//{	
//	int i;
//	double min_y = LARGENUM;
//	double max_y = -LARGENUM;
//	
//	for (i=0; i < flp->n_units; i++) {
//		/* core is that part of the chip excluding the l2 and rim	*/
//		if (strstr(flp->units[i].name, l2_label) != flp->units[i].name &&
//			strstr(flp->units[i].name, RIM_PREFIX) != flp->units[i].name) {
//			if (flp->units[i].bottomy < min_y)
//				min_y = flp->units[i].bottomy;
//			if (flp->units[i].bottomy + flp->units[i].height > max_y)
//				max_y = flp->units[i].bottomy + flp->units[i].height;
//		}		
//	}
//
//	return (max_y - min_y);
//}

//double get_total_area(PDN_flp_t *flp)
//{
//	int i;
//	double area = 0.0;
//	for(i=0; i < flp->n_units; i++)
//		area += flp->units[i].width * flp->units[i].height;
//	return area;	
//}

//double get_core_area(PDN_flp_t *flp, char *l2_label)
//{
//	int i;
//	double area = 0.0;
//	for(i=0; i < flp->n_units; i++)
//		if (strstr(flp->units[i].name, l2_label) != flp->units[i].name &&
//			strstr(flp->units[i].name, RIM_PREFIX) != flp->units[i].name)
//			area += flp->units[i].width * flp->units[i].height;
//	return area;		
//}

/* excluding the dead blocks	*/
//double get_core_occupied_area(PDN_flp_t *flp, char *l2_label)
//{
//	int i, num;
//	double dead_area = 0.0;
//	for(i=0; i < flp->n_units; i++) {
//		/* 
//		 * there can be a max of n-1 dead blocks where n is the
//		 * number of non-dead blocks (since each cut, vertical
//		 * or horizontal, can correspond to a maximum of one
//		 * dead block
//		 */
//		if ((sscanf(flp->units[i].name, DEAD_PREFIX"%d", &num) == 1) &&
//			(num < (flp->n_units-1) / 2))
//			dead_area += flp->units[i].width * flp->units[i].height;
//	}		
//	return get_core_area(flp, l2_label) - dead_area;	
//}

//double get_manhattan_dist(PDN_flp_t *flp, int i, int j)
//{
//	double x1 = flp->units[i].leftx + flp->units[i].width / 2.0;
//	double y1 = flp->units[i].bottomy + flp->units[i].height / 2.0;
//	double x2 = flp->units[j].leftx + flp->units[j].width / 2.0;
//	double y2 = flp->units[j].bottomy + flp->units[j].height / 2.0;
//	return (fabs(x2-x1) + fabs(y2-y1));
//}
