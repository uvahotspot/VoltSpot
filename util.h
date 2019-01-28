#ifndef __UTIL_H
#define __UTIL_H

#include <stdio.h>
#include "slu_ddefs.h"

#define MAX(x,y)		(((x)>(y))?(x):(y))
#define MIN(x,y)		(((x)<(y))?(x):(y))
#define MAX3(a,b,c)		MAX(MAX(a,b),c)
#define MIN3(a,b,c)		MIN(MIN(a,b),c)
#define MID3(a,b,c)		((MIN(a,b)<(c))?(MIN(MAX(a,b),c)):(MAX(MIN(a,b),c)))
#define MAX4(a,b,c,d)	MAX(MAX(MAX(a,b),c),d)
#define MIN4(a,b,c,d)	MIN(MIN(MIN(a,b),c),d)
#define DELTA			1.0e-6
#define LARGENUM		1.0e100
#define LARGEINT		99999999
#define NULLFILE		"(null)"
#define VAPER		    1.256637e-6

#define PI	3.1416

#define TRUE	  		1
#define	FALSE	  		0

#define MAXIMUM			1
#define	MINIMUM			0

#define RAND_SEED		1500450271
#define SEED_UPPER	1000000000
#define SEED_LOWER	100000000

#define STR_SIZE		512
#define LINE_SIZE		65536
#define MAX_ENTRIES		512

int eq(double x, double y);
int le(double x, double y);
int ge(double x, double y);
void fatal (char *s);
void warning (char *s);
void swap_ival (int *a, int *b);
void swap_dval (double *a, double *b);
int tolerant_ceil(double val);
int tolerant_floor(double val);

/* vector routines	*/
double 	*dvector(int n);
void free_dvector(double *v);
void dump_dvector(double *v, int n);
void copy_dvector (double *dst, double *src, int n);
void zero_dvector (double *v, int n);
void abs_dvector (double *v, int n);
void negL_dvector (double *v, int n);

double **dvector_2D(int n, int d);
void free_dvector_2D(double **v, int d);

int *ivector(int n);
void free_ivector(int *v);
void dump_ivector(int *v, int n);
void copy_ivector (int *dst, int *src, int n);
void zero_ivector (int *v, int n);
void abs_ivector (int *v, int n);

int **ivector_2D(int n, int d);
void free_ivector_2D(int **v, int d);

/* matrix routines - Thanks to Greg Link
 * from Penn State University for the 
 * memory allocators/deallocators
 */
double **dmatrix(int nr, int nc);
void free_dmatrix(double **m);
void dump_dmatrix(double **m, int nr, int nc);
void dump_dmatrix_file(double **m, int nr, int nc, char *file);
void dump_dmatrix_skip0_file(double **m, int nr, int nc, char *file);
void copy_dmatrix(double **dst, double **src, int nr, int nc);
void zero_dmatrix(double **m, int nr, int nc);
void resize_dmatrix(double **m, int nr, int nc);
/* mirror the lower triangle to make 'm' fully symmetric	*/
void mirror_dmatrix(double **m, int n);

int **imatrix(int nr, int nc);
void free_imatrix(int **m);
void dump_imatrix(int **m, int nr, int nc);
void dump_imatrix_file(int **m, int nr, int nc, char *file);
void dump_imatrix_skip0_file(int **m, int nr, int nc, char *file);
void copy_imatrix(int **dst, int **src, int nr, int nc);
void resize_imatrix(int **m, int nr, int nc);
int  sum_imatrix(int **m, int nr, int nc);

/* Matrix analysis*/
double sum_dmatrix(double **m, int nr, int nc);
double ave_dmatrix(double **m, int nr, int nc);
double max_dmatrix(double **m, int nr, int nc);
double max_dmatrix_pos(double **m, int nr, int nc, int* i, int *j);
double min_dmatrix(double **m, int nr, int nc);
double min_dmatrix_pos(double **m, int nr, int nc, int* i, int *j);
double var_dmatrix(double **m, int nr, int nc, double ave);
int    above_threshold_dmatrix(double **m, int nr, int nc, double th);

/* skip0 only deals with non-zero values */
double ave_dmatrix_skip0(double **m, int nr, int nc);
double max_dmatrix_skip0(double **m, int nr, int nc);
double min_dmatrix_skip0(double **m, int nr, int nc);
double var_dmatrix_skip0(double **m, int nr, int nc, double ave);

/* vector analysis*/
double sum_dvector(double *m, int n);
double ave_dvector(double *m, int n);
double max_dvector(double *m, int n);
double min_dvector(double *m, int n);
double var_dvector(double *m, int n, double ave);
double max_dvector_skip0(double *m, int n);
int    above_threshold_dvector(double *m, int n, double th);

/* vector operation */
void add_dvector(double *dst, double *rs1, double *rs2, int len); // dst = rs1+rs2
void sub_dvector(double *dst, double *rs1, double *rs2, int len); // dst = rs1-rs2
void mul_val_dvector(double *m, double val, int len); // m = m*val
void sin_dvector(double *dst, double *amp, double time, double freq, int len); // dst = amp*sin(2*pi*f*t)

/* routines for 3-d matrix with tail	*/
/* allocate 3-d matrix with 'nr' rows, 'nc' cols, 
 * 'nl' layers	and a tail of 'xtra' elements 
 */
double ***dcuboid_tail(int nr, int nc, int nl, int xtra);
/* destructor	*/
void free_dcuboid(double ***m);

/* initialize random number generator	*/
void init_rand(void);
/* random number within the range [0, max-1]	*/
int rand_upto(int max);
/* random number in the range [0, 1)	*/
double rand_fraction(void);

/* a table of name value pairs	*/
typedef struct str_pair_st
{
	char name[STR_SIZE];
	char value[STR_SIZE];
}str_pair;
/* 
 * reads tab-separated name-value pairs from file into
 * a table of size max_entries and returns the number 
 * of entries read successfully
 */
int read_str_pairs(str_pair *table, int max_entries, char *file);
/* same as above but from command line instead of a file	*/
int parse_cmdline(str_pair *table, int max_entries, int argc, char **argv);
/* append the table onto a file	*/
void dump_str_pairs(str_pair *table, int size, char *file, char *prefix);
/* table lookup	for a name */
int get_str_index(str_pair *table, int size, char *str);
/* 
 * remove duplicate names in the table - the entries later 
 * in the table are discarded. returns the new size of the
 * table
 */
int str_pairs_remove_duplicates(str_pair *table, int size);
/* 
 * binary search a sorted double array 'arr' of size 'n'. if found,
 * the 'loc' pointer has the address of 'ele' and the return 
 * value is TRUE. otherwise, the return value is FALSE and 'loc' 
 * points to the 'should have been' location
 */
int bsearch_double(double *arr, int n, double ele, double **loc);
/* 
 * binary search and insert an element into a partially sorted
 * double array if not already present. returns FALSE if present
 */
int bsearch_insert_double(double *arr, int n, double ele);

/* 
 * Full search an arbitary array, if found, return index
 * return -1 if not found
 */
int full_search_int(int *arr, int size, int ele);

/* 
 * population count of an 8-bit integer - using pointers from 
 * http://aggregate.org/MAGIC/
 */
unsigned int ones8(register unsigned char n);
/* 
 * find the number of non-empty, non-comment lines
 * in a file open for reading
 */
int count_significant_lines(FILE *fp);
/* if even, return 1, odd 0 */
int even_or_odd(int x);
void print_time();

void int_to_str(int input, int width, char *strOut);

/* Sparse Matrix operations */
void SparseMatrix_mul_Vector(SuperMatrix *A, double *rhs);
void Vector_add_Vector(int m, double *v1, double *v2);
void SparseMatrix_add_Iden(SuperMatrix *A, int mul);
void SparseMatrix_dump_diag(SuperMatrix *A);
void SparseMatrix_mul_SingleNum(SuperMatrix *A, double num);
void SparseMatrix_add(SuperMatrix *A, SuperMatrix *B, int mul);
void Print_CompCol_Matrix_detail(char *what, SuperMatrix *A);

int coo2csc(int size, int nnz, 
            int *cooX, int *cooY, double *cooV,
            int *cscRowInd, int *cscColPtr, double *cscV);
int c2c_cmp( const void *a , const void *b);

#endif
