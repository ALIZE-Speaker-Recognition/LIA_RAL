#ifndef SVDUTIL_H
#define SVDUTIL_H

#include "svdlib.h"

#define SAFE_FREE(a) {if (a) {free(a); a = NULL;}}

/* Allocates an array of longs. */
extern long *svd_longArray(long size, char empty, char *name);
/* Allocates an array of doubles. */
extern double *svd_doubleArray(long size, char empty, char *name);

extern void svd_debug(char *fmt, ...);
extern void svd_error(char *fmt, ...);
extern void svd_fatalError(char *fmt, ...);
extern FILE *svd_fatalReadFile(char *filename);
extern FILE *svd_readFile(char *fileName);
extern FILE *svd_writeFile(char *fileName, char append);
extern void svd_closeFile(FILE *file);

extern char svd_readBinInt(FILE *file, int *val);
extern char svd_readBinFloat(FILE *file, float *val);
extern char svd_writeBinInt(FILE *file, int x);
extern char svd_writeBinFloat(FILE *file, float r);

/************************************************************** 
 * returns |a| if b is positive; else fsign returns -|a|      *
 **************************************************************/ 
extern double svd_fsign(double a, double b);

/************************************************************** 
 * returns the larger of two double precision numbers         *
 **************************************************************/ 
extern double svd_dmax(double a, double b);

/************************************************************** 
 * returns the smaller of two double precision numbers        *
 **************************************************************/ 
extern double svd_dmin(double a, double b);

/************************************************************** 
 * returns the larger of two integers                         *
 **************************************************************/ 
extern long svd_imax(long a, long b);

/************************************************************** 
 * returns the smaller of two integers                        *
 **************************************************************/ 
extern long svd_imin(long a, long b);

/************************************************************** 
 * Function scales a vector by a constant.     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
extern void svd_dscal(long n, double da, double *dx, long incx);

/************************************************************** 
 * function scales a vector by a constant.	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
extern void svd_datx(long n, double da, double *dx, long incx, double *dy, long incy);

/************************************************************** 
 * Function copies a vector x to a vector y	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
extern void svd_dcopy(long n, double *dx, long incx, double *dy, long incy);

/************************************************************** 
 * Function forms the dot product of two vectors.      	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
extern double svd_ddot(long n, double *dx, long incx, double *dy, long incy);

/************************************************************** 
 * Constant times a vector plus a vector     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
extern void svd_daxpy (long n, double da, double *dx, long incx, double *dy, long incy);

/********************************************************************* 
 * Function sorts array1 and array2 into increasing order for array1 *
 *********************************************************************/
extern void svd_dsort2(long igap, long n, double *array1, double *array2);

/************************************************************** 
 * Function interchanges two vectors		     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
extern void svd_dswap(long n, double *dx, long incx, double *dy, long incy);

/***************************************************************** 
 * Function finds the index of element having max. absolute value*
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/ 
extern long svd_idamax(long n, double *dx, long incx);

/**************************************************************
 * multiplication of matrix B by vector x, where B = A'A,     *
 * and A is nrow by ncol (nrow >> ncol). Hence, B is of order *
 * n = ncol (y stores product vector).		              *
 **************************************************************/
extern void svd_opb(SMat A, double *x, double *y, double *temp);

/***********************************************************
 * multiplication of matrix A by vector x, where A is 	   *
 * nrow by ncol (nrow >> ncol).  y stores product vector.  *
 ***********************************************************/
extern void svd_opa(SMat A, double *x, double *y);

/***********************************************************************
 *                                                                     *
 *				random2()                              *
 *                        (double precision)                           *
 ***********************************************************************/
extern double svd_random2(long *iy);

/************************************************************** 
 *							      *
 * Function finds sqrt(a^2 + b^2) without overflow or         *
 * destructive underflow.				      *
 *							      *
 **************************************************************/ 
extern double svd_pythag(double a, double b);


#endif /* SVDUTIL_H */
