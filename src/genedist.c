/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998, 2001  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include "R_ext/Arith.h"
#include "R_ext/Error.h"
#include "R_ext/Applic.h" /* machar */

static int intcompare(const void *p1, const void *p2)
{
  int i = *((int *)p1);
  int j = *((int *)p2);
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}

static double mm_euclidean(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;

    count= 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = (x[i1] - x[i2]);
	    dist += dev * dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return sqrt(dist);
}

static double mm_maximum(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist = -DBL_MAX;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = fabs(x[i1] - x[i2]);
	    if(dev > dist)
		dist = dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    return dist;
}

static double mm_manhattan(double *x, int nr, int nc, int i1, int i2)
{
    double dist;
    int count, j;

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dist += fabs(x[i1] - x[i2]);
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return dist;
}

static double xmin = 0.0;

static double mm_canberra(double *x, int nr, int nc, int i1, int i2)
{
    double dist, sum, diff;
    int count, j;

    if(xmin == 0.0) {
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	double eps, epsneg,  xmax;
	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp,
	       &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax);
    }

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    sum = fabs(x[i1] + x[i2]);
	    diff = fabs(x[i1] - x[i2]);
	    if (sum > xmin || diff > xmin) {
		dist += diff/sum;
		count++;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return dist;
}

static double mm_dist_binary(double *x, int nr, int nc, int i1, int i2)
{
    int total, count, dist;
    int j;

    total = 0;
    count = 0;
    dist = 0;

    for(j = 0 ; j < nc ; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    if(x[i1] || x[i2]){
		count++;
		if( ! (x[i1] && x[i2]) ) dist++;
	    }
	    total++;
	}
	i1 += nr;
	i2 += nr;
    }

    if(total == 0) return NA_REAL;
    if(count == 0) return 0;
    return (double) dist / count;
}

enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, BINARY };
/* == 1,2,..., defined by order in the R function dist */

void mm_distance(double *x, int *nr, int *nc, double *d, int *iRow, int *numResults, int *method)
{
    int  i;
    size_t size;
    double *tmp;
    double (*distfun)(double*, int, int, int, int) = NULL;

    /* Check to insure that the number of requested results is not */
    /* actually greater then the number of rows being acquired, if */
    /* so, just use the number of rows */
    if (*numResults > *nr) {
      warning("Number of results selected is greater than number of rows, using the number of rows instead\n");
      numResults = nr;
    }

    /* Allocate memory via R */
    tmp = (double *)R_alloc(*nr, sizeof(double));

    switch(*method) {
    case EUCLIDEAN:
	distfun = mm_euclidean;
	break;
    case MAXIMUM:
	distfun = mm_maximum;
	break;
    case MANHATTAN:
	distfun = mm_manhattan;
	break;
    case CANBERRA:
	distfun = mm_canberra;
	break;
    case BINARY:
	distfun = mm_dist_binary;
	break;
    default:
	error("distance(): invalid distance");
    }

    /* Get the distances and feed them into the temp array */
    for(i = 0 ; i < (*nr) ; i++) {
      tmp[i] = distfun(x, *nr, *nc, (*iRow)-1, i);
    }
    /* Run a sort on the temp array */
    qsort((void *)tmp, *nr, sizeof(double), intcompare);


    /* Move over the numResults first cells in the sorted array */
    size = *numResults * sizeof(double);
    d = memmove(d, tmp, size);
}
