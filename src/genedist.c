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

typedef struct {
    int geneNum;
    double geneDist;
} gene_t;


void detectTies(int geneNum, int nResults, int nRows, gene_t *data) {
    /* Will scan through the first nResults+1 distances in the */
    /* data array, and if it detects any ties, will flag a R */
    /* warning */
    int i; /* Loop indices */
    
    /* If nResults == nRows, do not exceed nResults - otherwise exceed it */
    /* by 1 in order to see if there were trailing ties */
    if (nResults == nRows) {
	nResults = nRows-1;
    }
    
    for (i = 1; i < nResults; i++) {
	if (data[i].geneDist == data[i+1].geneDist) {
	    warning("There are distance ties in the data for gene %d\n",geneNum);
	    break;
	}
    }
}

static int distCompare(const void *p1, const void *p2)
{
    const gene_t *i = p1;
    const gene_t *j = p2;

    if (!R_FINITE(i->geneDist ))
      return(1);
    if (!R_FINITE(j->geneDist))
      return(-1);

    if (i->geneDist > j->geneDist) 
	return (1);
    if (i->geneDist < j->geneDist) 
	return (-1);
    return (0);
    
}

static double mm_correlation(double *x, int nr, int nc, int i1, int i2) {
    int i; /* Loop index */
    int a,b; /* Used as array indices for i1 and i2 */
    double xAvg, yAvg; /* Averages of the i1 and i2 rows */
    double upTot = 0; /* Upper summation */
    double botTotL, botTotR; /* The lower two summations */
    double botVal; /* Bottom value for Rho */
    double Rho, dist;

    botTotL = botTotR = 0;
    xAvg = yAvg = 0;
    a = i1;
    b = i2;
    
    /* Calculate the averages for the i1 and i2 rows */
    for (i = 0; i < nc; i++) {
	if (R_FINITE(x[a])) {
	    xAvg += x[a];
	}
	if (R_FINITE(x[b])) {
	    yAvg += x[b];
	}
	a += nr;
	b += nr;
    }
    xAvg /= (double)nc;
    yAvg /= (double)nc;
    /* Reset a & b */
    a = i1; b = i2;
    
    /* Build up the three summations in the equation */
    for (i = 0; i < nc; i++) {
	if (R_FINITE(x[a]) && R_FINITE(x[b])) {
	    upTot += ((x[a] - xAvg) * (x[b] - yAvg));
	    botTotL += pow((x[a] - xAvg),2);
	    botTotR += pow((x[b] - yAvg),2);
	}
	a += nr;
	b += nr;    
    }

    /* Compute Rho & Distance (1 - R) */
    botVal = sqrt((botTotL * botTotR));
    Rho = upTot / botVal;
    dist = 1 - Rho;
    
    return(dist);
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

enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, CORRELATION, BINARY };
/* == 1,2,..., defined by order in the R function dist */

void mm_distance(double *x, int *nr, int *nc, int *g, double *d, 
		 int *iRow, int *nInterest, int *nResults, int *method) {
    
    /*
      x -> Data Array
      nr -> Number of rows in X
      nc -> number of columns in X
      g -> The nResults closest genes to the genes of interest
      d -> The distances of the genes from g, 1 to 1 mapping
      iRow -> rows of X that we are interested in
      nInterest -> Number of elements in iRow
      nResults -> The top X results to pass back
      method -> which distance method to use
    */
    
    int  i,j, k;  /* Loop indices */
    int baseIndex; /* Used to index data arrays */
    gene_t *tmp; /* Temporary array to hold the distance data */
    double (*distfun)(double*, int, int, int, int) = NULL;

    /* Sanity check the nResults vs. number of rows in the data */
    if (*nResults > *nr) {
	warning("Number of results selected is greater than number of rows, using the number of rows instead\n");
	nResults = nr;
    }
    
    /* Size of tmp == *nr, as each gene we're interested in will generate */
    /* nr number of distance points */
    tmp = (gene_t *)R_alloc(*nr, sizeof(gene_t));
    
    /* Determine which distance function to use */
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
    case CORRELATION:
	distfun = mm_correlation;
	break;
    case BINARY:
	distfun = mm_dist_binary;
	break;
    default:
	error("distance(): invalid distance");
    }
/**/    
    for (j = 0; j < *nInterest; j++) {  
	/* Get the distances for this gene, store in tmp array */

	for(i = 0 ; i < (*nr) ; i++) {
	    tmp[i].geneNum = i; 
	    tmp[i].geneDist = distfun(x, *nr, *nc, iRow[j]-1, i);       
	}
	
	/* Run a sort on the temp array */
	qsort(tmp, *nr, sizeof(gene_t), distCompare);    

	/* Detect any ties */
	detectTies(iRow[j], *nResults, *nr, tmp); 

	/* Copy the 1<->nResults data points into the final array */
	baseIndex = *nResults * j;
	for (k = 1; k <= *nResults; k++) {
	    g[baseIndex + (k-1)] = tmp[k].geneNum; 
	    d[baseIndex + (k-1)] = tmp[k].geneDist; 
	}
    }
}



