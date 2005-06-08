/*
 * Copyright W. Huber 2005, all rights reserved
 */
 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 

#include <stdlib.h>

/* #define DEBUG */

char errmsg[256];

/*-----------------------------------------------------------------
  which=0:  t-test by row
  which=1:  t-test by column
-----------------------------------------------------------------*/
void rowcolttests_c(double *x, int *fac, int nr, int nc, int no, int nt, 
                    int which, int nrgrp, double *statistic, double *dm, 
                    double *df) {

    int i, j, grp, dof;
    double z;

    /* Currently we only provide for one- and two-sample t-tests (nrgrp=1 or
       2), but it should be possible to generalize this code to more samples
       (F-test) without too many changes */

    int n[2];   
    double mean[2];
    double* s[2];
    double* ss[2];

    if(nrgrp>2)
	error("Please do not use 'nrgrp' >2 with 'rowcolttests'");

    /* allocate and initialize storage for intermediate quantities
       (namely first and second moments for each group) */
    for(grp=0; grp<nrgrp; grp++) {
	s[grp]  = (double*) R_alloc(nt, sizeof(double));
	ss[grp] = (double*) R_alloc(nt, sizeof(double));
	n[grp]  = 0;
	for(i=0; i<nt; i++)  
	    s[grp][i] = ss[grp][i] = 0;
    }

    /* To determine first and second moments, we work through the 
       large matrix x in the order in which it is in memory -
       this may speed up things considerably */
    switch(which) {
	case 0:  /* by row */
	    for(i=0; i<nr; i++) {
		for(j=0; j<nc; j++) {
		    grp = fac[j];
		    if(!R_IsNA(grp)) {
			z = x[i+nr*j];
			s[grp][i]  += z;
			ss[grp][i] += z*z;
		    }
		} /* for j */
	    } /* for i */
	    break;
	case 1:  /* by column */
	    for(i=0; i<nr; i++) {
		grp = fac[i];
		if(!R_IsNA(grp)) {
		    for(j=0; j<nc; j++) {
			z = x[i+nr*j];
			s[grp][j]  += z;
			ss[grp][j] += z*z;
		    } /* for j */
		} /* if */ 
	    } /* for i */
	    break;
	default:
	    error("Bummer!");
    }

    /* determine group sizes */
    for(i=0; i<no; i++) {
	grp = fac[i];
	if(!R_IsNA(grp))
	    n[grp]++;
    }

    /* calculate t statistic */
    for(i=0; i<nt; i++) {
        z = 0;
        for(grp=0; grp<nrgrp; grp++) {
	    mean[grp] = s[grp][i]/n[grp]; /* mean of group 'grp' */
	    z        += ss[grp][i] - mean[grp]*mean[grp]*n[grp];
	}
        switch(nrgrp) {
	case 1:
	    *df = dof = n[0]-1;
	    statistic[i] = mean[0] / sqrt(z/(dof*n[0]));
	    dm[i]        = mean[0];
            break;
	case 2:
	    *df = dof = n[0]+n[1]-2;
	    statistic[i] = (mean[0]-mean[1]) / sqrt(z*(n[0]+n[1])/(dof*n[0]*n[1]));
	    dm[i]        = mean[0]-mean[1];
            break;
	default:
	    error("Bummer!");
	} /* switch */
    }

    return;
} 

/*-----------------------------------------------------------------

   R interface 
   x :    matrix
   fac:   int with values 0 and 1, defining the two groups.
   which: int. For 0, do the tests along the rows, for 1, 
          along the columns 
------------------------------------------------------------------*/
SEXP rowcolttests(SEXP _x, SEXP _fac, SEXP _nrgrp, SEXP _which) 
{
  SEXP dimx;  /* dimensions of x */
  SEXP res, namesres;      /* return value: a list */
  SEXP statistic, dm, df;  /* list elements for constructing 
                              the return value */

  double *x;
  int *fac;
  int i, which, nrgrp;
  int nr;  /* number of rows     */
  int nc;  /* number of columns  */
  int no;  /* number of objects  */
  int nt;  /* number of tests    */

  /* check input argument x */
  PROTECT(dimx = getAttrib(_x, R_DimSymbol));
  if((!isReal(_x)) | isNull(dimx) | (LENGTH(dimx)!=2))
      error("Invalid argument 'x': must be a real matrix."); 
  x   = REAL(_x);
  nr  = INTEGER(dimx)[0];
  nc  = INTEGER(dimx)[1];
  UNPROTECT(1);          
  /* done with dimx */

  /* check input argument which */
  if(!isInteger(_which) || length(_which)!=1) 
      error("'which' must be integer of length 1.");
  which = INTEGER(_which)[0];

  /* check input argument nrgrp */
  if(!isInteger(_nrgrp) || length(_nrgrp)!=1) 
      error("'nrgrp' must be integer of length 1.");
  nrgrp = INTEGER(_nrgrp)[0];

  /* check input argument fac */
  if(!isInteger(_fac))
      error("'fac' must be an integer.");
  switch(which) {
      case 0: 
	  if(length(_fac)!=nc) {
	      sprintf(errmsg, "length(fac)=%d, ncol(x)=%d, should be the same.",
		      length(_fac), nc);
	      error(errmsg);
	  }
          no = nc;
          nt = nr;
	  break;
      case 1:
	  if(length(_fac)!=nr) {
	      sprintf(errmsg, "length(fac)=%d, nrow(x)=%d, should be the same.",
		      length(_fac), nr);
	      error(errmsg);
	  }
          no = nr;
          nt = nc;
	  break;
      default:
	  error("'which' must be 0 or 1.");
  }
  
  fac = INTEGER(_fac);
  for(i=0; i<no; i++)
      if(! (R_IsNA(fac[i]) || ((fac[i]>=0)&&(fac[i]<nrgrp))) )
	  error("Elements of 'fac' must be >=0 and < 'nrgrp'.");
  /* done with fac */

  PROTECT(statistic = allocVector(REALSXP, nt));
  PROTECT(dm        = allocVector(REALSXP, nt));
  PROTECT(df        = allocVector(REALSXP, 1));
  /* currently df is trivial, since we only do the normal t-test with
     equal variances; at some point, var.equal=FALSE might be implemented */

  /* Do it! */
  rowcolttests_c(x, fac, nr, nc, no, nt, which, nrgrp, REAL(statistic), REAL(dm), REAL(df));

  /* return value: a list with two elements, statistic and df */
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, statistic);
  SET_VECTOR_ELT(res, 1, dm);
  SET_VECTOR_ELT(res, 2, df);

  PROTECT(namesres = allocVector(STRSXP, 3));
  SET_STRING_ELT(namesres, 0, mkChar("statistic"));
  SET_STRING_ELT(namesres, 1, mkChar("dm"));
  SET_STRING_ELT(namesres, 2, mkChar("df"));
  setAttrib(res, R_NamesSymbol, namesres);

  UNPROTECT(5); /* done with res, namesres, statistic, dm, df */
  return(res);
}

