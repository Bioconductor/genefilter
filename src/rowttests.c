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
-----------------------------------------------------------------*/
void rowcolttests_c(double *x, int *fac, int nr, int nc, int no, int nt, 
                    int which, double *statistic, double *df) {

    int i, j, grp, k, dof;
    double z, m0, m1;
    int n[2];
    double* s[2];
    double* ss[2];

    s[0]  = (double*) R_alloc(nt, sizeof(double));
    s[1]  = (double*) R_alloc(nt, sizeof(double));
    ss[0] = (double*) R_alloc(nt, sizeof(double));
    ss[1] = (double*) R_alloc(nt, sizeof(double));

    /* initialize */
    n[0] = n[1] = 0;
    for(i=0; i<nt; i++)  
	s[0][i] = s[1][i] = ss[0][i] = ss[1][i] = 0;

    /* To determine first and second moments, we work through the 
       large matrix x in the order in which it is in memory; as 
       the price for this, we have if statements in the inner 
       loop - but hopefully these stay within the CPU. */
    for(i=0; i<nr; i++) {
	for(j=0; j<nc; j++) {
            grp = (which==0) ? fac[j] : fac[i];
            k   = (which==0) ?   i    :    j  ;
	    if(!R_IsNA(grp)) {
	      z = x[i+nr*j];
	      s[grp][k]  += z;
	      ss[grp][k] += z*z;
	    }
	} /* for j */
    } /* for k */

    /* determine group sizes */
    for(i=0; i<no; i++) {
	grp = fac[i];
	if(!R_IsNA(grp))
	    n[grp]++;
    }

    /* calculate t statistic */
    *df = dof = n[0]+n[1]-2;
    for(i=0; i<nt; i++) {
        m0 = s[0][i]/n[0]; /* mean of group 0 */
        m1 = s[1][i]/n[1]; /* mean of group 1 */
        z  = ss[0][i] - m0*m0*n[0] + ss[1][i] - m1*m1*n[1];
	statistic[i] = (m1-m0)/sqrt(z*(n[0]+n[1])/(dof*n[0]*n[1]));
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
SEXP rowcolttests(SEXP _x, SEXP _fac, SEXP _which) 
{
  SEXP dimx;  /* dimensions of x */
  SEXP res, namesres;    /* return value: a list */
  SEXP statistic, df;    /* list elements for constructing 
                            the return value */

  double *x;
  int *fac;
  int i, which;
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
      if(!(R_IsNA(fac[i])||(fac[i]==0)||(fac[i]==1)))
	  error("Elements of 'fac' must be 0 or 1.");
  /* done with fac */

  PROTECT(statistic = allocVector(REALSXP, nt));
  PROTECT(df        = allocVector(REALSXP, 1));
  /* currently df is trivial, since we only do the normal t-test with
     equal variances; at some point, var.equal=FALSE might be implemented */

  /* Do it! */
  rowcolttests_c(x, fac, nr, nc, no, nt, which, REAL(statistic), REAL(df));

  /* return value: a list with two elements, statistic and df */
  PROTECT(res = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, statistic);
  SET_VECTOR_ELT(res, 1, df);

  PROTECT(namesres = allocVector(STRSXP, 2));
  SET_STRING_ELT(namesres, 0, mkChar("statistic"));
  SET_STRING_ELT(namesres, 1, mkChar("df"));
  setAttrib(res, R_NamesSymbol, namesres);

  UNPROTECT(4); /* done with res, namesres, statistic, df */
  return(res);
}

