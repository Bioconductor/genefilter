/*
 * F. Hahne  10/24/2006
 */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 

#include <stdlib.h>

/*-----------------------------------------------------------------
internal c function for calculation of pAUCs
-----------------------------------------------------------------*/

void pAUC_c(double *spec, double *sens, double *area, double *auc, double *p,
	     int columns, int rows) {

  int i, j, k, d;
  double *x, *y;
  double a, ta, tmp, lim, xsum ,ysum;

  x   = (double *) R_alloc(columns+1, sizeof(double));
  y   = (double *) R_alloc(columns+1, sizeof(double));


  /* this computes pAUC for roc curve in row k*/
  printf("Computing area under the curve...\n");
  for(k=0; k<rows; k++){   /* iterate over rows (genes) */
    for(i=k*columns,d=0; i<k*columns+columns; i++,d++){
      x[d] = 1 - spec[i];
      y[d] = sens[i];
      xsum += x[d];
      ysum += y[d];
    }/* for i,d */
    /*rotate 180° if necessary*/
    if(xsum > ysum){
      for(i=k*columns,d=0; i<k*columns+columns; i++,d++){
	spec[i] = 1 - sens[i];
        sens[i] = x[d];
	x[d] = 1-spec[i];
	y[d] = sens[i];
      }/* for i,d */
    }
    d--;

    /* reverse order if necessary */
    if(x[0] > x[d]){    
      for(i=0, j=d; i<=d/2; i++, j--){
        tmp=x[i]; x[i]=x[j]; x[j]=tmp;
        tmp=y[i]; y[i]=y[j]; y[j]=tmp;
      } 
      for(i=0; i<=d; i++){
      }
    } 
    x[columns]=1;
    y[columns]=y[columns-1];
  

    /* compute area by trapezoidal rule*/
    lim = x[0] < (*p) ? x[0] : *p; /*right border of first segment*/
    a = (lim*y[0])/2; /*area of 1. segement (from x1=0 to x2=lim)*/
    i=1;
    while(x[i] < (*p)){
      a += ((x[i]-x[i-1])*(y[i]-y[i-1])/2) + ((x[i]-x[i-1])*y[i-1]);
      i++;
    }
   
    if(i > 2) /*last segment (from xn to p)*/
      a += (((*p)-x[i-1])*(y[i]-y[i-1])/2) + (((*p)-x[i-1])*y[i-1]);
    ta = a;
    /*compute full AUC and flip curve if necessary*/ 
    if((*p) < 1){
      ta =  ta += ((x[i]-(*p))*(y[i]-y[i-1])/2) + ((x[i]-(*p))*y[i-1]);
      i++;
      while(i < columns+1 && x[i] < 1){
	ta =  ta += ((x[i]-x[i-1])*(y[i]-y[i-1])/2) + ((x[i]-x[i-1])*y[i-1]);
	i++;
      }
      ta =  ta += ((1-x[i-1])*(1-y[i-1])/2) + ((1-x[i-1])*y[i-1]);
    }else{
      d=1;
    }
    if((*p)==1 && ta < 0.5){ /*rotate 180° if area < 0.5*/
      a = (*p) - a;
      ta = 1-ta;
    }
    if(a>1)
      error("Internal error");
    area[k] = a;
    auc[k] = ta;
  }
}







/*-----------------------------------------------------------------
   interface to R with arguments:
     spec :    matrix of numerics (specificity)
     sens:   matrix of numerics (sensitivity)
     p:        numeric in 0<p<1, limit to integrate pAUC to 
------------------------------------------------------------------*/

SEXP pAUC(SEXP _spec, SEXP _sens, SEXP _p)
{ 
  SEXP res, namesres;      /* return value: a list */
  SEXP area;  /* list element for constructing 
                              the return value */
  SEXP auc;  /* list element for constructing 
                              the return value */
  SEXP dimSpec; /* dimensions for spec and sens matrices */
  SEXP dimSens;

  double *spec;
  double *sens;
  double *p;
  int rows, columns;  /* dimensions of spec and sens  */
  int i;
 

  /* check input argument spec */
  PROTECT(dimSpec = getAttrib(_spec, R_DimSymbol));
  if((!isReal(_spec)) | isNull(dimSpec) | (LENGTH(dimSpec)!=2))
      error("Invalid argument 'spec': must be a real matrix."); 
  spec   = REAL(_spec);
  rows = INTEGER(dimSpec)[1];
  columns = INTEGER(dimSpec)[0];
  UNPROTECT(1);          
  /* done with spec */

  /* check input argument sens */
  PROTECT(dimSens = getAttrib(_sens, R_DimSymbol));
  if((!isReal(_sens)) | isNull(dimSens) | (LENGTH(dimSens)!=2))
      error("Invalid argument 'sens': must be a real matrix."); 
  sens   = REAL(_sens);
  if(rows != INTEGER(dimSens)[1] | columns != INTEGER(dimSens)[0])
    error("'spec' and 'sens' must be matrices with equal dimensions");
  UNPROTECT(1);          
  /* done with sens */
   
  /* check input argument p */
  if(!isReal(_p) || length(_p)!=1) 
    error("'p' must be numeric.");
  p = REAL(_p);
  if(((*p)<0)||((*p)>1))
    error("'p' must be between 0 and 1.");
 /* done with p */


  /* allocate memory for return values */
  PROTECT(area = allocVector(REALSXP, columns));
  PROTECT(auc = allocVector(REALSXP, columns));

  /* Do it! */
  pAUC_c(spec, sens, REAL(area), REAL(auc), p, rows, columns);

  /* return value: a list with elements spec sens and area */
  PROTECT(res = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, area);
  SET_VECTOR_ELT(res, 1, auc);


  PROTECT(namesres = allocVector(STRSXP, 2));
  SET_STRING_ELT(namesres, 0, mkChar("pAUC"));
  SET_STRING_ELT(namesres, 1, mkChar("AUC"));
  setAttrib(res, R_NamesSymbol, namesres);

  UNPROTECT(4); /* done with res, namesres, pAUC, auc */
  return(res);
}
