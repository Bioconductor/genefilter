/*
 * F. Hahne  10/26/2005
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

void pAUC_c(double *data, int nrd, int ncd, 
            double *cutp, int ncc, 
            int *truth, 
	    double *spec, double *sens, double *area, double *p) {

  int i, j, k, pred, d, rsum, csum, rcount, ccount;
  double *x, *y;
  double a, tmp, lim;

  x   = (double *) R_alloc(ncc+1, sizeof(double));
  y   = (double *) R_alloc(ncc+1, sizeof(double));

  /* this code computes roc for a given n * n matrix at given 
     cut points */
  for(k=0; k<nrd; k++){   /* iterate over rows (genes) */
    for(i=k; i<ncc*nrd; i+=nrd){   /* iterate over cut points */
      rsum = csum = rcount = ccount = 0;  
      for(j=k, d=0; j<nrd*ncd; j+=nrd, d++){   /* iterate over columns (samples) */
		  pred = (data[j] > cutp[i]) ? 1 : 0;
		  if(truth[d] == 1){
			 rsum += pred;
			 rcount++;
		  }
		  else{
			 csum+=(1-pred);
			 ccount++;
		  }
      }   /* for j (columns)*/
      sens[i] = (double)rsum/rcount;
      spec[i] = (double)csum/ccount;
    }   /* for i (cutpoints)*/
    

    /* this computes pAUC for roc curve in row k*/
    for(i=k,d=0; i<ncc*nrd; i+=nrd,d++){
      x[d] = 1 - spec[i];
      y[d] = sens[i];
    }/* for i,d */


    if(x[0] > x[d]){   /* reverse order if necessary */ 
      for(i=0, j=d-1; i<(d/2); i++, j--){
        tmp=x[i]; x[i]=x[j]; x[j]=tmp;
        tmp=y[i]; y[i]=y[j]; y[j]=tmp;
      } 
    }

    /* fill x and y to span the whole interval [0,1] */
    x[ncc+1]=1;
    y[ncc+1]=y[ncc];

    /* compute area by trapezoidal rule*/
    lim = x[0] < (*p) ? x[0] : *p;
    a = (lim*y[0])/2;
    i=1;
    while(x[i] < (*p)){
      a += ((x[i]-x[i-1])*(y[i]-y[i-1])/2)+((x[i]-x[i-1])*y[i-1]);
      i++;
    }
    if(i > 2)
      a += (((*p)-x[i-1])*(y[i]-y[i-1])/2)+(((*p)-x[i-1])*y[i-1]);
    area[k] = a;   
  }
}


/*-----------------------------------------------------------------
   interface to R with arguments:
     data :    matrix of numerics
     cutpts:   matrix with treshholds for ROC curve calculation
     truth:    int with values 0 and 1, defining the real classification
     p:        numeric in 0<p<1, limit to integrate pAUC to 
------------------------------------------------------------------*/

SEXP pAUC(SEXP _data, SEXP _cutpts, SEXP _truth, SEXP _p)
{ 
  SEXP dimData;  /* dimensions of data */
  SEXP dimCutpts;  /* dimensions of cutpts */
  SEXP res, namesres;      /* return value: a list */
  SEXP spec, sens, area;  /* list elements for constructing 
                              the return value */
  SEXP dim; /* dimensions for spec and sens matrices in return value */

  double *data;
  double *cutp;
  int *truth;
  double *p;
  int nrd, ncd;  /* dimensions of data    */
  int nrc, ncc;  /* dimensions of cutpts  */
  int i;
 

  /* check input argument data */
  PROTECT(dimData = getAttrib(_data, R_DimSymbol));
  if((!isReal(_data)) | isNull(dimData) | (LENGTH(dimData)!=2))
      error("Invalid argument 'data': must be a real matrix."); 
  data   = REAL(_data);
  nrd = INTEGER(dimData)[0];
  ncd = INTEGER(dimData)[1];
  UNPROTECT(1);          
  /* done with dimData */

  PROTECT(dimCutpts = getAttrib(_cutpts, R_DimSymbol));
  if((!isReal(_data)) | isNull(dimCutpts) | (LENGTH(dimCutpts)!=2))
      error("Invalid argument 'cutpts': must be a real matrix."); 
  cutp   = REAL(_cutpts);
  nrc  = INTEGER(dimCutpts)[0];
  ncc  = INTEGER(dimCutpts)[1];
  UNPROTECT(1);          
  if(nrc!=nrd)
      error("nrc and nrd must be the same."); 
  /* done with dimCutpts */

  /* check input argument truth */
  if(!isInteger(_truth))
      error("'truth' must be an integer.");
  if(length(_truth)!=ncd) {
    error("length(truth) and ncol(data) should be the same.");
  }
  truth = INTEGER(_truth);
  for(i=0; i<ncd; i++)
      if(! (R_IsNA(truth[i]) || ((truth[i]>=0)&&(truth[i]<=1))) )
	  error("Elements of 'truth' must be 0 or 1.");
  /* done with truth */

  /* check input argument p */
  if(!isReal(_p) || length(_p)!=1) 
    error("'p' must be numeric.");
  p = REAL(_p);
  if(((*p)<0)||((*p)>1))
    error("'p' must be between 0 and 1.");
 /* done with p */


  /* allocate memory for return values */
  PROTECT(spec = allocVector(REALSXP, nrd*ncc));
  PROTECT(sens = allocVector(REALSXP, nrd*ncc));
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = nrd;
  INTEGER(dim)[1] = ncc;
  SET_DIM(spec, dim); 
  SET_DIM(sens, dim); 

  PROTECT(area = allocVector(REALSXP, nrd));

  /* Do it! */
  /* note nrc is the same as nrd */
  pAUC_c(data, nrd, ncd, cutp, ncc, truth, REAL(spec), REAL(sens), REAL(area), p);

  /* return value: a list with  elements spec sens and pAUC */
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, spec);
  SET_VECTOR_ELT(res, 1, sens);
  SET_VECTOR_ELT(res, 2, area);


  PROTECT(namesres = allocVector(STRSXP, 3));
  SET_STRING_ELT(namesres, 0, mkChar("spec"));
  SET_STRING_ELT(namesres, 1, mkChar("sens"));
  SET_STRING_ELT(namesres, 2, mkChar("pAUC"));
  setAttrib(res, R_NamesSymbol, namesres);

  UNPROTECT(6); /* done with res, namesres, spec, sens, dim, pAUC */
  return(res);
}
