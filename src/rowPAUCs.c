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
  internal c function for calculation of ROC curves and pAUCs
  -----------------------------------------------------------------*/

void ROCpAUC_c(double *data, int nrd, int ncd, double *cutp, int ncc,
	       int *truth, double *spec, double *sens, double *area,
	       double *auc, double *p, int flip) {

    int i, j, k, pred, d, rsum, csum, rcount, ccount;
    double *x, *y;
    double a, ta, tmp, lim, xsum, ysum;

    x   = (double *) R_alloc(ncc+1, sizeof(double));
    y   = (double *) R_alloc(ncc+1, sizeof(double));

    /* this code computes roc for a given n * n matrix at given
       cut points */
    //printf("Computing ROC curves for %d rows at %d cutpoints ...\n", nrd, ncc);
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
	xsum = ysum = 0;
	for(i=k,d=0; i<ncc*nrd; i+=nrd,d++){
	    x[d] = 1 - spec[i];
	    y[d] = sens[i];
	    xsum += x[d];
	    ysum += y[d];
	}/* for i,d */

	/*rotate 180° if necessary*/
	if(flip && xsum > ysum){
	    for(i=k,d=0; i<ncc*nrd; i+=nrd,d++){
		spec[i] = 1 - sens[i];
		sens[i] = x[d];
		x[d] = 1-spec[i];
		y[d] = sens[i];
	    }/* for i,d */
	}
	d--;

	/* reverse order if necessary */
	if(x[0] > x[d]){
	    for(i=0, j=d; i<=(d+1)/2; i++, j--){
		tmp=x[i]; x[i]=x[j]; x[j]=tmp;
		tmp=y[i]; y[i]=y[j]; y[j]=tmp;
	    }
	}
	x[ncc] = 1;
	y[ncc] = y[ncc-1];

	/* compute area by trapezoidal rule*/
	lim = x[0] < (*p) ? x[0] : *p; /*right border of first segment*/
	a = (lim*y[0])/2; /*area of 1. segement (from x1=0 to x2=lim)*/
	i=1;
	while(x[i] < (*p)){
	    a += ((x[i]-x[i-1])*(y[i]-y[i-1])/2) + ((x[i]-x[i-1])*y[i-1]);
	    i++;
	}

	if(i > 2){ /*last segment (from xn to p)*/
	    a += (((*p)-x[i-1])*(y[i]-y[i-1])/2) + (((*p)-x[i-1])*y[i-1]);
	}
	ta = a;
	/*compute full AUC and flip curve if necessary*/
	if((*p) < 1){
	    ta += ((x[i]-(*p))*(y[i]-y[i-1])/2) + ((x[i]-(*p))*y[i-1]);
	    i++;
	    while(i < ncc+1 && x[i] < 1){
		ta += ((x[i]-x[i-1])*(y[i]-y[i-1])/2) + ((x[i]-x[i-1])*y[i-1]);
		i++;
	    }
	    ta += ((1-x[i-1])*(1-y[i-1])/2) + ((1-x[i-1])*y[i-1]);
	}
	if(flip && (*p)==1 && ta < 0.5){ /*rotate 180° if area < 0.5*/
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
  data :    matrix of numerics
  cutpts:   matrix with treshholds for ROC curve calculation
  truth:    int with values 0 and 1, defining the real classification
  p:        numeric in 0<p<1, limit to integrate pAUC to
  ------------------------------------------------------------------*/

SEXP ROCpAUC(SEXP _data, SEXP _cutpts, SEXP _truth, SEXP _p, SEXP _flip)
{
    SEXP dimData;  /* dimensions of data */
    SEXP dimCutpts;  /* dimensions of cutpts */
    SEXP res, namesres;      /* return value: a list */
    SEXP spec, sens, area, auc;  /* list elements for constructing
				    the return value */
    SEXP dim; /* dimensions for spec and sens matrices in return value */

    double *data;
    double *cutp;
    int *truth;
    double *p;
    int flip;
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

    /* check input argument flip */
    if(!isInteger(_flip))
	error("'flip' must be an integer.");
    flip = (int)INTEGER(_flip)[0];
    /* done with flip */

    /* allocate memory for return values */
    PROTECT(spec = allocVector(REALSXP, nrd*ncc));
    PROTECT(sens = allocVector(REALSXP, nrd*ncc));
    PROTECT(dim = allocVector(INTSXP, 2));
    INTEGER(dim)[0] = nrd;
    INTEGER(dim)[1] = ncc;
    SET_DIM(spec, dim);
    SET_DIM(sens, dim);

    PROTECT(area = allocVector(REALSXP, nrd));
    PROTECT(auc = allocVector(REALSXP, nrd));

    /* Do it! */
    /* note nrc is the same as nrd */
    ROCpAUC_c(data, nrd, ncd, cutp, ncc, truth, REAL(spec), REAL(sens),
	      REAL(area), REAL(auc), p, flip);

    /* return value: a list with  elements spec sens and pAUC */
    PROTECT(res = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(res, 0, spec);
    SET_VECTOR_ELT(res, 1, sens);
    SET_VECTOR_ELT(res, 2, area);
    SET_VECTOR_ELT(res, 3, auc);


    PROTECT(namesres = allocVector(STRSXP, 4));
    SET_STRING_ELT(namesres, 0, mkChar("spec"));
    SET_STRING_ELT(namesres, 1, mkChar("sens"));
    SET_STRING_ELT(namesres, 2, mkChar("pAUC"));
    SET_STRING_ELT(namesres, 3, mkChar("AUC"));
    setAttrib(res, R_NamesSymbol, namesres);

    UNPROTECT(7); /* done with res, namesres, spec, sens, dim, pAUC */
    return(res);
}
