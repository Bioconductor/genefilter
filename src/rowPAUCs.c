#include <R.h>
#include <Rdefines.h>

void pAUC(double *data, int *nd, double *cutp, int *nc, int *truth, 
	  double *spec, double *sens, double *area, double *b) {

  int i, j, k, pred, d, rsum, csum, rcount, ccount, lx;
  double x[*nc];
  double y[*nc];
  double m;
  double h;
  double yb;

  /* this code computes roc for a given n * n matrix at given 
     cut points */
  for(k=0; k<nd[0]; k++){   /* iterate over rows (genes) */
    for(i=k; i<*nc*nd[0]; i+=nd[0]){   /* iterate over cut points */
      rsum = csum = rcount = ccount = 0;  
      for(j=k, d=0; j<nd[0]*nd[1]; j+=nd[0], d++){   /* iterate over 
							columns (samples) */
		  pred = data[j] > cutp[i] ? 1 : 0;
		  if(truth[d] == 1){
			 rsum += pred;
			 rcount++;
		  }
		  else{
			 csum+=(1-pred);
			 ccount++;
		  }
      }   /* for j */
      sens[i] = (double)rsum/rcount;
      spec[i] = (double)csum/ccount;
    }   /* for i */
    

    /* this computes pAUC for roc curve in row k*/
    for(i=k,d=0; i<*nc*nd[0]; i+=nd[0],d++){
      x[d] = 1 - spec[i];
      y[d] = sens[i];
    }
    lx = d;
 
    if(x[0] > x[lx]){   /* reverse order if necessary */ 
      double xdummy[lx];  
      double ydummy[lx];
      for(i=0; i<lx; i++){
	xdummy[i] = x[i];
	ydummy[i] = y[i];
      }
      for(i=1; i<lx+1; i++){
	x[i-1] = xdummy[lx-i];
	y[i-1] = ydummy[lx-i];
      }
    }

    double xin[lx];
    double yin[lx];
    d = 0;
    for(i=0; i<lx; i++){  /* subset for iteration */
      if(x[i] <= *b){
	xin[d] = x[i];
	yin[d] = y[i];
	d++;
      }
    }
    lx = d-1;  /* extend subset to include boundaries */
    double xall[lx+2];
    double yall[lx+2];
    xall[0] = 0;
    yall[0] = 0;
    for(i=0; i<=lx; i++){
      xall[i+1] = xin[i];
      yall[i+1] = yin[i];
    }
    m = (*b-x[lx]) / (x[lx+1] - x[lx]);  /* interpolation at right boundary */
    yb = y[lx] + (y[lx+1] - y[lx]) * m;
    xall[lx+2] = *b;
    yall[lx+2] = yb;
    lx = lx+2;
    for(i=0; i<lx; i++){   /* compute area */
      h = xall[i+1] - xall[i];
      area[k] += h * (yall[i+1] + yall[i]);
    }
    area[k] *= 0.5;   
  }
}

