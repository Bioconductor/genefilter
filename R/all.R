#copyright 2001 R. Gentleman

#FILTER FUNCTIONS

allNA <- function(x) !all(is.na(x))

anyNA <- function(x) !any(is.na(x))

kltA <- function(k, A=100, na.rm = TRUE) {
  function(x) {
      if(na.rm)
	x <- x[!is.na(x)]
      sum( x > A ) > k
  }
}

maxA <- function(A=75, na.rm=TRUE) {
    function(x) {max(x, na.rm=na.rm) >= A }
}


poverA <-  function(A=100, p = .05 ,na.rm = TRUE) {
  function(x) {
      if(na.rm)
	 x<-x[!is.na(x)]
      sum( x > A )/length(x) > p
  }
}

#coefficient of variation
cv <- function(a=1, b=Inf, na.rm=TRUE) {
  function(x) {
	sdx <- sd(x, na.rm=na.rm)
        if(sdx == 0 ) return(FALSE)
	val <- mean(x, na.rm=na.rm)/sdx
        if(val < a ) return(FALSE)
        if(val > b ) return(FALSE)
        return(TRUE)
    }
}

sdom <- function(a=1, b=Inf, na.rm=TRUE) {
    function(x) {
        	sdx <- sd(x, na.rm=na.rm)
        if(sdx == 0 ) return(FALSE)
	val <- sdx/mean(x, na.rm=na.rm)
        if(val < a ) return(FALSE)
        if(val > b ) return(FALSE)
        return(TRUE)
            }
}

Anova <- function(cov, p=0.05, na.rm=TRUE)
{
    function(x) {
        if( na.rm ) {
            drop <- is.na(x)
            x <- x[!drop]
            cov <- cov[!drop]
        }
        m1 <- lm(x~cov)
        m1s <- summary(m1)
        fstat <- 1 - pf(m1s$fstat[1], m1s$fstat[2], m1s$fstat[3])
        if( fstat < p )
            return(TRUE)
        return(FALSE)
    }
}

coxfilter <- function(surt, cens, p)
 {
	function(x) {
	   srvd <- coxph(Surv(surt,cens)~x)
	   ltest <- -2*(srvd$loglik[1] - srvd$loglik[2])
           pv <- 1 - pchisq(ltest, 1)
           if( pv < p )
             return(TRUE)
           return(FALSE)
       }
}


# normalize within rows

standardize <- function(x, na.rm=TRUE) {
    sdx<- sd(x, na.rm=na.rm)
    if(sdx == 0 )
        x <- 0
    else
        x <- (x-mean(x, na.rm=na.rm))/sdx
    return(x)
}

# Apply type functions

genefilter <- function(expr, flist) {
      for(filt in flist) {
          tmp <- apply(expr, 1, filt)
          expr <- expr[tmp,]
      }
      return(expr)
  }


genefilter2 <- function(expr, flist)
     apply(expr, 1, flist)

filterfun <- function(...) {
     flist <- list(...)
     f <- function( x ) {
         for( fun in flist )
             if( ! fun(x) )
                 return(FALSE)
         return(TRUE)
     }
     class(f) <- "filterfun"
     return(f)
 }
