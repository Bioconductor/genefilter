#copyright 2001 R. Gentleman

#FILTER FUNCTIONS -- some trivial changes

allNA <- function(x) !all(is.na(x))

anyNA <- function(x) !any(is.na(x))

kOverA <- function(k, A=100, na.rm = TRUE) {
  function(x) {
      if(na.rm)
	x <- x[!is.na(x)]
      sum( x > A ) > k
  }
}

maxA <- function(A=75, na.rm=TRUE) {
    function(x) {max(x, na.rm=na.rm) >= A }
}


pOverA <-  function(A=100, p = .05 ,na.rm = TRUE) {
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
        if(is.na(sdx) || sdx == 0 ) return(FALSE)
	val <- mean(x, na.rm=na.rm)/sdx
        if(val < a ) return(FALSE)
        if(val > b ) return(FALSE)
        return(TRUE)
    }
}

sdom <- function(a=1, b=Inf, na.rm=TRUE) {
    function(x) {
        	sdx <- sd(x, na.rm=na.rm)
        if(is.na(sdx) || sdx == 0 ) return(FALSE)
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
        m1s <- summary.lm(m1) # VC, 15XII01 -- problem w/generic summary in S4?
        fstat <- 1 - pf(m1s$fstat[1], m1s$fstat[2], m1s$fstat[3])
        if( fstat < p )
            return(TRUE)
        return(FALSE)
    }
}

coxfilter <- function(surt, cens, p) {
   autoload("coxph", "survival")
   function(x) {
       srvd <- try(coxph(Surv(surt,cens)~x))
       if( inherits(srvd, "try-error") )
           return(FALSE)
       ltest <- -2*(srvd$loglik[1] - srvd$loglik[2])
       pv <- 1 - pchisq(ltest, 1)
       if( pv < p )
           return(TRUE)
       return(FALSE)
   }
}

ttest <- function(m, p=0.05, na.rm=TRUE) {
    if( length(m) == 1)
        function(x) {
            n <- length(x)
            sub1 <- x[1:m]
            sub2 <- x[(m+1):(m+n)]
            if(na.rm) {
                drop <- is.na(x)
                sub1 <- sub1[!drop[1:m]]
                sub2 <- sub2[!drop[(m+1):(m+n)]]
            }
            t.test(sub1, sub2 )$p.value < p
        }
    else
        function(x) {
            if(na.rm) {
                drop <- is.na(x) | is.na(m)
                x<- x[!drop]
                m<- m[!drop]
            }
            t.test(x~m)$p.value < p
        }
  }


##a filter based on gaps

gapFilter <- function(Gap, IQR, Prop, na.rm=TRUE, neg.rm=TRUE) {
  function(x) {
     if(na.rm) x <- x[!is.na(x)]
     if(neg.rm) x <- x[x>0]
     lenx <- length(x)
     if( lenx < 4 || lenx < Prop+1 )
       return(FALSE)
     srtd <- sort(x)
     lq <- lenx*.25
     uq <- lenx*.75
     if( (srtd[uq] - srtd[lq]) > IQR )
        return(TRUE)
     if(Prop < 1)
        bot <- lenx*Prop
     else
        bot <- Prop
     top <- lenx - bot
     lag1 <- srtd[2:lenx]-srtd[1:(lenx-1)]
     if( max(lag1[bot:top]) > Gap )
       return(TRUE)
    return(FALSE)
  }
}


# Apply type functions


genefilter <- function(expr, flist)
     apply(expr, 1, flist)

filterfun <- function(...) {
     flist <- list(...)
 #let the user supply a list
     if( length(flist) == 1 && is.list(flist[[1]]) )
         flist <- flist[[1]]
     f <- function( x ) {
         for( fun in flist ) {
             fval <- fun(x)
             if( is.na(fval) || ! fval )
                 return(FALSE)
             return(TRUE)
         }
     }
     class(f) <- "filterfun"
     return(f)
 }





