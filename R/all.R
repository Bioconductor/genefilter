"allNA" <- function(x) !all(is.na(x))

"fundiff" <- function(x, fun, na.rm=TRUE) {
	drop <- is.na(x)
	if( na.rm) {
	   x <- x[!drop]
           var <- var[!drop]
        }
        gps <- split(x, var)
        if(length(gps)>2)
           stop("wrong number of groups")
        if(length(gps)<2)
	   return(NA)
        mns <- lapply(gps, fun$fun, fun$args)
        if( abs(mns[[1]]-mns[[2]]) < fun$accept )
          return(FALSE)
        TRUE
   }

"genefilter" <-
function(expr, flist) {
      for(filt in flist) {
          tmp <- apply(expr, 1, filt)
          expr <- expr[tmp,]
      }
      return(expr)
  }

"kltA" <- function(k, A=100, na.rm = TRUE) {
  function(x) {
      if(na.rm)
	x <- x[!is.na(x)]
      sum( x > A ) > k
  }
}

"maxA" <-
function(A=75, na.rm=TRUE) {
    function(x) {max(x, na.rm=na.rm) >= A }
}

"mk2fun" <-
function(fun2g, var2g) {
    fundiff <- function(x, na.rm=TRUE) {
	drop <- is.na(x)
	if( na.rm) {
	   x <- x[!drop]
           var2g <- var2g[!drop]
        }
        gps <- split(x, var2g)
        if(length(gps)>2)
           stop("wrong number of groups")
        if(length(gps)<2)
	   return(NA)
        mns <- lapply(gps, fun2g$fun, fun2g$args)
        if( abs(mns[[1]]-mns[[2]]) < fun2g$accept )
          return(FALSE)
        TRUE
     }
     return(fundiff)
   }

"poverA" <-  function(A=100, p = .05 ,na.rm = TRUE) {
  function(x) {
      if(na.rm)
	 x<-x[!is.na(x)]
      sum( x > A )/length(x) > p
  }
}

"sd" <- function(x, na.rm=T) {
  sqrt(var(x, na.rm=na.rm))
}

#coefficient of variation
"cv" <- function(a=1, b=Inf, na.rm=TRUE) {
  function(x) {
	sdx <- sd(x, na.rm=na.rm)
        if(sdx == 0 ) return(FALSE)
	val <- mean(x, na.rm=na.rm)/sdx
        if(val < a ) return(FALSE)
        if(val > b ) return(FALSE)
        return(TRUE)
    }
}

anova <- function(cov, p=0.05, na.rm=TRUE)
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


# normalize within rows

 stdize <- function(x, na.rm=TRUE) {
	sdx<- sd(x, na.rm=na.rm)
        if(sdx == 0 ) 
             x <- 0
        else 
           x <- (x-mean(x, na.rm=na.rm))/sdx
        return(x)
  }
          
