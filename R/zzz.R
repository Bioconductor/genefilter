
.First.lib <- function(lib, pkgname, where) 
  {
  library.dynam("genefilter", pkgname, lib)
  require(methods)
  require(Biobase)
    if(missing(where)) {
        where <- match(paste("package:", pkgname, sep=""), search())
        if(is.na(where)) {
            warning(paste("Not a package name: ",pkgname))
            return()
        }
        where <- pos.to.env(where)
    }
  .initFinder(where)
  cacheMetaData(as.environment(where))
  }

