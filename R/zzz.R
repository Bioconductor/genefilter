
.First.lib <- function(lib, pkgname, where) 
  {
  library.dynam("genefilter", pkgname, lib)
    if(missing(where)) {
        where <- match(paste("package:", pkgname, sep=""), search())
        if(is.na(where)) {
            warning(paste("Not a package name: ",pkgname))
            return()
        }
        where <- pos.to.env(where)
    }
  cacheMetaData(as.environment(where))
  if(.Platform$OS.type == "windows" && require("Biobase") && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("genefilter")
   }
 }

