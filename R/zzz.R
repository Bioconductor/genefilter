.onLoad <- function(lib, pkgname) {
    require("methods", quietly=TRUE) || stop("methods package not found")
    
    if(.Platform$OS.type == "windows" && require("Biobase") && interactive()
        && .Platform$GUI ==  "Rgui"){
      addVigs2WinMenu("genefilter")
   }
 }

