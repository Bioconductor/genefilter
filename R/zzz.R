.onLoad <- function(lib, pkgname) {
    if(.Platform$OS.type == "windows" && interactive() && .Platform$GUI ==  "Rgui"){
      addVigs2WinMenu("genefilter")
   }
 }

.onUnload <- function( libpath ) {
  library.dynam.unload( "genefilter", libpath )
}
