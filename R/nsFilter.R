setMethod("nsFilter", "ExpressionSet",
          function(eset,
                   require.entrez=TRUE,
                   require.symbol=TRUE,
                   require.GOBP=FALSE,
                   require.GOCC=FALSE,
                   require.GOMF=FALSE,
                   remove.dupEntrez=TRUE,
                   var.func=IQR, var.cutoff=0.5,
                   var.filter=TRUE)
          {
              if (!is.function(var.func))
                stop("'var.func' must be a function")

              annPkg <- annotation(eset)
              if (nchar(annPkg) == 0)
                stop("'eset' must have a valid annotation slot")
              getAnnEnv <- function(map) {
                  get(paste(annPkg, map, sep=""))
              }

              nfeat <- function(eset) length(featureNames(eset))
              filter.log <- new.env(parent=emptyenv())

              requireID <- function(eset, map) {
                  IDs <- mget(featureNames(eset), envir=getAnnEnv(map))
                  haveID <- names(IDs)[sapply(IDs, function(x) !is.na(x))]
                  logvar <- paste("numRemoved", map, sep=".")
                  assign(logvar, nfeat(eset) - length(haveID), envir=filter.log)
                  eset[haveID, ]
              }

              if (require.entrez) {
                  eset <- requireID(eset, "ENTREZID")
              }

              if (require.symbol) {
                  eset <- requireID(eset, "SYMBOL")
              }

              filterGO <- function(eset, ontology) {
                  haveGo <- sapply(mget(featureNames(eset), getAnnEnv("GO")),
                                   function(x) {
                                       if (length(x) == 1 && is.na(x))
                                         FALSE
                                       else {
                                           onts <- sapply(x, function(z) z$Ontology)
                                           ontology %in% onts
                                       }
                                   })
                  logvar <- paste("numNoGO", ontology, sep=".")
                  assign(logvar, sum(!haveGo), envir=filter.log)
                  eset[haveGo, ]
              }

              if (require.GOBP) {
                  eset <- filterGO(eset, "BP")
              }

              if (require.GOCC) {
                  eset <- filterGO(eset, "CC")
              }

              if (require.GOMF) {
                  eset <- filterGO(eset, "MF")
              }

              if (var.filter) {
                  esetIqr <- apply(exprs(eset), 1, var.func)
                  selected <- esetIqr > var.cutoff
                  eset <- eset[selected, ]
                  logvar <- "numLowVar"
                  assign(logvar, sum(!selected), filter.log)
              }

              if (remove.dupEntrez) {
                  ## Reduce to unique probe <--> gene mapping here by keeping largest IQR
                  ## We will want "unique genes" in the non-specific filtered gene
                  ## set.
                  numNsWithDups <- nfeat(eset)
                  nsFilteredIqr <- esetIqr[selected]
                  uniqGenes <- findLargest(featureNames(eset), nsFilteredIqr,
                                           annotation(eset))
                  eset <- eset[uniqGenes, ]
                  logvar <- "numDupsRemoved"
                  assign(logvar, numNsWithDups - nfeat(eset), envir=filter.log)
              }
              numSelected <- length(featureNames(eset))
              list(eset=eset, filter.log=as.list(filter.log))
          })
