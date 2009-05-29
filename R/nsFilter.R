##RG introduces two new functions, varFilter that does what nsFilter
##was supposed to, but never did, and featureFilter that does the only
##useful stuff that nsFilter does
rowIQRs <- function(eSet) {
  numSamp <- ncol(eSet)
  lowQ <- rowQ(eSet, floor(0.25 * numSamp))
  upQ <- rowQ(eSet, ceiling(0.75 * numSamp))
  upQ - lowQ
}

varFilter <- function(eset, var.func=IQR, var.cutoff=0.5,filterByQuantile=TRUE
)
{
    if (deparse(substitute(var.func)) == "IQR") {
        vars <- rowIQRs(eset)
    } else {
        vars <- apply(exprs(eset), 1, var.func)
    }
    if (filterByQuantile) {
        if ( 0 < var.cutoff && var.cutoff < 1 ) {
            quant = quantile(vars, probs = var.cutoff)
            selected = !is.na(vars) & vars > quant
        } else stop("Cutoff Quantile has to be between 0 and 1.")
    } else {
        selected <- !is.na(vars) & vars > var.cutoff
    }
    eset <- eset[selected, ]
}

featureFilter <- function(eset, require.entrez=TRUE,
                   require.GOBP=FALSE, require.GOCC=FALSE,
                   require.GOMF=FALSE, remove.dupEntrez=TRUE,
                   feature.exclude="^AFFX") {

    annChip <- annotation(eset)
    if (nchar(annChip) == 0) stop("'eset' must have a valid annotation slot")
    getAnnEnv <- function(map) { getAnnMap(map=map, chip=annChip) }
    
    nfeat <- function(eset) length(featureNames(eset))
    requireID <- function(eset, map) {
        IDs <- mget(featureNames(eset), envir=getAnnEnv(map))
        haveID <- names(IDs)[sapply(IDs, function(x) !is.na(x))]
         eset[haveID, ]
    }
    if (require.entrez)
        eset <- requireID(eset, "ENTREZID")

    filterGO <- function(eset, ontology) {
        haveGo <- sapply(mget(featureNames(eset), getAnnEnv("GO")),
                         function(x) {
                             if (length(x) == 1 && is.na(x))
                                 FALSE
                             else {
                                 onts <- subListExtract(x, "Ontology", simplify=TRUE)
                                 ontology %in% onts
                             }
                         })
        eset[haveGo, ]
    }
    
    if (require.GOBP) 
        eset <- filterGO(eset, "BP")
    if (require.GOCC) 
        eset <- filterGO(eset, "CC")
    if (require.GOMF) 
        eset <- filterGO(eset, "MF")

    if (length(feature.exclude)) {
        fnms <- featureNames(eset)
        badIdx <- integer(0)
        for (pat in feature.exclude) {
            if (nchar(pat) == 0)
                next
            badIdx <- c(grep(pat, fnms), badIdx)
        }
        if (length(badIdx)) {
            badIdx <- unique(badIdx)
            eset <- eset[-badIdx, ]
        }
    }
    if (remove.dupEntrez ) {
     ## Reduce to unique probe <--> gene mapping here by keeping largest IQR
     ## We will want "unique genes" in the non-specific filtered gene
     ## set.
        uniqGenes <- findLargest(featureNames(eset), rowIQRs(eset),
                                 annotation(eset))
        eset <- eset[uniqGenes, ]
    }
    eset
}


setMethod("nsFilter", "ExpressionSet",
          function(eset,
                   require.entrez=TRUE,
                   require.GOBP=FALSE,
                   require.GOCC=FALSE,
                   require.GOMF=FALSE,
                   remove.dupEntrez=TRUE,
                   var.func=IQR, var.cutoff=0.5,
                   var.filter=TRUE,
                   filterByQuantile=TRUE,
                   feature.exclude="^AFFX", ...)
          {
              if (!is.function(var.func))
                stop("'var.func' must be a function")

              annChip <- annotation(eset)
              if (nchar(annChip) == 0)
                stop("'eset' must have a valid annotation slot")
              getAnnEnv <- function(map) {
                  getAnnMap(map=map, chip=annChip)
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

              filterGO <- function(eset, ontology) {
                  haveGo <- sapply(mget(featureNames(eset), getAnnEnv("GO")),
                                   function(x) {
                                       if (length(x) == 1 && is.na(x))
                                         FALSE
                                       else {
                                           onts <- subListExtract(x, "Ontology", simplify=TRUE)
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

              if (length(feature.exclude)) {
                  fnms <- featureNames(eset)
                  badIdx <- integer(0)
                  for (pat in feature.exclude) {
                      if (nchar(pat) == 0)
                        next
                      badIdx <- c(grep(pat, fnms), badIdx)
                  }
                  if (length(badIdx)) {
                      badIdx <- unique(badIdx)
                      eset <- eset[-badIdx, ]
                      logvar <- "feature.exclude"
                      assign(logvar, length(badIdx), filter.log)
                  }
              }


              if (remove.dupEntrez) {
                  ## Reduce to unique probe <--> gene mapping here by keeping largest IQR
                  ## We will want "unique genes" in the non-specific filtered gene
                  ## set.
                  if (deparse(substitute(var.func)) == "IQR") {
                      esetIqr <- rowIQRs(exprs(eset))
                  } else {
                      esetIqr <- apply(exprs(eset), 1, var.func)
                  }
                  numNsWithDups <- nfeat(eset)
                  uniqGenes <- findLargest(featureNames(eset), esetIqr,
                                           annotation(eset))
                  eset <- eset[uniqGenes, ]
                  logvar <- "numDupsRemoved"
                  assign(logvar, numNsWithDups - nfeat(eset), envir=filter.log)
              }


              if (var.filter) {
                  if (deparse(substitute(var.func)) == "IQR") {
                      esetIqr <- rowIQRs(exprs(eset))
                  } else {
                      esetIqr <- apply(exprs(eset), 1, var.func)
                  }
                  ##note this was not happening in the first
                  ##version - despite the documentation
                  if (filterByQuantile) {
                      if ( 0 < var.cutoff && var.cutoff < 1 ) {
                          var.cutoff = quantile(esetIqr, var.cutoff)
                      } else stop("Cutoff Quantile has to be between 0 and 1.")
                  }
                  selected <- esetIqr > var.cutoff
                  eset <- eset[selected, ]
                  logvar <- "numLowVar"
                  assign(logvar, sum(!selected), filter.log)
              } 

              numSelected <- length(featureNames(eset))
              list(eset=eset, filter.log=as.list(filter.log))
          })
