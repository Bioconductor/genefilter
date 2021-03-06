# This widget allows users to pick filters in the order they are going
# to be used to filer genes and set the parameters for
# each filter.
#
# Copyright 2003, J. Zhang. All rights reserved.
#

eSetFilter <- function(eSet){
    #
    # 2020-12-31 Wolfgang Huber: marking this as deprecated due to dubious coding style
    # (leading to many notes in R CMD check), lack of maintenance, apparently minuscule
    # update (e.g. I did not find a single mention of this function with
    # https://support.bioconductor.org/post/search/?query=esetfilter today) and
    # availability of better alternatives.
    #
    .Deprecated(new = "iSEE", old = "eSetFilter")
    
    require("tkWidgets", character.only = TRUE) ||
    stop("eSetFilter requires the tkWidgets ",
         "package. Please have it installed")

    descList <- getFuncDesc()

    buildGUI <- function(){
        END <<- FALSE

        selectedNames <- NULL
        filterWithArgs <- list()

        setFilter <- function(){
            currentFilter <- as.character(tkget(filters,
                                             (tkcurselection(filters))))
            args <- setESetArgs(currentFilter)
            if(!is.null(args)){
                expression <- paste(currentFilter, "(",
                    paste(names(args), args, sep = "=", collapse = ","),
                    ")", sep = "")
                filterWithArgs[[currentFilter]] <<- eval(parse(text =
                                                               expression))
                selectedNames <<- unique(c(selectedNames, currentFilter))
                writeList(pickedF, selectedNames)
                tkconfigure(selectBut, state = "disabled")
            }
        }
        cancel <- function(){
            tkdestroy(base)
        }
        finish <- function(){
            END <<- TRUE
            tkdestroy(base)
        }
        viewFilter <- function(){
            currentFilter <- as.character(tkget(filters,
                                             (tkcurselection(filters))))
            tkconfigure(description, state = "normal")
            writeText(description, descList[[currentFilter]])
            tkconfigure(description, state = "disabled")
            tkconfigure(selectBut, state = "normal")
        }
        pickedSel <- function(){
            tkconfigure(remBut, state = "normal")
        }
        remove <- function(){
             filter <- as.character(tkget(pickedF,
                                             (tkcurselection(pickedF))))
             selectedNames <<- setdiff(selectedNames, filter)
             writeList(pickedF, selectedNames)
             tkconfigure(remBut, state = "disabled")

        }

        base <- tktoplevel()
        tktitle(base) <- "BioC Filter Master"
        # Pack the top frame with a brief description
        introText <- tktext(base, width = 30, height = 4, wrap = "word")
        text <- paste("Bioconductor's gene filtering functons are",
                      "listed below. Select one from the list to view the",
                      "description and formal arguments for each filter.",
                      "A filter can be selected to the set of filters",
                      "for filtering genes using the select button.")
        writeText(introText, text)
        tkconfigure(introText, state = "disabled")
        tkpack(introText, expand = FALSE, fill  = "both", padx = 5)
        # Pack a frame with a list box for selected filters and
        # buttons manipulate the selected filters
        infoFrame <- tkframe(base)
        filterFrame <- tkframe(infoFrame)
        tkpack(tklabel(filterFrame, text = "Filters"), expand = FALSE,
               fill = "x")
        listFrame <- tkframe(filterFrame)
        filters <- makeViewer(listFrame, vHeight = 10, vWidth = 12,
                                vScroll = TRUE, hScroll = TRUE,
                                what = "list")
        tkbind(filters, "<B1-ButtonRelease>", viewFilter)
        tkbind(filters, "<Double-Button-1>", setFilter)
        writeList(filters, getFilterNames())
        tkpack(listFrame, expand = TRUE, fill = "both")
        selectBut <- tkbutton(filterFrame, text = "Select",
                              command = setFilter, state = "disabled")
        tkpack(selectBut, expand = FALSE, fill = "x")
        tkpack(filterFrame, side = "left", expand = FALSE, fill = "both")
        descFrame <- tkframe(infoFrame)
        tkpack(tklabel(descFrame, text = "Description"), expand = FALSE,
               fill = "x")
        dListFrame <- tkframe(descFrame)
        description <- makeViewer(dListFrame, vHeight = 10, vWidth = 30,
                                vScroll = TRUE, hScroll = TRUE,
                                what = "text")
        tkconfigure(description, wrap = "word", state = "disabled")
        tkpack(dListFrame, expand = TRUE, fill = "both")
        tkpack(descFrame, side = "left", expand = TRUE, fill = "both")
        selFrame <- tkframe(infoFrame)
        tkpack(tklabel(selFrame, text = "Selected"),
               expand = FALSE, fill = "x")
        selFFrame <- tkframe(selFrame)
        pickedF <- makeViewer(selFFrame, vHeight = 10, vWidth = 12,
                                vScroll = TRUE, hScroll = TRUE,
                                what = "list")
        tkbind(pickedF, "<B1-ButtonRelease>", pickedSel)
        tkbind(pickedF, "<Double-Button-1>", remove)
        tkpack(selFFrame, expand = TRUE, fill = "both")
        remBut <- tkbutton(selFrame, text = "Remove", command = remove,
                           state = "disabled")
        tkpack(remBut, expand = FALSE, fill = "x")
        tkpack(selFrame, expand = FALSE, fill = "both")
        tkpack(infoFrame, expand = TRUE, fill = "both", padx = 5)
        # Pack the bottom frame with cancel and finish buttons
        endFrame <- tkframe(base)
        cancelBut <- tkbutton(endFrame, width = 8, text = "Cancel",
                           command = cancel)
        tkpack(cancelBut, side = "left", expand = TRUE, fill = "x",
               padx = 10)
        finishBut <- tkbutton(endFrame, width = 8, text = "finish",
                           command = finish)
        tkpack(finishBut, side = "left", expand = TRUE, fill = "x",
               padx = 10)
        tkpack(endFrame, expand = FALSE, fill = "x", pady = 5)

        showESet(eSet)
        tkwait.window(base)

        if(END){
            tempList <- list()
            for(i in  selectedNames){
                tempList[[i]] <- filterWithArgs[[i]]
            }
            return(tempList)
        }else{
            return(NULL)
        }
    }

    filters <- buildGUI()
    if(!is.null(filters)){
        filters <- filterfun(unlist(filters))
        return(genefilter(exprs(eSet), filters))
    }else{
        return(NULL)
    }
}

getFilterNames <- function(){
    .Deprecated()
    return(sort(c("Anova", "coxfilter", "cv", "gapFilter", "kOverA",
    "maxA", "pOverA", "ttest")))
}

getFuncDesc <- function(lib = "genefilter", funcs = getFilterNames()){
    .Deprecated()
    descList <- list()

    lines <- getRdAsText(lib)
    for(i in funcs){
        rd <- lines[grep(paste("\\\\name\\{", i, "\\}", sep = ""), lines)]
        desc <- parseDesc(rd)
        args <- parseArgs(rd)
        if(length(args) > 0){
            temp <- "\n\nArguments:"
            for(j in names(args)){
                temp <- c(temp, paste(j, "-", args[[j]]))
            }
            args <- paste(temp, sep = "", collapse = "\n")
        }
        descList[[i]] <- paste(desc, args, sep = "", collapse = "")
    }
    return(descList)
}

getRdAsText <- function(lib){
    .Deprecated()
    fileName <- gzfile(file.path(.path.package(lib), "man",
                          paste(lib, ".Rd.gz", sep = "")), open = "r")
    lines <- readLines(fileName)

    lines <- paste(lines, sep = "", collapse = " ")
    lines <- unlist(strsplit(lines, "\\\\eof"))
    return(lines)
}

parseDesc <- function(text){
    .Deprecated()
    descRExp <- ".*\\\\description\\{(.*)\\}.*\\\\usage\\{.*"
    text <- gsub(descRExp, "\\1", text)
    text <- gsub("(\\\\[a-zA-Z]*\\{|\\})", "", text)
    return(text)
}

parseArgs <- function(text){
    .Deprecated()
    argsList <- list()
    text <- gsub(".*\\\\arguments\\{(.*)\\}.*\\\\details\\{.*", "\\1", text)
    text <- gsub(".*\\\\arguments\\{(.*)\\}.*\\\\value\\{.*", "\\1", text)
    text <- unlist(strsplit(text, "\\\\item\\{"))
    text <- gsub("(\\\\[a-zA-Z]*\\{|\\})", "", text)
    for(i in text){
        i <- unlist(strsplit(i, "\\{"))
        if(length(i) > 1){
            argsList[[i[1]]] <- i[2]
        }
    }
    return(argsList)
}

showESet <- function(eSet){
    .Deprecated()
    end <- function(){
        tkdestroy(base)
    }
    if(!is(eSet, "eSet")){
        stop()
    }
    colNRow <- dim(exprs(eSet))
    vl <- varLabels(eSet)
    text <- c(paste("Genes: ", colNRow[1]),
              paste("Samples: ", colNRow[2], sep = ""),
              "Variable labels:",
              paste(names(vl), ": ", vl[1:length(vl)], sep = ""))

    base <- tktoplevel()
    tktitle(base) <- "BioC ExpressionSet viewer"
    dataDescFrame <- tkframe(base)
    data <- makeViewer(dataDescFrame, vHeight = 10, vWidth = 25,
                       vScroll = TRUE, hScroll = TRUE,
                       what = "list")
    writeList(data, text)
    tkpack(dataDescFrame, expand = TRUE, fill = "both")
    endBut <- tkbutton(base, text = "Finish", command = end)
    tkpack(endBut, expand = FALSE, fill = "x", pady = 5)
}

setESetArgs <- function(filter){
    .Deprecated()
    on.exit(tkdestroy(base))

    cancel <- function(){
        tkdestroy(base)
    }
    end <- function(){
        END <<- TRUE
        tkdestroy(base)
    }
    END <- FALSE
    argsVar <- list()
    desc <- list()
    entries <- list()
    ftFun <- list()

    args <- getRdAsText("genefilter")
    args <- args[grep(paste("\\\\name\\{", filter, "\\}", sep = ""), args)]
    args <- parseArgs(args)
    argValues <- formals(filter)

    base <- tktoplevel()
    tktitle(base) <- "BioC Filter Argument input"

    tkgrid(tklabel(base, text = "Arguments"),
           tklabel(base, text = "Descriptions"),
           tklabel(base, text = "Values"))
    for(i in names(args)){
        argsVar[[i]] <- tclVar(as.character(argValues[[i]]))
        tempFrame <- tkframe(base)
        desc[[i]] <- makeViewer(tempFrame, vHeight = 3, vWidth = 55,
                                vScroll = FALSE, hScroll = FALSE,
                                what = "text")
        writeText(desc[[i]], args[[i]])
        tkconfigure(desc[[i]], wrap = "word", state = "disabled")
        entries[[i]] <- tkentry(base, textvariable = argsVar[[i]],
                                    width = 10)
        tkgrid(tklabel(base, text = i), tempFrame, entries[[i]])
        if(any(as.character(argValues[[i]]) == c("FALSE", "TRUE"))){
             ftFun[[i]] <- function(){}
             body <- list(as.name("{"),
                          substitute(eval(if(tclvalue(argsVar[[j]]) ==
                                                      "TRUE"){
                              writeList(entries[[j]], "FALSE")}else{
                                  writeList(entries[[j]], "TRUE")}),
                                     list(j = i)))
             body(ftFun[[i]]) <- as.call(body)
             tkbind(entries[[i]],"<B1-ButtonRelease>", ftFun[[i]])
        }
        tkgrid.configure(tempFrame, sticky = "eswn")
    }

    butFrame <- tkframe(base)
    canBut <- tkbutton(butFrame, text = "cancel", width = 8,
                       command = cancel)
    endBut <- tkbutton(butFrame, text = "Finish", width = 8,
                       comman = end)
    tkpack(canBut, side = "left", expand = FALSE, fill = "x")
    tkpack(endBut, side = "left", expand = FALSE, fill = "x")
    tkgrid(butFrame, columnspan = 3)

    tkwait.window(base)
    if(END){
        for(i in names(argValues)){
            argValues[[i]] <- tkWidgets:::formatArg(tclvalue(argsVar[[i]]))
        }
        return(argValues)
    }else{
        return(NULL)
    }
}

isESet <- function(eSet){
    .Deprecated()
    if(missing(eSet) ||
       (!is(eSet, "ExpressionSet"))) {
        tkmessageBox(title = "Input Error",
                     message = paste("filterMaster has to take",
                     "an object of class ExpressionSet"),
                     icon = "warning",
                     type = "ok")
        return(FALSE)
    }else{
        return(TRUE)
    }
}
