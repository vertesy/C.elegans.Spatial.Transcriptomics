######################################################################################################
# source("~/Github_repos/Spermatogenesis/Functions.SP.r")





multi.stripchart.genes.per.cluster.SP <- function(genes = HE_genes, sc_ =sc, use='ndata', color=sc@fcol, cols_ = 2, rows_ = 3, ...,
                                                  pname_=paste0(substitute(genes), ".stripcharts.per.cluster.multiple.genes"), Plot2NewFile=T) {
  if(Plot2NewFile) pdfA4plot_on(pname = pname_, cols = cols_, rows = rows_, ...)
  if(Plot2NewFile) print(11111)
  for (i in 1:length(genes)) {
    g = genes[i]
    r = unlist( if (use=='ndata') { sc@ndata[g, ]-.1   } else {sc@expdata[g, ] })

    PerGene = split(x = r, f = sc@cluster$kpart)
    pname = id2name(g)
    wstripchart(PerGene, col =color, plotname = pname, jitter = .4, colorbyColumn = T, pch=20, pchcex = .5, ylab="Transcripts per cell", savefile = F, ...)
  }
  if(Plot2NewFile) pdfA4plot_off()
}

make96wp <- function(vector){
  plate = t(matrix(data = vector, nrow = 12))
  rownames(plate) = LETTERS[1:8]
  colnames(plate) = 1:12
  return(plate)
}

make384wp <- function(vector){
  plate = t(matrix(data = vector, nrow = 24))
  rownames(plate) = LETTERS[1:16]
  colnames(plate) = 1:24
  return(plate)
}

wplot_platelayout <- function(plate, restart_graphics = T, coll=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)){
  if (restart_graphics) { try(dev.off(), silent = T) }
  ttl =deparse(substitute(plate))
  pheatmap::pheatmap(plate, cluster_rows = F, cluster_cols = F, main = ttl, color = coll)
  wplot_save_this(plotname = ttl,h = 3, w=5)
  # if (restart_graphics) { try(dev.off(), silent = T) }
}
# wplot_platelayout(sp21_gateSortedPlate , coll = sc@fcol[1:5])

### color tnse maps ---------------------------------------------------------------------------------
# Original author: Mauro Muraro, modified

hlgene<-function(geneIDs, folder_suffix =""){
  g <- rownames(sc@ndata)[grep(geneIDs,rownames(sc@ndata))]
  dir = paste0(OutDir,"/tSNE_log2", folder_suffix); try(dir.create(dir), silent = T); setwd(dir);
  if ( length(g) == 0 ) next
  fname <- paste("tsne_map",n=geneIDs,"log",sep="_")
  plot.2.file(fname)
  plotexptsne(sc,g,n=geneIDs,logsc=TRUE)
  try(dev.off(), silent = T)
  dir = paste0(OutDir,"/tSNE",folder_suffix); try(dir.create(dir), silent = T); setwd(dir);
  fname <- paste("tsne_map",n=geneIDs,sep="_")
  plot.2.file(fname)
  plotexptsne(sc,g,n=geneIDs)
  try(dev.off(), silent = T)
  setwd(OutDir)
}


### diff genes functions ---------------------------------------------------------------------------------


# Define dummy function 4 normalization ----------------------------
plot.Distributions.4.norm <- function(NormExprData = NN, Genes2plot = SelectedGenes, Libraries = LibraryComposition, pname_ = "Normalized.Expression.Density.Distributions.per.lib.before.normalization") { # REFERS TO OBJECTS IN THE MEM SPACE
  try.dev.off()
  stopifnot( NCOL(NormExprData) == length(Libraries) ) # Number of cells in Libraries!= in NormExprData
  pdfA4plot_on(pname = pname_, cols = 4)
  i=1
  for (i in 1:length(Genes2plot) ) {
    g = Genes2plot[i]
    GenesExpr = as.named.vector(NormExprData[ g, , drop=F], WhichDimNames = 2)+1
    median.Total = median(GenesExpr)

    ExprPerLib = split(GenesExpr, Libraries) # split per library
    median.PerLib = unlapply(ExprPerLib, median)
    j=1
    for (j in 1:length(ExprPerLib) ) {
      plot(density(ExprPerLib[[j]]), main=names(ExprPerLib)[j], ylab=g)
      abline(v=median.Total, lty=2, col="red")
      lines(density(ExprPerLib[[j]]), main=names(ExprPerLib)[j], ylab=g)
      abline(v=median.PerLib[j], lty=2, col="grey33")
    } #for
    printEveryN(i=i,N = 20)
  } #for
  pdfA4plot_off()
}


PlotBatchEff.dummy <- function(df = df_rowMeans.log2, idx_rem = idx_removed, TTL = "Genes with extreme protocol bias are discarded") {
  CN = colnames(df)
  LLB  = p0( "log10(av. mRNA in " , CN, ")")
  CCC = c(3,1)
  names(CCC) = paste0(c("Normalized", "Removed"), ": ", percentage_formatter(table(idx_rem)/length(idx_removed)))
  table(idx_rem)
  plot(df, col=idx_rem, pch=18, cex=.33, xlab =LLB[1], ylab =LLB[2], main = TTL)
  abline(0,1)
  points(df[ !idx_rem, ], pch=18, cex=.33, col=3)
  wlegend(CCC, cex=.75, poz = 1)
}



calc.rel.expr <- function(expr.mat = mvexpdata, genes = GO_0007283_spermatogenesis, Factor=100,id.to.name=F) {
  iprint(length(genes), "genes are checked.")
  Subset.mat = getRows(expr.mat, genes)
  Subset.Total.expr = colSums(Subset.mat, na.rm = T)
  iprint("Range: ", (iround(range(Subset.Total.expr))))
  Total.expr = colSums(expr.mat, na.rm = T)
  return(Factor* (Subset.Total.expr/Total.expr))
}


# LinearizeNestedList Source: Akhil S Bhel via @mrdwab at https://gist.github.com/mrdwab/4205477
LinearizeNestedList <- function(NList, LinearizeDataFrames=FALSE,
                                NameSep="/", ForceNames=FALSE) {
  # LinearizeNestedList:
  #
  # https://sites.google.com/site/akhilsbehl/geekspace/
  #         articles/r/linearize_nested_lists_in_r
  #
  # Akhil S Bhel
  #
  # Implements a recursive algorithm to linearize nested lists upto any
  # arbitrary level of nesting (limited by R's allowance for recursion-depth).
  # By linearization, it is meant to bring all list branches emanating from
  # any nth-nested trunk upto the top-level trunk s.t. the return value is a
  # simple non-nested list having all branches emanating from this top-level
  # branch.
  #
  # Since dataframes are essentially lists a boolean option is provided to
  # switch on/off the linearization of dataframes. This has been found
  # desirable in the author's experience.
  #
  # Also, one'd typically want to preserve names in the lists in a way as to
  # clearly denote the association of any list element to it's nth-level
  # history. As such we provide a clean and simple method of preserving names
  # information of list elements. The names at any level of nesting are
  # appended to the names of all preceding trunks using the `NameSep` option
  # string as the seperator. The default `/` has been chosen to mimic the unix
  # tradition of filesystem hierarchies. The default behavior works with
  # existing names at any n-th level trunk, if found; otherwise, coerces simple
  # numeric names corresponding to the position of a list element on the
  # nth-trunk. Note, however, that this naming pattern does not ensure unique
  # names for all elements in the resulting list. If the nested lists had
  # non-unique names in a trunk the same would be reflected in the final list.
  # Also, note that the function does not at all handle cases where `some`
  # names are missing and some are not.
  #
  # Clearly, preserving the n-level hierarchy of branches in the element names
  # may lead to names that are too long. Often, only the depth of a list
  # element may only be important. To deal with this possibility a boolean
  # option called `ForceNames` has been provided. ForceNames shall drop all
  # original names in the lists and coerce simple numeric names which simply
  # indicate the position of an element at the nth-level trunk as well as all
  # preceding trunk numbers.
  #
  # Returns:
  # LinearList: Named list.
  #
  # Sanity checks:
  #
  stopifnot(is.character(NameSep), length(NameSep) == 1)
  stopifnot(is.logical(LinearizeDataFrames), length(LinearizeDataFrames) == 1)
  stopifnot(is.logical(ForceNames), length(ForceNames) == 1)
  if (! is.list(NList)) return(NList)
  #
  # If no names on the top-level list coerce names. Recursion shall handle
  # naming at all levels.
  #
  if (is.null(names(NList)) | ForceNames == TRUE)
    names(NList) <- as.character(1:length(NList))
  #
  # If simply a dataframe deal promptly.
  #
  if (is.data.frame(NList) & LinearizeDataFrames == FALSE)
    return(NList)
  if (is.data.frame(NList) & LinearizeDataFrames == TRUE)
    return(as.list(NList))
  #
  # Book-keeping code to employ a while loop.
  #
  A <- 1
  B <- length(NList)
  #
  # We use a while loop to deal with the fact that the length of the nested
  # list grows dynamically in the process of linearization.
  #
  while (A <= B) {
    Element <- NList[[A]]
    EName <- names(NList)[A]
    if (is.list(Element)) {
      #
      # Before and After to keep track of the status of the top-level trunk
      # below and above the current element.
      #
      if (A == 1) {
        Before <- NULL
      } else {
        Before <- NList[1:(A - 1)]
      }
      if (A == B) {
        After <- NULL
      } else {
        After <- NList[(A + 1):B]
      }
      #
      # Treat dataframes specially.
      #
      if (is.data.frame(Element)) {
        if (LinearizeDataFrames == TRUE) {
          #
          # `Jump` takes care of how much the list shall grow in this step.
          #
          Jump <- length(Element)
          NList[[A]] <- NULL
          #
          # Generate or coerce names as need be.
          #
          if (is.null(names(Element)) | ForceNames == TRUE)
            names(Element) <- as.character(1:length(Element))
          #
          # Just throw back as list since dataframes have no nesting.
          #
          Element <- as.list(Element)
          #
          # Update names
          #
          names(Element) <- paste(EName, names(Element), sep=NameSep)
          #
          # Plug the branch back into the top-level trunk.
          #
          NList <- c(Before, Element, After)
        }
        Jump <- 1
      } else {
        NList[[A]] <- NULL
        #
        # Go recursive! :)
        #
        if (is.null(names(Element)) | ForceNames == TRUE)
          names(Element) <- as.character(1:length(Element))
        Element <- LinearizeNestedList(Element, LinearizeDataFrames,
                                       NameSep, ForceNames)
        names(Element) <- paste(EName, names(Element), sep=NameSep)
        Jump <- length(Element)
        NList <- c(Before, Element, After)
      }
    } else {
      Jump <- 1
    }
    #
    # Update book-keeping variables.
    #
    A <- A + Jump
    B <- length(NList)
  }
  return(NList)
}



plotdiffgenesnb.nopoints <- function(x, pthr=.05, padj = TRUE, lthr = 1, mthr=-Inf, geothr=NULL, lthr_line = TRUE,
                                     show_names = TRUE, genenameCex=.75, label.nr = 50, label.entriched.only=F,
                                     Aname = NULL, Bname = NULL, short.axisname=T, axnPozX =0, axnPozY =c(.25, -.25), axCex=1, ...){
  y <- as.data.frame(x$res)
  if ( is.null(Aname) ) Aname <- "baseMeanA"
  if ( is.null(Bname) ) Bname <- "baseMeanB"

  YLB = if (short.axisname) "log2(Fold change)" else { p0("log2( (mRNA[", Aname, "] + mRNA[", Bname, "])/2 )") }
  XLB = if (short.axisname) "log2(Average transcript count)" else { p0("log2( mRNA[", Bname, "]) - log2(mRNA[", Aname, "])") }
  XLM = range(log2(y$baseMean))
  YLM = range(y$log2FoldChange)

  plotData = (y[ ,c("baseMean",  "log2FoldChange")])
  plotData$"baseMean" = log2(plotData$"baseMean")

  if (is.null(pthr)) {
    plot(plotData, pch = 20, cex=0,..., col="grey", xlab = XLB, ylab = YLB)
  } else if ( !is.null(pthr) ){
    p.used <- if ( padj ) y$padj else y$pval
    names(p.used) = rownames(y)
    f <- (p.used < pthr )
    par(bg=NA)
    plot(plotData[!f, ], pch = 20, cex=0,..., col="grey",xlim=XLM, ylim=YLM, xlab = XLB, ylab = YLB) # plot only the insignificant points > smaller file
    # points(plotData[f, ], col="red", pch = 20)
  }
  abline(0, 0)
  if(lthr_line){abline(h = c(-lthr,lthr), lty = 3, col="grey")}

  if (short.axisname) {text(axnPozX, axnPozY[2], labels = p0("Enriched in ", Aname), cex=axCex);  text(axnPozX, axnPozY[1], labels =  p0("Enriched in ", Bname), cex=axCex)}
  if ( !is.null(lthr) )   f <- f & abs( y$'log2FoldChange' ) > lthr
  if ( !is.null(geothr) ) f <- f & abs( log2(y$'foldChange_GeoMean' )) > log2(geothr)
  if ( !is.null(mthr) )   f <- f & log2(y$'baseMean') > mthr
  if ( show_names ){
    sorter = as.named.vector(log2(y[f,"baseMean", drop=F]) * y[f,"log2FoldChange"])
    NZ = if (label.nr) {
      if (label.entriched.only) names(tail(sort(sorter), label.nr)) else names(trail(sort(sorter), label.nr))
    } else rownames(y)[f]
    if ( sum(f) > 0 ) text(log2(y[NZ, "baseMean"]), y[NZ, "log2FoldChange"], labels = id2name(NZ), cex=genenameCex)
  }
  legend("bottomright", legend = p0("p < ", pthr ), bty="n", text.col = 2)
}


plotdiffgenesnb.onlypoints <- function(x, pthr=.05, padj = TRUE, lthr = 1, mthr=-Inf, geothr=NULL, lthr_line = TRUE,
                                       show_names = TRUE, genenameCex=.75, label.nr = 50, label.entriched.only=F,
                                       Aname = NULL, Bname = NULL, short.axisname=T, axnPozX =0, axnPozY =c(.25, -.25), axCex=1, ...){
  y <- as.data.frame(x$res)
  if ( is.null(Aname) ) Aname <- "baseMeanA"
  if ( is.null(Bname) ) Bname <- "baseMeanB"

  YLB = if (short.axisname) "log2(Fold change)" else { p0("log2( (mRNA[", Aname, "] + mRNA[", Bname, "])/2 )") }
  XLB = if (short.axisname) "log2(Average transcript count)" else { p0("log2( mRNA[", Bname, "]) - log2(mRNA[", Aname, "])") }
  XLM = range(log2(y$baseMean))
  YLM = range(y$log2FoldChange)

  plotData = (y[ ,c("baseMean",  "log2FoldChange")])
  plotData$"baseMean" = log2(plotData$"baseMean")

  if (is.null(pthr)) {
    plot(plotData, pch = 20, ..., col="grey", xlab = XLB, ylab = YLB)
  } else if ( !is.null(pthr) ){
    p.used <- if ( padj ) y$padj else y$pval
    names(p.used) = rownames(y)
    f <- (p.used < pthr )
    plot(plotData[!f, ], pch = 20, ..., col="grey",xlim=XLM, ylim=YLM, xlab = "", ylab = "", main = "", axes = F) #
    points(plotData[f, ], col="red", pch = 20)
  }
  # abline(0, 0)
  # if(lthr_line){abline(h = c(-lthr,lthr), lty = 3, col="grey")}

  # if (short.axisname) {text(axnPozX, axnPozY[2], labels = p0("Enriched in ", Aname), cex=axCex);  text(axnPozX, axnPozY[1], labels =  p0("Enriched in ", Bname), cex=axCex)}
  if ( !is.null(lthr) )   f <- f & abs( y$'log2FoldChange' ) > lthr
  if ( !is.null(geothr) ) f <- f & abs( log2(y$'foldChange_GeoMean' )) > log2(geothr)
  if ( !is.null(mthr) )   f <- f & log2(y$'baseMean') > mthr
  # if ( show_names ){
  #   sorter = as.named.vector(log2(y[f,"baseMean", drop=F]) * y[f,"log2FoldChange"])
  #   NZ = if (label.nr) {
  #     if (label.entriched.only) names(tail(sort(sorter), label.nr)) else names(trail(sort(sorter), label.nr))
  #   } else rownames(y)[f]
  #   if ( sum(f) > 0 ) text(log2(y[NZ, "baseMean"]), y[NZ, "log2FoldChange"], labels = id2name(NZ), cex=genenameCex)
  # }
  # legend("bottomright", legend = p0("p < ", pthr ), bty="n", text.col = 2)
}





multi.timecourse.colors <- function(genes = HE_genes, sc_ =sc, log10_=F, use='ndata', time =m$c$pstime # plot gene expression in Pseudo time usin
                                    , avOver=50, vline=F, xlb = "Pseudotime", returnData=F, RColorBrewerSet_ = "Paired"
                                    , AllGenes=F, POZ=1, ZSCORE=F, CEXX=.75, Plot2NewFile =T, addTitle=T, addSubTitle=T
                                    , pname_=substitute(genes), rows_ = 3, cols_ = (rows_-1), ...) {
  if(Plot2NewFile) try(dev.off(), silent = T)

  TIME =names(sort(time))
  xData = t( if (use=='ndata') { sc_@ndata[genes, TIME]-.1   } else { sc_@expdata[genes, TIME] })
  stopif( min(dim(xData)) ==0 )

  pname_ = paste('time', pname_, use, nameiftrue(AllGenes), nameiftrue(ZSCORE), nameiftrue(log10_), flag.name_value(avOver), sep = '.')
  kpart = sc_@cluster$kpart
  if (ZSCORE) {
    xData= apply(xData, 2, scale)
    avOver = 2*avOver} #if
  YLB = if (log10_) { "log10(Normalized Transcripts)" } else if (ZSCORE) { "Gene Expression (z-score)"} else if (use != 'ndata') { "Raw transcripts" } else { "Transcripts (normalized)"}

  MovingAverages = apply(xData, 2, movingAve, oneSide = avOver/2)
  if (log10_)  {MovingAverages = log10(MovingAverages+0.1)}

  SUBB = if (addSubTitle) { p0("Rolling average window:", avOver) } else ""
  if (AllGenes) { # try.dev.off()
    XX = id2name(genes)

    XX = paste0(id2name(genes), " (",id2chr(genes),")")
    ccc = wcolorize(XX, RColorBrewerSet =RColorBrewerSet_, ReturnCategoriesToo = T)

    plot(rowMax(MovingAverages), type='n', lwd=3, ylim=range(MovingAverages, na.rm = T)
         , panel.first=grid()
         , main=paste0(substitute(genes), " Expression"), ylab=YLB, xlab=xlb, sub=SUBB)
    for (i in 1:length(genes)) {
      g = genes[i]; print(g)
      lines(MovingAverages[, g], lwd=3, col=ccc$vec[i])
    } #for
    assign("plotnameLastPlot", pname_, envir = .GlobalEnv)
    if(CEXX) wlegend(ccc$categ, poz = POZ, cex=CEXX ) # , OverwritePrevPDF= Plot2NewFile
  } else {
    if(Plot2NewFile) pdfA4plot_on(pname = pname_, cols = cols_, rows = rows_)
    for (i in 1:length(genes)) {
      g = genes[i]; print(g)
      XD = jitter(xData[ , i], amount = .1)
      pchx=(20+kpart[names(XD)])
      pname = if (addTitle) { id2titlecaseitalic.sp(g, suffix = "expression") } else {""}
      plot(XD, cex=.75, pch=pchx, bg=sc_@fcol[kpart[names(XD)]], lwd=.5
           , panel.first=grid()

           , main=pname, ylab=YLB, xlab="Pseudotime", sub=SUBB)

      if (vline) abline(v=vline, lty=3)
      wlegend.label(id2chr(g))
      text(labels = "MSCI", x = 725, y = 0.95*max(XD) )
      lines(MovingAverages[, g] , lwd=3 )

    } #for
    if(Plot2NewFile) pdfA4plot_off()
  }
  if(Plot2NewFile) assign("plotnameLastPlot", pname_, envir = .GlobalEnv)
  if (returnData) return(MovingAverages)
} #fun


