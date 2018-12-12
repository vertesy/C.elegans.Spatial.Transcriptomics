 ######################################################################
# zz.Functions.CeleTomo.R
######################################################################
# source ('~/Github_repos/CeleTomo/zz.Functions.CeleTomo.R')


# Functions ------------------------

clone2gene <- function(CloneNames, RenameDuplicated =T) {
  if (!exists("metadata_Genes")) print('metadata_Genes object not found. Load it: metadata_Genes = read.simple.tsv("~/Google_Drive/Avano/CT/metadata_Genes/Genes_to_Transcripts_c_elegans.tsv")' )
  names(CloneNames) =CloneNames

  gNames = metadata_Genes[,"GeneName"]
  names(gNames)= metadata_Genes[,"CloneName"]

  NotFound = setdiff(CloneNames, metadata_Genes$"CloneName")
  if(l(NotFound)) llprint("! ",l(NotFound), " of ",l(CloneNames), "(", percentage_formatter(l(NotFound) /l(CloneNames)), ") clone names are not found. Clone names are kept for these: ", NotFound)

  Found = intersect(CloneNames, metadata_Genes$"CloneName"); l(Found)
  CloneNames[Found] = gNames[Found]

  dupl =any.duplicated(CloneNames)
  if(dupl) {
    iprint(dupl, " returned gene names are duplicated, suffixed to unique by default. These are: ", CloneNames[which(duplicated(CloneNames))])
    if (RenameDuplicated) { CloneNames= make.unique(CloneNames) }
  }
  iprint(l(CloneNames), " gene names are returned")
  CloneNames # these are now replaced with gene names
}


# source: http://www.cookbook-r.com/Graphs/Histogram_and_density_plot/
plot.multi.dens <- function(s) {
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s)) {
    junk.x = c(junk.x, density(s[[i]])$x)
    junk.y = c(junk.y, density(s[[i]])$y)
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  plot(density(s[[1]]), xlim = xr, ylim = yr, main = "")
  for(i in 1:length(s)) {
    lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
  }
}



wbarplot_CT <-function (variable, ..., col = "gold1", sub = F, plotname = substitute(variable), main = substitute(variable), tilted_text = F, ylimits = NULL,
                     hline = F, vline = F, filtercol = 1, lty = 1, lwd = 2, lcol = 2,
                     errorbar = F, upper = 0, lower = upper, arrow_width = 0.1, arrow_lwd = 1,
                     w = 7, h = 7, incrBottMarginBy = 0, mdlink = F) {

  fname = kollapse(plotname, ".barplot")
  if (incrBottMarginBy) { .ParMarDefault <- par("mar"); 	par(mar=c(par("mar")[1]+incrBottMarginBy, par("mar")[2:4]) ) } 	# Tune the margin
  cexNsize = 0.8/abs(log10(length(variable)))
  cexNsize = min(cexNsize, 1)
  if (sub == T) {	subtitle = paste("mean:", iround(mean(variable, na.rm = T)), "CV:", percentage_formatter(cv(variable)))	} else if (sub == F) { subtitle = "" } else { subtitle = sub }
  if (hline & filtercol == 1) { col = (variable >= hline) + 2	}
  if (hline & filtercol == -1) { col = (variable < hline) + 2	}
  if (errorbar & is.null(ylimits)) {	ylimits = range(c(0, (variable + upper + abs(0.1 * variable)), variable - lower - abs(0.1 * variable)), na.rm = T) } # else {	ylimits = range(0, variable)	}
  if (tilted_text) {	xlb = NA	}	else {		xlb = names(variable)	}

  x = barplot(variable, ylim = ylimits, ..., names.arg = xlb, main = main, sub = subtitle, col = col, las = 2, cex.names = cexNsize)
  if (hline) { abline(h = hline, lty = lty, lwd = lwd, col = lcol)	}
  if (vline[1]) { abline(v = x[vline], lty = lty, lwd = lwd, col = lcol)	}
  if (errorbar) {	arrows(x, variable + upper, x, variable - lower, angle = 90, code = 3, length = arrow_width, lwd = arrow_lwd, ...)	}
  if (tilted_text) {
    text(x = x - 0.25, y = 0, labels = names(variable), xpd = TRUE, srt = 45, cex = cexNsize, adj = c(1,3))
  }

  # dev.copy2pdf(file = FnP_parser(fname, "pdf"), width = w, height = h, title = paste0(basename(fname), " by ", scriptname))
  if (incrBottMarginBy) { par("mar" = .ParMarDefault )}
  assign("plotnameLastPlot", fname, envir = .GlobalEnv)
  if (mdlink) { MarkDown_Img_Logger_PDF_and_PNG(fname_wo_ext = fname)	}
}


# Temporary until updating Markd package

continue_logging_markdown <-function (scriptname) {
  if (exists("OutDir")) {	path = OutDir } else {	path = getwd(); iprint("OutDir not defined !!! Saving in working directory.") }
  path_of_report <- kollapse(path, "/", scriptname, ".log.md", print = F)
  iprint("Writing report in:", path_of_report)
  assign("path_of_report", path_of_report, envir = .GlobalEnv)

  BackupDir = kollapse(OutDir, "/", substr(scriptname, 1, (nchar(scriptname) - 2)), format(Sys.time(), "%Y_%m_%d-%Hh"), print = F)
  if (!exists(BackupDir)) {
    dir.create(BackupDir)
    assign("BackupDir", BackupDir, envir = .GlobalEnv)
  }
}

peakFinder <- function(rleResults, MinConsecutiveLength) {    # find values that are 1 (True) in 'MinConsecutiveLength' consecutive sections (Run Length Encoding)
  any(rleResults$lengths[rleResults$values==1] > MinConsecutiveLength)  }


annot_col.create.old <- function(df, annot_vec) { # Auxiliary function for pheatmap. Prepares the 2 variables needed for "annotation_col" and "annotation_colors" in pheatmap
  stopifnot( l(annot_vec) == dim(df)[2] )
  print(substitute(annot_vec))
  df = as.data.frame(annot_vec)
  df[,1] = as.character(df[,1])
  assign(x = "annot", value = df, envir = .GlobalEnv)
  assign(x = "annot_col", value = list("annot_vec" = val2col(annot_vec[!duplicated(annot_vec)])), envir = .GlobalEnv)
  print("annot and annot_col variables are created. Use: pheatmap(..., annotation_col = annot, annotation_colors = annot_col)")
}

# similar.genes <- function(gene = "flp-11", Expression, number = 0, corThr =.9, dist_ = "pearson", plotResults=T){
#   stopifnot(gene %in% rownames(Expression))
#   stopifnot(dist_ %in% c("pearson", "spearman", "euclidean","manhattan")  )
#   target = as.named.vector(Expression[gene, ], WhichDimNames = 2)
#
#   if (dist_ == "pearson" | dist_ == "spearman") {
#     corWt = as.named.vector(  cor(x = t(Expression), y =target, method =  dist_,use = "na.or.complete" )      );
#     if(number) {simGenes = tail (sort(corWt), n = number)
#     } else if(corThr) {    simGenes = corWt[corWt >= corThr]}
#     mth = "correlation"
#   }
#   if (dist_ == "euclidean" | dist_ == "manhattan") {
#     distances = signif(apply(t(Expression), 2, function(x) dist(rbind(target, x ), method = dist_)))
#     if(number) {simGenes = head(sort(distances), n = number)
#     } else if(corThr) {    simGenes = distances[distances <= corThr]}
#     mth = "distance in gene expression"
#   }
#   print(paste("The ",number,"most similar genes based on", dist_, mth,"are: ", paste(names(simGenes), collapse = ", ")))
#
#   if (plotResults) {
#     try(dev.off, silent = T)
#     pheatmap::pheatmap(Expression[names(simGenes), ], cluster_cols = F)
#   }
#   return(simGenes)
# }



similar.genes <- function(gene = "flp-11", Expression, number = 0, corThr =.9, dist_ = "pearson", plotResults=T){
  stopifnot(gene %in% rownames(Expression))
  stopifnot(dist_ %in% c("pearson", "spearman", "euclidean","manhattan")  )
  target = as.named.vector(Expression[gene, ], WhichDimNames = 2)

  if (dist_ == "pearson" | dist_ == "spearman") {
    corWt = as.named.vector(  cor(x = t(Expression), y =target, method =  dist_,use = "na.or.complete" )      );
    if(number) {simGenes = tail (sort(corWt), n = number)
    } else if(corThr) {    simGenes = corWt[corWt >= corThr]; number=length(simGenes)}
    mth = "correlation"
  }
  if (dist_ == "euclidean" | dist_ == "manhattan") {
    distances = signif(apply(t(Expression), 2, function(x) dist(rbind(target, x ), method = dist_)))
    if(number) {simGenes = head(sort(distances), n = number)
    } else if(corThr) {    simGenes = distances[distances <= corThr]; number=length(simGenes)}
    mth = "distance in gene expression"
  }
  print(paste("The ",number,"most similar genes based on", dist_, mth,"are: ", paste(names(simGenes), collapse = ", ")))

  if (plotResults) {
    try(dev.off, silent = T)
    DAT = Expression[names(simGenes), ]
    DAT = log10(DAT+1)
    MinCor = if(!number) corThr else min(iround(simGenes))
    plotnameLastPlot = p0("The ", number," most similar genes to ", gene, " (",dist_,"; thr. ",MinCor,")")
    pheatmap::pheatmap(DAT, cluster_cols = F, main = plotnameLastPlot)
    assign("plotnameLastPlot", value = plotnameLastPlot, envir = .GlobalEnv)
  }
  return(simGenes)
}


# reproducible.peak.genes <- function(gene = "flp-11", Expression_ls, number = 10, dist_ = "pearson", plotResults=T){
#   # stopifnot(gene %in% rownames(Expression))
#   stopifnot(gene %in% sort(unique(unlapply(Expression_ls, rownames ))))
#   stopifnot(dist_ %in% c("pearson", "spearman", "euclidean","manhattan")  )
#
#   target = as.named.vector(Expression[gene, ], WhichDimNames = 2)
#
#   if (dist_ == "pearson" | dist_ == "spearman") {
#     corWt = as.named.vector(cor(x = t(Expression), y =target, method =  dist_) );
#     simGenes = tail (sort(corWt), n = number)
#     mth = "correlation"
#   }
#
#   if (dist_ == "euclidean" | dist_ == "manhattan") {
#     distances = signif(apply(t(Expression), 2, function(x) dist(rbind(target, x ), method = dist_)))
#     simGenes = head(sort(distances), n=number)
#     mth = "distance in gene expression"
#   }
#   print(paste("The ",number,"most similar genes based on", dist_, mth,"are: ", paste(names(simGenes), collapse = ", ")))
#
#   if (plotResults) {
#     try(dev.off, silent = T)
#     pheatmap::pheatmap(Expression[names(simGenes), ], cluster_cols = F)
#   }
#   return(simGenes)
# }

# Functions for worm alignment ---------------------

geneExists <- function(gene, df)   {gene %in% rownames(df)} # geneExists(gene = "cwp-1", df = zData[["e1_01"]])

geneExists.inAll <- function(gene, ls_df)   {for(i in 1:l(ls_df)){   print( geneExists(gene, ls_df[[i]]) ) } }

wstepbarplot <- function(numvec, type="s"){ # use s or S
  plot(numvec, type=type)
}

corWgene <- function(gene, df, dist_ = "pearson", thr=.95){
  stopifnot(geneExists(gene,df))
  target = as.named.vector(df[gene, ], WhichDimNames = 2)
  corWithgene = as.named.vector(cor(x = t(df), y =target, method =  dist_) );
  idx = which_names(corWithgene>=thr)
  return(sort(corWithgene[idx]))
}

corWgene.MultiWorm.plotCor <- function(gene, Expression_ls, dist_ = "pearson", thr=.9){
  pname= kollapse("Correlations.with.", gene,".",dist_)
  pdfA4plot_on(pname =  pname)
  for (i in 1:l(Expression_ls)) {
    w = names(Expression_ls)[i]; print(w)
    crz = corWgene(gene = gene,df = Expression_ls[[w]], thr = 0)
    pass = (crz >= thr)
    whist(variable = crz, plotname = w,vline = thr, breaks = 50, col = getsexcolor(w), savefile = F, xlim=c(0,1)
          , main = w, sub=paste0("Cor >= ",thr,": ", pc_TRUE(pass, NumberAndPC=T)), xlb = paste0(dist_, " correlation"))
  }
  pdfA4plot_off()
} # corWgene.MultiWorm.plotCor(gene = gene,Expression_ls = Expression_ls)



# Alignment ----------------------------------------------------------------------------------------------------

maxSlice <- function(gene = "flp-11", Expression = zData[[1]], printit=T, UseSliceNames=F ){ # dont UseSliceNames!
  if (!require("IRanges")) { print("Please install IRanges: install.packages('IRanges')") }
  mx= IRanges::which.max(Expression[gene,])
  poz =if(UseSliceNames) names(mx) else mx
  if (printit) iprint("Max slice is", poz,", value:",mx)
  return(poz)
}

GeneExprOffSet.Max <- function(gene = "cwp-1", Expression_ls = zData[c( 'e1_01', 'e1_02', 'e1_03', 'e1_04')], printit=T){
  Pass = which_names(unlapply(Expression_ls, function(x) geneExists(gene, x)))
  if(l(Pass)<l(Expression_ls)) iprint("Gene found only in",l(Pass),"worms:", Pass)
  unlapply(Expression_ls[Pass], function(x) maxSlice(gene, x, printit = T))
}

MatrixAlignmentCorrection <- function(OffSetVector = GeneExprOffSet.Max("fip-2", zData), Expression_ls = zData){
  for(i in 1:l(OffSetVector) ) {
    NrOfAddedCols = max(OffSetVector)-OffSetVector[i]
    if(NrOfAddedCols){
      PlaceHolder = rep(NaN, NROW(Expression_ls[[i]])*NrOfAddedCols)
      PlaceHolder_mat = matrix(data = PlaceHolder, ncol = NrOfAddedCols)
      colnames(PlaceHolder_mat) = paste0("offs", 1:NrOfAddedCols)
      Expression_ls[[i]] = cbind(PlaceHolder_mat, Expression_ls[[i]])
    }
  }
  return(Expression_ls)
}


# ----------------------------------------------------------------------------------------------------

getsexcolor <- function(wormnames)  {x =SexColors[sex[wormnames]]; names(x)=wormnames; return(x)}

  # gene="ccb-1"#GENE
geneExpression.Comparison <- function(gene = "flp-11", data= zData)  {
  print(gene)
  ls = lapply(data, function (x) getRows(mat = x, rownamez = gene, silent = T) )
  if(min(unlapply(ls, l))==0 | min(unlapply(ls, nrow))==0){
    ls= fixEmptyListElements(ls)
  }
  list2df(ls)
}

fixEmptyListElements <- function(ls)  {
  fix_these = which(unlapply(ls, NROW) ==0)
  for(i in 1:l(fix_these)){
    worm = fix_these[i]
    ls[[worm]] = rep(NaN, 10)
  }
  return(ls)
}

pheatmap.noclust <- function(...) pheatmap::pheatmap(..., cluster_cols = F, cluster_rows = F)

image.heatmap.worm <- function(data, xlab = "Slices", ylab="", main = "a",axes=F, ...) {
  x <- (1:nrow(data))
  y <- (1:ncol(data))
  image(y, x, t(data), col=HeatMapCol_BGR(24), axes=FALSE, main = main, xlab =xlab, ylab=ylab)
  axis(1, at = 1:ncol(data), labels=colnames(data),las=2,tick=FALSE)
  axis(side=2, at = 1:nrow(data), labels=rownames(data),las=1,tick=FALSE)
}


list2df <- function(ls, byRow=T)  { # Filling up from the left hand
  Lz = unlapply(ls, l)
  m = matrix(data = NA, nrow = l(ls), ncol = max(Lz)); dim(m)
  if(l(names(ls))) rownames(m)= names(ls)
  for(i in 1:l(ls)){
    m[i, 1:Lz[i] ]   = unlist(ls[[i]])
  }
  if(!byRow) {m = t(m)}
  return(m)
}



findGene.CeleTomo <- function(gene="cyp-43A1")  {
  allgenz = unique(unlist(lapply(sfData, rownames)))
  sort(grep(gene, allgenz, value = T))
}
#
# clone2gene <- function(CloneNames, KeepNotFoundOriginalNames =T, RenameDuplicated =T) {
#   if (!exists("metadata_Genes")) print('metadata_Genes object not found. Load it: metadata_Genes = read.simple.tsv("~/Google_Drive/Avano/CT/metadata/Genes_to_Transcripts_c_elegans.tsv")' )
#   NotFound = setdiff(CloneNames, metadata_Genes$CloneName)
#   if(l(NotFound)) iprint("! ",l(NotFound), " of ",l(CloneNames), "(", percentage_formatter(l(NotFound) /l(CloneNames)), ") clone names are not found: ", head(NotFound), "... Clone names are kept for these.")
#
#   ind = which(unique(metadata_Genes$CloneName) %in% CloneNames)
#   ind_of_input = which(CloneNames %in% metadata_Genes$CloneName)
#   if(KeepNotFoundOriginalNames) {
#     Gene_Names = CloneNames; names(Gene_Names) = CloneNames
#     Gene_Names[ind_of_input] = metadata_Genes$GeneName[ind]
#     names(Gene_Names)[ind_of_input] = metadata_Genes$CloneName[ind]
#   } else {
#     Gene_Names = metadata_Genes$GeneName[ind]
#     names(Gene_Names) = metadata_Genes$CloneName[ind]
#   }
#   iprint(sum(duplicated(Gene_Names)), " returned gene names are duplicated, suffixed to unique by default. These are: ", Gene_Names[which(duplicated(Gene_Names))])
#   iprint(l(Gene_Names), " gene names are returned")
#   if (RenameDuplicated) { Gene_Names = make.names(Gene_Names, unique = T); CloneNames = substr(CloneNames, 2, 1000)  }
#   Gene_Names
# }

# "Might not work as well"
gene2clone <- function(Gene_Names, KeepNotFoundOriginalNames =T, RenameDuplicated =T) {
  if (!exists("metadata_Genes")) print('metadata_Genes object not found. Load it: metadata_Genes = read.simple.tsv("~/Google_Drive/Avano/CT/metadata/Genes_to_Transcripts_c_elegans.tsv")' )
  NotFound = setdiff(Gene_Names, metadata_Genes$GeneName)
  if(l(NotFound)) iprint("! ",l(NotFound), " of ",l(Gene_Names), "(", percentage_formatter(l(NotFound) /l(Gene_Names)), ")  gene names are not found: ", head(NotFound), "... Gene names are kept for these.")

  ind = which(unique(metadata_Genes$GeneName) %in% Gene_Names)
  ind_of_input = which(Gene_Names %in% metadata_Genes$GeneName)
  if(KeepNotFoundOriginalNames) {
    CloneNames = Gene_Names; names (CloneNames) = Gene_Names
    CloneNames[ind_of_input] = metadata_Genes$CloneName[ind]
    names(CloneNames)[ind_of_input] = metadata_Genes$GeneName[ind]
  } else {
    CloneNames = metadata_Genes$CloneName[ind]
    names(CloneNames) = metadata_Genes$GeneName[ind]
  }
  iprint(sum(duplicated(CloneNames)), " returned clone names are duplicated, suffixed to unique by default. These are: ", Gene_Names[which(duplicated(Gene_Names))])
  iprint(l(CloneNames), " clone names are returned")
  if (RenameDuplicated) { CloneNames = make.names(CloneNames, unique = T); CloneNames = substr(CloneNames, 2, 1000)  }
  CloneNames
}



prefix.colnames <- function(df, prefix = ".1") {
  colnames(df) = p0(prefix, colnames(df))
  return(df)
}

suffix.colnames <- function(df, suffix = ".1") {
  colnames(df) = p0(colnames(df), suffix)
  return(df)
}

suffix.rownames <- function(df, suffix = ".1") {
  rownames(df) = p0(rownames(df), suffix)
  return(df)
}

md.tableWriter.DF.w.dimnames <- function (df, FullPath = path_of_report, percentify = F, title_of_table = NA, print2screen=F, WriteOut =F) {
  if (is.na(title_of_table)) {    t = paste0(substitute(df), collapse = " ")
  } else {                        t = title_of_table  }

  title_of_table = paste("\n#### ", t)
  if (variable.or.path.exists(FullPath)) {  write(title_of_table, FullPath, append = T)

    h = paste(colnames(df), collapse = " \t| ")
    h = paste("\n| |", h, " |", collapse = "")
    ncolz = dim(df)[2] + 1
    nrows = dim(df)[1]
    rn = rownames(df)
    sep = kollapse(rep("| ---", ncolz), " |", print = F)

    write(h, FullPath, append = T)
    write(sep, FullPath, append = T)
    for (r in 1:nrows) {
      if (is.numeric(unlist(df[r, ]))) {
        print(22)
        b = iround(df[r, ])
        if (percentify) {  b = percentage_formatter(b)  }
      } else {
        print(11)
        b = df[r, ] }
      b = paste(b, collapse = " \t| ")
      b = paste("|", rn[r], "\t|", b, " |", collapse = "")
      write(b, FullPath, append = T)
    }
  } else { print("NOT LOGGED: Log path and filename is not defined in FullPath")  }
  if (WriteOut) { write.simple.tsv(df, ManualName = p0(substitute(df),".tsv")) }
  # if (print2screen) { print(b) }
  "It was not working above"
}
