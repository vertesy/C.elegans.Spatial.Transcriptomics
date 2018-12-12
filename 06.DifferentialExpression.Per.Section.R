######################################################################
# 06.DifferentialExpression.Per.Section.R
######################################################################
# source ('~/Github_repos/CeleTomo/06.DifferentialExpression.Per.Section.R')
try(dev.off(), silent = T)

# Functions ------------------------
# library(calibrate); # install.packages("calibrate")
require(DESeq2)
source ('~/Github_repos/CeleTomo/zz.Functions.CeleTomo.R')


# Setup ------------------------
OutDir = "~/Google_Drive/Avano/CT/Analysis/06.DifferentialExpression.Per.Section"
setup_MarkdownReports(OutDir = OutDir, scriptname = "06.DifferentialExpression.Per.Section.R", title = "Differential Expression Between Sections", append = F)

# Metadata ------------------------


# Parameters ------------------------

# DEseq
thr_padj =    0.05
thr_log2fc =  1


# Prepare ------------------------

i=1
for (i in Nr_HQ_sexual_worms:1) {
  w = HQ_sexual_worms[i];print(w)

  ANN = SectionAnnotation[[w]]
  Sectionz = names(ANN); l(Sectionz)
  SecTypes = unique(ANN)
  NrSecTypes = l(SecTypes)
  coldata = matrix.fromNames(colname_vec = SecTypes, rowname_vec = Sectionz)
  for (j in 1:NrSecTypes ) {  coldata[ , SecTypes[j]] = as.character(as.numeric(SecTypes[j] == ANN))} #for

  # Do ------------------------

  Filt.Res.DEseq =  Res.DEseq.per.Section = DEseq.per.Section =  list.fromNames(SecTypes)
  if (sex_HQ[w] ==3)  {
    DEseq.per.Section[["Head"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Head ) )
    DEseq.per.Section[["Germ1"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Germ1 ) )
    DEseq.per.Section[["Spermatheca"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Spermatheca ) )
    DEseq.per.Section[["Germ2"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Germ2 ) )
    DEseq.per.Section[["Vulva"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Vulva ) )
    DEseq.per.Section[["Germ3"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Germ3 ) )

    Res.DEseq.per.Section[["Head"]]		 <-	as.data.frame(results( DEseq.per.Section[["Head"]] ,		 contrast=c("Head","1","0")) )
    Res.DEseq.per.Section[["Germ1"]]		 <-	as.data.frame(results( DEseq.per.Section[["Germ1"]] ,		 contrast=c("Germ1","1","0")) )
    Res.DEseq.per.Section[["Spermatheca"]] <-	as.data.frame(results( DEseq.per.Section[["Spermatheca"]] ,contrast=c("Spermatheca","1","0")) )
    Res.DEseq.per.Section[["Germ2"]]		 <-	as.data.frame(results( DEseq.per.Section[["Germ2"]] ,		 contrast=c("Germ2","1","0")) )
    Res.DEseq.per.Section[["Vulva"]]		 <-	as.data.frame(results( DEseq.per.Section[["Vulva"]] ,		 contrast=c("Vulva","1","0")) )
    Res.DEseq.per.Section[["Germ3"]]		 <-	as.data.frame(results( DEseq.per.Section[["Germ3"]] ,		 contrast=c("Germ3","1","0")) )

  } else if (sex_HQ[w] == 2)  { # if males
    DEseq.per.Section[["Head1"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Head1 ) )
    DEseq.per.Section[["Head2"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Head2 ) )
    DEseq.per.Section[["Germ1"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Germ1 ) )
    DEseq.per.Section[["Germ2"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Germ2 ) )
    DEseq.per.Section[["Tail"]] <- DESeq(DESeqDataSetFromMatrix("countData" = floor(sfData[[w]][,Sectionz]), "colData" = coldata, "design" = ~ Tail ) )

    Res.DEseq.per.Section[["Head1"]]		 <-	as.data.frame(results( DEseq.per.Section[["Head1"]] ,		 contrast=c("Head1","1","0")) )
    Res.DEseq.per.Section[["Head2"]] <-	as.data.frame(results( DEseq.per.Section[["Head2"]] ,contrast=c("Head2","1","0")) )
    Res.DEseq.per.Section[["Germ1"]]		 <-	as.data.frame(results( DEseq.per.Section[["Germ1"]] ,		 contrast=c("Germ1","1","0")) )
    Res.DEseq.per.Section[["Germ2"]]		 <-	as.data.frame(results( DEseq.per.Section[["Germ2"]] ,		 contrast=c("Germ2","1","0")) )
    Res.DEseq.per.Section[["Tail"]]		 <-	as.data.frame(results( DEseq.per.Section[["Tail"]] ,		 contrast=c("Tail","1","0")) )

  } #if


  # Filtering and writing out ---------------------------------------------
  llprint("### Filtering")
  for (k in 1:NrSecTypes ) {
    tsvName = kollapse("DiffExp.", w,".", SecTypes[k],".VS.rest.tsv")
    Filt.Res.DEseq[[w]]  = filter_DESeq(DESeq_results = Res.DEseq.per.Section[[k]], thr_log2fc_ =thr_log2fc, thr_padj_=thr_padj)

    # Add links -----------
    if (nrow(Filt.Res.DEseq[[w]])) {
      Filt.Res.DEseq[[w]]  = cbind(Filt.Res.DEseq[[w]], "link_wormbase" = link_wormbase(rownames(Filt.Res.DEseq[[w]]), writeOut = F, Open = F), "link_String" = link_String(rownames(Filt.Res.DEseq[[w]]), organism ="elegans", writeOut = F, Open = F) )
    } #if
    write.simple.tsv(Filt.Res.DEseq[[w]], ManualName = tsvName)

  } #for

  # Plotting ---------------------------------------------
  for (m in 1:NrSecTypes ) {
    pname = kollapse("plot.DiffExp.", w,".", SecTypes[m],".VS.rest")
    ylbb = kollapse("log2(",SecTypes[m],"/rest)")
    plotMA(prepare4plotMA(Res.DEseq.per.Section[[m]], thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=SecTypes[m], ylab=ylbb); wplot_save_this(plotname = pname, mdlink = T)

    pname = kollapse("volcano.DiffExp.", w,".", SecTypes[m],".VS.rest")
    wplot_Volcano(Res.DEseq.per.Section[[m]],   thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname)
  } #for

} # for each sexual worm



