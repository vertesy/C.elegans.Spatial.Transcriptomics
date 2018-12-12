######################################################################
# ReproduciblePeakGenes
######################################################################
# source ('/Users/abelvertesy/Github_repos/CeleTomo/01e.Reproducible.Peak.Genes.R')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try (source ('~/Github_repos/TheCorvinas/R/CodeAndRoll.R'),silent= F)
# try (source("~/Github_repos/Spermatogenesis/Mapping/functions_sp.r") , silent= F)
# source ('~/Github_repos/TheCorvinas/R/DatabaseLinke.r')

wlegend <-function(fill_ = NA, poz=4, legend, bty = "n", ..., w_=7, h_=w_, OverwritePrevPDF = trueUnlessSpec("b.save.wplots")) { # Add a legend, and save the plot immediately
  fNames = names(fill_)
  LF = length(fill_)
  LN = length(fNames)
  stopif( ( LN != LF & missing(legend) ), message = "The color vector (fill_) has less names than entries / the variable 'legend' is not provided.")
  # stopif( ( LF  != length(legend)), message = "Fill and legend are not equally long.")
  legend = if( LN == LF & missing(legend) ) fNames else legend
  pozz = translate(poz, oldvalues = 1:4, newvalues = c("topleft", "topright", "bottomright", "bottomleft"))
  legend(x=pozz, legend=legend, fill=fill_, ..., bty=bty)
  if (OverwritePrevPDF) {   wplot_save_this(plotname = plotnameLastPlot, w= w_, h = h_)  }
}


md.import <- function(from.file, to.file = path_of_report) {
  linez = readLines(from.file)
  if (!exists("path_of_report")) { print("Log path and filename is not defined in path_of_report") } else iprint(length(linez), "lines from",basename(from.file) ,"are concatenated to:", basename(path_of_report))
  for(LogEntry in linez) {
    write(LogEntry, path_of_report, append = T)
  }
}

# Setup ------------------------
setup_MarkdownReports(OutDir = p0(ProjectDir, "/Reproducible.Peak.Genes"), scriptname = "01e.Reproducible.Peak.Genes.R", b.mdlink = T, b.usepng = T, b.png4Github = F)

md.import("/Users/abelvertesy/Google_Drive/Avano/CT/Notes/Reproducible.Peak.Genes.md")


# Metadata ------------------------
mc=NULL
WormCol_Sex_HQ = WormCol_Sex[worms_HQ]

# Parameters ------------------------
percent =T
# Read In ------------------------
# lapply(zData)


# QC ------------------------
NrPlateauGenes = unlapply(PlateauGenes, l)
wbarplot(NrPlateauGenes, col=WormCol_Sex_HQ)

GeneCounts.VS.NrPlateauGenes = cbind("TotalGeneCount"=TotalGeneCount[worms_HQ], NrPlateauGenes) # X=cbind(TotalTrCount, NrPlateauGenes)
wplot(GeneCounts.VS.NrPlateauGenes, col=WormCol_Sex[worms_HQ], cex=2)

TranscripCounts.VS.NrPlateauGenes = cbind("TotalTrCount"=TotalTrCount[worms_HQ], NrPlateauGenes) # X=cbind(TotalTrCount, NrPlateauGenes)
wplot(TranscripCounts.VS.NrPlateauGenes, col=WormCol_Sex[worms_HQ], cex=2, savefile = F)
wlegend(fill_ = SexColors, 3, legend = SexLabels) #

# ------------------------
# Plateaugenes------------------------
llprint("## Reproducibility of Plateaugenes - In how many worms do we see the peaking genes?")

### Venn.Diagrams ------------------------
Venn.Diagrams = T
if (Venn.Diagrams) {
  OutDir
  create_set_SubDir("Venn.Diagrams.4")
  for (i in 1:l(HQ_per_Sex) ) {
    SX = HQ_per_Sex[[i]]
    # for (percent in T:F )  wvenn(PlateauGenes[SX],  plotname = kollapse("Reproducible.Plateaugenes.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent)
    for (percent in T:F )  wvenn(ls_gene.names[SX],  plotname = kollapse("Reproducible.HE.genes.in.sfData.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent, fontfamily =rep("Arial", 1),  main.fontfamily = "Arial", sub.fontfamily ="Arial"  )
    # for (percent in T:F )  wvenn(ls_gene.names[SX],  plotname = kollapse("Reproducible.HE.genes.in.sfData.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent,  main.fontfamily = "Arial", sub.fontfamily ="Arial", cex =0)
    # for (percent in F:T )    GR <-venn.diagram(x = PlateauGenes[SX], filename = NULL,  main = kollapse("Reproducible.Plateaugenes.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1]); grid.draw(GR)
  }
  create_set_Original_OutDir()
} #if

### Summary.Barplots ------------------------
Summary.Barplots = T
if (Summary.Barplots) { try.dev.off()
  llprint("## Summary Barplots - How many times are plateau genes seen?")
  GeneOccurence = list.fromNames(SexLabels)
  pdfA4plot_on("Multiple.Occurence.Plateaugenes", rows = 4)
  for (i in 1:l(HQ_per_Sex) ) {
    SX = HQ_per_Sex[[i]]
    GeneOccurence[[i]] = table(unlist(PlateauGenes[SX]))
    Occurences = sortbyitsnames(table(GeneOccurence[[i]]), decreasing = T)
    wbarplot(Occurences, col=SexColors[i], main =  p0("Occurences in ", SexLabels[i]))
    PC = percentage_formatter(Occurences/sum(Occurences))
    barplot_label(Occurences, labels = PC, bottom = T)
  } #for

  for (i in 1:l(HQ_per_Sex) ) {
    SX = HQ_per_Sex[[i]]
    Occurences = sortbyitsnames(table(table(unlist(PlateauGenes[SX]))), decreasing = T)
    OccurenceCumSum = cumsum(Occurences/sum(Occurences))
    wbarplot(OccurenceCumSum, col=SexColors[i], main =  p0("Cumulative Sum Coverage", SexLabels[i]))
    PC = percentage_formatter(OccurenceCumSum)
    barplot_label(OccurenceCumSum, labels = PC, TopOffset = .1)
  } #for
  pdfA4plot_off()
} #if
OutDir
### Reproducibility by Relative Position ------------------------

Relative.Position = T
if (Relative.Position) { try.dev.off()
  llprint("## Reproducibility of Relative Position of plateau genes seen multiple times.")


  RelativePositionOfPlateauGenes = GenesSeenMultip = list.fromNames(worms_HQ)
  for (i in 1:NrHQworms ) {
    w= worms_HQ[i]; iprint(w)
    PGw = PlateauGenes[[w]]
    ZD = zData[[w]]
    MaxPozz = apply(ZD[ PGw, ], 1, which.max )
    MaxSection = colnames(ZD)[MaxPozz]
    RelPozW = RelativePosition.nData[[w]][MaxSection]
    names(RelPozW) = PGw
    RelativePositionOfPlateauGenes[[w]] = RelPozW
  } #for

  try.dev.off()
  pdfA4plot_on("PlateauGenes_Slice_w_Max")
  for (i in 1:NrHQworms ) {
    w= worms_HQ[i]; iprint(w)
    RelP.W = RelativePositionOfPlateauGenes[[w]]
    PlateauGenes_Slice_w_Max = table(RelP.W)
    names(PlateauGenes_Slice_w_Max) = percentage_formatter( as.numeric(names(PlateauGenes_Slice_w_Max)) )
    # wbarplot(PlateauGenes_Slice_w_Max, col=WormCol_Sex_HQ[i],main = p0(w) )

    SEX = sex_HQ[i]
    GOC = GeneOccurence[[ SEX ]]
    NRWperSex = table(sex_HQ)[SEX]
    GenesSeenMultip[[w]] = which_names(GOC[GOC > (NRWperSex-1) ])
    # GenesSeenMultip[[w]] = which_names(GOC[GOC > 1 ])

    CATEG = (RelativePosition.nData[[WRM]])
    CATEG = sort(unique(RelP.W))
    PG.Slice_w_Max = rbind(
      'FreqNonRepr' = table_fixed_categories(na.omit.strip(RelP.W[ setdiff(names(RelP.W),GenesSeenMultip[[w]]) ]) , categories_vec =CATEG),
      'FreqRepr' = table_fixed_categories(na.omit.strip(RelP.W[    GenesSeenMultip[[w]]]), categories_vec =CATEG)
    )
    # View(PG.Slice_w_Max)
    CCC = c("Not in all worms" = "darkgrey", "Reproducible" = as.character(WormCol_Sex_HQ[i]))

    wbarplot(PG.Slice_w_Max, col = CCC, main=w)
    wlegend(fill_ = CCC, poz = 1, cex=.75 , title="Genes")
  } #for
  pdfA4plot_off()


  llprint("### Reproducibility Peaking Position within Plateaugenes")

  ReprPltPos = list.fromNames(worms_HQ)
  for (i in 1:NrHQworms ) {
    w= worms_HQ[i]; iprint(w)
    ReprPltPos[[w]] = RelativePositionOfPlateauGenes[[w]][ GenesSeenMultip[[w]] ]
  }
  RelPozMatrixPerSex = list.fromNames(SexLabels)
  for (i in 1:l(HQ_per_Sex) ) RelPozMatrixPerSex[[i]] = list2df(ReprPltPos[ HQ_per_Sex[[i]] ])

  pdfA4plot_on("Most Plateau Genes are peaking in the same region", rows = 4)
  for (i in 1:l(HQ_per_Sex) ) {   SX = HQ_per_Sex[[i]]
    Variation.in.Maximum.Expr.SD = 100*apply(RelPozMatrixPerSex[[i]],1,sd)
    whist(Variation.in.Maximum.Expr.SD, col = SexColors[i], main= SexLabels[i], xlb = "Standard deviation of peaking position (% length)")
  }
  for (i in 1:l(HQ_per_Sex) ) {   SX = HQ_per_Sex[[i]]
    Variation.in.Maximum.Expr.SEM = 100*rowSEM(RelPozMatrixPerSex[[i]])
    whist(Variation.in.Maximum.Expr.SEM, col = SexColors[i], main= SexLabels[i], xlb = "Standard Error of the Mean of peaking position (% length)")
  }
  for (i in 1:l(HQ_per_Sex) ) { SX = HQ_per_Sex[[i]]
    Variation.in.Maximum.Expr.CV = sort(rowCV(RelPozMatrixPerSex[[i]]))
    whist(Variation.in.Maximum.Expr.CV, col = SexColors[i], main= SexLabels[i], xlb = "CV of peaking position (% length)")
  }
  pdfA4plot_off()
} #if

# High-Z-Score Peak Genes ------------------------
setup_MarkdownReports(OutDir = p0(OutDirOrig,"Z-score-high"),scriptname = "01e.Reproducible.Peak.Genes.R",title = "High-Z-Score Peak Genes" )


Reproducible.Peak.Genes = T
if (Reproducible.Peak.Genes) {
  zPass = zPassNames = list.fromNames(worms_HQ)
  llprint("## Filtering genes by min. Z-score")
  zData.rowMax = lapply(zData[worms_HQ], rowMax)

  pdfA4plot_on("Z.score.selection")
  for (i in 1:NrHQworms ) {
    w= worms_HQ[i];w
    ZM = zData.rowMax[[w]]
    whist(ZM, breaks = 50, vline = mc$minZ, col=WormCol_Sex_HQ[w], main=p0("Max Z-scores in ", w))
    zPass[[w]] = filter_HP(ZM, mc$minZ)
    zPassNames[[w]] = which_names(zPass[[w]] )
  } #for

  Genes.Above.Zscore.X = unlapply(zPassNames, l)
  wbarplot(Genes.Above.Zscore.X, col=WormCol_Sex_HQ)
  X = cbind(Genes.Above.Zscore.X, NrPlateauGenes)
  wplot(X)
  NrGenes = unlapply(zPass, l)
  NrGenes.vs.NrPlateauGenes = cbind(NrGenes, NrPlateauGenes)
  wplot(NrGenes.vs.NrPlateauGenes, col = WormCol_Sex_HQ, cex=2)

  NrGenes.vs.Z.genes = cbind(NrGenes, Genes.Above.Zscore.X)
  wplot(NrGenes.vs.Z.genes, col = WormCol_Sex_HQ, cex=2)
  pdfA4plot_off()

  for (i in 1:l(HQ_per_Sex) ) {
    SX = HQ_per_Sex[[i]]
    percent =T
    wvenn(zPassNames[SX],  plotname = kollapse("Reproducible.Peak.Genes.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent)
    # for (percent in T:F )  wvenn(zPassNames[SX],  plotname = kollapse("Reproducible.Peak.Genes.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent)
  } #for
} #if


# ------------------------



# ------------------------
Z_score.VS.Expression = T

mc$minTotExpr = 100
mc$minZ=4

PASS_Z = PASS_Z_names = list.fromNames(worms_HQ)
if (Z_score.VS.Expression) {
  nData.rowSum =  lapply(nData, rowSums)

  llprint("## Genes with high Z-score are typically lowly expressed")
  pdfA4plot_on("Z_score.VS.Expression")
  for (i in 1:NrHQworms ) {
    w= worms_HQ[i]; iprint(w)
    Expr.VS.Z = cbind(
      "Total Expression log10" = log10(nData.rowSum[[i]]),
      "Maximum Z-score" = zData.rowMax[[i]]      )

    THR = quantile(as.numeric(nData.rowSum[[i]]),probs = seq(0, 1, 0.1))[2]
    PASS_Z[[w]] = (zData.rowMax[[i]] > mc$minZ & nData.rowSum[[i]] > THR)
    PASS_Z_names[[w]] =  which_names(PASS_Z[[w]])
    wplot(Expr.VS.Z, cex=.5, plotname = w, col=PASS_Z[[w]]+2)
  } #for
  pdfA4plot_off()

  lapply(X, idim)
} #if

for (i in 1:l(HQ_per_Sex) ) {
  SX = HQ_per_Sex[[i]]
  percent =T
  wvenn(PASS_Z_names[SX],  plotname = kollapse("Reproducible.Peak.2x.Filtered.Genes.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent)
  # for (percent in T:F )  wvenn(PASS_Z_names[SX],  plotname = kollapse("Reproducible.Peak.2x.Filtered.Genes.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent)
} #for

# ------------------------
setwd(OutDirOrig)


