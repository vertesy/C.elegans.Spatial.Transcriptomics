######################################################################
# 01.QC_C.elegans.Spatial.Transcriptomics.R
######################################################################
# source ('~/Github_repos/C.elegans.Spatial.Transcriptomics/01.QC_C.elegans.Spatial.Transcriptomics.R')
rm(list=ls(all.names = TRUE));
try(dev.off(), silent = T)

# Functions ------------------------
require(MarkdownReportsDev)
require(gtools)

try (source ('~/Github_repos/TheCorvinas/R/CodeAndRoll.R'),silent= F)
try (source("~/Github_repos/C.elegans.Spatial.Transcriptomics/Functions.SP.r") , silent= F)
source ('~/Github_repos/TheCorvinas/R/DatabaseLinke.r')
source ('~/Github_repos/C.elegans.Spatial.Transcriptomics/zz.Functions.C.elegans.Spatial.Transcriptomics.R')
source("~/Github_repos/TheCorvinas/R/RNA_seq_specific_functions.r")
# Setup ------------------------

ProjectDir = "~/Google_Drive/Avano/CT/Analysis4"; try(dir.create(ProjectDir), silent = T)
OutDir = p0(ProjectDir,"/QC")
setup_MarkdownReports(OutDir = OutDir, scriptname = "01.QC_C.elegans.Spatial.Transcriptomics.R",
                      title = "Quality Control of C. elegans TomoSeq experiments.")

InputDir = "~/Google_Drive/Avano/CT/count_tables_new/elegans/Transcripts/"
Suffix = ".TranscriptCounts.tsv.gz"

Files = list.files(InputDir, pattern = Suffix)
Worms = stringr::str_split_fixed(Files, pattern = "\\.",n=2)[,1]
sex = as.numeric(substr(Worms,2,2))+1; names(sex) =Worms
SexColors = gplots::rich.colors(6, plot=T)[c(3,5,6)]; names(SexColors) =1:3
WormCol_Sex = SexColors[sex]; names(WormCol_Sex) =Worms

NrExp = l(Worms)

# Used for cleaned and normalized data files
WormSlicesFolder = "~/Google_Drive/Avano/CT/count_tables_new/elegans/WormSlices"; dir.create(WormSlicesFolder, showWarnings = F)
NormSeqDepthFolder = "~/Google_Drive/Avano/CT/count_tables_new/elegans/NormWormSlices"; dir.create(NormSeqDepthFolder, showWarnings = F)
Z_ScoreNormFolder =  "~/Google_Drive/Avano/CT/count_tables_new/elegans/ZScoreNormalized"; dir.create(Z_ScoreNormFolder, showWarnings = F)

# Metadata ------------------------
metadata_Genes = read.simple.tsv("~/Google_Drive/Avano/CT/Metadata/Genes/Genes_to_Transcripts_c_elegans.tsv.gz")
GeneSymbols = metadata_Genes$GeneName
names(GeneSymbols) = metadata_Genes$CloneName

# Parameters ------------------------
ReadIn = T
WriteOut_sf_and_zData= F

SectionMinTrCount = 1500
SectionMinGeneCount = 50
MinTotalTrCount = 1e6
MinTotalGeneCount = 1e4
MinExpr4WormLength = 20

# Keep these genes
MinExpr = 10
MinSlices = 1

SexLabels = c("Asex", "Male", "Herma"); names(SexLabels) =1:3
SexColors2 = wcolorize(3:1, show = T, RColorBrewerSet = "Set1")
SexColors2Labeled = SexColors2; names(SexColors2Labeled) =SexLabels

# Read In ------------------------

if (ReadIn) {
  if (!exists("ExpData")) {
    ExpData = list.fromNames(Worms)
    for (i in 1:NrExp) {
      infname = kollapse(InputDir, Worms[i], Suffix, print = F)
      transcripts = read.simple.tsv(infname)
      matx = na.replace(data.matrix(transcripts), replace = 0); idim(ExpData[[i]])
      idx.nonERCC = grepv('^ERCC-',rownames(matx),invert = T)
      ExpData[[i]] = matx[idx.nonERCC,]
      } # for
  } # !exists("ExpData")
}


# Calculate ------------------------
TotalTrCount = TotalGeneCount = ls.is.detected =NULL
# ExpData = ExpData_bac
ExpData_bac = ExpData

setup_MarkdownReports(OutDir = OutDir, scriptname = "01.QC_C.elegans.Spatial.Transcriptomics.R - NameConversion",
                      title = "Clone to Gene-name conversion in C. elegans TomoSeq experiments.")

CloneAndGene_names =list.fromNames(Worms)
for (i in 1:NrExp) {
  rnn = rownames(ExpData[[i]])
  rnx = GeneSymbols[rnn ]

  CloneAndGene_names[[i]] = rownames(ExpData[[i]])
  rownames(ExpData[[i]]) = names(CloneAndGene_names[[i]]) = clone2gene(rownames(ExpData[[i]]), RenameDuplicated =F) # rename to from clone-  to gene-name, where possible

  TotalTrCount[i] = sum(ExpData[[i]], na.rm = T)
    ls.is.detected[[i]] = rowSums(ExpData[[i]], na.rm = T) >= 1
  TotalGeneCount[i] = sum( ls.is.detected[[i]], na.rm = T)
}

names(ExpData) = names(TotalTrCount)= names(TotalGeneCount) = Worms
Dimensions = t(as.data.frame(lapply(ExpData, dim))); colnames(Dimensions) = c("Genes", "Sections"); Dimensions

continue_logging_markdown("01.QC_C.elegans.Spatial.Transcriptomics.R")


# QC ----------------------------------------------------------------------
llprint("## Overall Sample Quality")
llprint("We see the following sample statistics:")
try(dev.off(), silent = T)


Total_TrCount = (sort(TotalTrCount, decreasing = T))
wbarplot(Total_TrCount, hline= (MinTotalTrCount), tilted_text = T, mdlink = T)

# MinTotalGeneCount = 1e4
Total_GeneCount = (sort(TotalGeneCount, decreasing = T))
wbarplot(Total_GeneCount, tilted_text = T, hline= (MinTotalGeneCount), mdlink = T)

HQ_tr  = filter_HP(TotalTrCount, threshold =  MinTotalTrCount, prepend = "From all experiments, ")
HQ_gene =filter_HP(TotalGeneCount, threshold =  MinTotalGeneCount, prepend = "From all experiments, ")

HQ = (HQ_tr & HQ_gene)
isHQ = translate(HQ, oldvalues = c(T,F), newvalues = c("HQ", "LQ"))
sex_HQ = sex[HQ]
HQ_per_Sex = splititsnames_byValues(sex_HQ)

llprint("Due to low overall gene count, we exlude experiments:", which_names(HQ==F))

worms_HQ = which_names(HQ==T)
NrHQworms = l(worms_HQ)
llprint("Worms used for analysis:", worms_HQ)
  asex_HQ =which_names(sex_HQ==1)
  males_HQ =which_names(sex_HQ==2)
  hermaphrodites_HQ = which_names(sex_HQ==3)

worms_HQ.proper.name = c( 'glp1_1', 'glp1_3', 'Male_1', 'Male_2', 'Male_3', 'Male_4', 'Herm_1','Herm_2','Herm_3','Herm_4'); names(worms_HQ.proper.name) = worms_HQ




Tr_vs_Gene_Counts = cbind(TotalGeneCount, TotalTrCount)

ylmm= c(0,1.1*max(Tr_vs_Gene_Counts))
CCC= SexColors2[sex[rownames(Tr_vs_Gene_Counts)]]
wplot(Tr_vs_Gene_Counts, pch= 18, col =HQ+2, ylim = ylmm, cex=2)
points(Tr_vs_Gene_Counts, pch= 21, bg=CCC)

abline(v=MinTotalGeneCount, lty=3)
abline(h=MinTotalTrCount, lty=3)

text(Tr_vs_Gene_Counts, labels = Worms, srt=-45)
wplot_save_this(plotname = plotnameLastPlot, mdlink = T)

GeneCount.per.Sex = T
if (GeneCount.per.Sex) {
  Total_GeneCount.per.Sex = split(Total_GeneCount[worms_HQ], sex_HQ)
  names(Total_GeneCount.per.Sex) = names(SexColors2Labeled)
  wstripchart_list(Total_GeneCount.per.Sex, coll = SexColors2Labeled, tilted_text = T, pchcex = 2)


  ls_isexpr.gene = lapply(ls.is.detected, which_names)
  names(ls.is.detected) = names(ls_isexpr.gene) = Worms

  TotalGeneCount.in.each.sex = c(
    l(unique(unlist(ls_isexpr.gene[asex_HQ] ) ) ),
    l(unique(unlist(ls_isexpr.gene[males_HQ] ) ) ),
    l(unique(unlist(ls_isexpr.gene[hermaphrodites_HQ] ) ) )
  )

  GeneCounts.per.sex = cbind(
    "Median" = round(unlapply(Total_GeneCount.per.Sex, median)),
    "Mean" = round(unlapply(Total_GeneCount.per.Sex, mean)),
    "Worms" = table(sex_HQ),
    "Total" = TotalGeneCount.in.each.sex
  )
  md.tableWriter.DF.w.dimnames(GeneCounts.per.sex, print2screen = T)

} #if

try.dev.off()
Fig.1c = T
if (Fig.1c) {
  Fig.1c.Tr_vs_Gene_HQ = as.data.frame(Tr_vs_Gene_Counts[worms_HQ, 2:1])
  SexSymbols = 15:17; names(SexSymbols)=1:3
  SX = sex_HQ[worms_HQ]
  wplot(Fig.1c.Tr_vs_Gene_HQ, SexColors2[SX], pch = SexSymbols[SX], cex=1.5, panel_first=grid())
  wlegend(SexColors2Labeled,3)
} #if

# Section QC ------------------------------------------------------------------------------------------------
llprint("## Filtering Sections")
source("~/Github_repos/C.elegans.Spatial.Transcriptomics/01b.QC.SliceBarplots.R")


# QC Plots ------------------------------------------------

subb = paste("CV = ", percentage_formatter(cv(NrOfUsefulSections[worms_HQ])))
wbarplot(NrOfUsefulSections[worms_HQ], sub = subb, col = WormCol_Sex[worms_HQ], tilted_text = T, mdlink = T)

subb = paste("CV = ", percentage_formatter(cv(WormLength[worms_HQ])))
wbarplot(WormLength[worms_HQ], sub = subb, col = WormCol_Sex[worms_HQ], tilted_text = T)

WormLengthBySex = split(WormLength[worms_HQ], f=sex[worms_HQ])

MeanSex = unlist(lapply(WormLengthBySex, mean))
MedianSex = unlist(lapply(WormLengthBySex, median))
CV_Sex = unlist(lapply(WormLengthBySex, cv))

Sex_SummarStat = cbind(MeanSex, MedianSex, CV_Sex)
write.simple.tsv(Sex_SummarStat); md.tableWriter.DF.w.dimnames(Sex_SummarStat)

subb = paste("Median length: ", paste(MedianSex, collapse = " | "))
wstripchart_list(WormLengthBySex, sub = subb, ylim= c(0, max(WormLength)), bg =SexColors, incrBottMarginBy = 1)


llprint("### Overall Quality seems to correlate with the number of HQ slices.")
try(dev.off(), silent = T)

corr_OverallQual_vs_SliceCount =  cbind.data.frame(NrOfUsefulSections[ worms_HQ], "TotalGeneCount log2(#)" = log2(TotalGeneCount[worms_HQ]))

wplot(corr_OverallQual_vs_SliceCount, col =WormCol_Sex[worms_HQ], cex =2)
text(corr_OverallQual_vs_SliceCount, labels = worms_HQ, srt=-45)
wplot_save_this(plotname = plotnameLastPlot, mdlink = T)

llprint("### Worm length seem to correlate less")

corr_OverallQual_vs_WormLength =  cbind.data.frame(WormLength[worms_HQ], "TotalGeneCount log2(#)" = log2(TotalGeneCount[worms_HQ]))
plot(corr_OverallQual_vs_WormLength, pch = 18, col =WormCol_Sex[worms_HQ], cex =2) # xlim =c(min(WormLength),max(WormLength)), ylim =c(16, 19)
text(corr_OverallQual_vs_WormLength, labels = worms_HQ, srt=-45)
wplot_save_this(plotname = plotnameLastPlot, mdlink = T)


# Filtered Data tables ----------------------------------------------------
llprint("## Lowly expressed genes are filtered out")
llprint("We keep genes that are expressed min.", MinExpr, " reads in at least", MinSlices, "slices. See sfData.")

GeneExprPerSample = sfData = nData = zData = list.fromNames(worms_HQ)

i=1
HE_Genes =  list.fromNames(worms_HQ)
for (i in 1:NrHQworms) {

  w = worms_HQ[i]
  WormSlices = Boundaries[w,1]:Boundaries[w,2]
  "We keep gene"
  HE_Genes[[w]] = rowSums(ExpData[[w]] >= MinExpr) >= MinSlices
  llprint("In worm",w ," - ", sum(HE_Genes[[w]]),"or",pc_TRUE(HE_Genes[[w]]), "of the transcripts have min.",MinExpr, "in at least", MinSlices,"sections.")
  HE_genes = which_names(HE_Genes[[i]])

  "WormSlices are all slices between the head and the tail regardless of thier quality"
  sfData[[w]] = iround(ExpData[[w]][ HE_genes, WormSlices ])
  GeneExprPerSample[[w]]  = rowSums(sfData[[w]])

  # z-score transformed Data tables ----------------------------------------------------
  zData[[w]]  = iround(t(scale(t(sfData[[w]]))))

  # Normalized Data tables ----------------------------------------------------
  # SeqDepthNormFactor is calculated on the raw unfiltered ExpData in '01b.QC.SliceBarplots.R'.
  nData[[w]] = sfData[[w]] * SeqDepthNormFactor[w]

}

# TPM transcript per million ----------------------------------------------------
tpmData = lapply(sfData, TPM_normalize, SUM = 1e6) # to a million per section
"Z-score transformation after TPM normalization reduces technical slice-to-slice effects, but erases the huge germline-to-non-germline difference in RNA-content."
tpm.zData = lapply(tpmData, function(x) iround(t(scale(t(x)))) )
tpmData = lapply(tpmData, iround) # round only after z-scores are calculated


try.dev.off()


# If Venn diagrams for Fig1.G ----------------------------------------------------
Fig1.G = T
if (Fig1.G) {
  ls_gene.names = lapply(sfData, rownames)
  create_set_SubDir( "Fig.1G.Venn.Diagrams", setDir = T)
  for (i in 1:l(HQ_per_Sex) ) {
    SX = HQ_per_Sex[[i]]
    for (percent in T:F )  wvenn(ls_gene.names[SX],  plotname = kollapse("Reproducible.HE.genes.in.sfData.", SexLabels[i], nameiftrue(percent,prefix="." )), print.mode=c('raw' , 'percent')[percent+1], mdlink = percent, fontfamily =rep("Arial", 1),  main.fontfamily = "Arial", sub.fontfamily ="Arial"  )
  }
  create_set_Original_OutDir()
} #if


# Write out ----------------------------------------------------
WriteOut_sf_and_zData=F

for (i in 1:NrHQworms) {
  if (WriteOut_sf_and_zData) {
    w = worms_HQ[i]
    fname = paste0("/Worm_", w, "_len_", WormLength[w], "_between_", Boundaries[w, 1],"_", Boundaries[w, 2],".tsv")
    write.simple.tsv(sfData[[i]], ManualName = paste0(WormSlicesFolder,  fname), gzip = T)
    write.simple.tsv(zData[[i]], ManualName = paste0(Z_ScoreNormFolder, fname), gzip = T)
    write.simple.tsv(nData[[i]], ManualName = paste0(NormSeqDepthFolder, fname), gzip = T)
  }

   # Write out with Clone Name ----------------------------------------------------
  sfData.cloneNames = sfData
  zData.cloneNames = zData
  nData.cloneNames = nData

  rownames(sfData.cloneNames[[i]]) = CloneAndGene_names[[i]][rownames(sfData[[i]])]
  rownames(zData.cloneNames[[i]]) = CloneAndGene_names[[i]][rownames(zData[[i]])]
  rownames(nData.cloneNames[[i]]) = CloneAndGene_names[[i]][rownames(nData[[i]])]

  if (WriteOut_sf_and_zData) {
    fname = paste0("/Worm_", w, "_len_", WormLength[w], "_between_", Boundaries[w, 1],"_", Boundaries[w, 2],".CLONE.NAMES.tsv")
    write.simple.tsv(sfData.cloneNames[[i]], ManualName = paste0(WormSlicesFolder,  fname), gzip = T)
    write.simple.tsv(zData.cloneNames[[i]], ManualName = paste0(Z_ScoreNormFolder, fname), gzip = T)
    write.simple.tsv(nData.cloneNames[[i]], ManualName = paste0(NormSeqDepthFolder, fname), gzip = T)
  }
}



# Region and relative position ------------------------

RelativePosition.nData = list.fromNames(worms_HQ)
for (i in 1:NrHQworms ) {
  WRM = worms_HQ[i];WRM
  RL =colnames(nData[[WRM]])
  ABSPOZ = as.numeric(substr(RL, 2, 3))
  ABSPOZ = 1+ABSPOZ-min(ABSPOZ) # offset
  X = iround(ABSPOZ/max(ABSPOZ))
  names(X) = RL
  RelativePosition.nData[[WRM]] = X
} #for

# ----------------------------------------------------

asexuals = which_names(sex==1)
males = which_names(sex==2)
hermaphrodites = which_names(sex==3)


GeneCorZscoreSpearman = GeneCorZscorePearson = list.fromNames(worms_HQ)
Zcorrelation = F
if (Zcorrelation) {
  "This generates 2.5 GB of files"
  Gcor = "/Users/abelvertesy/Google_Drive/Avano/CT/GeneCorrelations/"
  for (i in 10:length(worms_HQ) ) {
    w = worms_HQ[8]; print(w)

    GeneCorZscoreSpearman[[i]] = iround(cor(t(zData.cloneNames[[i]]), method = "spearman") )
    fname = paste0(Gcor, "GeneCorZscoreSpearman.", w, "_len_", WormLength[w], "_between_", Boundaries[w, 1],"_", Boundaries[w, 2],".CLONE.NAMES.tsv")
    write.simple.tsv(GeneCorZscoreSpearman[[i]], ManualName = fname, gzip = T)

    GeneCorZscorePearson[[i]]  = iround(cor(t(zData.cloneNames[[i]]), method = "pearson") )
    fname = paste0(Gcor, "GeneCorZscorePearson.", w, "_len_", WormLength[w], "_between_", Boundaries[w, 1],"_", Boundaries[w, 2],".CLONE.NAMES.tsv")
    write.simple.tsv(GeneCorZscorePearson[[i]], ManualName = fname, gzip = T)

  } #for
} #if

save.image("CT.2018.08.21.Rdata")




