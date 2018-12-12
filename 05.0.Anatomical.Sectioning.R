######################################################################
# 05.0.Anatomical.Sectioning.R
######################################################################
# source ('~/Github_repos/CeleTomo/05.0.Anatomical.Sectioning.R')
try(dev.off(), silent = T)

# Functions ------------------------
source ('~/Github_repos/CeleTomo/zz.Functions.CeleTomo.R')
require(pheatmap)
require(DESeq2)




# Setup ------------------------
setup_MarkdownReports(OutDir = p0(ProjectDir,"/Anatomical.Sectioning"), scriptname = "05.0.Anatomical.Sectioning.R",
                      title = "Anatomical Sectioning of Worms", append = F)

ConsistentColorscales = T
if (ConsistentColorscales) {
  UnionColors = wcolorize(c(  '2Head1', '1Head2', '4Germ1','4Repr.Tract', '6Sperm', '3Tail',  '0Vulva'), set = 'rich',show = T, randomize = F)
  names(UnionColors) = substr(names(UnionColors) ,2,100)
  UnionColors.male = UnionColors.herm =UnionColors
  names(UnionColors.male) = translate(names(UnionColors.herm), 'Repr.Tract', 'Germ2' )
} #if

# Metadata ------------------------
AllGeneClusteringSectioning = F


HQ_male_worms =  which_names(sex_HQ == 2)
HQ_herma_worms =  which_names(sex_HQ == 3)

HQ_sexual_worms = c(HQ_male_worms, HQ_herma_worms)
Nr_HQ_sexual_worms = l(HQ_sexual_worms)
sfData_perSection = nData_perSection = list.fromNames(HQ_sexual_worms)


# Correlations for all ------------------------------------------------------------------------------------------------------------
corDataPearson =  list.fromNames(worms_HQ);worms_HQ
try(dev.off(), silent = T)
i=8


HQ_Slices_above100 = list.fromNames(worms_HQ)
for (i in 1:NrHQworms ) { w = worms_HQ[i]
  HQ_Slices_above100[[w]] = which_names(GeneCountsRaw[[w]][HQ_Slices[[w]]] >= 100)
} #for



for (i in NrHQworms:1) {
  w = worms_HQ[i]; print(w)
  ds = select.rows.and.columns(df=nData[[w]][ , HQ_Slices_above100[[w]]], RowIDs = HE_genes, ColIDs = HQ_Slices[[w]] )
  # ds = select.rows.and.columns(df=tmpData[[w]], RowIDs = HE_genes, ColIDs = HQ_Slices[[w]] )
  corDataPearson[[w]] = cor(ds, method = "pearson")
  # pheatmap(corDataPearson[[w]], col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F, main = w)
  # pheatmap(corDataPearson[[w]], col = colorRamps::matlab.like(50), cutree_cols = 3, main = w)
  # wplot_save_this(plotname = p0("Sectioned.heatmap.", w))
}


# Plot for all ------------------------------------------------------------------------------------------------------------
SliceAnnot = list.fromNames(worms_HQ)
SectionAnnotation = lapply(HQ_Slices,vec.fromNames)

# Males
source("~/Github_repos/CeleTomo/05.1.Sectioning.Male.2017.12.12.R")

# Hermaphrodites
source("~/Github_repos/CeleTomo/05.2.Sectioning.Herma.2017.12.12.R")

# Multi.Male.Correlation ------------------------
Multi.Male.Correlation = F
if (Multi.Male.Correlation) { source ('~/Github_repos/CeleTomo/05.3.Multi.Male.Correlation.R') } #if

# Multi.Herma.Correlation ------------------------
Multi.Herma.Correlation = F
if (Multi.Herma.Correlation) { source ('~/Github_repos/CeleTomo/05.4.Multi.Herma.Correlation.R') } #if

# Section.Length.Comparison ------------------------
Section.Length.Comparison = F
if (Section.Length.Comparison) { source ('~/Github_repos/CeleTomo/03.z.AnatomicalSections.Length.Comparison.Old.R') } #if



# Add up and normalise read counts ------------------------
AddUp = T
if (AddUp) {
  for (i in Nr_HQ_sexual_worms:1) {
    w = HQ_sexual_worms[i];print(w)

    ExpressionPerSectionDF = matrix.fromNames(rowname_vec = rownames(sfData[[w]]), colname_vec = unique(SectionAnnotation[[w]]), fill = NA)

    ls_per_Section = colsplit(sfData[[w]][,HQ_Slices[[w]]], f = SectionAnnotation[[w]])
    RowSumzLS = lapply(lapply(ls_per_Section, as.data.frame), rowSums) # You need the as.data.frame beciase vulva is a single vector, and RowSums does not work on it
    ExpressionPerSectionDF = list2df(RowSumzLS) # Add up columns in each section (=each list element)
    sfData_perSection[[w]] = ExpressionPerSectionDF
  }
} #if


# Calculate normalised data ----------------------------------------
nData_perSection = lapply(sfData_perSection, TPM_normalize)
lapply(nData_perSection, colSums)

# write out data ----------------------------------------
WriteOut = T
if (WriteOut) {
  for (i in 1:Nr_HQ_sexual_worms ) {
    tsvName = p0(" s", names(sfData_perSection)[i], ".tsv")
    write.simple.tsv(sfData_perSection[[i]], ManualName = tsvName, gzip = T)

    tsvName = p0(tsvName, ".norm.tsv")
    write.simple.tsv(nData_perSection[[i]], ManualName = tsvName, gzip = T)
  } #for
} #if



# --------------------------------------------------------------------------------

Fig.2AB = F
if (Fig.2AB) {

  heatmapData =  list.fromNames(worms_HQ)
  for (i in NrHQworms:1) { try.dev.off()
    try(dev.off(), silent = T)
    w = worms_HQ[i]; print(w)

    d = heatmapData[[w]] = tpm.zData[[w]][ , HQ_Slices[[w]] ];    idim(heatmapData[[w]])
    d= clip.outliers(d, showhist = T) # signal clipping; clip the exterme 1-1% of the z-score distribution

    mat_breaks <- quantile_breaks(d, n = 51);mat_breaks
    ccc = colorRamps::matlab.like(length(mat_breaks) - 1)

    ANN = cbind('Gene count' = GeneCounts.PerSection.HQ[[w]][HQ_Slices[[w]]]    )
    if (sex_HQ[w] > 1) {          ANN = cbind(ANN, "glh-1" = d["glh-1",], "spe-11" = d["spe-11",]) }
    if (sex_HQ[w] == 2) {         ANN = cbind(ANN, "clec-207" = d["clec-207",])
    } else if (sex_HQ[w] == 3) {  ANN = cbind(ANN, "pes-8" = d["pes-8",]) }

    annot_col.create.pheatmap.df(d, annot_df_per_column = ANN)

    plotnameLastPlot = p0( "Z-score Gene Expression (HE) ", w)
    zz=pheatmap(d, silent = F,col = colorRamps::matlab.like(50), cluster_cols = F, show_rownames = F, treeheight_row = 0, # , cutree_rows = l(gapz)
                annotation_col = annot,# annotation_colors = annot_col,
                main = plotnameLastPlot, labels_row = F)
    wplot_save_this(plotnameLastPlot)

  } # for
} # if


# --------------------------------------------------------------------------------
