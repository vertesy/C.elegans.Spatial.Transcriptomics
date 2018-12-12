######################################################################
# 02.CeleTomo_SampleComparison.R
######################################################################
# source ('~/Github_repos/CeleTomo/02.CeleTomo_SampleComparison.R')
try(dev.off(), silent = T)

# Functions ------------------------

# Setup ------------------------
OutDir = "~/Google_Drive/Avano/CT/Analysis/SampleComparisons"
setup_MarkdownReports(OutDir = OutDir, scriptname = "02.CeleTomo_SampleComparison.R", title = "Sample Comaparison for C. elegans TomoSeq experiments.", append = F)

PairwisePlots =F
thr_log2fc =1

# Set to equal length / common genes only ------------------------------------
GeneNames = lapply(GeneExprPerSample, names)

AllGenes = Reduce(union, GeneNames); l(AllGenes)
GeneExprPerSample_ = matrix(data = 0, nrow =  l(AllGenes), ncol = NrHQworms)
rownames(GeneExprPerSample_) = AllGenes;
colnames(GeneExprPerSample_) = worms_HQ


# Add up slices and unify them in a single data frame ------------------------
i=1
for (i in 1:NrHQworms) {
  w = worms_HQ[i]
  x = rowSums(sfData[[i]])  # Add up slices
  GeneExprPerSample_[ names(x), i] = x
}

TrPerWorm=as.data.frame(GeneExprPerSample_); dim(TrPerWorm); dimnames(TrPerWorm)
whist(log2(unlist(TrPerWorm)))

MinWorms = 2
MinExprPerWorm = 500
RowIndex = (rowSums(TrPerWorm>=MinExprPerWorm, na.rm = T)>= MinWorms);
llprint("We keep", sum(RowIndex), "genes above", MinExprPerWorm, "transcripts in at least",MinWorms," wroms")
TrPerWor = TrPerWorm[ RowIndex,]

# Pairwise Scatterplots ------------------------

if (PairwisePlots){
  # continue_logging_markdown(scriptname = "01.QC_CeleTomo.R")
  llprint("## Pairwise Worm-to-Worm Gene Expression Correlations reflect gender differences")
  # linear ------------------------------------------------------------------------------------------------------
  par("pch" =".")
  cormethod = "pearson"
  pname = "Worm-to-Worm Pearson Correlation in Linear Expression Space"

  pairs(TrPerWorm, lower.panel=panel.smooth, upper.panel=panel.cor, main = pname) # , cex.labels = .5
  wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)

  cormethod = "spearman"
  pname = "Worm-to-Worm Spearman Correlation in Linear Expression Space"
  pairs(TrPerWorm, lower.panel=panel.smooth, upper.panel=panel.cor, main = pname) # , cex.labels = .5
  wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)

  # log2 ------------------------------------------------------------------------------------------------------
  TrPerWorm_log2 =  log2(TrPerWorm+1)

  cormethod = "pearson"
  pname = "Worm-to-Worm Pearson Correlation in Log2 Expression Space"
  pairs(TrPerWorm_log2, lower.panel=panel.smooth, upper.panel=panel.cor, main = pname) # , cex.labels = .5
  wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)


  cormethod = "spearman"
  pname = "Worm-to-Worm Spearman Correlation in Log2 Expression Space"
  pairs(TrPerWorm_log2, lower.panel=panel.smooth, upper.panel=panel.cor, main = pname) # , cex.labels = .5
  wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)
}



# Parameters ------------------------
MinWormExpr = 64 # 2^4=16 transcripts in the whole worm
MinWorms = 2 # seen above MinWormExpr in at least MinWorms

cormethod_HM =  "spearman"
cormethod_HM =  "pearson"

# Read In ------------------------


index = na.omit.strip(rowSums(TrPerWorm >= MinWormExpr)>=MinWorms)
llprint("We keep", sum(index), "genes for hclustering worms. Criteria: min",MinWormExpr,"observed transcripts in at least",MinWorms,"worms.")

TrCountperHQworm = colSums(TrPerWorm, na.rm = T)
wbarplot(TrCountperHQworm, mdlink = T)



# Do ---pheatmap::---------------------
llprint("## Sample Clustering")
llprint("Clustering seem to work best on spearman correlation of Z-score normalized raw transcript counts(no median normalization of samples).")
TrPerWorm_HE = TrPerWorm[ index,]

TrPerWorm_norm = median_normalize(TrPerWorm_HE) # Median Normalize Samples
barplot(colSums(TrPerWorm_norm))

idim(TrPerWorm_HE)
TrPerWorm_norm_Zscore = iround(t(scale(t(TrPerWorm_norm))), digitz = 2); # Z-score normaoze genes
TrPerWorm_HE_Zscore = iround(t(scale(t(TrPerWorm_HE))), digitz = 2); # Z-score normaoze genes

try(dev.off(), silent = T)
# pheatmap(TrPerWorm_norm_Zscore)

pname = kollapse("Worms clusters by sex using ",cormethod_HM," correlation of total expression")
Wcor = cor(TrPerWorm_HE_Zscore, method =  cormethod_HM)
Wcor = cor(TrPerWorm_norm, method =  cormethod_HM)
# write.simple.tsv(TrPerWorm_norm)

SexLabels
SEX = SexLabels[sex_HQ[colnames(Wcor)]]
names(SEX) = colnames(Wcor)
annot_col.create.pheatmap.vec(data = Wcor, annot_vec =SEX )
annot_col$SEX = SexColors2Labeled
pheatmap::pheatmap(Wcor, main = pname, cutree_rows = 3, cutree_cols = 3, annotation_col = annot, annotation_colors = annot_col)
wplot_save_this(plotname = pname, mdlink = T)

