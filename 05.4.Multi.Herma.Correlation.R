######################################################################
# 05.4.Multi.Herma.Correlation.R
######################################################################
# source ('~/Github_repos/CeleTomo/05.4.Multi.Herma.Correlation.R')
# rm(list=ls(all.names = TRUE));
try(dev.off(), silent = T)


# Functions ------------------------
# Setup ------------------------
# Metadata ------------------------
# Parameters ------------------------
cormeth = "pearson"
CUTT = 4
cmm = "ward.D"

# Plot ------------------------------------------------------------------------------------------------------

# AnnotationForCoClustering ------------------------------------------------------------------------------------------------------------

# Co-Correlations for Males ------------------------------------------------------------------------------------------------------------
n_e2_05 = nData$e2_05[ ,HQ_Slices_above100[["e2_05"]] ]
n_e2_06 = nData$e2_06[ ,HQ_Slices_above100[["e2_06"]] ]
n_e2_07 = nData$e2_07[ ,HQ_Slices_above100[["e2_07"]] ]
n_e2_08 = nData$e2_08[ ,HQ_Slices_above100[["e2_08"]] ]
colnames(n_e2_05)

n_e2_05 = prefix.colnames(n_e2_05, prefix = "e2_05.")
n_e2_06 = prefix.colnames(n_e2_06, prefix = "e2_06.")
n_e2_07 = prefix.colnames(n_e2_07, prefix = "e2_07.")
n_e2_08 = prefix.colnames(n_e2_08, prefix = "e2_08.")


# Merge and calculate correlation ------------------------------------------------------------------------------------------------------------
x1=merge_numeric_df_by_rn(n_e2_05, n_e2_06); idim(x1)
x2=merge_numeric_df_by_rn(n_e2_07, n_e2_08); idim(x2)
x=merge_numeric_df_by_rn(x1, x2); idim(x)
HermaCor = cor(x, method = cormeth)
colnames(HermaCor)


# Annotation ------------------------------------------------------------------------------------------------------------


All.Herma.Annot = cbind(
  "Regions" = unlist(SectionAnnotatation.Herma),
  "Hermaphrodite" = as.factor.numeric(substr(names(unlist(SectionAnnotatation.Herma)),1,5))
)

annot_col.create.pheatmap.df( data = HermaCor, annot_df_per_column = All.Herma.Annot)
nn = names(annot_col$Hermaphrodite)
annot_col$Hermaphrodite = terrain.colors(4)
names(annot_col$Hermaphrodite) =nn

annot_col$Regions = annot_col$Manual.Region = UnionColors.herm[names(annot_col$Regions)  ]

# Heatmap ------------------------------------------------------------------------------------------------------------



pname = ppp("Fig.S3.B.All.Herma.Cor", cormeth, cmm, CUTT, idate())
pheatmap(HermaCor, col = colorRamps::matlab.like(50), cutree_rows = CUTT, cutree_cols = CUTT, main = pname, annotation_col = annot, clustering_method = cmm, annotation_colors = annot_col)
wplot_save_this(plotname = pname, w = 18, h=18)



