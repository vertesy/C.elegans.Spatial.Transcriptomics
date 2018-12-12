######################################################################
# 05.3.Multi.Male.Correlation.R
######################################################################
# source ('~/Github_repos/CeleTomo/05.3.Multi.Male.Correlation.R')
# rm(list=ls(all.names = TRUE));
try(dev.off(), silent = T)


# Functions ------------------------
# Setup ------------------------
# Metadata ------------------------
# Parameters ------------------------
cormeth = "pearson"
CUTT = 6
cmm = "ward.D"

# Plot ------------------------------------------------------------------------------------------------------



# Co-Correlations for Males ------------------------------------------------------------------------------------------------------------
n_e1_01 = nData$e1_01[ , HQ_Slices_above100[["e1_01"]] ]
n_e1_02 = nData$e1_02[ , HQ_Slices_above100[["e1_02"]] ]
n_e1_03 = nData$e1_03[ , HQ_Slices_above100[["e1_03"]] ]
n_e1_04 = nData$e1_04[ , HQ_Slices_above100[["e1_04"]] ]

n_e1_01 = prefix.colnames(n_e1_01, prefix = "e1_01.")
n_e1_02 = prefix.colnames(n_e1_02, prefix = "e1_02.")
n_e1_03 = prefix.colnames(n_e1_03, prefix = "e1_03.")
n_e1_04 = prefix.colnames(n_e1_04, prefix = "e1_04.")



# Merge and calculate correlation ------------------------------------------------------------------------------------------------------------

x1=merge_numeric_df_by_rn(n_e1_01, n_e1_02); idim(x1)
x2=merge_numeric_df_by_rn(n_e1_03, n_e1_04); idim(x2)
x=merge_numeric_df_by_rn(x1, x2); idim(x)
MaleCor = cor(x, method = cormeth)
colnames(MaleCor)

# Annotation ------------------------------------------------------------------------------------------------------------
All.Male.Annot = cbind(
  "Regions" = unlist(SectionAnnotatation.Male),
  "Male" = as.factor.numeric(substr(names(unlist(SectionAnnotatation.Male)),1,5))
)

annot_col.create.pheatmap.df( data = MaleCor, annot_df_per_column = All.Male.Annot)
nn = names(annot_col$Male)
annot_col$Male = terrain.colors(4)
names(annot_col$Male) =nn

annot_col$Regions = annot_col$Manual.Region = UnionColors.male[names(annot_col$Regions)  ]

# Heatmap ------------------------------------------------------------------------------------------------------------
try.dev.off()
pname = ppp("Fig.3C.All.Male.Cor", cormeth, cmm, CUTT, idate())

pheatmap(MaleCor, col = colorRamps::matlab.like(50), cutree_rows = CUTT, cutree_cols = CUTT, clustering_method = cmm, main = pname, annotation_col = annot, annotation_colors = annot_col)
wplot_save_this(plotname = pname, w = 18, h=18)




# create_set_Original_OutDir(NewOutDir = ProjectDir)
# save.image("ct.2017.12.22.RData")
