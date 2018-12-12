######################################################################
# 05.1.Sectioning.Male.R
######################################################################
# source("~/Github_repos/CeleTomo/05.1.Sectioning.Male.2017.12.12.R")
try(dev.off(), silent = T)


# Functions ------------------------

# Setup ------------------------
SilentIntermediateHeatmaps = T
Fig.2.AB.panels = T
Fig.S3.panels =T
Fig.3.panels = T
Subclusters.of.the.Head = T
Show.the.Zscore.of.Markers = T


# Metadata ------------------------
SectionAnnotatation.Male =  list.fromNames(HQ_male_worms)
MaleSections.Final = c("Head1", "Head2", "Germ1", "Germ2", "Tail" )
MaleSections.Intermediate = c("HeadTail", "Head")
Sections.ls = list.fromNames( c(MaleSections.Final, MaleSections.Intermediate) )

SectionSizes.Male = matrix.fromNames(rowname_vec = MaleSections.Final, colname_vec = HQ_male_worms)


ManualMarkerGenes.Male = c( 'gst-41', 'inx-6', 'cwp-1', 'flp-1', 'ceh-37', 'glh-1', 'msp-3', 'nas-19', 'ins-31', 'clec-207', 'cwp-1')
names(ManualMarkerGenes.Male) = p0("R", c(1:11))




# Go ------------------------
i=1
for (i in 1:l(HQ_male_worms)) {
  w = HQ_male_worms[i]; print(w)
  CorMat = corDataPearson[[w]]
  zz = pheatmap(CorMat, col = colorRamps::matlab.like(50), cutree_cols = 3, clustering_method = "ward.D", main = w,  silent = SilentIntermediateHeatmaps)

  hzz = hclust.getClusterID.col(zz)
  # BELOW NOT NECESSARY AS 2017.12.12
  # newNames = substr(names(hzz),start = 2, stop = 4)
  # names(hzz) = newNames
  ls_IDs = splititsnames_byValues(hzz); ls_IDs
  idx.HeadTail = which.max(unlapply(ls_IDs, l)) # Largest cluster is Head & Tail
  Sections.ls$"HeadTail" = ls_IDs[[idx.HeadTail]]

  idx.Germs = setdiff(names(ls_IDs), idx.HeadTail) # Others are the 2 germ cluster
  germ.sect.begin = unlapply(ls_IDs[idx.Germs], min)
  germ.sect.end = unlapply(ls_IDs[idx.Germs], max)
  idx.Germ1 = which.min(germ.sect.begin) # Small section-ID is anterior
  idx.Germ2 = which.max(germ.sect.begin) # Small section-ID is anterior
  # idx.Germ2 = which.max(germ.sect.end) # Large section-ID is posterior

  Sections.ls$"Germ1" = ls_IDs[idx.Germs][[idx.Germ1]]
  Sections.ls$"Germ2" = ls_IDs[idx.Germs][[idx.Germ2]]

  if (w == "e1_01") {
    "In worm 1 slices 60 and 65 cluster together with cluster with  Germ1, although individual correlation values show clear resemblance to germ 2. To resolve we cluster all germ regions separetely into 2 clusters."
    xxx = unlist(ls_IDs[2:3])
    zz = pheatmap(CorMat[xxx,xxx], col = colorRamps::matlab.like(50), cutree_cols = 2, clustering_method = "ward.D", main = w,  silent = SilentIntermediateHeatmaps)
    Germ_male1 = splititsnames_byValues(hclust.getClusterID.col(zz,k = 2))
    Germ1 = which.min(lapply(Germ_male1, min))
    Germ2 = which.max(lapply(Germ_male1, min))
    Sections.ls$"Germ1" = Germ_male1[[Germ1]]
    Sections.ls$"Germ2" = Germ_male1[[Germ2]]
  } #if

  idx.Tail = which(Sections.ls$"HeadTail" > max(germ.sect.end)) # we identified the head region anterior to Germ1.
  idx.Head = which(Sections.ls$"HeadTail" < min(germ.sect.end)) # we identified the Tail region posterior to Germ2.

  if (Subclusters.of.the.Head) {
    SHX = Sections.ls$"Head" = Sections.ls$"HeadTail"[idx.Head]
    Sections.ls$"Tail" = Sections.ls$"HeadTail"[idx.Tail]
    zhh = pheatmap(CorMat[SHX,SHX], col = colorRamps::matlab.like(50), cutree_cols = 2,cutree_rows = 2,  main = w,  silent = SilentIntermediateHeatmaps)
    hzz = hclust.getClusterID.col(zhh, k = 2)
    ls_head_IDs = splititsnames_byValues(hzz)
    head.sect.starts = unlapply(ls_head_IDs, min)
    idx.Head1 = which.min(head.sect.starts)
    idx.Head2 = which.max(head.sect.starts)

    Sections.ls$"Head1" = ls_head_IDs[[idx.Head1]]
    Sections.ls$"Head2" = ls_head_IDs[[idx.Head2]]
  }

  SectionRanges.ls = lapply(Sections.ls, range) # calulcate the ranges
  Sections.ls.Final = Sections.ls[MaleSections.Final] # Final categories
  SectionRanges.ls.Final = SectionRanges.ls[MaleSections.Final]
  SectionSizes.Male[ , w] = unlapply(Sections.ls.Final, l) # Stats

  # Plot Sectioning Results ------------------------------------------------------------

  idx.Slices = unlist(Sections.ls.Final)
  SectionAnnotatation.Male[[w]] = list.2.replicated.name.vec(Sections.ls.Final)
  COLANN = as.data.frame(list(
    "Sections" = SectionAnnotatation.Male[[w]],
    "TranscriptCounts" = iround((TranscriptCountsRaw[[w]][idx.Slices])), # FROM RAW DATA (expData)
    "GeneCounts" = iround((GeneCountsRaw[[w]][idx.Slices])),
    "Manual.Region" = SectionAnnotatation.Male[[w]]
  ))

  annot_col.create.pheatmap.df(CorMat, annot_df_per_column =  COLANN)
  annot_col.fix.numeric(ListOfColnames = names(annot_col)[2:3])
  gapz = cumsum(SectionSizes.Male[ , w])
  annot_col$Sections = annot_col$Manual.Region = UnionColors.male[names(annot_col$Sections)  ]

  Show.the.Zscore.of.Markers=F
  MarkerExpressions = if (Show.the.Zscore.of.Markers) {    t(zData[[w]][ManualMarkerGenes.Male, HQ_Slices_above100[[w]] ])
  } else {                        t(nData[[w]][ManualMarkerGenes.Male, HQ_Slices_above100[[w]] ])  } #if
  annot_row.create.pheatmap.df(CorMat, MarkerExpressions)

  if (Fig.S3.panels) {
    pname =p0("Data.S3.", worms_HQ.proper.name[w],".clustering.based.sequential.sectioning")
    zz=pheatmap(CorMat, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F, main = pname,
                gaps_col =  gapz, gaps_row =  gapz, annotation_col = annot, annotation_row = annot_rows[,1:11], annotation_colors = annot_col)
    wplot_save_this(pname, w=20, h=17)
  }

  if (Fig.3.panels) {
    pname =p0("Fig.3A.", worms_HQ.proper.name[w],".clustering.based.sequential.sectioning.Main.Figure")
    zz=pheatmap(CorMat, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F, main = pname,
                gaps_col =  gapz, gaps_row =  gapz, annotation_col = annot, annotation_colors = annot_col)
    wplot_save_this(pname, w=10, h=8)
  } #if

  if (Fig.2.AB.panels) { try.dev.off()
    pname =p0("Fig.2B.",worms_HQ.proper.name[w],".clustering.based.sequential.sectioning.genes.tpm.Zscore")

    lapply(tpm.zData, idim)
    d =  tpm.zData[[w]][, HQ_Slices_above100[[w]]]; idim(d)

    d = clip.outliers(d, showhist = T) # signal clipping; clip the exterme 1-1% of the z-score distribution
    d = remove.na.rows(d); range(d)
    zz=pheatmap(d, col = colorRamps::matlab.like(50), cluster_cols = F, show_rownames = F, treeheight_row = 0, # , cutree_rows = l(gapz)
                annotation_col = annot, annotation_colors = annot_col, # gaps_col =  gapz,
                main = p0( "Z-score Gene Expression ", w), labels_row = F)
    wplot_save_this(pname)
  } #if
} # for

md.tableWriter.DF.w.dimnames(SectionSizes.Male, WriteOut = T)

GeneCountsRaw$e2_07

# save.image("ct.2017.12.13.RData")
