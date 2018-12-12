######################################################################
# 05.2.Sectioning.Herma.2017.12.12.R
######################################################################
# source("~/Github_repos/CeleTomo/05.2.Sectioning.Herma.2017.12.12.R")
try(dev.off(), silent = T)



# Functions ------------------------
findBoundary <- function(names_vec_of_categ =SectionAnnotation[[w]] ) {
  firstOfEach = rle(names_vec_of_categ)$values # first elements with corresponding names
  namefirst = flip_value2name(firstOfEach) # get the names
  which(names(names_vec_of_categ) %in% namefirst) # get the index
}

# Setup ------------------------
Fig.2.AB.panels = T
Fig.3.panels = T
Fig.S3.panels =T

SilentIntermediateHeatmaps = T
Show.the.Zscore.of.Markers = F # show normalized readcounts instead(nData)

HQ_herma_worms =  worms_HQ[7:10]

# SpTheca = "fkh-6"
SpermGene = "msp-3"
VulvaGene = "pes-8"

findGene.CeleTomo("fkh")


ManualMarkerGenes.Herma = c( 'pgp-14', 'flp-3', 'flp-1', 'ceh-37', 'glh-1', 'msp-3', 'fkh-6', 'pes-8', 'sdz-24')
names(ManualMarkerGenes.Herma) = p0("R", c(1:8,12))


# Metadata ------------------------
HermaSections.Final = c("Head1", "Head2", "Germ1", "Sperm", "Repr.Tract", "Vulva", "Tail")
HermaSections.Intermediate = c("HeadGerm", "TailGerm", "Germ1.ant")
Sections.ls = list.fromNames( c(HermaSections.Final, HermaSections.Intermediate) )
poz.Vulva = poz.SpermPoz = SectionAnnotatation.Herma = list.fromNames(HQ_herma_worms)

SectionSizes.Herma = matrix.fromNames(rowname_vec = HermaSections.Final, colname_vec = HQ_herma_worms)

# Find the Vulva and the Sperm ------------------------
try.dev.off()
pdfA4plot_on(pname = paste("Fig.S3a.Vulva.and.Sperm.Positions", VulvaGene, SpermGene, sep = '.'), cols = 2)
i=2
for (i in 1:l(HQ_herma_worms)) {
  w = HQ_herma_worms[i]; print(w)

  # Find the Vulva ------------------------
  VulvaGeneExpression = zData[[w]][VulvaGene, HQ_Slices_above100[[w]]]
  poz.Vulva[[w]] = names(which.max(VulvaGeneExpression))
  CCC = (HQ_Slices_above100[[w]] %in% poz.Vulva[[w]]) +2
  wbarplot(VulvaGeneExpression, main=p0("Vulva in ", w, " (",VulvaGene,")"), ylab = p0(VulvaGene, " expression (z-Score)"), col = CCC)

  # Find and the Sperm ------------------------
  colnames(zData[[w]])
  SpermGeneExpression = zData[[w]][SpermGene, HQ_Slices_above100[[w]]]
  idx.ant.or.post = (HQ_Slices_above100[[w]] > poz.Vulva[[w]])
  ls.SpermGeneExpression = split(SpermGeneExpression, idx.ant.or.post)
  SpermPoz = unlapply(ls.SpermGeneExpression, which.max)
  poz.SpermPoz[[w]] = names(SpermPoz) = substrRight(names(SpermPoz),2)

  CCC = (HQ_Slices_above100[[w]] %in% poz.SpermPoz[[w]]) +2
  wbarplot(SpermGeneExpression, main=p0("Sperm in ", w, " (",SpermGene,")"), ylab = p0(SpermGene, " expression (z-Score)"), col = CCC)

}# for
pdfA4plot_off()



# Reiterate ------------------------
i=2
for (i in l(HQ_herma_worms):1) {
  w = HQ_herma_worms[i]; print(w)

  Sections.ls$'Sperm' = poz.SpermPoz[[w]]
  Sections.ls$'Vulva'       = poz.Vulva[[w]]

  # Find sections defined by marker genes ------------------------


  idx.anterior        = HQ_Slices_above100[[w]] < min(Sections.ls$'Sperm')
  idx.poster          = HQ_Slices_above100[[w]] > max(Sections.ls$'Sperm')
  idx.notVulva        = (HQ_Slices_above100[[w]] != Sections.ls$'Vulva')
  idx.betweenSperm  = (HQ_Slices_above100[[w]] > min(Sections.ls$'Sperm') & HQ_Slices_above100[[w]] < max(Sections.ls$'Sperm') )

  idx.anterior        = HQ_Slices_above100[[w]] < min(Sections.ls$'Sperm')
  idx.poster          = HQ_Slices_above100[[w]] > max(Sections.ls$'Sperm')
  idx.notVulva        = (HQ_Slices_above100[[w]] != Sections.ls$'Vulva')
  idx.betweenSperm  = (HQ_Slices_above100[[w]] > min(Sections.ls$'Sperm') & HQ_Slices_above100[[w]] < max(Sections.ls$'Sperm') )

  Sections.ls$'HeadGerm'	= HQ_Slices_above100[[w]][idx.anterior]
  Sections.ls$'TailGerm'	= HQ_Slices_above100[[w]][idx.poster]
  Sections.ls$'Repr.Tract'		  = HQ_Slices_above100[[w]][idx.notVulva & idx.betweenSperm]


  # Find remaining sections by subclustering ------------------------
  HGS = Sections.ls$'HeadGerm'
  CorMat =corDataPearson[[w]]

  # Subclustering of anteriror part ------------------------
  try.dev.off()
  CUTinto=3
  zhh = pheatmap(CorMat[HGS,HGS], col = colorRamps::matlab.like(50), cutree_cols = CUTinto,cutree_rows = CUTinto,  main = w,  silent = SilentIntermediateHeatmaps)
  hzz = hclust.getClusterID.col(zhh, k = CUTinto)

  ls_head_IDs = splititsnames_byValues(hzz)
  head.sect.starts = names(sort(unlapply(ls_head_IDs, min), decreasing = F)) # sort the most anterior

  Sections.ls$'Head1'	= ls_head_IDs[[head.sect.starts[1]]]
  Sections.ls$'Head2'	= ls_head_IDs[[head.sect.starts[2]]]
  Sections.ls$'Germ1.ant'	= ls_head_IDs[[head.sect.starts[3]]]

  # Subclustering of posteriror part ------------------------
  # if (w ==  "e2_08" ) { # e2_08 is missing the tail
  #   Sections.ls$'Germ1'	= c(Sections.ls$'Germ1.ant',  Sections.ls$'TailGerm')
  #   Sections.ls$'Tail'	= NULL
  # } else if  (w != "e2_08" ) { # e2_08 is missing the tail
    try.dev.off()
    TGS = Sections.ls$'TailGerm'
    ztt = pheatmap(CorMat[TGS,TGS], col = colorRamps::matlab.like(50), cutree_cols = 2,cutree_rows = 2,  main = w,  silent = SilentIntermediateHeatmaps)
    tzz = hclust.getClusterID.col(ztt, k = 2)

    ls_tail_IDs = splititsnames_byValues(tzz)
    tail.sect.starts = names(sort(unlapply(ls_tail_IDs, min) ) )
    idx.germ1 = tail.sect.starts[1]
    idx.tail = tail.sect.starts[2]

    Sections.ls$'Tail'	= ls_tail_IDs[[idx.tail]]
    Sections.ls$'Germ1'	= c(Sections.ls$'Germ1.ant', ls_tail_IDs[[idx.germ1]])
  # }
  SectionRanges.ls = lapply(Sections.ls, range) # calulcate the ranges
  Sections.ls.Final = Sections.ls[HermaSections.Final] # Final categories
  SectionRanges.ls.Final = SectionRanges.ls[HermaSections.Final]
  SectionSizes.Herma[ , w] = unlapply(Sections.ls.Final, l) # Stats

  # Plot Sectioning Results ------------------------------------------------------------


  SectionAnnotatation.Herma[[w]] = sortbyitsnames(list.2.replicated.name.vec(Sections.ls.Final))
  idx.Slices = unlist(Sections.ls.Final)
  COLANN = as.data.frame(list(
    # "Sections" = gsub('[0-9]+', '', SectionAnnotatation.Herma[[w]])  ,
    "Sections" = SectionAnnotatation.Herma[[w]],
    "TranscriptCounts" = iround((TranscriptCountsRaw[[w]][idx.Slices])), # FROM RAW DATA (expData)
    "GeneCounts" = iround((GeneCountsRaw[[w]][idx.Slices])),
    "Manual.Region" = SectionAnnotatation.Herma[[w]]
  ))
  annot_col.create.pheatmap.df(CorMat, annot_df_per_column =  COLANN)
  annot_col.fix.numeric(ListOfColnames = names(annot_col)[2:3])
  annot_col$Sections = annot_col$Manual.Region = UnionColors.herm[names(annot_col$Sections)  ]

  gapz = cumsum(rle(as.character(COLANN$'Sections'))$'lengths')

  MarkerExpressions = if (Show.the.Zscore.of.Markers) {    t(zData[[w]][ManualMarkerGenes.Herma, HQ_Slices_above100[[w]] ])
                      } else {                        t(nData[[w]][ManualMarkerGenes.Herma, HQ_Slices_above100[[w]] ])  } #if
  annot_row.create.pheatmap.df(CorMat, MarkerExpressions)


  if (Fig.S3.panels) {
    pname =p0("Data.S3.", worms_HQ.proper.name[w],".clustering.based.sequential.sectioning")
    zz=pheatmap(CorMat, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F, main = pname,
                gaps_col =  gapz, gaps_row =  gapz, annotation_col = annot, annotation_row = annot_rows, annotation_colors = c(annot_col))
    wplot_save_this(pname, w=20, h=17)
  } #if



  if (Fig.3.panels) {

    pname =p0("Fig.3B.",worms_HQ.proper.name[w],".clustering.based.sequential.sectioning.Main.Figure")
    zz=pheatmap(CorMat, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F, main = pname,
                gaps_col =  gapz, gaps_row =  gapz, annotation_col = annot, annotation_colors = annot_col)
    wplot_save_this(pname, w=10, h=8)
  } #if

  if (Fig.2.AB.panels) { try.dev.off()

    # pname =p0(w,".clustering.based.sequential.sectioning.genes.Zscore")
    # d =  zData[[w]][, HQ_Slices_above100[[w]]]; idim(d)

    pname =p0("Fig.2A.", worms_HQ.proper.name[w],".clustering.based.sequential.sectioning.genes.tpm.Zscore")
    d =  tpm.zData[[w]][, HQ_Slices_above100[[w]]]; idim(d)

    d = remove.na.rows(d); range(d)
    d= clip.outliers(d) # signal clipping; clip the exterme 1-1% of the z-score distribution

    zz=pheatmap(d, col = colorRamps::matlab.like(50), cluster_cols = F, show_rownames = F, treeheight_row = 0, # , cutree_rows = l(gapz)
                annotation_col = annot, annotation_colors = annot_col, # gaps_col =  gapz,
                main = p0( "Z-score Gene Expression ", w), labels_row = F)
    wplot_save_this(pname)

    } #if

} # for

md.tableWriter.DF.w.dimnames(SectionSizes.Herma, WriteOut = T)
