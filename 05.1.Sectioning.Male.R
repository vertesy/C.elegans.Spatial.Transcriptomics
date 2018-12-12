######################################################################
# 05.1.Sectioning.Male.R
######################################################################
# source("~/Github_repos/CeleTomo/05.1.Sectioning.Male.R")
try(dev.off(), silent = T)


# Functions ------------------------

# Setup ------------------------
AllGeneClusteringSectioning = F

# Metadata ------------------------


# Go ------------------------

for (i in l(HQ_male_worms):1) {
  w = HQ_male_worms[i]; print(w)
  CorMat = corDataPearson[[w]]
  zz = pheatmap(CorMat, col = colorRamps::matlab.like(50), cutree_cols = 3, main = w)

  hzz = hclust.getClusterID.col(zz)
  newNames = substr(names(hzz),start = 2, stop = 4)
  names(hzz) = newNames
  ls_IDs = splititsnames_byValues(hzz)
  nzz  = names(sort(unlapply(ls_IDs, l))[1:2]) # Select the 2 smallest clusters, they are the germ clusters in all males
  RangesGerm = lapply(ls_IDs[nzz], range)
  RangesGerm = lapply(RangesGerm, as.numeric)
  names(RangesGerm) =  2:3

  HeadCLuster = c(min(newNames), min(unlist(RangesGerm))-1) # WTF???
  RangesGerm$`1` = HeadCLuster

  TailCLuster = c(as.numeric(max(unlist(RangesGerm)))+1, max(newNames) )
  RangesGerm$`4` = TailCLuster

  RangesGerm = sortbyitsnames(RangesGerm)# reorder
  RangesGerm = lapply(RangesGerm, as.numeric); RangesGerm

  Subclusters.of.the.Head = T
  if (Subclusters.of.the.Head) {
    pname = p0("Subclusters.of.the.Head.", w)
    rzz = RangesGerm[[1]]
    SubSlices = p0("X",rzz[1]:rzz[2])
    SubSlices = intersect(colnames(CorMat), SubSlices) # Only exisitng slices

    zhh = pheatmap(CorMat[ SubSlices, SubSlices ], col = colorRamps::matlab.like(50), cutree_rows = 2, cutree_cols = 2,
                   main = pname)
    wplot_save_this(pname)

    hzz = hclust.getClusterID.col(zhh, k = 2)

    ls_IDs = splititsnames_byValues(hzz)
    HeadClusters = lapply(ls_IDs, range)

  } #if

  CombineBoundaries = T
  if (CombineBoundaries) {

    GermClusters = RangesGerm[2:4]
    GermClusters  =lapply(GermClusters, function(x) p0("X",x)) # Change to properslice name

    AllBoundaries = c(HeadClusters, GermClusters)
    names(AllBoundaries) = c("head1", "head2", "germ1", "germ2", "tail")
    AllBoundaries_df = t(list2df(AllBoundaries))

    gapz = which(colnames(CorMat) %in% AllBoundaries_df[ ,2])
    zz=pheatmap(CorMat, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F, gaps_col =  gapz)


    ColAnnot = T
    if (ColAnnot) {
      Sizez = gapz - c(0,gapz[-l(gapz)])

      SectionAnnotation[[w]] = VecAnn = rep(names(AllBoundaries), Sizez)
      names(VecAnn) = HQ_Slices[[w]]

      annot_col.create.pheatmap.vec(data = CorMat, annot_vec = VecAnn, annot_names = "Sections")
      SliceAnnot[[w]] = annot

      pname =p0(w,".clustering.based.sequential.sectioning")
      zz=pheatmap(CorMat, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F,
                  gaps_col =  gapz, annotation_col = annot, annotation_colors = annot_col,
                  main = pname)
      wplot_save_this(pname)

      if (AllGeneClusteringSectioning) { try.dev.off()
        pname =p0(w,".clustering.based.sequential.sectioning.genes.Zscore")

        d =  zData[[w]][, HQ_Slices[[w]]]; idim(d)
        zz=pheatmap(d, col = colorRamps::matlab.like(50), cluster_cols = F, cutree_rows = 6, show_rownames = F,
                    gaps_col =  gapz, annotation_col = annot, annotation_colors = annot_col,
                    main = p0( "Z-score Gene Expression ", w), labels_row = F)
        wplot_save_this(pname)
      } #if

    } #if ColAnnot

  } #if

} # for

