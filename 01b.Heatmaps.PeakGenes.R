######################################################################
# 01b.Heatmaps.PeakGenes.R
######################################################################
# source ('~/Github_repos/CeleTomo/01b.Heatmaps.PeakGenes.R')
# rm(list=ls(all.names = TRUE));
try(dev.off(), silent = T)

# Functions ------------------------
require(pheatmap)

# Setup ------------------------
OutDir = "~/Google_Drive/Avano/CT/Analysis/Heatmaps"
setup_MarkdownReports(OutDir = OutDir, scriptname = "01b.Heatmaps.PeakGenes.R",
                      title = "Global Gene Expression Patterns" )

# Metadata ------------------------
# Z_ScoreNormFolder = "~/Google_Drive/Avano/CT/count_tables/elegans/ZScoreNormalized"

# Parameters ------------------------


PlotPlateauGenes = T
MinExprPlat = 15
MinConsecutiveLength = 2
SectionMinGeneCount_4heatmap = 50

# pheatmap settings
MissingSections = NULL
GeneClusters =    NULL

# Select genes with a broader expression pattern to identify head / body / tail
# Plot ------------------------------------------------------------------------------------------------------
if (PlotPlateauGenes) {
  heatmapData = PlateauGenes = as.list(rep(NA,NrHQworms))
  worms_HQ
  names(PlateauGenes) = names(heatmapData) = worms_HQ


i=1
  for (i in NrHQworms:1) {
    try(dev.off(), silent = T)
    w = worms_HQ[i]; print(w)

    isAboveThr= (nData[[w]] >= MinExprPlat)
    r = apply(isAboveThr, MARGIN = 1, rle) ## count consecutive values; summarize it in a table for each gene (row)
    PlateauGenes[[w]] = which_names(unlapply(r, peakFinder, MinConsecutiveLength = 2)) # Peakfinder checks if there is at least 'MinConsecutiveLength' TRUE values in the rle results.
    NrPlGenes = l(PlateauGenes[[w]])
    llprint("In worm",w ," - ", NrPlGenes,"or",pc_TRUE( NrPlGenes / NROW(nData[[w]])), "of the transcripts have min.",MinExprPlat, "in at least", MinConsecutiveLength,"sections.")

    heatmapData[[w]] = tpm.zData[[w]][ PlateauGenes[[w]], HQ_Slices[[w]] ];    idim(heatmapData[[w]])



    d = iround(t(scale(t(heatmapData[[w]]))))
    mat_breaks <- quantile_breaks(d, n = 51);mat_breaks
    ccc = colorRamps::matlab.like(length(mat_breaks) - 1)

    # Show gaps on heatmap where sections are missing
    if (l(MissingSections)) {
      x= as.numeric(substr(HQ_Slices[[w]], 2,10))
      ind = x[-1]-x[-l(x)]-1
      MissingSections = which(ind==1) }


    # Annotation column showing complexity of the slice
    GenesExpr = GeneCountsRaw[[w]][HQ_Slices[[w]]]; idim(GenesExpr)


    annot_col.create.pheatmap.vec(d, annot_vec = GenesExpr, annot_names = "GenesExpr")

    # Calculate correlation to next slice
    calcCorrelations = T
    if (calcCorrelations) {
      pdfA4plot_on(pname = p0("Correlations2nextSlice.",w))
      sigD = apply(d,2,sign)
      corz=c("Pearson" , "Pearson_log", "Spearman")
      Correlations_ls = list.fromNames(corz)
      names(Correlations_ls)
      sigD = apply(d,2,sign)
      Correlations_ls = list( cor(d, method = "pearson"),
                              cor(log10(abs(d))*sigD, method = "pearson"),
                              cor(d, method = "spearman")  )

      for (k in 1:length(corz) ) {
        dist2next_cor = c(zero.omit(as.numeric(diag(Correlations_ls[[k]][-1,]))), NA);
        names(dist2next_cor) = colnames(d)
        annot[,corz[k]] =dist2next_cor
        wbarplot(dist2next_cor, main =paste(w, corz[k]), col=i , savefile = F)
      } #for
      for (k in 1:length(corz) ) {
        xx = Correlations_ls[[k]]
        # diag(xx)  = NaN
        image(xx, col = colorRamps::matlab.like(50), axes=F)
      }
      pdfA4plot_off()

      k=3
      zz=pheatmap(Correlations_ls[[1]], col = colorRamps::matlab.like(50), cutree_rows = k, cluster_cols = F)
      bz = hclust.ClusterSeparatingLines.row(zz)
      zz=pheatmap(Correlations_ls[[1]][1:bz[1], 1:bz[1]], col = colorRamps::matlab.like(50), cutree_rows = 2, cluster_cols = F)

      bz[1] =  bz[1]+1
      bz[2] =  bz[2]+1
      zz=pheatmap(Correlations_ls[[1]][bz[1]:bz[2], bz[1]:bz[2]], col = colorRamps::matlab.like(50), cutree_rows = k, cluster_cols = F)
      pheatmap(Correlations_ls[[1]], col = colorRamps::matlab.like(50), cluster_rows =F)
      # source ('~/Github_repos/CeleTomo/03c.Sectioning.R')

    } #if

  } # for

  for (i in NrHQworms:1) {
    try(dev.off(), silent = T)
    w = worms_HQ[i]; print(w)

    r = apply(nData[[w]] >= MinExprPlat, MARGIN = 1, rle) ## count consecutive values; summarize it in a table for each gene (row)
    UbiqGeme = NULL
    for (j in 1:l(r)) {
      UbiqGeme[j] = !(l(HQ_Slices[[w]]) <= max(r[[j]]$lengths))
      # stopif(l(HQ_Slices[[w]]) < max(r[[j]]$lengths)) # DOES not work because we work with sfData
    } ; names(UbiqGeme) = names(r)

    y =unlist(lapply(r, peakFinder, MinConsecutiveLength))
    PlateauGenes[[w]] = which_names(y)
    llprint("In worm",w ," - ", l(PlateauGenes[[w]]),"or",pc_TRUE(y), "of the transcripts have min.",MinExprPlat, "in at least", MinConsecutiveLength,"sections.")

    # interestingGenes = intersect(which_names(UbiqGeme), PlateauGenes[[w]])
    heatmapData[[w]] = iround(median_normalize(sfData[[w]][  ,HQ_Slices[[w]] ]));    idim(heatmapData[[w]])
    interestingGenes = intersect(which_names(UbiqGeme), rownames(heatmapData[[w]]))
    heatmapData[[w]] = heatmapData[[w]][interestingGenes,]

    d = iround(t(scale(t(heatmapData[[w]]))))
    mat_breaks <- quantile_breaks(d, n = 51);mat_breaks
    ccc = colorRamps::matlab.like(length(mat_breaks) - 1)

    # Show gaps on heatmap where sections are missing
    if (l(MissingSections)) {
      x= as.numeric(substr(HQ_Slices[[w]], 2,10))
      ind = x[-1]-x[-l(x)]-1
      MissingSections = which(ind==1) }


    # Annotation column showing complexity of the slice
    GenesExpr = GeneCountsRaw[[w]][HQ_Slices[[w]]]; idim(GenesExpr)
    # annot = as.data.frame(GenesExpr)

    annot_col.create.pheatmap.vec(d, annot_vec = GenesExpr, annot_names = "GenesExpr")
    # annot_col.create.old(d, annot_vec = GenesExpr)

    UnsectionedHeatmaps = F
    if (UnsectionedHeatmaps) {
      GeneClusters=6
      NrRowClusters = NULL
      try(dev.off(), silent = T)
      x =pheatmap::pheatmap(d, cluster_cols = F, main = w, show_rownames = F, color = ccc,
                            gaps_col = MissingSections, cutree_rows = GeneClusters, kmeans_k = GeneClusters,
                            annotation_col = annot, annotation_colors = annot_col, annotation_legend = FALSE)#, cutree_rows = NrRowClusters
      plotnameLastPlot = kollapse("Grouped.PlateauGenes.",w,".",MinExprPlat,".reads.per.slices_Sections.above.",SectionMinGeneCount_4heatmap)
      wplot_save_this()

      fname = kollapse(w,".genes.displayed.Z.score.tsv")
      write.simple.tsv(d, ManualName = fname, gzip = T)
      RX = rownames(d)

      fname = kollapse(w,".genes.displayed.RelativeExpression.tsv")
      write.simple.tsv(nData[[w]][RX , ], ManualName = fname, gzip = T)

      x =pheatmap::pheatmap(d, cluster_cols = F, main = w, show_rownames = F, color = ccc,
                            gaps_col = MissingSections,
                            annotation_col = annot, annotation_colors = annot_col, annotation_legend = FALSE)#, cutree_rows = NrRowClusters
      plotnameLastPlot = kollapse(w,"_heatmap_expr",MinExprPlat,".genes_",SectionMinGeneCount_4heatmap)
      wplot_save_this()

    } #if
  } # for

}

