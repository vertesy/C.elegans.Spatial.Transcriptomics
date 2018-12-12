######################################################################
# 05.2.Sectioning.Herma.R
######################################################################
# source("~/Github_repos/CeleTomo/05.2.Sectioning.Herma.R")
try(dev.off(), silent = T)



# Functions ------------------------
findBoundary <- function(names_vec_of_categ =SectionAnnotation[[w]] ) {
  firstOfEach = rle(names_vec_of_categ)$values # first elements with corresponding names
  namefirst = flip_value2name(firstOfEach) # get the names
  Boundaries_ = which(names(names_vec_of_categ) %in% namefirst) # get the index
  Boundaries_
}

# Setup ------------------------
AllGeneClusteringSectioning = F

SpThecaGene = "msp-53"
minZ.spTheca =2
VulvaGene = "pes-8"
minZ.Vulva =2


# Metadata ------------------------

# Go ------------------------

try.dev.off()
SelectThreshold = T
ls_SpermathecaSectionIDs = list.fromNames(HQ_herma_worms)
if (SelectThreshold) {
  pdfA4plot_on(pname = p0("SelectThreshold.", SpThecaGene), cols = 3)
  for (i in 1:l(HQ_herma_worms)) {
    w = HQ_herma_worms[i]; print(w)
    wbarplot(zData[[w]][SpThecaGene,], main=w, hline = minZ.spTheca, savefile = F)
    wbarplot(as.numeric(nData[[w]][SpThecaGene,]), main=w, hline = 100, savefile = F)
    wbarplot(as.numeric(tmpData[[w]][SpThecaGene,]), main=w, hline = 500, savefile = F)
    # Select slices where
    # ls_SpermathecaSectionIDs[[w]] = which_names(tmpData[[w]][SpThecaGene,] > 500 )

    # Has to be manually set for e2_05
    ls_SpermathecaSectionIDs[[w]] = if (w == "e2_05") c("X42", "X56") else which_names(zData[[w]][SpThecaGene,] >= minZ.spTheca )

    SectionAnnotation[[w]][ ls_SpermathecaSectionIDs[[w]] ] = "Spermatheca"
  } # for

  pdfA4plot_off()

  ls_VulvaSection = list.fromNames(HQ_herma_worms)
  pdfA4plot_on(pname = p0("SelectThreshold.", VulvaGene), cols = 3)
  HQ_herma_worms =  worms_HQ[7:10]
  i=1
  for (i in 1:l(HQ_herma_worms)) {
    w = HQ_herma_worms[i]; print(w)
    VulvaGeneExpression = zData[[w]][VulvaGene,]
    wbarplot(VulvaGeneExpression, main=w, ylab = p0("zScore.",VulvaGene), hline = 1, savefile = F)
    ls_VulvaSection[[w]] = names(which.max(VulvaGeneExpression))
    SectionAnnotation[[w]][ls_VulvaSection[[w]]] = "Vulva"
  }# for
  pdfA4plot_off()
} #if

# SelectSpermathecaMarkerGene = F
# if (SelectSpermathecaMarkerGene) {
#   w = "e2_08"
#   mxx= apply(zData[[w]],1,which.max)
#   colnames(zData[[w]])[33]
#
#
#   HQ_Slices[[w]]
#   GenesInSpermatheca = mxx[mxx==33]
#
#   pheatmap(zData[[w]][names(GenesInSpermatheca),], cluster_cols = F)
#
#   ndata_w = nData[[w]][names(GenesInSpermatheca),]
#   pheatmap(ndata_w, cluster_cols = F, fontsize_row = 2, fontsize_col = 5)
#   wplot_save_this(plotname ="heatmap.GenesInSpermatheca" )
# } #if



i=2
for (i in 1:l(HQ_herma_worms)) {
  w = HQ_herma_worms[i]; print(w)
  CorMat_herma =corDataPearson[[w]]
  spx = ls_SpermathecaSectionIDs[[w]]
  idx_non_spTheca = setdiff(colnames(CorMat_herma), spx)
  # Remove Spermatheca
  ds_non_spTheca = CorMat_herma[idx_non_spTheca, idx_non_spTheca]

  # Separate front, mid and tail sections
  idx_head = idx_non_spTheca[idx_non_spTheca<min(spx)]
  idx_tail = idx_non_spTheca[idx_non_spTheca>max(spx)]
  idx_mid =  idx_non_spTheca[idx_non_spTheca<max(spx) & idx_non_spTheca>min(spx)]

  ds_head = ds_non_spTheca[ idx_head, idx_head]; idim(ds_head)
  ds_tail = ds_non_spTheca[ idx_tail, idx_tail]; idim(ds_tail)
  ds_mid = ds_non_spTheca[ idx_mid, idx_mid]; idim(ds_mid)
  breaksList = seq(0, 1, by = .02)

  pheatmap(CorMat_herma, cluster_cols = F, cluster_rows = F, col = colorRamps::matlab.like(50))

  # HEAD AND GERM ----------------
  print("Head and Germ")
  idim(ds_head)
  pheatmap(ds_head, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F, breaks = breaksList, cutree_cols = 2 )
  zz = pheatmap(ds_head, col = colorRamps::matlab.like(50), breaks = breaksList, cutree_cols = 3 )
  headClIDs = hclust.getClusterID.col(zz, k = 3)

  LastSlice = max(idx_head) # find the tail-most sections clusterID - its a germ cluster
  GermClustInHead = headClIDs[LastSlice]


  idxGermClean = which_names(headClIDs == GermClustInHead)
  SectionAnnotation[[w]][idxGermClean] <- "Germ1" # name germ

  idxHeadClean = which_names(headClIDs != GermClustInHead)
  SectionAnnotation[[w]][idxHeadClean] <- "Head" # overwrite the head sections




  # MID-GERM AND VULVA ----------------
  pheatmap(ds_mid, col = colorRamps::matlab.like(50), breaks = breaksList)
  # Add annotation to non-vulva midsection
  SectionAnnotation[[w]][idx_mid][ (SectionAnnotation[[w]][idx_mid]) == "0"] <- "Germ2"



  # TAIL  ----------------


  if (w !=  "e2_08" ) { # e2_08 is missing the tail
    zz = pheatmap(ds_tail, col = colorRamps::matlab.like(50), breaks = breaksList, cutree_cols = 2, cutree_rows = 2)
    TailClusers = hclust.getClusterID.col(zz, k = 2)

    TailCluserID = TailClusers[max(names(TailClusers))] # find the tail-most sections clusterID
    idxRealTail = names(TailClusers[TailClusers == TailCluserID])
    SectionAnnotation[[w]][idxRealTail] <- "Tail"
    idxRealTail = names(TailClusers[TailClusers != TailCluserID])
    SectionAnnotation[[w]][idxRealTail] <- "Germ3"
    SectionAnnotation


  } else if  (w == "e2_08" ) { # e2_08 is missing the tail
    pheatmap(ds_tail, col = colorRamps::matlab.like(50), breaks = breaksList)
    # Add annotation to tail sections
    SectionAnnotation[[w]][idx_tail] <- "Germ3"

  }#if

  unique(SectionAnnotation[[w]])
  rle(SectionAnnotation[[w]])$values

  gapz= findBoundary(SectionAnnotation[[w]])

  annot_col.create.pheatmap.vec(data = CorMat_herma, annot_vec = SectionAnnotation[[w]], annot_names = "Sections")

  annot$"msp-53" = zData[[w]][SpThecaGene,rownames(annot)]
  annot$"pes-8" = zData[[w]][VulvaGene,rownames(annot)]

  pname =p0(w,".clustering.based.sequential.sectioning")
  zz=pheatmap(CorMat_herma, col = colorRamps::matlab.like(50), cluster_cols = F, cluster_rows = F,
              gaps_col =  gapz, gaps_row = gapz, annotation_col = annot, annotation_colors = annot_col,
              main = pname )
  wplot_save_this(pname, w=15,h=15)

  AllGeneClusteringSectioning = F
  if (AllGeneClusteringSectioning) { try.dev.off()
    pname =p0(w,".clustering.based.sequential.sectioning.genes")
    zz=pheatmap(d, col = colorRamps::matlab.like(50), cluster_cols = F, cutree_rows = 6,
                gaps_col =  gapz, annotation_col = annot, annotation_colors = annot_col,
                main = pname, labels_row = F)
    wplot_save_this(pname)
  } #if






} # for



FigS.HermaSectioning =T
if (FigS.HermaSectioning) {
  pdfA4plot_on(pname = "FigS.HermaSectioning", cols = 4)
  for (i in 1:l(HQ_herma_worms)) {
    w = HQ_herma_worms[i]; print(w)

    wbarplot(zData[[w]][SpThecaGene,], main=w, hline = minZ.spTheca, savefile = F)
    img_data = select.rows.and.columns(df =zData[[w]], RowIDs = OtherSpThecaGenes)
    image.heatmap.worm(img_data, main = "Spermatheca Marker Gene Z-score expression")

    wbarplot(VulvaGeneExpression, main=w, ylab = p0("zScore.",VulvaGene), hline = 1, savefile = F)
    img_data = select.rows.and.columns(df =zData[[w]], RowIDs = OtherVulvaGenes)
    image.heatmap.worm(img_data, main = "Vulva Marker Gene Z-score expression")
  } # for
  pdfA4plot_off()
} #if



# save.image("YEAH.Rdata")
# load("~/Google_Drive/Avano/CT/Analysis/Sectioning/YEAH.Rdata")
# ~/Google_Drive/Avano/CT/Analysis/Sectioning
