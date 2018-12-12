######################################################################
# 03b.DifferentialExpression.HeadAndBody.R
######################################################################
# source ('~/Github_repos/CeleTomo/03b.DifferentialExpression.HeadAndBody.R')
try(dev.off(), silent = T)

# Functions ------------------------
library(calibrate); # install.packages("calibrate")
source ('~/Github_repos/CeleTomo/zz.Functions.CeleTomo.R')
require(DESeq2)

# Setup ------------------------
OutDir = "~/Google_Drive/Avano/CT/Analysis/DifferentialExpression.Head"
setup_MarkdownReports(OutDir = OutDir, scriptname = "03.DifferentialExpression.R", title = "Whole Worm Differential Expression Analysis", append = F)

# Metadata ------------------------
AnatomicalSections = read.simple.tsv("~/Google_Drive/Avano/CT/Metadata/Anatomy_Manual/AnatomicalSections.tsv")
attach_w_rownames(AnatomicalSections)

# Parameters ------------------------
DiffEx_Unique = T
DiffEx_Pairwise = T

# DEseq ------------------------
thr_padj_head =    0.05
thr_log2fc_head =   log2(3)

# Prepare ------------------------

Sex <- substr(worms_HQ, 2,2); names(Sex ) <- worms_HQ
wormNr<-substr(worms_HQ, 4,5)

# load("~/Downloads/_Annabel.RData")

# DiffEx_Head ------------------------
DiffEx_Head = T
if (DiffEx_Head) {
  llprint("## Male and Hermaphrodite Head comparison")

  nonGLPworms = which_names(sex_HQ>1)
  NrOfnonGLPworms = l(nonGLPworms)
  # Malez = which_names(sex==2)
  # Hermaz = which_names(sex==3)

  Heads = HeadSums = list.fromNames(nonGLPworms)
  for (w in nonGLPworms) {  print(w)
    headendW = which(colnames(sfData[[w]]) == paste0("X",head_end[w]))
    Heads[[w]] = sfData[[w]][ ,1:headendW]
    HeadSums[[w]] = rowSums(Heads[[w]])
  }
  AllHeadGeneNames = sort(unique(unlist(lapply(HeadSums, names))))
  AllHeadGenes = matrix(data=0, nrow = l(AllHeadGeneNames), ncol = NrOfnonGLPworms )
  rownames(AllHeadGenes) = AllHeadGeneNames; colnames(AllHeadGenes) = nonGLPworms

  for (w in nonGLPworms) {  print(w)
    AllHeadGenes[names(HeadSums[[w]]),w]   = HeadSums[[w]]
  }
  wbarplot(colSums(AllHeadGenes))

  AllHeadGenes = median_normalize(AllHeadGenes)

  coldata_ = coldata[ nonGLPworms,]

  DE_isMale <- DESeq(DESeqDataSetFromMatrix("countData" = floor(AllHeadGenes), "colData" = coldata_,    "design" = ~ isMale))

  res_1vs2 <- as.data.frame(results(DE_isMale, contrast=c("isMale","1","0") ) )
  pname = "Differerntial Expression Male VS Herma Heads";
  # Plot results, with accepted hits higlighted -----------
  plotMA(prepare4plotMA(res_1vs2, thr_padj_ = thr_padj_head, thr_log2fc_ = thr_log2fc_head), main=pname); wplot_save_this(plotname = pname, mdlink = T)
  # wplot_Volcano(res_1vs2, thr_padj_ = thr_padj_head, thr_log2fc_ = thr_log2fc_head, pname_ = pname)

  # Select top hits -----------
   llprint("### Filtering")
   Hits_1vs2  = filter_DESeq(DESeq_results = res_1vs2, thr_log2fc_ = thr_log2fc_head, thr_padj_ = thr_padj_head)

   # Add links -----------
   Hits_1vs2  = cbind(Hits_1vs2, "link_wormbase" = link_wormbase(rownames(Hits_1vs2), writeOut = F, Open = F), "link_String" = link_String(rownames(Hits_1vs2), organism ="elegans", writeOut = F, Open = F) )
   write.simple.tsv(Hits_1vs2)
}



xx= t(scale(t(AllHeadGenes)))[rownames(Hits_1vs2),]
dev.off()
pheatmap::pheatmap(xx, cluster_cols = F, fontsize_row = 10)
wplot_save_this(plotname = "Heatmap.Heads.1vs2", h = 20)


AllHeadGenes["her-1",]
