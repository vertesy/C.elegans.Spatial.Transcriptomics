######################################################################
# 03.DifferentialExpression.R
######################################################################
# source ('~/Github_repos/CeleTomo/03.DifferentialExpression.R')
try(dev.off(), silent = T)

# Functions ------------------------
library(calibrate); # install.packages("calibrate")
require(DESeq2)
source ('~/Github_repos/CeleTomo/zz.Functions.CeleTomo.R')


# Setup ------------------------
OutDir = "~/Google_Drive/Avano/CT/Analysis/DifferentialExpression"
setup_MarkdownReports(OutDir = OutDir, scriptname = "03.DifferentialExpression.R", title = "Whole Worm Differential Expression Analysis", append = F)

# Metadata ------------------------


# Parameters ------------------------
DiffEx_Unique = T
DiffEx_Pairwise = T

# DEseq
thr_padj =    0.01
thr_log2fc =  3

# Prepare ------------------------


# Sex <- substr(worms_HQ, 2,2); names(Sex ) <- worms_HQ
wormNr<-substr(worms_HQ, 4,5)


# Do ------------------------
Sex = sex_HQ -1
coldata <-as.data.frame(cbind(Sex,
          "isAsex" = as.character(as.numeric(Sex==0)),
          "isMale" = as.character(as.numeric(Sex==1)),
          "isHerma" = as.character(as.numeric(Sex==2) )))
dim(TrPerWorm_norm)
dim(coldata)
# DiffEx_Unique ------------------------
if (DiffEx_Unique) {
  llprint("## One sex compared to the two others")
  DE_isAsex <- DESeq(DESeqDataSetFromMatrix("countData" = floor(TrPerWorm_norm), "colData" = coldata,    "design" = ~ isAsex))
  DE_isMale <- DESeq(DESeqDataSetFromMatrix("countData" = floor(TrPerWorm_norm), "colData" = coldata,    "design" = ~ isMale))
  DE_isHemra <- DESeq(DESeqDataSetFromMatrix("countData" = floor(TrPerWorm_norm), "colData" = coldata,    "design" = ~ isHerma))

  res_isAsex  <- as.data.frame(results(DE_isAsex, contrast=c("isAsex","1","0")) )
  res_isMale  <- as.data.frame(results(DE_isMale, contrast=c("isMale","1","0")) )
  res_isHerma <- as.data.frame(results(DE_isHemra, contrast=c("isHerma","1","0")) )

 # Plot results, with accepted hits higlighted -----------
  pname = "Differerntial Expression isAsex";
  plotMA(prepare4plotMA(res_isAsex, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=pname); wplot_save_this(plotname = pname, mdlink = T)
  wplot_Volcano(res_isAsex, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname)

  pname = "Differerntial Expression isMale";
  plotMA(prepare4plotMA(res_isMale, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=pname); wplot_save_this(plotname = pname, mdlink = T)
  wplot_Volcano(res_isMale, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname)

  pname = "Differerntial Expression isHerma";
  plotMA(prepare4plotMA(res_isHerma, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=pname); wplot_save_this(plotname = pname, mdlink = T)
  wplot_Volcano(res_isHerma, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname)

 # Select top hits -----------
 llprint("### Filtering")

 Hits_isAsex  = filter_DESeq(DESeq_results = res_isAsex)
 Hits_isMale  = filter_DESeq(DESeq_results = res_isMale)
 Hits_isHerma  = filter_DESeq(DESeq_results = res_isHerma)

 # Add links -----------
 Hits_isAsex  = cbind(Hits_isAsex, "link_wormbase" = link_wormbase(rownames(Hits_isAsex), writeOut = F, Open = F), "link_String" = link_String(rownames(Hits_isAsex), organism ="elegans", writeOut = F, Open = F) )
 Hits_isMale  = cbind(Hits_isMale, "link_wormbase" = link_wormbase(rownames(Hits_isMale), writeOut = F, Open = F), "link_String" = link_String(rownames(Hits_isMale), organism ="elegans", writeOut = F, Open = F) )
 Hits_isHerma  = cbind(Hits_isHerma, "link_wormbase" = link_wormbase(rownames(Hits_isHerma), writeOut = F, Open = F), "link_String" = link_String(rownames(Hits_isHerma), organism ="elegans", writeOut = F, Open = F) )

 write.simple.tsv(Hits_isAsex)
 write.simple.tsv(Hits_isMale)
 write.simple.tsv(Hits_isHerma)

}

# DiffEx_Pairwise ------------------------
if (DiffEx_Pairwise) {
  llprint("## Pariwise comparison of sexes")
  DE_sex <- DESeq(DESeqDataSetFromMatrix("countData" = floor(TrPerWorm_norm), "colData" = coldata, "design" = ~ Sex))

  res_0vs1  <- as.data.frame(results(DE_sex, contrast=c("Sex","0","1") ) )
  res_0vs2  <- as.data.frame(results(DE_sex, contrast=c("Sex","0","2") ) )
  res_1vs2 <- as.data.frame(results(DE_sex, contrast=c("Sex","1","2") ) )

  # Plot results, with accepted hits higlighted -----------
  pname = "Differerntial Expression 0 vs 1";
  plotMA(prepare4plotMA(res_0vs1, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=pname); wplot_save_this(plotname = pname, mdlink = T)
  wplot_Volcano(res_0vs1, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname)

  pname = "Differerntial Expression 0 vs 2";
  plotMA(prepare4plotMA(res_0vs2, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=pname); wplot_save_this(plotname = pname, mdlink = T)
  wplot_Volcano(res_0vs2, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname)

  pname = "Differerntial Expression 1 vs 2";
  plotMA(prepare4plotMA(res_1vs2, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc), main=pname); wplot_save_this(plotname = pname, mdlink = T)
  wplot_Volcano(res_1vs2, thr_padj_ = thr_padj, thr_log2fc_ = thr_log2fc, pname_ = pname)

 # Select top hits -----------
 llprint("### Filtering")
 Hits_0vs1  = filter_DESeq(DESeq_results = res_0vs1)
 Hits_0vs2  = filter_DESeq(DESeq_results = res_0vs2)
 Hits_1vs2  = filter_DESeq(DESeq_results = res_1vs2)

 # Add links -----------
 Hits_0vs1  = cbind(Hits_0vs1, "link_wormbase" = link_wormbase(rownames(Hits_0vs1), writeOut = F, Open = F), "link_String" = link_String(rownames(Hits_0vs1), organism ="elegans", writeOut = F, Open = F) )
 Hits_0vs2  = cbind(Hits_0vs2, "link_wormbase" = link_wormbase(rownames(Hits_0vs2), writeOut = F, Open = F), "link_String" = link_String(rownames(Hits_0vs2), organism ="elegans", writeOut = F, Open = F) )
 Hits_1vs2  = cbind(Hits_1vs2, "link_wormbase" = link_wormbase(rownames(Hits_1vs2), writeOut = F, Open = F), "link_String" = link_String(rownames(Hits_1vs2), organism ="elegans", writeOut = F, Open = F) )

 write.simple.tsv(Hits_0vs1)
 write.simple.tsv(Hits_0vs2)
 write.simple.tsv(Hits_1vs2)

}

# whist(fromClipboard.as_num_vec())
## Control results
g = "ssq-2"
wstripchart(split(TrPerWorm_norm[ g,], Sex), col = SexColors, colorbyColumn = T, plotname = g)

