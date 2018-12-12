######################################################################
# 04.Peak_Genes.R
######################################################################
# source ('~/Github_repos/CeleTomo/02.CeleTomo_SampleComparison.R')
try(dev.off(), silent = T)

# Functions ------------------------



# Setup ------------------------
OutDir = "~/Google_Drive/Avano/CT/Analysis/04.Peak_Genes.R"
setup_MarkdownReports(OutDir = OutDir, scriptname = "04.Peak_Genes.R", title = "Finding Peaking genes in C. elegans TomoSeq experiments.")


# Metadata ------------------------

# Parameters ------------------------

# Read In ------------------------


#  --------------------
llprint("## Peak Gene- and Transcript counts are uncorrelated")

MaxZscores = lapply(zData, rowMax)


Zhi_genes= NULL
thr_x = 4
for (i in 1:NrHQworms) {
  Zhi_genes[[i]] = which_names(MaxZscores[[i]] >=thr_x)
}
names(Zhi_genes) = names(MaxZscores) = worms_HQ
NrOfPeakGenes = unlapply(Zhi_genes, l)
wbarplot(NrOfPeakGenes, sub = kollapse("Genes with a Z-score above: ",thr_x), col = WormCol_Sex)

AllReprPeakGenes = Reduce(intersect, Zhi_genes); l(AllReprPeakGenes)
TotalTrCount_VS_NrOfPeakGenes  = cbind("TotalTrCount" = TotalTrCount[worms_HQ], NrOfPeakGenes)
wplot(TotalTrCount_VS_NrOfPeakGenes, col=WormCol_Sex, cex=2, mdlink = T)
wlegend(fill_ = SexColors, poz = "topleft")

# Venn Diagram, Sex comparison -------------------------------------------------------
ReprPeakGenes_perSex =NULL
for (i in 1:3) {
  ReprPeakGenes_perSex[[i]] = Reduce(intersect, Zhi_genes[which_names(sex[worms_HQ] == i)]); iprint(l(ReprPeakGenes_perSex[[i]]))
}
names(ReprPeakGenes_perSex) = SexLabels
wvenn(ReprPeakGenes_perSex, fill = SexColors)

NrOfReproduciblePeakGenesPerSex =  c(unlapply(ReprPeakGenes_perSex, l), "Across All" = l(AllReprPeakGenes))
wbarplot(NrOfReproduciblePeakGenesPerSex, col = c(SexColors,"grey"),tilted_text = T)


# PeakGenesAll ------------------------

try.dev.off()
pdfA4plot_on(pname = "PeakGenesAll", rows = 1, cols = 1)
for (i in 1:NrHQworms) {
  pkz = zData[[i]][AllReprPeakGenes, ]
  pheatmap::pheatmap(pkz, cluster_cols = F)
}
pdfA4plot_off()


plot(zData[[i]][gene, ], type = "s")

x=(1:3)

plot(x, type = "S")

x <- rnorm(1000)
y <- hist(x)
plot(y$breaks,
     c(y$counts,0)
     ,type="s",col="blue")

# Do ------------------------
c("spe44", "spe9", "spe11")

# Do ------------------------

Males = which_names(sex_HQ==2)
x =similar.genes(gene = "cwp-1", zData[[Males[1]]], number = 25)


# Do ------------------------



lapply(zData, dim)
