######################################################################
# 03.Sample_Alignment.R
######################################################################
# source ('~/Github_repos/CeleTomo/03.Sample_Alignment.R')
try(dev.off(), silent = T)

# Functions ------------------------
# save.image("~/Downloads/_Annabel.RData")
# load("~/Downloads/_Annabel.RData")


# Setup ------------------------
OutDir = "~/Google_Drive/Avano/CT/Analysis/SampleAlignment"
setup_MarkdownReports(OutDir = OutDir, scriptname = "03.Sample_Alignment.R", title = "Gene Expression Maps and Alignment", append = F)

wormcol = rainbow(NrHQworms)

# Metadata ------------------------
ExpressionMarkers = read.simple.tsv("~/Google_Drive/Avano/CT/metadata/Genes/ExpressionMarkers.tsv")
attach_w_rownames(ExpressionMarkers)
ExpressionMarkerNames = rownames(ExpressionMarkers)
ExpressionMarkerNames = c( 'egg-1', 'egl-20', 'myo-1', 'flp-11', 'flp-6')
toClipboard(ExpressionMarkerNames)

metadata_Genes = read.simple.tsv("~/Google_Drive/Avano/CT/metadata/Genes/Genes_to_Transcripts_c_elegans.tsv.gz")

UseRiks = T
if (UseRiks) {
  SecGenes = read.simple.tsv("~/Google_Drive/Avano/CT/Metadata/Anatomy.Genes/Sections.genes.Abel.tsv")
  SecGeneNames =rownames(SecGenes)
  SecGenes_Function = SecGenes$Function; names(SecGenes_Function) =SecGeneNames
  setdiff(SecGeneNames, GeneSymbols)
} #if

MarkerGenes = clone2gene(ExpressionMarkerNames)
if (UseRiks) {
  MarkerGenes = union(MarkerGenes, SecGeneNames)
  Expr = c(Expr, SecGenes_Function)
  } #if
md.tableWriter.VEC.w.names(MarkerGenes)
NrMarkerGenes = l(MarkerGenes);NrMarkerGenes

# Do ------------------------

g = ExpressionMarkerNames[j]


MarkerEx = MarkerEx_log2 = list.fromNames(worms_HQ)
for (i in 1:NrHQworms) {
  w = worms_HQ[i]
  # expr.FOUND = getRows(ExpData.norm[[w]], rownamez = MarkerGenes)
  expr.FOUND = getRows(ExpData[[w]], rownamez = MarkerGenes)
  MarkerEx[[i]]       = expr.FOUND
  MarkerEx_log2[[i]]  = log2(expr.FOUND+2)-1
} ; names(MarkerEx) = names(MarkerEx_log2) = worms_HQ


for (j in 1:l(MarkerGenes)) {
  g = MarkerGenes[j]
  # gg =names(g)
  pname = kollapse("Expression of ", g, " (",Expr[g],")")
  pdfA4plot_on(pname = pname, rows = 5, 2, title = pname)
  for (i in 1:NrHQworms) {
    w = worms_HQ[i]
    idx_genes.found = rownames(MarkerEx[[i]])
    if(g %in% idx_genes.found){
      expr = as.numeric((MarkerEx_log2[[i]][ g,]))
      expr[is.na(expr)] <- -1
      names(expr) = colnames(MarkerEx_log2[[i]])
      mx = if (max(expr) >= 5) max(expr) else 5
      ylm =  ceiling(c(0,mx))
      barplot(expr, col = wormcol[i], main = w, las=2, cex.names = .5, ylab = "log2(#mRNA)", ylim = ylm)
    }

  } # i
  pdfA4plot_off()
} # j

# Do ------------------------

Pozz = matrix.fromNames(rowname_vec = MarkerGenes, colname_vec = worms_HQ, fill = NaN)

for (j in 1:l(MarkerGenes)) {
  g = MarkerGenes[j]
  # gg = names(g)
  for (i in 1:NrHQworms) {
    w = worms_HQ[i]
    idx_genes.found = rownames(MarkerEx[[i]])
    if(g %in% idx_genes.found){

      expr = as.numeric((MarkerEx_log2[[i]][ g,]))
      poz = which.max(expr)
      Pozz[j, i]= if (!l(poz)) NaN else poz
    } #if
  } # i
} # j
Pozz

PositionalVarianceInMarkerEpression =cbind("cv" = apply(Pozz, 1, cv), "var" = apply(Pozz, 1, var))
md.tableWriter.DF.w.dimnames(PositionalVarianceInMarkerEpression)



# Correlation ------------------------
llprint("### Find other markers by correlation")

e2_06=sfData[["e2_06"]]
egl20 = as.numeric(e2_06[ clone2gene("egl-20"),])
barplot(egl20)
Corr_W_egl20 = as.named.vector(cor(x = t(e2_06), y = egl20))
Corr_W_egl20

CorrGenes = which_names(Corr_W_egl20 >= .9)
llprint("We find", l(CorrGenes), "genes correlating with egl-20")
index =which(rowMax(e2_06[CorrGenes,]) >= .5*max(egl20))
Hits = CorrGenes[index]
llprint("These are: ", Hits)
llprint("Or: ", (Hits))

yy=e2_06[CorrGenes[index],]
xx= as.list.df.by.row(yy)
lapply(xx,l)
plot.multi.dens(xx)

barplot(as.numeric(yy[1,]))

# Do ------------------------
llprint("## Q")




wStarts = colMeans(Pozz[c("flp-11-pharynx", "myo-1-pharynx"), ])
wEnds = Pozz[ "egl-20-tail", ]
realLength = wEnds - wStarts
WormLength

percentage_formatter(cv(realLength))
percentage_formatter(cv(WormLength))

