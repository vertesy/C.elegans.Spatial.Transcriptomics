######################################################################
# 01.QC_CeleTomo.R
######################################################################
# source ('~/Github_repos/CeleTomo/01.QC_CeleTomo.R')
require(MarkdownReportsDev)

# Setup ------------------------
NrOfUsefulSections = NULL
WormSections = HQ_Slices = ls.GeneCounts.PerSection.Raw = list.fromNames(worms_HQ)

# Metadata ------------------------
# Parameters ------------------------


# Plot Transcripts ------------------------------------------------------------------------------------------------
NrOfUsefulSections = NULL
ylb =  subscript_in_plots(quantity = "Transcripts per section", subscr = 10)
TranscriptCountsRaw = GeneCountsRaw = list.fromNames(Worms)

try(dev.off(), silent = T)
pdfA4plot_on(pname = "TranscriptCountsRaw_log10", rows = 5,cols = 2)
i=13
for (i in 1:NrExp) {
  w = Worms[i]
  TranscriptCountsRaw[[w]] = colSums(ExpData[[w]])
  TranscriptCountsRaw_log10 = log10(TranscriptCountsRaw[[w]]+1)
  x = filter_HP(TranscriptCountsRaw[[w]], SectionMinTrCount, passequal = T, prepend = paste0("In ", w, ", "))
  NrOfUsefulSections[i]  = sum(x, na.rm = T)

  pname = paste(w, isHQ[i])
  subb =  kollapse(NrOfUsefulSections[i], " sections are above: ", SectionMinTrCount, " transcripts.")
  if (HQ[i]) {
    wbarplot_CT(TranscriptCountsRaw_log10, hline = log10(SectionMinTrCount), ylimits = c(0,6), main = pname, sub = subb, ylab=ylb)
  } else {
    wbarplot_CT(TranscriptCountsRaw_log10, hline=log10(SectionMinTrCount), filtercol = 0,ylimits = c(0,6), main = pname, sub = subb, ylab=ylb)
  }

}
pdfA4plot_off()
names(NrOfUsefulSections) = Worms

# Plot Transcripts LINEAR------------------------------------------------------------------------------------------------
NrOfUsefulSections = NULL
ylb =  "Transcripts per section"
try(dev.off(), silent = T)
pdfA4plot_on(pname = "TranscriptCountsRaw_LINEAR", rows = 5,cols = 2)
for (i in 1:NrExp) {
  w = Worms[i]
  x = filter_HP(TranscriptCountsRaw[[w]], SectionMinTrCount, passequal = T, prepend = paste0("In ", w, ", "))

  pname = paste(w, isHQ[i])
  subb =  kollapse(sum(x, na.rm = T), " sections are above: ", SectionMinTrCount, " transcripts.")
  if (HQ[i]) {
    wbarplot_CT(TranscriptCountsRaw[[w]], col = (x+2), main = pname, sub = subb, ylab=ylb) # , ylimits = c(0,2^17)
  } else {
    wbarplot_CT(TranscriptCountsRaw[[w]], main = pname, sub = subb, ylab=ylb)
  }
}
pdfA4plot_off()

# Plot gene Counts ------------------------------------------------------------------------------------------------


Boundaries = matrix(data = NA, nrow = NrExp, ncol = 2)
rownames(Boundaries) = Worms; colnames(Boundaries) = c("Start", "End")
WormLength = 1:NrExp; names(WormLength) = Worms

ylb = subscript_in_plots(quantity = paste("Genes with >", MinExpr4WormLength, " transcripts"), subscr = 10)

SeqDepth = unlapply(ExpData, sum, na.rm=T)
SeqDepthNormFactor = (1/ (SeqDepth /1e7))
wbarplot(SeqDepthNormFactor, col = WormCol_Sex)
wbarplot(SeqDepth, col = WormCol_Sex)

NrOfUsefulSections = vec.fromNames(Worms)
try(dev.off(), silent = T)
pdfA4plot_on(pname = kollapse("GeneCountsRaw_log10_minExpr_",MinExpr4WormLength ), rows = 5,cols = 2)
i=4
for (i in 1:NrExp) {
  w = Worms[i]; print(w)

  GeneCountsRaw[[w]] = colSums((SeqDepthNormFactor[w] * ExpData[[w]]) > MinExpr4WormLength)

  SlicesAbove = filter_HP(GeneCountsRaw[[w]], SectionMinGeneCount, passequal = T, prepend = paste0("In ", w, ", "))
  pass = which(SlicesAbove)

  pass_no_lonely =  which_names(diff(pass)==1) # filter out single HQ slices not followed by another HQ slice
  Bound= c(pass_no_lonely[1], pass_no_lonely[l(pass_no_lonely)])
  Bound= which( names(SlicesAbove) %in% Bound) # get numeric position
  Bound[1] =Bound[1]-1 # correct for the fact diff() gives back the next element

  ccc  = rep(2,l(SlicesAbove)); ccc[Bound[1]:Bound[2]] = 3
  ccc[!SlicesAbove]= 2 # intermediate no-pass sections

  WormSections[[w]] = SlicesAbove[Bound[1]:Bound[2]] # SlicesAbove is a boolean, TRUE if 50+ genes per section
  HQ_Slices[[w]]=  which_names(WormSections[[w]])

  Boundaries[w, ] =c(Bound[1],Bound[2])
  WormLength[w] = Bound[2] - Bound[1] +1

  GeneCountsRaw_log10 = log10(GeneCountsRaw[[w]]+1)
  pname = paste(w, isHQ[i])
  subb = kollapse(WormLength[w], " sections are above ",SectionMinGeneCount," genes between: ", Boundaries[i,1], " and ", Boundaries[i,2])
  if (HQ[i]) {
    wbarplot_CT(GeneCountsRaw_log10, hline = log10(SectionMinGeneCount), ylimits = c(0,4), main = pname, sub = subb, ylab=ylb, filtercol = 0, col = ccc)
  } else {
    wbarplot_CT(GeneCountsRaw_log10, hline=log10(SectionMinGeneCount), filtercol = 0,ylimits = c(0,4), main = pname, sub = subb, ylab=ylb)
  }
  NrOfUsefulSections[i]  = sum(SlicesAbove)
}
pdfA4plot_off()



# Plot gene Counts LINEAR ------------------------------------------------------------------------------------------------
ylb = paste("Number of genes with >", MinExpr4WormLength, " transcripts")
# ls.GeneCountsRaw[[w]] = list.fromNames(names(ExpData))

try(dev.off(), silent = T)
pdfA4plot_on(pname = kollapse("GeneCountsRaw_LINEAR_minExpr_",MinExpr4WormLength ), rows = 5,cols = 2)
i=10
for (i in 1:NrExp) {
  w = Worms[i]

  CCC = (names(GeneCountsRaw[[w]]) %in% HQ_Slices[[w]]) +2
  pname = paste(w, isHQ[i])
  subb = kollapse(WormLength[w], " sections are selected between between: ", Boundaries[i,1], " and ", Boundaries[i,2])
  if (HQ[i]) {
    wbarplot_CT(GeneCountsRaw[[w]], col = CCC, main = pname, sub = subb, ylab=ylb)
  } else {
    wbarplot_CT(GeneCountsRaw[[w]], main = pname, sub = subb, ylab=ylb)
  }
  ls.GeneCounts.PerSection.Raw[[w]] = GeneCountsRaw[[w]][HQ_Slices[[w]] ]
}
pdfA4plot_off()


Fig.1D = T
if (Fig.1D) {
  try.dev.off()
  GeneCounts.PerSection.HQ = ls.GeneCounts.PerSection.Raw[worms_HQ]
  wstripchart(GeneCounts.PerSection.HQ, tilted_text = T, colorbyColumn = T, col = SexColors2[sex[names(ls.GeneCounts.PerSection.Raw[worms_HQ])]])

  LS_GC_perSex_pooled = list(
    "asex" = unlist(GeneCounts.PerSection.HQ[asex_HQ]),
    "males" = unlist(GeneCounts.PerSection.HQ[males_HQ]),
    "hermaphrodites" = unlist(GeneCounts.PerSection.HQ[hermaphrodites_HQ])
  )

  Genes.Per.Section.perSex.pooled = cbind(
    "median" = round(unlapply(LS_GC_perSex_pooled, median, na.rm=T)),
    "mean" = round(unlapply(LS_GC_perSex_pooled, mean, na.rm=T))
  )
  md.tableWriter.DF.w.dimnames(Genes.Per.Section.perSex.pooled, WriteOut = T)
} #if



