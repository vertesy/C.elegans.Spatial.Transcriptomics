######################################################################
# 05.z.AnatomicalSections.Length.Comparison.Old.R
######################################################################
# source ('~/Github_repos/CeleTomo/05.z.AnatomicalSections.Length.Comparison.Old.R')

try(dev.off(), silent = T)

# Functions ------------------------

# Setup ------------------------

OutDir = p0(OutDir_orig,"/AnatomicalSections.Length")
setup_MarkdownReports(OutDir = OutDir, scriptname = "05.z.AnatomicalSections.Length.Comparison.Old.R", title = "Statistics on manually curated anatomical sections", append = F)


# Metadata ------------------------
AnatomicalSections = read.simple.tsv("~/Google_Drive/Avano/CT/Metadata/Anatomy_Manual/AnatomicalSections.tsv")
attach_w_rownames(AnatomicalSections)

# Do ------------------------
head	 = c("worm_start", "head_end")
male1	 = c("male1_start", "male1_end")
male2	 = c("male2_start", "male2_end")
male3	 = c("male3_start", "male3_end")
tail	 = c("tail_start", "worm_end")
midbody	 = c("head_end", "tail_start")

Length_head	 = head_end	- worm_start
Length_male1	 = male1_end	- male1_start
Length_male2	 = male2_end	- male2_start
Length_male3	 = male3_end	- male3_start
Length_tail	 = worm_end	- tail_start
Length_midbody	 = tail_start	- head_end

ls_Length_head 	= split(Length_head	, f= sex[worms_HQ])
ls_Length_male1 	= split(Length_male1	, f= sex[worms_HQ])
ls_Length_male2 	= split(Length_male2	, f= sex[worms_HQ])
ls_Length_male3 	= split(Length_male3	, f= sex[worms_HQ])
ls_Length_tail 	= split(Length_tail	, f= sex[worms_HQ])
ls_Length_midbody 	= split(Length_midbody	, f= sex[worms_HQ])

names(ls_Length_head) = names(ls_Length_male1) = names(ls_Length_male2) = names(ls_Length_male3) = names(ls_Length_tail) = names(ls_Length_midbody) = c("glp1", "Male", "Herma")

wstripchart_list( ls_Length_head, bg = SexColors)
wstripchart_list( ls_Length_male1, bg = SexColors)
wstripchart_list( ls_Length_male2, bg = SexColors)
wstripchart_list( ls_Length_male3, bg = SexColors)
wstripchart_list( ls_Length_tail, bg = SexColors)
wstripchart_list( ls_Length_midbody, bg = SexColors)


OutDir = OutDir_orig