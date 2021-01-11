# this script will be used to generate 7-way venn to show overlap of proteins 
# among different categories
# make sure to change the working directory and include input file kf_combined_prop.csv

library(dplyr) 
library(venn)

setwd("~/Dropbox/Killifish-Collaboration/Systems_Aggregate_Paper/20201006_Code-Data/")
fpath <- ""
fpathout <- "Overlap"

tissues <- list("Brain", "Gut", "Heart", "Liver", "Muscle", "Skin", "Testis")
tissue_colors <- list("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594", "#3288bd")
tissue_proteins <- vector("list", length(tissues))
names(tissue_proteins) <- tissues

######## SECTION 1: FIND THE OVERLAP IN IDENTIFIED PROTEINS AMONG DIFFERENT TISSUES ########
# start with TL or AGG data
fdata <- "kf_combined_prop.csv"
df_table <- read.csv(fdata)
# find the protein identified in each tissue and assign them to tissue_proteins
# also create the label for them
# make separate plot for TL and AGG
for (st in list("TL", "AGG")){
  for (i in seq(1:length(tissues))){
    tissue <- tissues[[i]]
    tissue_proteins[[i]] <- df_table[(df_table$Tissue == tissue) & (df_table$Sample.Type == st), c("N..furzeri.Protein.Id")]}
  fout = paste(fpathout, paste(st,"_venn7.pdf", sep=""), sep="/")
  # plot the 7-way venn digram
  pdf(file = fout, width=100, height=100, pointsize = 500)
  venn(tissue_proteins, ilab=TRUE, zcolor=tissue_colors, opacity = 0.6)
  dev.off()}
##################################### END OF SECTION 1 #####################################


######## SECTION 2: FIND THE OVERLAP IN PROTEINS WITH INCREASED EXPRESSION OR AGGREGATES  ########
# filter with q-value <=0.05
fdata <- "kf_combined_prop.csv"
df_table <- read.csv(fdata)
cutoff <- 0.05
# find the protein identified in each tissue and assign them to tissue_proteins
# also create the label for them
for (st in list("TL", "AGG")){
  for (agecomp in list("OvY", "TvO")){
    if (grepl("T", agecomp)) n <- 6 else n <-7
    for (i in seq(1:n)){
      tissue <- tissues[[i]]
      tissue_proteins[[i]] <- df_table[ (df_table$Tissue == tissue) & (df_table[[paste(agecomp,"_logFC", sep="")]] >= 0) & (df_table[[paste(agecomp,"_pval", sep="")]] <= cutoff) & (df_table$Sample.Type == st), c("N..furzeri.Protein.Id")]}
      fout = paste(fpathout, paste(c(st, agecomp,"SigPos_venn7.pdf"), collapse = "_"), sep="/")
      # plot the 7-way venn digram
      pdf(file = fout, width=100, height=100, pointsize = 500)
      venn(tissue_proteins[1:n], ilab=TRUE, zcolor=tissue_colors[1:n], opacity = 0.6)
      dev.off()}}
# ##################################### END OF SECTION 2 #####################################

######## SECTION 3: FIND THE OVERLAP IN PROTEINS WITH INCREASED AGGREGATION PROPENSITY  ########
# filter with q-value <=0.05
fdata <- "kf_combined_prop.csv"
df_table <- read.csv(fdata)
cutoff <- 0.05
# find the protein identified in each tissue and assign them to tissue_proteins
# also create the label for them
for (agecomp in list("OvY", "TvO")){
  if (grepl("T", agecomp)) n <- 6 else n <-7 
for (i in seq(1:n)){
  tissue <- tissues[[i]]
  tissue_proteins[[i]] <- df_table[ (df_table$Tissue == tissue) & (df_table[[paste(agecomp,"_prop_logFC", sep="")]] >= 0) & (df_table[[paste(agecomp,"_prop_pval", sep="")]] <= cutoff) & (df_table$Sample.Type == "AGG"), c("N..furzeri.Protein.Id")]}
  fout = paste(fpathout, paste(c("Prop_", agecomp,"_SigPos_venn7.pdf"), collapse = ""), sep="/")
# plot the 7-way venn digram
pdf(file = fout, width=100, height=100, pointsize = 500)
venn(tissue_proteins[1:n], ilab=TRUE, zcolor=tissue_colors[1:n], opacity = 0.6)
dev.off()}
# ##################################### END OF SECTION 3 #####################################