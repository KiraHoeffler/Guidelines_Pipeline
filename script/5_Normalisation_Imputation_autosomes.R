
### OUTLIER SELECTION IN PCA AT THE END OF THE SCRIPT ###

################################################################################
# SETUP
################################################################################

print("7 - NORMALISATION & IMPUTATION - AUTOSOMES")
print("Setting up and loading libraries")

# LOAD PATHS
source("2_Paths.R")

if (array.type == "EPICv1"){
  annofile_path <- file.path(resources,"EPICv1_manifest.csv") #downloaded Illumina annotation file
}
if (array.type == "EPICv2"){
  annofile_path <- file.path(resources,"EPICv2_manifest.csv")
  replic_path <- file.path(resources,"HandleReplicProbes.RData")
}
if(array.type == "450K"){
  annofile_path <- file.path(resources,"450k_manifest.csv")
}



################################################################################
# LOAD LIBERIES
################################################################################

suppressMessages({
  
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  #library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
  #library(shiny, lib.loc = "Z:/Bioconductorpackages/319")
  
  # LOAD LIBRARIES
  library(minfi)
  library(limma)
  library(ENmix)
  library(wateRmelon)
  library(ggplot2)
  library(ChAMP)
  library(readxl)
  library(writexl)
  library(grid)
  library(gridExtra)
  if (tissue_type == "saliva"){
    library(EpiDISH)
  }
  
  
  if (array.type == "EPICv2"){
    array <- "IlluminaHumanMethylationEPICv2" #array
    anno_type <- "20a1.hg38" #select type of annotation
    library(IlluminaHumanMethylationEPICv2manifest)
    library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  }
  
  if (array.type == "EPICv1"){
    array <- "IlluminaHumanMethylationEPIC" #array
    anno_type <- "ilm10b4.hg19" #select type of annotation
    library(IlluminaHumanMethylationEPICmanifest)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }
  if (array.type == "450K"){
    array <- "IlluminaHumanMethylation450k" #array
    anno_type <- "ilmn12.hg19" #select type of annotation
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
})


################################################################################
# IMPORT & ADAPT
################################################################################

# SET PATHS
setwd(dir_gen)

# LOAD INTENSITIES & DETECTION P VALUES
load('output/RData/intensities_filtered_1e16.RData')

# ANNOTATion
annotation <- read.csv(annofile_path, skip = 7)

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx"))

# IMPORT THEME FOR GGPLOTS
load("output/RData/theme.RData")

################################################################################
# QUANTILE-NORMALISATION
################################################################################

print("Step 1/5: Quantile normalisation of intensities")

# NORMALISE
TypeII.Green = normalizeQuantiles(TypeII.Green)
TypeII.Red = normalizeQuantiles(TypeII.Red)
TypeI.Green.M = normalizeQuantiles(TypeI.Green.M)
TypeI.Green.U = normalizeQuantiles(TypeI.Green.U)
TypeI.Red.M = normalizeQuantiles(TypeI.Red.M)
TypeI.Red.U = normalizeQuantiles(TypeI.Red.U)
save(TypeII.Green, TypeII.Red, TypeI.Green.M, TypeI.Green.U, TypeI.Red.M, TypeI.Red.U, file = "output/RData/intensities_normalised.RData")

# CALCULATE BETA VALUES & EXPORT
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))

rm(TypeII.Green,TypeII.Red,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U,TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
#save(beta, file='output/RData/beta_QN_1e16.RData',compress=F)

# VISUALISE BETA DISTRIBUTION
pdf("output/Figures/Distributions/beta_distribution_filtered_QN_1e16.pdf", width = 7, height = 5)
densityPlot(beta, main = "Filtered and QN beta values")
dev.off()






################################################################################
# PREPARE FOR BMIQ
################################################################################

print("Step 2/5: BMIQ normalisation of beta values (takes longer to run)")

# CHANGE ORDER BETA
beta <- beta[order(rownames(beta)), ]

# GET SAMPLE AND PROBE IDs
sampleID <- colnames(beta)

if (array.type == "EPICv2"){
  anno_red <- annotation[which(annotation$IlmnID %in% rownames(beta)), c("IlmnID",  "Infinium_Design")]
  rownames(anno_red) <- anno_red$IlmnID
  anno_red <- anno_red[order(rownames(anno_red)), ]
  probe_types <- anno_red$Infinium_Design
}

if (array.type == "EPICv1"){
  anno_red <- annotation[which(annotation$IlmnID %in% rownames(beta)), c("IlmnID",  "Infinium_Design_Type")]
  rownames(anno_red) <- anno_red$IlmnID
  anno_red <- anno_red[order(rownames(anno_red)), ]
  anno_red$Infinium_Design_Type[which(anno_red$Infinium_Design_Type == "I")] <- 1
  anno_red$Infinium_Design_Type[which(anno_red$Infinium_Design_Type == "II")] <- 2
  probe_types <- as.numeric(anno_red$Infinium_Design_Type)
}

if (array.type == "450K"){
  anno_red <- annotation[which(annotation$IlmnID %in% rownames(beta)), c("IlmnID",  "Infinium_Design_Type")]
  rownames(anno_red) <- anno_red$IlmnID
  anno_red <- anno_red[order(rownames(anno_red)), ]
  anno_red$Infinium_Design_Type[which(anno_red$Infinium_Design_Type == "I")] <- 1
  anno_red$Infinium_Design_Type[which(anno_red$Infinium_Design_Type == "II")] <- 2
  probe_types <- as.numeric(anno_red$Infinium_Design_Type)
}


################################################################################
# BMIQ
################################################################################

# CALCULATE BMIQ SAMPLE BY SAMPLE AND ADD THE NORMALISED BETA VALUES TO BMIQ_comb
BMIQ_betas <- matrix(data = NA, nrow = nrow(beta), ncol = ncol(beta))
colnames(BMIQ_betas) <- colnames(beta)
rownames(BMIQ_betas) <- rownames(beta)

suppressMessages({
  for (i in c(1:ncol(beta))){
    sample <- sampleID[i]
    beta_redX <- as.matrix(beta[, i])
    BMIQ_resultX <- BMIQ(beta_redX, design.v = probe_types, nfit = 50000, sampleID = sample, plots = FALSE, pri = FALSE)
    BMIQbetaX <- BMIQ_resultX$nbeta
    BMIQ_betas[, i] <- BMIQbetaX
  }
})


# VISUALISE BETA DISTRIBUTION
pdf("output/Figures/Distributions/beta_distribution_filtered_QN_BMIQ_1e16.pdf", width = 7, height = 5)
densityPlot(BMIQ_betas, main = "Filtered and normalised beta values")
dev.off()




################################################################################
# DISTRIBUTION NA VALUES PER CPG FOR BETA VALUES
################################################################################

print("Step 3/5: impute missing values.")

# COUNT NAs PER CPG
num_na <- as.data.frame(rowSums(is.na(BMIQ_betas)))
colnames(num_na) <- "nr_NA"
write_xlsx(num_na, path = "output/Tables/Nr_NAs_per_CpGs.xlsx")


################################################################################
# IMPUTATION OF BETA VALUES
################################################################################

suppressMessages({
  invisible(capture.output({
    imputed_beta <- champ.impute(beta = BMIQ_betas, pd = samplesheet, method = "KNN")
  }))
})

beta_imp <- imputed_beta$beta

# EXPORT
if (array.type == "EPICv2"){
  save(beta_imp, file = "output/RData/beta_filtered_QN_BMIQ_imputed_1e16.RData", compress=F)
}
if (array.type == "EPICv1" | array.type == "450K"){
  beta <- beta_imp
  save(beta, file = "output/RData/beta_final_autosomes.RData", compress=F)
}



################################################################################
# CALCULATE M VALUES
################################################################################

if (array.type == "EPICv1" | array.type == "450K"){
  # REPLACE ALL ZEROS IN BETA VALUE WITH LOWEST VALUE
  replacement_zero <- min(beta_imp[which(beta_imp > 0 & !is.na(beta_imp))]) 
  beta_temp <- beta_imp
  beta_temp[beta_temp == 0] <- replacement_zero
  
  # CALCULATE M VALUES 
  Mvalues <- B2M(beta_temp)
  rm(beta_temp)
  save(Mvalues, file='output/RData/Mvalues_final_autosomes.RData', compress=F)
}


#################################################################################
# REPLICATED PROBES (ONLY EPICv2)
#################################################################################

print("Extra step EPICv2: Handle replicated probes.")

if (array.type == "EPICv2"){
  load(replic_path)
  
  beta <- as.matrix(HandleReplicProbes(beta_imp))
  
  # SAVE FINAL BETA VALUE DF
  save(beta, file = "output/RData/beta_final_autosomes.RData")
  
  # REPLACE ALL ZEROS IN BETA VALUE WITH LOWEST VALUE
  replacement_zero <- min(beta[which(beta > 0 & !is.na(beta))])
  beta_temp <- beta
  beta_temp[beta_temp == 0] <- replacement_zero
  
  # CALCULATE M VALUES
  Mvalues <- B2M(beta_temp)
  rm(beta_temp)
  
  #EXPORT M VALUES
  save(Mvalues, file='output/RData/Mvalues_final_autosomes.RData')
}


################################################################################
# VISUALISE DISTRIBUTION
################################################################################

print("Step 4/5: Visualise final distributions.")

pdf("output/Figures/Distributions/Beta_distribution_final_autosomes.pdf", width = 7, height = 5)
densityPlot(beta, main = "Final beta values autosomes")
dev.off()

# VISUALISE M VALUE DISTRIBUTION
pdf("output/Figures/Distributions/Mvalue_distribution_final_autosomes.pdf", width = 7, height = 5)
densityPlot(Mvalues, main = "Final M values autosomes")
dev.off()


################################################################################
# PCA
################################################################################

print("Step 5/5: Run PCA.")

# FILTER SAMPLESHEET
filtered_samplesheet <- samplesheet[samplesheet$Basename %in% colnames(Mvalues), ] # filter on individuals that were not excluded
filtered_samplesheet <- filtered_samplesheet[order(filtered_samplesheet$Basename, decreasing = FALSE), ] # match order of samplesheet to betas df


# CONVERT CLASSES OF SAMPLESHEET TO FACTORS
filtered_samplesheet$AMP_Plate <- as.factor(filtered_samplesheet$AMP_Plate)
filtered_samplesheet$Sentrix_ID <- as.factor(filtered_samplesheet$Sentrix_ID)
filtered_samplesheet$Sentrix_Position <- as.factor(filtered_samplesheet$Sentrix_Position)
filtered_samplesheet$Sex <- as.factor(filtered_samplesheet$Sex)
filtered_samplesheet$Case_Control <- as.factor(filtered_samplesheet$Case_Control)
filtered_samplesheet$Ancestry <- as.factor(filtered_samplesheet$Ancestry)


#PCA & PLOTS
Mvalues <- Mvalues[, colnames(Mvalues) %in% filtered_samplesheet$Basename]
PCA_Mvalues <- prcomp(t(Mvalues)) # PCA 
pca_scores <- data.frame(PCA_Mvalues$x[,])

# SAVE PCA SCORES
save(pca_scores, file = "output/RData/PCA_scores_Mvalues_final.RData")

#PLOT
pc_plot <- qplot(x=PC1, y=PC2, data=pca_scores, ) + th  + theme(legend.position = "none") 
pc_plot
ggsave("output/Figures/PCA_FINAL_Mvalues.pdf", pc_plot, units = "cm", device = "pdf", width = 8, height = 8, dpi = 400)

# Diagnosis / CaseControl
pc12_Diagnosis <- qplot(x=PC1, y=PC2, data=pca_scores, colour = filtered_samplesheet$Case_Control) + th  + 
  theme(legend.position = "none")
pc23_Diagnosis <- qplot(x=PC2, y=PC3, data=pca_scores, colour = filtered_samplesheet$Case_Control) + th  + 
  theme(legend.position = "none")
pc34_Diagnosis <- qplot(x=PC3, y=PC4, data=pca_scores, colour = filtered_samplesheet$Case_Control) + th  + 
  labs(color = "CaseControl")
y_label <- textGrob("CaseControl", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Diagnosis <- grid.arrange(y_label, pc12_Diagnosis, pc23_Diagnosis, pc34_Diagnosis, nrow = 1, widths = c(1, 10, 10, 18))

#Sex
pc12_Sex <- qplot(x=PC1, y=PC2, data=pca_scores, colour = filtered_samplesheet$Sex) + th  + 
  theme(legend.position = "none")
pc23_Sex <- qplot(x=PC2, y=PC3, data=pca_scores, colour = filtered_samplesheet$Sex) + th  + 
  theme(legend.position = "none")
pc34_Sex <- qplot(x=PC3, y=PC4, data=pca_scores, colour = filtered_samplesheet$Sex) + th  + 
  labs(color = "Sex")
pc34_Sex
y_label <- textGrob("Sex", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Sex <- grid.arrange(y_label, pc12_Sex, pc23_Sex, pc34_Sex, nrow = 1, widths = c(1, 10, 10, 15))

#Ancestry
pc12_Anc_classif <- qplot(x=PC1, y=PC2, data=pca_scores, colour = filtered_samplesheet$Ancestry) + th  + 
  theme(legend.position = "none")
pc23_Anc_classif <- qplot(x=PC2, y=PC3, data=pca_scores, colour = filtered_samplesheet$Ancestry) + th  + 
  theme(legend.position = "none")
pc34_Anc_classif <- qplot(x=PC3, y=PC4, data=pca_scores, colour = filtered_samplesheet$Ancestry) + th  + 
  labs(color = "Ancestry")
y_label <- textGrob("Ancestry", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Anc_classif <- grid.arrange(y_label, pc12_Anc_classif, pc23_Anc_classif, pc34_Anc_classif, nrow = 1, widths = c(1, 10, 10, 20))

#AMP plate
pc12_AMPplate <- qplot(x=PC1, y=PC2, data=pca_scores, colour = filtered_samplesheet$AMP_Plate) + th  + 
  theme(legend.position = "none")
pc23_AMPplate <- qplot(x=PC2, y=PC3, data=pca_scores, colour = filtered_samplesheet$AMP_Plate) + th  + 
  theme(legend.position = "none")
pc34_AMPplate <- qplot(x=PC3, y=PC4, data=pca_scores, colour = filtered_samplesheet$AMP_Plate) + th  + 
  theme(legend.position = "none")
y_label <- textGrob("AMP plate", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_AMPplate <- grid.arrange(y_label, pc12_AMPplate, pc23_AMPplate, pc34_AMPplate, nrow = 1, widths = c(1, 10, 10, 10))

#Sentrix position
pc12_Sentrix_Position <- qplot(x=PC1, y=PC2, data=pca_scores, colour = filtered_samplesheet$Sentrix_Position) + th  + 
  theme(legend.position = "none")
pc23_Sentrix_Position <- qplot(x=PC2, y=PC3, data=pca_scores, colour = filtered_samplesheet$Sentrix_Position) + th  + 
  theme(legend.position = "none")
pc34_Sentrix_Position <- qplot(x=PC3, y=PC4, data=pca_scores, colour = filtered_samplesheet$Sentrix_Position) + th  + 
  theme(legend.position = "none")
y_label <- textGrob("Sentrix_Position", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Sentrixpos <- grid.arrange(y_label, pc12_Sentrix_Position, pc23_Sentrix_Position, pc34_Sentrix_Position, nrow = 1, widths = c(1, 10, 10, 10))

#Sentrix ID
pc12_Sentrix_ID <- qplot(x=PC1, y=PC2, data=pca_scores, colour = filtered_samplesheet$Sentrix_ID) + th  + 
  theme(legend.position = "none")
pc23_Sentrix_ID <- qplot(x=PC2, y=PC3, data=pca_scores, colour = filtered_samplesheet$Sentrix_ID) + th  + 
  theme(legend.position = "none")
pc34_Sentrix_ID <- qplot(x=PC3, y=PC4, data=pca_scores, colour = filtered_samplesheet$Sentrix_ID) + th  +
  theme(legend.position = "none")
y_label <- textGrob("Sentrix_ID", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_SentrixID <- grid.arrange(y_label, pc12_Sentrix_ID, pc23_Sentrix_ID, pc34_Sentrix_ID, nrow = 1, widths = c(1, 10, 10, 10))

#COMBINE & EXPORT
comb_individual <- grid.arrange(comb_Diagnosis, comb_Sex, comb_Anc_classif, comb_AMPplate, comb_Sentrixpos, comb_SentrixID, ncol = 1)
ggsave("output/Figures/PCA_final_Groups.pdf", comb_individual, units = "cm", device = "pdf", width = 20, height = 40, dpi = 400)



# if there is an outlier in the final PCA plot, select a PCA value that separate the outlier from the other samples. 
# Run the code below to identify the outlier:

# ADAPT AND RUN THIS CODE (REMOVE THE # AT THE START OF THE NEXT LINE), CHANGE WHICH PC SEPARATES THE OUTLIER AFTER THE $ SYMBOL AND THE NUMBER
#outlier <- rownames(pca_scores[which(pca_scores$PC2 > 200), ])
#print(outlier)


################################################################################
# CELL TYPE PROPORTIONS SALIVA WITH HEPIDISH (ZHENG, 2018)
################################################################################

if (tissue_type == "saliva"){
  print("Saliva-specific extra step: Estimate cell type proportions.")
  
  # LOAD HEPIDISH CELL TYPE PROPORTIONS
  data("centEpiFibFatIC.m")
  data("centBloodSub.m")
  centEpiFibIC.m_wofat <- centEpiFibFatIC.m[, -3] # remove fat cells (was used for breast tissue)
  
  # CALCULATE CELL TYPE PROPORTIONS
  cell_type_prop_hepi <- hepidish(beta.m = beta, ref1.m = centEpiFibIC.m_wofat, ref2.m = centBloodSub.m, h.CT.idx =  3)
  
  # MAKE OVERVIEW TABLE
  cell_type_prop_hepi <- as.data.frame(cell_type_prop_hepi)
  
  # ADD COLUMN FOR COMBINED IMMUNE CELLS
  cell_type_prop_hepi$comb_ICs <- cell_type_prop_hepi$B + cell_type_prop_hepi$NK + cell_type_prop_hepi$CD4T + cell_type_prop_hepi$CD8T +
    cell_type_prop_hepi$Mono + cell_type_prop_hepi$Neutro + cell_type_prop_hepi$Eosino
  
  # EXPORT CELL TYPE PROPORTIONS
  save(cell_type_prop_hepi, file = "output/RData/Cell_type_proportions.RData")
  
  
  se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))) # define function standard error
  summary_table <- data.frame(
    mean = sapply(cell_type_prop_hepi, mean),
    SD = sapply(cell_type_prop_hepi, sd),
    SE = sapply(cell_type_prop_hepi, se)
  )
  summary_table$Cell_type <- rownames(summary_table)
  
  #EXPORT
  write_xlsx(summary_table, path = "output/Tables/Cell_type_proportions_saliva.xlsx")
  
  samplesheet$Epi <- NA
  samplesheet$Fib <- NA
  samplesheet$B <- NA
  samplesheet$NK <- NA
  samplesheet$CD4T <- NA
  samplesheet$CD8T <- NA
  samplesheet$Mono <- NA
  samplesheet$Neutro <- NA
  samplesheet$Eosino <- NA
  samplesheet$comb_ICs <- NA
  
  for (i in c(1:nrow(cell_type_prop_hepi))){
    Basename <- rownames(cell_type_prop_hepi)[i]
    Epi <- cell_type_prop_hepi$Epi[i]
    Fib <- cell_type_prop_hepi$Fib[i]
    B <- cell_type_prop_hepi$B[i]
    NK <- cell_type_prop_hepi$NK[i]
    CD4T <- cell_type_prop_hepi$CD4T[i]
    CD8T <- cell_type_prop_hepi$CD8T[i]
    Mono <- cell_type_prop_hepi$Mono[i]
    Neutro <- cell_type_prop_hepi$Neutro[i]
    Eosino <- cell_type_prop_hepi$Eosino[i]
    comb_ICs <- cell_type_prop_hepi$comb_ICs[i]
    
    samplesheet$Epi[which(samplesheet$Basename == Basename)] <- Epi
    samplesheet$Fib[which(samplesheet$Basename == Basename)] <- Fib
    samplesheet$B[which(samplesheet$Basename == Basename)] <- B
    samplesheet$NK[which(samplesheet$Basename == Basename)] <- NK
    samplesheet$CD4T[which(samplesheet$Basename == Basename)] <- CD4T
    samplesheet$CD8T[which(samplesheet$Basename == Basename)] <- CD8T
    samplesheet$Mono[which(samplesheet$Basename == Basename)] <- Mono
    samplesheet$Neutro[which(samplesheet$Basename == Basename)] <- Neutro
    samplesheet$Eosino[which(samplesheet$Basename == Basename)] <- Eosino
    samplesheet$comb_ICs[which(samplesheet$Basename == Basename)] <- comb_ICs
  }
  write_xlsx(samplesheet, path = "output/Tables/Samplesheet.xlsx")
}

samplesheet_final <- samplesheet[which(is.na(samplesheet$to_exclude1e16)), ]
write_xlsx(samplesheet_final, path = "output/Tables/Samplesheet_after_QC.xlsx")

rm(list = ls());gc()


