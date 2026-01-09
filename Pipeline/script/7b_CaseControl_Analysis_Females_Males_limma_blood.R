
################################################################################
# SETUP
################################################################################

# LOAD PATHS
source("2_Paths.R")

# LAMBDA FUNCTION
calculate_lambda <- function(pvalues){
  chisq <- qchisq(1-pvalues, 1)
  lambda <- median(chisq)/qchisq(0.5, 1)
  return(lambda)
}
  


################################################################################
# LOAD LIBRARIES & MAKE FOLDERS TO EXPORT
################################################################################

suppressMessages({
  
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  #library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
  
  # LOAD LIBRARIES
  library(readxl)
  library(qqman)
  library(ggplot2)
  library(parallel)
  library(doParallel)
  library(writexl)
  library(car)
  library(sva)
  library(limma)
  library(grid)
  library(gridExtra)
})


# MAKE FOLDERS
if (!dir.exists("output/Figures/Sex_adjusted")) dir.create("output/Figures/Sex_adjusted")
if (!dir.exists("output/Tables/Sex_adjusted")) dir.create("output/Tables/Sex_adjusted")
if (!dir.exists("output/RData/Sex_adjusted")) dir.create("output/RData/Sex_adjusted")


################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

# LOAD SAMPLESHEET & EUROPEAN INFO
samplesheet <- as.data.frame(read_excel("output/Tables/samplesheet_after_QC.xlsx"))

# LOAD M VALUES & BETA VALUES
load("output/RData/Mvalues_final_autosomes.RData")
load("output/RData/beta_final_autosomes.RData")

# CTRL PROBES
load("output/RData/Positive_ctrlprobe_intensities.RData") #ctrl probes

# CELL TYPE PROPORTIONS
load("output/RData/Cell_type_proportions.RData") 

# SNPs FOR ANCESTRY PCs
load("output/RData/comb_SNPs.RData")

# GENOTYPING DATA
if (available_genotyping_data == "yes"){
  load(genotyping_path)
}

# IMPORT THEME FOR GGPLOTS
load("output/RData/theme.RData")


################################################################################
# GENERAL
################################################################################

# SELECT SAMPLES
sample_names <- intersect(samplesheet$Basename, colnames(Mvalues))
samples <- samplesheet[which(samplesheet$Basename %in% sample_names), ]
sample_names <- samples$Basename

print(paste0("Nr of individuals included: ", length(samples$Basename))) #798
write_xlsx(samples, path = "output/Tables/Sex_adjusted/Final_samplesheet_CaseCtrl.xlsx")

Mvalues <- Mvalues[, samples$Basename]
beta <- beta[, samples$Basename]

# VARIABLES
Age <- scale(as.numeric(samples$Age))[, 1]
smoking <- scale(as.numeric(Mvalues["cg05575921",]))[, 1]
CaseCtrl <- factor(samples$Case_Control, levels = c("Ctrl", "Case"))
Sex <- factor(samples$Sex)

# CELL TYPE PCs
cellcounts_df <- cellcounts$prop[, c(1:2, 10)]
cellcounts_df = cellcounts_df[sample_names, ]

cell_PCs = prcomp(cellcounts_df)
cell_PCs = as.data.frame(cell_PCs$x)
colnames(cell_PCs) <- paste0("cell_", colnames(cell_PCs))
cell_PCs <- sapply(cell_PCs, scale)


# CTRL PROBE PCs
ctrl <- ctrl[, sample_names]
pca <- prcomp(na.omit(t(ctrl))) # run pca
ctrlprobe_PCAscores <- as.data.frame(pca$x) #extract PCA scores
colnames(ctrlprobe_PCAscores) <- paste0("Ctrl_", colnames(ctrlprobe_PCAscores))
ctrlprobe_PCAscores <- sapply(ctrlprobe_PCAscores, scale)

#### ANCESTRY PCs ###
# PCA
if (available_genotyping_data == "yes"){
  
  if (!sample_names %in% colnames(genotype_matrix)){
    stop("PIPELINE STOPPED. Not all sample names are available in colnames(genotype_matrix). Please either ensure overlapping sample names (Basename in samplesheet) or use available_genotyping_data == no in 2_Paths.R")
  }
  
  genotype_matrix_red <- genotype_matrix[, sample_names]
  pc <- prcomp(t(genotype_matrix_red))
}

if (available_genotyping_data == "no"){
  comb_SNPs_red <- comb_SNPs[, sample_names]
  pc <- prcomp(t(comb_SNPs_red))
}

anc_PCs <- as.data.frame(pc$x)
colnames(anc_PCs) <- paste0("anc_", colnames(anc_PCs))
anc_PCs <- sapply(anc_PCs, scale)

# MAKE DF WITH VARIABLES FOR MODEL
variables_df <- data.frame(CaseCtrl = CaseCtrl, Sex = Sex, Age = Age, smoking = smoking)
variables_df <- cbind(variables_df, cell_PCs[,1:5], ctrlprobe_PCAscores[,1:15], anc_PCs[,1:10])               


# TEST MULTI-COLINEARITY
temp_variables_df <- variables_df
temp_variables_df$Mvalues <- Mvalues[1, ]

if (nrow(samples) >= 40){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ CaseCtrl + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + 
                                    Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                    Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                                    anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 + 
                                    anc_PC7 + anc_PC8 + anc_PC9 + anc_PC10, data = temp_variables_df)))
  
}

if (nrow(samples) >= 25 & nrow(samples) < 40){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ CaseCtrl + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + 
                                    Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                    anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5, data = temp_variables_df)))
  
}

if (nrow(samples) < 25){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ CaseCtrl + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +
                                    anc_PC1 + anc_PC2, data = temp_variables_df)))
  
}

vif_obj$Variable <- rownames(vif_obj)
colnames(vif_obj) <- c("vif", "variable")
write_xlsx(vif_obj, path = "output/Tables/Sex_adjusted/VIF_CaseCtrl_Sex_adjusted.xlsx")
rm(temp_variables_df)


if (nrow(samples) >= 40){
  high_vif <- vif_obj[10:21, ]
  high_vif <- high_vif$variable[which(high_vif$vif > 5)]
}

if (nrow(samples) >= 25 & nrow(samples) < 40){
  high_vif <- vif_obj[10:16, ]
  high_vif <- high_vif$variable[which(high_vif$vif > 5)]
}

if (nrow(samples) < 25){
  high_vif <- vif_obj[10:11, ]
  high_vif <- high_vif$variable[which(high_vif$vif > 5)]
}

variables_df <- variables_df[, !names(variables_df) %in% c(high_vif)]


################################################################################
# RUN MODEL
################################################################################

# SELECT M VALUES
Mvalues <- Mvalues[row.names(Mvalues) != "cg05575921", ] # remove cpg that is used to adjust for smoking


# TABLE TO SAVE ALL LAMBDA VALUES
lambda_table <- data.frame(Model = c("2Ctrl_1Anc", "5Ctrl_1Anc", "10Ctrl_1Anc", "15Ctrl_1Anc",
                                     "2Ctrl_2Anc", "5Ctrl_2Anc", "10Ctrl_2Anc", "15Ctrl_2Anc",
                                     "2Ctrl_5Anc", "5Ctrl_5Anc", "10Ctrl_5Anc", "15Ctrl_5Anc",
                                     "2Ctrl_10Anc", "5Ctrl_10Anc", "10Ctrl_10Anc", "15Ctrl_10Anc"), lambda = NA)


### 2 Ctrl PCs, 1 Anc PCs ###
model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + anc_PC1")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_2CtrlPCs_1AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[1,2] <- calculate_lambda(DMPs$P.Value)



### 5 Ctrl PCs, 1 Anc PCs ###
model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2  + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_5CtrlPCs_1AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[2,2] <- calculate_lambda(DMPs$P.Value)



### 10 Ctrl PCs, 1 Anc PCs ###
model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +  Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_10CtrlPCs_1AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[3,2] <- calculate_lambda(DMPs$P.Value)



### 15 Ctrl PCs, 1 Anc PCs ###
model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 + anc_PC1")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_15CtrlPCs_1AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[4,2] <- calculate_lambda(DMPs$P.Value)


#########################################################################################################


### 2 Ctrl PCs, 2 Anc PCs ###
model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + anc_PC1 + anc_PC2")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[5,2] <- calculate_lambda(DMPs$P.Value)



### 5 Ctrl PCs, 2 Anc PCs ###
model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1 + anc_PC2")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_5CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[6,2] <- calculate_lambda(DMPs$P.Value)




if (nrow(samples) > 25){
  ### 10 Ctrl PCs, 2 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1 + anc_PC2")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_10CtrlPCs_2AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[7,2] <- calculate_lambda(DMPs$P.Value)
}


if (nrow(samples) >= 35){
  ### 15 Ctrl PCs, 2 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 + anc_PC1 + anc_PC2")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_15CtrlPCs_2AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[8,2] <- calculate_lambda(DMPs$P.Value)
}


#########################################################################################################




if (nrow(samples) >= 25){
  ### 2 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + Ctrl_PC1 + Ctrl_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +  anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_2CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[9,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 5 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_5CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[10,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 10 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1  + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_10CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[11,2] <- calculate_lambda(DMPs$P.Value)
  
}

if (nrow(samples) >= 35){
  ### 15 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_15CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[12,2] <- calculate_lambda(DMPs$P.Value)
}


#########################################################################################################





if (nrow(samples) >= 35){
  ### 2 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + Ctrl_PC1 + Ctrl_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 +anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_2CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_2CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[13,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 5 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 +anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_5CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_5CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[14,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 10 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1  + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5  + anc_PC6 + anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_10CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_10CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[15,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 15 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ CaseCtrl + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +  anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5  + anc_PC6 +anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_15CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted/QQplot_CaseCtrl_Sex_adjusted_15CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[16,2] <- calculate_lambda(DMPs$P.Value)
  
}


# EXPORT LAMBDA TABLE
write_xlsx(lambda_table, path = "output/Tables/Sex_adjusted/Lambda_CaseCtrl_Sex_adjusted_limma.xlsx")




################################################################################
# PCA
################################################################################


# CONVERT CLASSES OF SAMPLESHEET TO FACTORS
samples$AMP_Plate <- as.factor(samples$AMP_Plate)
samples$Sentrix_ID <- as.factor(samples$Sentrix_ID)
samples$Sentrix_Position <- as.factor(samples$Sentrix_Position)
samples$Sex <- as.factor(samples$Sex)
samples$Case_Control <- as.factor(samples$Case_Control)
samples$Ancestry <- as.factor(samples$Ancestry)


#PCA & PLOTS
Mvalues <- Mvalues[, samples$Basename]
PCA_Mvalues <- prcomp(t(Mvalues)) # PCA 
pca_scores <- data.frame(PCA_Mvalues$x[,])

# SAVE PCA SCORES
save(pca_scores, file = "output/RData/PCA_scores_Mvalues_final_beforeExclusion.RData")

#PLOT
pc_plot <- qplot(x=PC3, y=PC4, data=pca_scores) + th + theme(legend.position = "none") 
pc_plot
ggsave("output/Figures/PCA_FINAL_Mvalues_BaselineOnly.pdf", pc_plot, units = "cm", device = "pdf", width = 8, height = 8, dpi = 400)
save(pc_plot, file = "output/Figures/PC.RData")




outliers_pc1 <- rownames(pca_scores[which((pca_scores$PC1 > (mean(pca_scores$PC1) + 4*sd(pca_scores$PC1))) |
                                            (pca_scores$PC1 < (mean(pca_scores$PC1) - 4*sd(pca_scores$PC1)))), ])
outliers_pc2 <- rownames(pca_scores[which((pca_scores$PC2 > (mean(pca_scores$PC2) + 4*sd(pca_scores$PC2))) |
                                            (pca_scores$PC2 < (mean(pca_scores$PC2) - 4*sd(pca_scores$PC2)))), ])
outliers_pc3 <- rownames(pca_scores[which((pca_scores$PC3 > (mean(pca_scores$PC3) + 4*sd(pca_scores$PC3))) |
                                            (pca_scores$PC3 < (mean(pca_scores$PC3) - 4*sd(pca_scores$PC3)))), ])
outliers_pc4 <- rownames(pca_scores[which((pca_scores$PC4 > (mean(pca_scores$PC4) + 4*sd(pca_scores$PC4))) |
                                            (pca_scores$PC4 < (mean(pca_scores$PC4) - 4*sd(pca_scores$PC4)))), ])
outliers <- unique(c(outliers_pc1, outliers_pc2, outliers_pc3, outliers_pc4))
save(outliers, file = "output/RData/Outliers.RData")


# Diagnosis / CaseControl
pc12_Diagnosis <- qplot(x=PC1, y=PC2, data=pca_scores, colour = samples$Case_Control) + th + 
  theme(legend.position = "none")
pc23_Diagnosis <- qplot(x=PC2, y=PC3, data=pca_scores, colour = samples$Case_Control) + th + 
  theme(legend.position = "none")
pc34_Diagnosis <- qplot(x=PC3, y=PC4, data=pca_scores, colour = samples$Case_Control) + th + 
  labs(color = "CaseControl")
y_label <- textGrob("CaseControl", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Diagnosis <- grid.arrange(y_label, pc12_Diagnosis, pc23_Diagnosis, pc34_Diagnosis, nrow = 1, widths = c(1, 10, 10, 18))
comb_Diagnosis
save(pc12_Diagnosis,pc23_Diagnosis,pc34_Diagnosis,y_label,comb_Diagnosis, file = "output/Figures/PC_Diagnosis.RData")

#Sex
pc12_Sex <- qplot(x=PC1, y=PC2, data=pca_scores, colour = samples$Sex) + th + 
  theme(legend.position = "none")
pc23_Sex <- qplot(x=PC2, y=PC3, data=pca_scores, colour = samples$Sex) + th + 
  theme(legend.position = "none")
pc34_Sex <- qplot(x=PC3, y=PC4, data=pca_scores, colour = samples$Sex) + th + 
  labs(color = "Sex")
pc34_Sex
y_label <- textGrob("Sex", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Sex <- grid.arrange(y_label, pc12_Sex, pc23_Sex, pc34_Sex, nrow = 1, widths = c(1, 10, 10, 15))
save(pc12_Sex,pc23_Sex,pc34_Sex,y_label,comb_Sex, file = "output/Figures/PC_Sex.RData")


#Ancestry
pc12_Anc_classif <- qplot(x=PC1, y=PC2, data=pca_scores, colour = samples$Ancestry) + th + 
  theme(legend.position = "none")
pc23_Anc_classif <- qplot(x=PC2, y=PC3, data=pca_scores, colour = samples$Ancestry) + th + 
  theme(legend.position = "none")
pc34_Anc_classif <- qplot(x=PC3, y=PC4, data=pca_scores, colour = samples$Ancestry) + th + 
  labs(color = "Ancestry")
y_label <- textGrob("Ancestry", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Anc_classif <- grid.arrange(y_label, pc12_Anc_classif, pc23_Anc_classif, pc34_Anc_classif, nrow = 1, widths = c(1, 10, 10, 20))
save(pc12_Anc_classif,pc23_Anc_classif,pc34_Anc_classif,y_label,comb_Anc_classif, file = "output/Figures/PC_Ancestry.RData")

#AMP plate
pc12_AMPplate <- qplot(x=PC1, y=PC2, data=pca_scores, colour = samples$AMP_Plate) + th + 
  theme(legend.position = "none")
pc23_AMPplate <- qplot(x=PC2, y=PC3, data=pca_scores, colour = samples$AMP_Plate) + th + 
  theme(legend.position = "none")
pc34_AMPplate <- qplot(x=PC3, y=PC4, data=pca_scores, colour = samples$AMP_Plate) + th + 
  theme(legend.position = "none")
y_label <- textGrob("AMP plate", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_AMPplate <- grid.arrange(y_label, pc12_AMPplate, pc23_AMPplate, pc34_AMPplate, nrow = 1, widths = c(1, 10, 10, 10))
save(pc12_AMPplate,pc23_AMPplate,pc34_AMPplate,y_label,comb_AMPplate, file = "output/Figures/PC_AMPplate.RData")

#Sentrix position
pc12_Sentrix_Position <- qplot(x=PC1, y=PC2, data=pca_scores, colour = samples$Sentrix_Position) + th + 
  theme(legend.position = "none")
pc23_Sentrix_Position <- qplot(x=PC2, y=PC3, data=pca_scores, colour = samples$Sentrix_Position) + th + 
  theme(legend.position = "none")
pc34_Sentrix_Position <- qplot(x=PC3, y=PC4, data=pca_scores, colour = samples$Sentrix_Position) + th + 
  theme(legend.position = "none")
y_label <- textGrob("Sentrix_Position", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_Sentrixpos <- grid.arrange(y_label, pc12_Sentrix_Position, pc23_Sentrix_Position, pc34_Sentrix_Position, nrow = 1, widths = c(1, 10, 10, 10))
save(pc12_Sentrix_Position,pc23_Sentrix_Position,pc34_Sentrix_Position,y_label,comb_Sentrixpos, file = "output/Figures/PC_Sentrixpos.RData")


#Sentrix ID
pc12_Sentrix_ID <- qplot(x=PC1, y=PC2, data=pca_scores, colour = samples$Sentrix_ID) + th + 
  theme(legend.position = "none")
pc23_Sentrix_ID <- qplot(x=PC2, y=PC3, data=pca_scores, colour = samples$Sentrix_ID) + th + 
  theme(legend.position = "none")
pc34_Sentrix_ID <- qplot(x=PC3, y=PC4, data=pca_scores, colour = samples$Sentrix_ID) + th +
  theme(legend.position = "none")
y_label <- textGrob("Sentrix_ID", rot = 90, gp = gpar(fontface = "bold", fontsize = 16))
comb_SentrixID <- grid.arrange(y_label, pc12_Sentrix_ID, pc23_Sentrix_ID, pc34_Sentrix_ID, nrow = 1, widths = c(1, 10, 10, 10))
save(pc12_Sentrix_ID,pc23_Sentrix_ID,pc34_Sentrix_ID,y_label,comb_SentrixID, file = "output/Figures/PC_SentrixID.RData")


#COMBINE & EXPORT
comb_individual <- grid.arrange(comb_Diagnosis, comb_Sex, comb_Anc_classif, comb_AMPplate, comb_Sentrixpos, comb_SentrixID, ncol = 1)
ggsave("output/Figures/PCA_final_Groups_BaselineOnly.pdf", comb_individual, units = "cm", device = "pdf", width = 20, height = 40, dpi = 400)




rm(list = ls());gc()
