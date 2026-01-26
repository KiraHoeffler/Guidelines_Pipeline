
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
if (!dir.exists("output/Figures/Sex_adjusted_contin")) dir.create("output/Figures/Sex_adjusted_contin")
if (!dir.exists("output/Tables/Sex_adjusted_contin")) dir.create("output/Tables/Sex_adjusted_contin")
if (!dir.exists("output/RData/Sex_adjusted_contin")) dir.create("output/RData/Sex_adjusted_contin")


################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

# LOAD SAMPLESHEET & EUROPEAN INFO
samplesheet <- as.data.frame(read_excel("output/Tables/samplesheet_after_QC.xlsx"))

# LOAD M VALUES
load("output/RData/Mvalues_final_autosomes.RData")

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
write_xlsx(samples, path = "output/Tables/Sex_adjusted_contin/Final_samplesheet_Continuous.xlsx")

Mvalues <- Mvalues[, samples$Basename]

# VARIABLES
Age <- scale(as.numeric(samples$Age))[, 1]
smoking <- scale(as.numeric(Mvalues["cg05575921",]))[, 1]
continPhenotype <- as.numeric(samples$pheno_cont)
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
variables_df <- data.frame(continPhenotype = continPhenotype, Sex = Sex, Age = Age, smoking = smoking)
variables_df <- cbind(variables_df, cell_PCs[,1:5], ctrlprobe_PCAscores[,1:15], anc_PCs[,1:10])               


# TEST MULTI-COLINEARITY
temp_variables_df <- variables_df
temp_variables_df$Mvalues <- Mvalues[1, ]

if (nrow(samples) >= 40){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ continPhenotype + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + 
                                    Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                    Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                                    anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 + 
                                    anc_PC7 + anc_PC8 + anc_PC9 + anc_PC10, data = temp_variables_df)))
  
}

if (nrow(samples) >= 25 & nrow(samples) < 40){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ continPhenotype + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + 
                                    Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                    anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5, data = temp_variables_df)))
  
}

if (nrow(samples) < 25){
  vif_obj <- as.data.frame(vif(lm(Mvalues ~ continPhenotype + Sex + Age + smoking + 
                                    cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +
                                    anc_PC1 + anc_PC2, data = temp_variables_df)))
  
}

vif_obj$Variable <- rownames(vif_obj)
colnames(vif_obj) <- c("vif", "variable")
write_xlsx(vif_obj, path = "output/Tables/Sex_adjusted_contin/VIF_continPhenotype_Sex_adjusted.xlsx")
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
# OUTLIER BOXPLOT
################################################################################

phenotype_boxplot <- ggplot(samples, aes(x="", y = pheno_cont)) + 
  geom_boxplot(fill = "grey") +
  labs(x=NULL, y = "continuous phenotype") +
  th
ggsave(plot = phenotype_boxplot,filename = "output/Figures/Sex_adjusted_contin/Boxplot_phenotype_continuous.pdf", width = 2, height = 3.5, dpi = 300)


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
model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + anc_PC1")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_2CtrlPCs_1AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[1,2] <- calculate_lambda(DMPs$P.Value)



### 5 Ctrl PCs, 1 Anc PCs ###
model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2  + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_5CtrlPCs_1AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[2,2] <- calculate_lambda(DMPs$P.Value)


if (nrow(samples) > 25){
  ### 10 Ctrl PCs, 1 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +  Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_10CtrlPCs_1AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[3,2] <- calculate_lambda(DMPs$P.Value)
}

if (nrow(samples) >= 35){
  ### 15 Ctrl PCs, 1 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 + anc_PC1")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_15CtrlPCs_1AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[4,2] <- calculate_lambda(DMPs$P.Value)
}

#########################################################################################################


### 2 Ctrl PCs, 2 Anc PCs ###
model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + anc_PC1 + anc_PC2")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[5,2] <- calculate_lambda(DMPs$P.Value)



### 5 Ctrl PCs, 2 Anc PCs ###
model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1 + anc_PC2")
model_formula <- update(model_formula, paste("~ . -", high_vif))

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_5CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[6,2] <- calculate_lambda(DMPs$P.Value)




if (nrow(samples) > 25){
  ### 10 Ctrl PCs, 2 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1 + anc_PC2")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_10CtrlPCs_2AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[7,2] <- calculate_lambda(DMPs$P.Value)
}


if (nrow(samples) >= 35){
  ### 15 Ctrl PCs, 2 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 + anc_PC1 + anc_PC2")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_15CtrlPCs_2AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[8,2] <- calculate_lambda(DMPs$P.Value)
}


#########################################################################################################




if (nrow(samples) >= 25){
  ### 2 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + Ctrl_PC1 + Ctrl_PC2 + cell_PC3 + cell_PC4 + cell_PC5 +  anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_2CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[9,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 5 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_5CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[10,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 10 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1  + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_10CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[11,2] <- calculate_lambda(DMPs$P.Value)
  
}

if (nrow(samples) >= 35){
  ### 15 Ctrl PCs, 5 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_15CtrlPCs_5AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[12,2] <- calculate_lambda(DMPs$P.Value)
}


#########################################################################################################





if (nrow(samples) >= 35){
  ### 2 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + Ctrl_PC1 + Ctrl_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 +anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_2CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_2CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[13,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 5 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 +anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_5CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_5CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[14,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 10 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + anc_PC1  + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5  + anc_PC6 + anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_10CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_10CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[15,2] <- calculate_lambda(DMPs$P.Value)
  
  
  
  ### 15 Ctrl PCs, 10 Anc PCs ###
  model_formula <- as.formula(" ~ continPhenotype + Sex + Age + smoking + cell_PC1 + cell_PC2 + cell_PC3 + cell_PC4 + cell_PC5 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 + Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +  anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5  + anc_PC6 +anc_PC7 +anc_PC8 +anc_PC9 +anc_PC10")
  model_formula <- update(model_formula, paste("~ . -", high_vif))
  
  design <- model.matrix(model_formula, data = variables_df)
  fit <- lmFit(Mvalues, design)
  fit <- eBayes(fit)
  
  DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
  DMPs <- DMPs[order(rownames(DMPs)), ]
  DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92
  
  # EXPORT
  #save(DMPs, file = "output/RData/Sex_adjusted_contin/DMPs_continPhenotype_Sex_adjusted_15CtrlPCs_10AncPCs_limma.RData")
  
  pdf(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
  qq(DMPs$P.Value)
  dev.off()
  
  png(file = "output/Figures/Sex_adjusted_contin/QQplot_continPhenotype_Sex_adjusted_15CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
  qq(DMPs$P.Value)
  dev.off()
  
  lambda_table[16,2] <- calculate_lambda(DMPs$P.Value)
  
}


# EXPORT LAMBDA TABLE
write_xlsx(lambda_table, path = "output/Tables/Sex_adjusted_contin/Lambda_continPhenotype_Sex_adjusted_limma.xlsx")




rm(list = ls());gc()
