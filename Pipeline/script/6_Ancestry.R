

################################################################################
# SETUP
################################################################################

print("6 - ANCESTRY")
print("Setting up and loading libraries")

# LOAD PATHS
source("2_Paths.R")


# ################################################################################
# # LOAD
# ################################################################################
# 
# suppressMessages({
#   library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
#   library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
#   library(shiny)
# })
#   

################################################################################
# IMPORT & ADAPT
################################################################################

print("Step 1/3: Import needed objects")

# SET PATHS
setwd(dir_gen)

#EpiAnceR (functions & folder)
source(paste0(resources, "Ancestry_PCs/Ancestry_functions.R")) 
path_SNP_folder <- paste0(resources, "Ancestry_PCs/")

#RGset
load("output/RData/RGset_bgcorrected.RData") #RGset

#Samplesheet
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx")) #samplesheet
samplesheet <- samplesheet[which(is.na(samplesheet$to_exclude1e16)), ]

# Detection p values
dp_folder <- "output/RData/detectionPvalue.RData"

#theme
load("output/RData/theme.RData")


################################################################################
# EXTRACT ANCESTRY INFORMATION
################################################################################

print("Step 2/3: Extract ancestry information")

if (tissue_type == "saliva"){
  comb_SNPs <- ancestry_info(background_corr_RGset = RGset, samplesheet = samplesheet, array_type = array.type, 
                             path_SNP_folder = path_SNP_folder, cell_cols = c("Epi", "Fib", "comb_ICs"), dp = dp_folder)
}

if (tissue_type == "blood"){
  comb_SNPs <- ancestry_info(background_corr_RGset = RGset, samplesheet = samplesheet, array_type = array.type, 
                             path_SNP_folder = path_SNP_folder, cell_cols = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), dp = dp_folder)
}

save(comb_SNPs, file = "output/RData/comb_SNPs.RData")


################################################################################
# RUN PCA
################################################################################

print("Step 3/3: Run PCA and visualize.")

### PCA ###

pc <- prcomp(t(comb_SNPs))
top10pc <- as.data.frame(pc$x[,1:10])
top10pc$Basename <- rownames(top10pc)

pc_plot <- qplot(x=PC1, y=PC2, data = top10pc, color = samplesheet$Ancestry, ) + th +
  theme(legend.background = element_rect(color = NA, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())
ggsave("output/Figures/PCA_Ancestry.pdf", pc_plot, units = "cm", device = "pdf", width = 13, height = 8, dpi = 400)


rm(list = ls());gc()


