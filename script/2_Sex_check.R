
################################################################################
# SETUP
################################################################################

print("2 - SEX CHECK")
print("Setting up and loading libraries")

# LOAD PATHS
source("2_Paths.R")

# SEX OUTLIER
# After looking at the plot for sex outliers, the values for sex outlier in line
# 92 can be adapted, using thresholds that separate the outliers from the male 
# and female clusters (see guideline)

################################################################################
# LOAD LIBRARIES, GENERAL GGPLOT THEME & DEFINE ARRAY
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

suppressMessages({
  # LOAD LIBRARIES
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  library(minfi)
  library(readxl)
  library(wateRmelon)
  library(writexl)
  library(ggplot2)
  
  # ARRAY-SPECIFIC
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


# LOAD THEME
load("output/RData/theme.RData")

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx"))

# IMPORT RGset
load("output/RData/RGset_bgcorrected.RData")

# MAKE SURE SAMPLES IN RGset AND SAMPLESHEET OVERLAP
samplesheet <- samplesheet[which(samplesheet$Basename %in% colnames(RGset)), ]

################################################################################
# SEX CHECK
################################################################################

print("Identify sex mismatches between reported and predicted sex and identify sex outliers")

# MAP METHYLATION DATA TO GENOME
RGset <- mapToGenome(RGset)

# GET SEX FROM DNA METHYLATION DATA
sex <- getSex(RGset)
save(sex, file="output/RData/PredictedSex.RData")
# estimation of sex is based on the median values of measurements on the X and Y
# chr. If yMed - xMed is less than cutoff, we predict a female, otherwise male.

# CHECK IF MISMATCH BETWEEN PREDICTED AND REPORTED SEX AND EXPORT
sex_reported <- data.frame(Basename = samplesheet$Basename, reportedSex = samplesheet$Sex)
sex_predicted <- data.frame(Basename = rownames(sex), predictedSex = sex$predictedSex)
sex_comparison <- merge(sex_reported, sex_predicted, by.x = "Basename", by.y = "Basename")
sex_comparison$match <- ifelse(sex_comparison$reportedSex == sex_comparison$predictedSex, "Y", "N")
mismatch <- sex_comparison[which(sex_comparison$match == "N"), ]
write_xlsx(mismatch, "output/Tables/Sex_mismatch.xlsx")

# CREATE SEX DF
sex_df <- data.frame(sex)
sex_df$Basename <- rownames(sex_df)


# SEX OUTLIER (TRESHOLDS SET BASED ON PLOT)
sex_outlier <- sex_df[which(sex_df$xMed < 10.6 & sex_df$xMed > 0 & sex_df$yMed < 10 & sex_df$yMed > 0), ]
write_xlsx(sex_outlier, path = "output/Tables/Sex_outlier.xlsx")

# SEX PLOT
sex_df$predictedSex <- ifelse(sex_df$Basename %in% mismatch$Basename, "mismatch", sex_df$predictedSex)
sex_df$predictedSex <- ifelse(sex_df$Basename %in% sex_outlier$Basename, "outlier", sex_df$predictedSex)
sex_df$predictedSex <- factor(sex_df$predictedSex, levels = c("F", "M", "mismatch", "outlier"))

sex_plot <- ggplot(sex_df) +
  geom_point(data = subset(sex_df, predictedSex %in% c("F", "M")),
             aes(x = xMed, y = yMed, color = predictedSex), size = 2, shape = 16) +
  geom_point(data = subset(sex_df, predictedSex %in% c("mismatch")),
             aes(x = xMed, y = yMed, color = predictedSex), size = 2, shape = 16) +
  geom_point(data = subset(sex_df, predictedSex %in% c("outlier")),
             aes(x = xMed, y = yMed, color = predictedSex), size = 2, shape = 16) +
  th +
  scale_color_manual(values = c("mismatch" = "#E119DA", "F" = "#6e1fa3", "M" = "#0abdc7", "outlier" = "darkgrey")) +
  theme(#legend.position = c(0.2,0.85),
        legend.background = element_rect(color = NA, fill = NA),
        legend.box.background = element_blank(),
        legend.position = "top", 
        legend.key = element_blank(),
        legend.title = element_blank(), 
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5, ),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_continuous(breaks = seq(floor(0), ceiling(15), by = 0.1)) +
  scale_y_continuous(breaks = seq(floor(0), ceiling(15), by = 0.5)) +
  labs(x= "Median signals X chr",
       y= "Median signals Y chr",
       color = "predicted sex")
sex_plot
ggsave("output/Figures/Sex_plot_mismatch.pdf", plot = sex_plot, device = "pdf", width =10, height = 7, dpi = 400)




# ADD MISMATCH AND OUTLIER TO SAMPLESHEET 
samplesheet$Sex_mismatch <- NA
samplesheet$Sex_outlier <- NA

for (i in c(1:nrow(samplesheet))){
  Basename <- samplesheet[i, "Basename"][1]
  if (Basename %in% sex_outlier$Basename){
    samplesheet[i, "Sex_outlier"] <- "Sex_outlier"
  }
  if (Basename %in% mismatch$Basename){
    samplesheet[i, "Sex_mismatch"] <- "Sex_mismatch"
  }
}

# TEST FOR SEX ANEUPLODIES
suppressWarnings(
  suppressMessages({
    beta <- getBeta(RGset)
    sextest <- wateRmelon::estimateSex(beta)
    sex_aneuploidies <- sextest[which(sextest$predicted_sex!= "Female" & sextest$predicted_sex!= "Male"),]
    sex_aneuploidies$Basename <- rownames(sex_aneuploidies)
  })
)
print(paste0("Individuals with sex aneuploidies: ", sex_aneuploidies$Basename, " (details in output/Tables/Sex_aneuploidy.xlsx)"))
samplesheet$Sex_outlier[which(samplesheet$Basename %in% sex_aneuploidies$Basename)] <- "Sex_aneuploidy"
write_xlsx(sex_aneuploidies, path = "output/Tables/Sex_aneuploidy.xlsx")


# EXPORT FINAL SAMPLESHEET
write_xlsx(samplesheet, "output/Tables/Samplesheet.xlsx")


###
rm(list = ls());gc()
