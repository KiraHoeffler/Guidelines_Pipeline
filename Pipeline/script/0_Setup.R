
################################################################################
# SET-UP
################################################################################

print("0 - SET-UP AND INITIAL TESTS")
source("2_Paths.R")

################################################################################
# TESTS
################################################################################

if (!tissue_type %in% c("blood", "saliva")){
  stop("Tissue type in 2_Paths.R needs to be either saliva or blood.")
}

if (!array.type %in% c("450K", "EPICv1", "EPICv2")){
  stop("array.type in 2_Paths.R needs to be either 450K, EPICv1, or EPICv2.")
}

if (!available_genotyping_data %in% c("yes", "no")){
  stop("available_genotyping_data in 2_Paths.R must be either yes or no.")
}


################## Packages are installed ######################################

needed.packages = c("BiocManager","readxl", "ggplot2", "writexl", "tidyr", "car", "doParallel", "parallel",
                    "stringr", "grid", "gridExtra", "factoextra", "plotly", "htmlwidgets", "qqman", "knitr",
                    "minfi","ENmix", "ChAMP", "limma", "sva", "wateRmelon")

if (tissue_type == "blood"){
  needed.packages <- c(needed.packages, "FlowSorted.Blood.EPIC")
}
if (tissue_type == "saliva"){
  needed.packages <- c(needed.packages, "EpiDISH") 
}

if (array.type == "EPICv2"){
  needed.packages = c(needed.packages, "devtools", "IlluminaHumanMethylationEPICv2manifest", "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
}

if (array.type == "EPICv1"){
  needed.packages = c(needed.packages, "IlluminaHumanMethylationEPICmanifest", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}

if (array.type == "450K"){
  needed.packages = c(needed.packages, "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19")
}

missing.packages = needed.packages[!(needed.packages %in% installed.packages()[,"Package"])]

if(length(missing.packages)>0){
  stop(paste0("The following packages are missing: ", missing.packages))
}



################## Check the samplesheet #######################################

samplesheet <- as.data.frame(readxl::read_xlsx(samplesheet_path))

missing_columns = setdiff(c("AMP_Plate","Sentrix_ID","Sentrix_Position","Sex","Age","Ancestry"),colnames(samplesheet))
if(length(missing_columns)>0){
  print("The following columns are missing or spelled incorrectly in the samplesheet:")
  print(missing_columns)
  stop("PIPELINE STOPPED: Please check that the samplesheet contains the columns listed in the Guide document")
}

missing_columns = setdiff(c("Case_Control"),colnames(samplesheet))
if(length(missing_columns)>0){
  print("The Case_Control column is missing or spelled incorrectly in the samplesheet")
}
if(length(setdiff(c("Case","Ctrl"),levels(factor(samplesheet$Case_Control))))>0){
  stop("PIPELINE STOPPED: Please check that the Case_Control column is coded as 'Case' and 'Ctrl'")
}
if(any(is.na(samplesheet$Case_Control))){
  stop("PIPELINE STOPPED: Please check that there is no NA in the Case_Control column")
}



if(length(setdiff(c("F","M"),levels(factor(samplesheet$Sex))))>0){
  stop("PIPELINE STOPPED: Please check that the Sex column is coded as 'M' and 'F'")
}




if(any(is.na(samplesheet$AMP_Plate))){
  stop("PIPELINE STOPPED: Please check that there is no NA in the AMP_Plate column")
}

if(any(is.na(samplesheet$Sentrix_ID))){
  stop("PIPELINE STOPPED: Please check that there is no NA in the Sentrix_ID column")
}

if(any(is.na(samplesheet$Sentrix_Position))){
  stop("PIPELINE STOPPED: Please check that there is no NA in the Sentrix_Position column")
}

if(any(is.na(samplesheet$Sex))){
  stop("PIPELINE STOPPED: Please check that there is no NA in the Sex column")
}

if(any(is.na(samplesheet$Age))){
  stop("PIPELINE STOPPED: Please check that there is no NA in the Age column")
}

if(any(is.na(samplesheet$Ancestry))){
  stop("PIPELINE STOPPED: Please check that there is no NA in the Ancestry column")
}



print("Passed all intitial tests.")


################################################################################
# SET UP FOLDER STRUCTURE FOR THE OUTPUT FILES
################################################################################

# OUTPUT FOLDER STRUCTURE
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/RData")) dir.create("output/RData")
if (!dir.exists("output/Figures")) dir.create("output/Figures")
if (!dir.exists("output/Tables")) dir.create("output/Tables")

print("A new output folder generated in working directory. QC and analysis output will be saved here.")

################################################################################
# DEFINE GENERAL THEME FOR THE FIGURES
################################################################################

suppressMessages({
  #LOAD PACKAGE
  library(ggplot2)
})


# DEFINE THEME
th <- theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(size = 14, color = "black", face = "bold"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
            axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
            axis.title.x = element_text(vjust = -0.5, size = 12, colour = "black"),
            axis.title.y = element_text(vjust = -0.5, size = 12, colour = "black", margin = margin(r = 10)),
            axis.title = element_text(face = "bold"),
            panel.border = element_blank(),
            legend.text = element_text(size = 12),
            legend.title = element_text(size =12),
            legend.key = element_blank())

save(th, file = "output/RData/theme.RData")


################################################################################
# COPY SAMPLESHEET
################################################################################

suppressMessages({
  #LOAD PACKAGE
  library(readxl)
  library(writexl)
  
})

samplesheet <- as.data.frame(read_xlsx(samplesheet_path))
samplesheet$Basename <- paste0(samplesheet$Sentrix_ID, "_", samplesheet$Sentrix_Position)
write_xlsx(samplesheet, path = "output/Tables/Samplesheet.xlsx")

print("Samplesheet copied to output/Tables/. During the QC process, information will be added to the samplesheet.")


rm(list = ls());gc()
