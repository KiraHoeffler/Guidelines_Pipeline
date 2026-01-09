
################################################################################
# SETUP
################################################################################

print("1 - IDAT TO RGSET")
print("Setting up and loading packages.")

# LOAD PATHS
source("2_Paths.R")

#beadcount threshold: percentage of samples with beadcounts < 3:
bc_threshold <- 0.05

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



################################################################################
# GET RGSET (BACKGROUND-CORRECTED)
################################################################################

print("Step 1/3: Read in IDAT files and save RGset (takes time)")
print("     Please check potential warning messages carefully.")


# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx"))

if(!("Basename" %in% colnames(samplesheet))){
  samplesheet$Basename = paste(samplesheet$Sentrix_ID,samplesheet$Sentrix_Position,sep="_")
}



##### SET IDAT FILE DIRECTORY #################
setwd(dir_idat)


# GET NAMES OF IDAT FILES
filenames <- sort(unique(sub("_(Grn|Red).*", "", dir(rec=T, pattern = '.*idat'))))

if (length(filenames) < 1){
  stop("PIPELINE STOPPED: No IDAT files found, check if the dir_idat directory is correct.")
}

# MAKE DATAFRAME WITH FILENAMES, SLIDE AND ARRAY
locations = data.frame(Basename = filenames, 
                       Slide = sub(".*/([^_]+)_.*","\\1" , filenames),
                       Array = sub(".*_","",filenames))

# IDENTIFY MISSING FILES (ON SAMPLESHEET BUT NOT IN FOLDER)
missing_files <- setdiff(samplesheet$Basename, gsub(".*/", "", locations$Basename))

if (length(missing_files)) {
  cat(missing_files, "\n")
  cat("WARNING: IDAT files are missing\n")
  cat("The missing samples will be ignored.")
}

# IDENTIFY UNKNOWN FILES (IN FOLDER BUT NOT IN SAMPLESHEET)
unknown_files <- setdiff(gsub(".*/", "", locations$Basename),samplesheet$Basename)
unknown_files <- sapply(unknown_files, function(x) locations$Basename[grep(x, locations$Basename)])

if (length(unknown_files)) {
  cat("WARNING: there are unknown idat files in folder\n")
  cat("these will be ignored:\n")
  cat(unknown_files,"\n")
  cat("please double-check if samplesheet is complete\n")
  filenames = setdiff(filenames,unknown_files)
}

# READ IN IDAT FILES TO CREATE RG SET
suppressMessages({
  RGset = read.metharray(filenames, T, T, T)
})

# SELECT THE RIGHT ANNOTATION FOR THE RGset
RGset@annotation <- c(array = array, annotation = anno_type)


##### SET WORKING DIRECTORY ###############
setwd(dir_gen)

# EXPORT RGset
save(RGset, file = "output/RData/RGset.RData", compress = F)
  
# ILLUMINA BACKGROUND SUBTRACTION
NegControls <- getControlAddress(RGset, controlType = "NEGATIVE")	
RGset <- bgcorrect.illumina(RGset) 

# EXPORT BACKGROUND-CORRECTED RGset
save(RGset, file = "output/RData/RGset_bgcorrected.RData", compress = F)	

# REMOVE WHAT IS NOT NEEDED ANYMORE
rm(locations, filenames, missing_files, NegControls, unknown_files)


################################################################################
# BEADCOUNTS
################################################################################

print("Step 2/3: Identify probes with a lot bead count.")

# GET BEAD COUNTS
bc <- beadcount(RGset)
bc_probe <- rowSums(is.na(bc))/ncol(bc) #NAs represent probes with beadcount < 3

# IDENTIY CpGs WITH LOW BEAD COUNTS
lowbeadcount <- bc_probe[bc_probe > bc_threshold]
CpGs_lowbeadcount = names(bc_probe[bc_probe > bc_threshold])

# EXPORT LOW BEAD VALUES
lowbeadcount_df <- data.frame(CpG = CpGs_lowbeadcount, bead_count = lowbeadcount)
nrow(lowbeadcount_df)
write_xlsx(lowbeadcount_df, "output/Tables/Low_beadcount.xlsx")
save(lowbeadcount_df, file = "output/RData/Low_beadcount.RData")

# REMOVE NOT NECESSARY
rm(bc, lowbeadcount_df, bc_probe, lowbeadcount, CpGs_lowbeadcount)


################################################################################
# CALCULATE DETECTION P VALUES
################################################################################

print("Step 3/3: Save detection p values (takes time).")
dp = minfi::detectionP(RGset)
save(dp, file = "output/RData/detectionPvalue.RData",compress=F)


rm(list = ls());gc()

