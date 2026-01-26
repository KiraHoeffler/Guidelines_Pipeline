
################################################################################
# SETUP
################################################################################

print("3 - RGSET TO INTENSITIES")
print("Setting up and loading libraries.")

# LOAD PATHS
source("2_Paths.R")

# PATH FOR BLOOD CELL PROPORTION REFERENCE
sorted_blood_path <- file.path(resources,"EPICv1_FlowSorted.Blood.EPIC.Robj")

if (array.type == "EPICv2"){
  annofile_path <- file.path(resources,"EPICv2_manifest.csv")
}

################################################################################
# LOAD LIBERIES
################################################################################

suppressMessages({
  # LOAD LIBRARIES
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  library(minfi)
  library(wateRmelon)
  library(readxl)
  library(writexl)
  if (tissue_type == "blood"){
    library(FlowSorted.Blood.EPIC)
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
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)


# LOAD RGSET & ANNOTATE
load("output/RData/RGset_bgcorrected.RData") # background-corrected
RGset@annotation <- c(array = array, annotation = anno_type) # set right annotation

# SAMPLESHEET
samplesheet <- as.data.frame(read_xlsx("output/Tables/Samplesheet.xlsx"))

# ANNOTATION
if (array.type == "EPICv2"){
  annotation <- read.csv(annofile_path, skip = 7)
}

load(file.path(resources,"HandleReplicProbes.RData"))


################################################################################
# SNP BETA VALUES
################################################################################

print("Step 1/6: extract beta values of rs/SNP probes.")

# GET BETA VALUES OF SNP PROBES AND ACCORDING SAMPLE NAME
snp_betas <- getSnpBeta(RGset)

# EXPORT SNP BETAS
save(snp_betas,  file = "output/RData/SNP_betas.RData",compress=F)


################################################################################
# BISULFITE CONVERSION
################################################################################

print("Step 2/6: Confirm successful bisulfite conversion.")

# samples with values under 90 are flagged, samples with values under 80 should be excluded
bisulfite_conv <- as.data.frame(bscon(RGset))
bisulfite_conv$Basename <- rownames(bisulfite_conv)
colnames(bisulfite_conv) <- c("percent_bisulfite_conv", "Basename")

# VISUALISE MEDIANS OF THE METH/UNMETH CHANNELS (cutoff 10.5)
pdf("output/Figures/Bisulfite_histogram.pdf", width = 7, height = 7)
hist(bisulfite_conv$percent_bisulfite_conv)
dev.off()

samplesheet$perc_bisulfiteConv <- NA
samplesheet$flag_bisulfiteConv <- NA
for (i in c(1:nrow(samplesheet))){
  Basename <- samplesheet$Basename[i]
  bis <- bisulfite_conv$percent_bisulfite_conv[which(bisulfite_conv$Basename == Basename)]
  
  if (!is.na(bis)){
    if (bis < 90){
      flag <- "under_90"
    }
    if (bis < 80){
      flag <- "under_80"
    }
    if (bis >= 90){
      flag <- NA
    }
  } else {
    flag <- "bisulfite_NA"
  }
  
  samplesheet$perc_bisulfiteConv[i] <- bis
  samplesheet$flag_bisulfiteConv[i] <- flag
}
rm(bisulfite_conv)



################################################################################
# MINFI INTENSITY QC
################################################################################

print("Step 3/6: Flag samples with low overall intensities.")

# CONVERT RED/GREEN CHANNEL INTO METHYLATION SIGNAL WITHOUT USING NORMALISATION
mset <- preprocessRaw(RGset) 

# GET CHIPWIDE MEDIANS OF THE METH AND UNMETH CHANNELS
qc <- getQC(mset)

# VISUALISE MEDIANS OF THE METH/UNMETH CHANNELS (cutoff 10.5)
pdf("output/Figures/MinfiIntensity_10.5.pdf", width = 7, height = 7)
plotQC(qc)
dev.off()

# MAKE DATAFRAME FROM QC
qc_df <- as.data.frame(qc)
qc_df$Basename <- rownames(qc_df)
rm(qc)

samplesheet$meth_signal <- NA
samplesheet$unmeth_signal <- NA

for (i in c(1:nrow(samplesheet))){
  Basename <- samplesheet$Basename[i]
  samplesheet$meth_signal[i] <- qc_df$mMed[which(qc_df$Basename == Basename)]
  samplesheet$unmeth_signal[i] <- qc_df$uMed[which(qc_df$Basename == Basename)]
}


samplesheet$med_intensity <- (samplesheet$meth_signal + samplesheet$unmeth_signal)/2
samplesheet$low_intensity <- ifelse(samplesheet$med_intensity < 10.5, "under 10.5", NA)


rm(qc_df)
rm(mset)


################################################################################
# INTENSITIES OF NORMAL PROBES
################################################################################

print("Step 4/6: Get intensities of normal probes from RGset.")


# TYPE II PROBES
TypeII.Name <- getProbeInfo(RGset, type = "II")$Name
TypeII.Green <- getGreen(RGset)[getProbeInfo(RGset, type = "II")$AddressA,]
TypeII.Red <- getRed(RGset)[getProbeInfo(RGset, type = "II")$AddressA,]
rownames(TypeII.Red) <- TypeII.Name
colnames(TypeII.Red) <- sampleNames(RGset)
rownames(TypeII.Green) <- TypeII.Name
colnames(TypeII.Green) <- sampleNames(RGset)


# TYPE I PROBES, SPLIT INTO GREEN AND RED CHANNELS
  # green
TypeI.Green.Name <- getProbeInfo(RGset, type = "I-Green")$Name
TypeI.Green.M <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressB,]
rownames(TypeI.Green.M) <- TypeI.Green.Name
colnames(TypeI.Green.M) <- sampleNames(RGset)
TypeI.Green.U <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressA,]
rownames(TypeI.Green.U) <- TypeI.Green.Name
colnames(TypeI.Green.U) <- sampleNames(RGset)

  # red
TypeI.Red.Name <- getProbeInfo(RGset, type = "I-Red")$Name
TypeI.Red.M <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressB,]
rownames(TypeI.Red.M) <- TypeI.Red.Name
colnames(TypeI.Red.M) <- sampleNames(RGset)
TypeI.Red.U <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressA,]
rownames(TypeI.Red.U) <- TypeI.Red.Name
colnames(TypeI.Red.U) <- sampleNames(RGset)

# EXPORT INTENSITIES OF NORMAL PROBES
save(TypeII.Red ,TypeII.Green ,TypeI.Red.M ,TypeI.Red.U ,TypeI.Green.M ,TypeI.Green.U , file="output/RData/intensities.RData",compress=F)

# REMOVE OBJECTS NOT NEEDED
rm(TypeII.Red ,TypeII.Green ,TypeI.Red.M ,TypeI.Red.U ,TypeI.Green.M ,TypeI.Green.U)

################################################################################
# INTENSITIES OF CONTROL PROBES
################################################################################

print("Step 5/6: Get intensities of positive control probes from RGset.")

# ALL CONTROL PROBES
controls=getProbeInfo(RGset, type = "Control")
types=unique(controls$Type)
types=types[types!='NEGATIVE']
ctrl.names=controls[controls$Type %in% types,'ExtendedType']
ctrl.address=controls[controls$Type %in% types,'Address']
ctrl.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(ctrl.names), dimnames = list(ctrl.names, sampleNames(RGset)))
ctrl.Green[ctrl.names,] <- getGreen(RGset)[ctrl.address,]
ctrl.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(ctrl.names), dimnames = list(ctrl.names, sampleNames(RGset)))
ctrl.Red[ctrl.names,] <- getRed(RGset)[ctrl.address,]
rownames(ctrl.Green)=paste(rownames(ctrl.Green), '.grn', sep='')
rownames(ctrl.Red)=paste(rownames(ctrl.Red), '.red', sep='')
ctrl = rbind(ctrl.Green, ctrl.Red)

# EXPORT CONTROL PROBE INTENSITIES
save(ctrl, file = "output/RData/Positive_ctrlprobe_intensities.RData",compress=F)


################################################################################
# BLOOD CELL TYPE PROPORTIONS
################################################################################

if (tissue_type == "blood"){
  print("Step 6/6: Calculate cell type proportions (takes longer for EPICv2).")
  
  load("output/RData/RGset.RData") #load RGset without background correction
  
    # LOAD SORTED BLOOD 
    load(sorted_blood_path)
    
    suppressMessages({
      # ESTIMATE CELL TYPE PROPORTIONS
      if (array.type == "450K"){# Needs update
        cellcounts = estimateCellCounts2(rgSet = RGset, referenceset = "FlowSorted.Blood.EPIC", meanPlot = FALSE)
      } 
      
      if(array.type == "EPICv1") {
        cellcounts = estimateCellCounts2(rgSet = RGset, referenceset = "FlowSorted.Blood.EPIC", meanPlot = FALSE)
        
      }
      if(array.type == "EPICv2") {##Needs update
        MSet <-preprocessNoob(RGset)
        Betas<-getBeta(MSet)
        Betas <- HandleReplicProbes(Betas)
        IDOLOptimizedCpGsBloodv2<- IDOLOptimizedCpGs[which(IDOLOptimizedCpGs%in%rownames(Betas))]
        identical(rownames(IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,]), IDOLOptimizedCpGsBloodv2)
        cellcounts <- projectCellType_CP(
          Betas[IDOLOptimizedCpGsBloodv2, ],
          IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,],
          contrastWBC = NULL, nonnegative = TRUE,
          lessThanOne = FALSE
        )
      }
    })
  
  save(cellcounts, file="output/RData/Cell_type_proportions.RData")
  
  if (array.type == "EPICv2"){
    cellcounts_df <- as.data.frame(cellcounts)
  } else {
    cellcounts_df = as.data.frame(cellcounts[[1]])
  }
  
  cellcounts_df$Basename <- rownames(cellcounts_df)
  samplesheet <- samplesheet[, !names(samplesheet) %in% c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")]
  samplesheet <- merge(samplesheet, cellcounts_df, by = "Basename", all.x = TRUE)
  
} else {
  print("Step 6/6: Calculate cell type proportions -> only blood, saliva will be estimated later")
}


write_xlsx(samplesheet, path = "output/Tables/Samplesheet.xlsx")


###
rm(list = ls());gc()
