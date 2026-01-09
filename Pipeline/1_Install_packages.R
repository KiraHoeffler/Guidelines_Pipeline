
################################################################################
# SETUP
################################################################################

# adapt your array type to one of "EPICv1", "EPICv2", "450K" before running
array.type = "EPICv2"


# using BiocManager and install_github might not work if you do not have access to the internet on your secure server
# if you get error messages, see our guide for detailed instructions on how to install packages manually
#
# please use these functions if you have to manually install packages
#     package_path <- "Path_to_unpacked_package_folder"
#     install.packages(package_path, repos = NULL, type = "source")

################################################################################
# INSTALL PACKAGES FOR THE QC AND THE ANALYSES
################################################################################


#################### INSTALL CRAN PACKAGES #####################################
# even if you have the cran packages already installed, it is probably good to update them.

needed.cran.packages = c("BiocManager","readxl", "ggplot2", "writexl", "tidyr", "car", "doParallel", "parallel",
                         "stringr", "grid", "gridExtra", "factoextra", "plotly", "htmlwidgets", "qqman", "knitr", "svglite")

if (array.type == "EPICv2"){
  needed.cran.packages = c(needed.cran.packages, "devtools")
}

install.packages(needed.cran.packages)


################### INSTALL BIOCONDUCTOR PACKAGES ##############################

if (array.type == "EPICv1"){
  needed.bioc.packages = c("minfi","wateRmelon","ENmix", "ChAMP", "limma", "sva", "FlowSorted.Blood.EPIC",
                           "IlluminaHumanMethylationEPICmanifest", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}

if (array.type == "EPICv2"){
  needed.bioc.packages = c("minfi","ENmix", "ChAMP", "limma", "sva", "FlowSorted.Blood.EPIC")
}

if (array.type == "450K"){
  needed.bioc.packages = c("minfi","wateRmelon","ENmix", "ChAMP", "limma", "sva", "FlowSorted.Blood.EPIC", 
                           "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19")
}


missing.bioc.packages = needed.bioc.packages[!(needed.bioc.packages %in% installed.packages()[,"Package"])]
print(paste0("missing bioconductor package: ", missing.bioc.packages))

if(length(missing.bioc.packages)>0){
  print(missing.bioc.packages)
  BiocManager::install(missing.bioc.packages)
  
}

############ INSTALL GITHUB PACKAGES (ONLY EPICv2!) ############################


if (array.type == "EPICv2"){
  library(devtools)
  
  needed.github.packages <- c("IlluminaHumanMethylationEPICv2manifest", "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "wateRmelon")
  missing.github.packages = needed.github.packages[!(needed.github.packages %in% installed.packages()[,"Package"])]
  print(paste0("missing github package: ", missing.github.packages))
  
  if ("IlluminaHumanMethylationEPICv2manifest" %in% missing.github.packages){
    install_github("jokergoo/IlluminaHumanMethylationEPICv2manifest")
  }
  if ("IlluminaHumanMethylationEPICv2anno.20a1.hg38" %in% missing.github.packages){
    install_github("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
  }
  if ("wateRmelon" %in% missing.github.packages){
    install_github("schalkwyk/wateRmelon")
  }
}

