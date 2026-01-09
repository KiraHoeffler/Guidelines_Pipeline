

# PLEASE READ ALL THE WARNING MESSAGES THAT ARE PRINTED WHEN RUNNING THE PIPELINE.

##########################################################################################################################
# SET PATHS
##########################################################################################################################

# THIS CODE AT THE START HAS TO BE RUN EVERY SINGLE TIME YOU RUN SCRIPTS OF THE PIPELINE
# please adapt the complete path to the 2_Paths.R script here

source("S:/Project/WP-epigenetics/15_Pipeline_GuidelinesPaper/2_Paths.R")
setwd(dir_gen)



##########################################################################################################################
# PART ONE - SET-UP
##########################################################################################################################

#RUN THE SETUP AND SOME BASIC TESTS
source("script/0_Setup.R")


##########################################################################################################################
# PART TWO - QC and Filtering pipeline I
##########################################################################################################################

# RUN THE PIPELINE SCRIPTS 1-2
source("script/1_IDAT_to_RGset.R")
source("script/2_Sex_check.R")

# CREATE A REPORT
source("2_Paths.R")
rmarkdown::render("script/Report_Part1.Rmd", 
                  output_file = paste0(dir_gen,"Report1_QC_Sex.html"),
                  params = list(dir_gen = dir_gen))

# IMPORTANT: 
# Check the QC report before continuing (see guideline)
#     -> in particular check if outliers have been identified in the sex plot
#     if they have not been identified correctly, adapt the 2_Sex_check.R script and re-run it (see guideline). 


##########################################################################################################################
# PART THREE - QC and Filtering pipeline II
##########################################################################################################################

# RUN THE PIPELINE SCRIPTS 3-4
source("script/3_RGset_to_Intensities.R")
source("script/4_Filtering_autosomes.R")

# CREATE A REPORT
source("2_Paths.R")
rmarkdown::render("script/Report_Part2.Rmd", 
                  output_file = paste0(dir_gen, "Report2_QC.html"),
                  params = list(dir_gen = dir_gen))


# IMPORTANT:
# Check the QC report before continuing (see guideline).
#      Specifically look at the density plot with the grey background at the 
#      end of the QC report to check that outliers have been removed.
#      if they have not been removed, adapt the respective filtering script and re-run it (see guideline). 



##########################################################################################################################
# PART FOUR - Normalization, imputation
##########################################################################################################################

# RUN THE PIPELINE SCRIPTS 5-6
source("script/5_Normalisation_Imputation_autosomes.R")
source("script/6_Ancestry.R")

# CREATE A REPORT
source("2_Paths.R")
rmarkdown::render("script/Report_Part3.Rmd", 
                  output_file = paste0(dir_gen, "Report3_Normalisation_Imputation_PCA.html"),
                  params = list(dir_gen = dir_gen)
)


# IMPORTANT:
# Check the QC report before continuing (see guideline).
#      Specifically look at the Outlier detection PCA plot to check that outliers have been removed.
#      if they have not been removed, follow the steps described in the guideline. 


##########################################################################################################################
# PART FIVE - Analysis
##########################################################################################################################

# RUN THE SCRIPTS PIPELINE SCRIPTS 7 AND CREATE A REPORT
source("script/7a_Run_analysis_limma.R")
