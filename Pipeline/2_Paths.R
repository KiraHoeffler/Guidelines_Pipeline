

#####################################################################################################
# Update below (no need to run the script, just adapt and save)
#####################################################################################################

# WHICH ARRAY VERSION DID YOU USE ("450K", "EPICv1", "EPICv2")
array.type <- "EPICv2"

# WHICH TISSUE TYPE ("saliva" or "blood")
tissue_type <- "saliva"

# PATH TO WORKING DIRECTORY (HAS TO CONTAIN THE "script" FOLDER):
dir_gen <- "S:/Project/WP-epigenetics/15_Pipeline_GuidelinesPaper/" #working directory

# PATH TO RESOURCES (PROVIDED BY US):
resources = "S:/Project/WP-epigenetics/02_Import/OCD_Pipeline/Resources/"

# PATH TO SAMPLESHEET (details on the samplesheet in our guide)
samplesheet_path <- "S:/Project/WP-epigenetics/15_Pipeline_GuidelinesPaper/samplesheet_PhaseI.xlsx"

# PATH TO THE IDAT FILES (RAW DNA METHYLATION DATA)
#dir_idat <- "S:/Project/WP-epigenetics/01_Raw_data/2023_April/2022-086-ILL_GSAUHB_METUHB_N=1672/2022-086-ILL_GSAUHB_METUHB_N=1672/UHB/"?
dir_idat <- "S:/Project/WP-epigenetics/15_Pipeline_GuidelinesPaper/Raw_test_data/"


# DO YOU HAVE GENOTYPING DATA FOR ALL SAMPLES IN THE SAMPLESHEET ("yes" or "no")
available_genotyping_data <- "no"

if (available_genotyping_data == "yes"){
  # add the path to the quality controlled genotyping data that can be used to calculate PCs
  # please have the data saved in a matrix called genotype_matrix as RData object (samples as columns, SNPs as rows.
  genotyping_path <- "add_complete_path"
}

########## # ADDITIONAL OUTLIERS (OPTIONAL) #################

# if you want to exclude additional samples in the filtering step for whatever reason, please put their Basename here, you can also add a reason, e.g. "genotype_mismatch"
# one reason for each outlier, so the number of outliers and reasons need to be the same

additional_outliers <- c()
reason_exclusion <- c()
