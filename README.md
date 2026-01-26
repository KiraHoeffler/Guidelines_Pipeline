# Guidelines_Pipeline

Here, we provide an example MWAS pipeline from raw IDAT files to a case-control analysis. 

Please read the Guide carefully before you start using this pipeline. 

Download the Resources folder from OneDrive: https://universityofbergen-my.sharepoint.com/:f:/g/personal/kira_hoeffler_uib_no/IgDy_tQ8RcM2SIRQi6KqKwisAbLa_75LqX5Kmpj2Zn37WKY?e=0Hj13G

Alternatively, download the folder here from Github and add the following objects/files manually:
- run the provided download_FlowSorted.Blood.EPIC_RGset script, save the RGset in the Resources folder
1. 450K: download HumanMethylation450 v1.2 Manifest File (CSV Format) from https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html, rename the csv file to 450K_manifest.csv and move it to the resources folder.
2. EPICv1: download Infinium MethylationEPIC v1.0 B5 Manifest File (CSV Format) from https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html, unpack the zip folder, rename the csv file to EPICv1_manifest.csv and move it to the resources folder.
3. EPICv2: download Infinium MethylationEPIC v2.0 Product Files from https://support.illumina.com/downloads/infinium-methylationepic-v2-0-product-files.html, unpack the folder, rename the EPIC-8v2-0_A2.csv file to EPICv2_manifest.csv and move it to the resources folder.

If only a subset of the QCed data is used, please add a filtering step in the analysis scripts and run an additional PCA on the filtered M values to check for outliers. 

Contact: kira.hoeffler@uib.no
