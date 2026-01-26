
# LOAD LIBRARY
suppressMessages({
  library(readxl)
})


if (tissue_type == "blood"){
  if (phenotype == "case_control"){
    print("Running sex-adjusted case-control analysis.")
    
    source("script/7b_CaseControl_Analysis_Females_Males_limma_blood.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4.Rmd", 
                      output_file = paste0(dir_gen,"Report4_CaseCtrl_Analysis_sex_adjusted.html"),
                      params = list(dir_gen = dir_gen)
    )
  }
  if (phenotype == "continuous"){
    print("Running sex-adjusted continuous analysis.")
    
    source("script/7d_Continuous_Analysis_Females_Males_limma_blood.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4_Continuous.Rmd", 
                      output_file = paste0(dir_gen,"Report4_Continuous_Analysis_sex_adjusted.html"),
                      params = list(dir_gen = dir_gen)
    )
  }
}

if (tissue_type == "saliva"){
  if (phenotype == "case_control"){
    print("Running sex-adjusted case-control analysis.")
    
    source("script/7c_CaseControl_Analysis_Females_Males_limma_saliva.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4.Rmd", 
                      output_file = paste0(dir_gen,"Report4_CaseCtrl_Analysis_sex_adjusted.html"),
                      params = list(dir_gen = dir_gen)
    )
  }
  if (phenotype == "continuous"){
    print("Running sex-adjusted continuous phenotype analysis.")
    
    source("script/7e_Continuous_Analysis_Females_Males_limma_saliva.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4_Continuous.Rmd", 
                      output_file = paste0(dir_gen,"Report4_Continuous_Analysis_sex_adjusted.html"),
                      params = list(dir_gen = dir_gen)
    )
  }
}
  
  
  