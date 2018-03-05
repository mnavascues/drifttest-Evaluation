#########################################################################################
#
# Script to evaluate the "Fc outlier method" modified from Goldringer and Bataillon 2002
#
# by A. Becheler, R. Vitalis & M. Navascués
# 
#########################################################################################

# This scripts uses SLiM (Messer 2013) for simulating the data
# http://messerlab.org/software/

# Run from command line examples:
# R --no-save --args seed4random  1234                 < main.R 
# R --no-save --args simID        "SelfAdapt_project"  < main.R 
# R --no-save --args Ne_only      "c(T)"               < main.R 
# R --no-save                                          < main.R

#---------------------------------
# FIXED SETTINGS AND DEPENDENCIES
#---------------------------------

# This script uses the following packages:
library(batch)
library(qqman)
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
library(qvalue)

# high penality for scientific notation
# (necessary to input large numbers in SLiM, which does not take scientific notation as input)
options("scipen"=999)

# functions to write/read SLiM input/output
source("src/fun/slim_tools.R")
source("src/fun/manipulate_output.R")
# functions to calculate F statistics and analyse data
source("src/fun/F_stats_tools.R")
# color blind palette
source("src/fun/ColorBlindPalette.R")



#----------------------------------------------------------------------
# PARAMETERS: DEFAULT VALUES + VALUES FROM COMMAND LINE CALL ARGUMENTS
#----------------------------------------------------------------------
source("src/step1_parameter_values.R")

#------------------------
# WRITE HEAD OF LOG FILE
#------------------------
source("src/step2_write_head_log.R")

#-----------
# SIMULATION
#-----------

# write SLiM input files and run SLiM for drift period
source("src/step3_simulation_drift.R")


source("src/step3.1_set_replica_loop.R")

while (replic <= number_of_replicates){

  print(paste(simID,": replicate",replic,"of",number_of_replicates))

  source("src/step3.2_settings_per_replica.R")
  
  # read SLiM output files from drift period
  # write SLiM input files and run SLiM for selection period
  source("src/step4_simulation_selection.R") 

  # read SLiM output files from selection period
  # sample individuals and loci
  source("src/step5_sample_output.R")

  # Analyse data (Fst outlier neutrality test)
  # 1. estimate Fst, Ne (all loci)
  # 2. test for outliers (locus by locus)
  source("src/step6_outlier_test.R")

}

# Calculate FPR and footprint of selection
source("src/step7_PositiveRate.R")






if (!file.exists(paste0("results/",simID,"/populationATequilibrium"))){
  file.copy( from=paste0("bin/",slim_init_drift),
             to=paste0("results/",simID,"/populationATequilibrium"))
}

# clean up files from simulation
if (remove_files){
  files2delete <- c(paste0("bin/",slim_out_drift), 
                    paste0("bin/",slim_out_drift2),
                    paste0("bin/",slim_out_selection), 
                    paste0("bin/",slim_init_selection),
                    paste0("bin/",slim_init_drift),
                    paste0("bin/",slim_init_drift2),
                    "populationATequilibrium",
                    slim_in_drift,
                    slim_in_drift2,
                    slim_in_selection,
                    slim_log_drift,
                    slim_log_drift2,
                    slim_log_selection)
  for (item in files2delete){
    system( paste("rm",item) )
  }
}


# test
# R --no-save --args sigma 0.0 sel_coef 0.0 simID "test" selection_mode "SV" selection_period_duration 20 seed4random 4444 replic 0 N 500 chr_num 2 genome_length 5e5 sample_size 50 sample_size_loci 5 remove_files "c(F)" < Main.R
# R --no-save --args sigma 0.0 sel_coef 0.0 simID "test" selection_mode "SV" selection_period_duration 20 seed4random 4444 replic 0 N 500 chr_num 2 genome_length 5e8 sample_size 50 sample_size_loci 1000 remove_files "c(F)" < Main.R
