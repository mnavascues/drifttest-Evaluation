# DEFAULT VALUES

# Set seed for random number generation
seed4random <- 1623656
# Simulation ID for file identification
simID <- "test"
# selfing rate
sigma <- 0.95 
# mutation rate per bp
u <- 1e-8
# recombination rate per bp
r <- 1e-8
# genome size
genome_length <- 5e8
# number of chromosomes
chr_num <- 2 # DO NOT CHANGE VALUE
# effective population size (Ne) (in number of diploid individuals)
N <- 500
# selection coeficient
sel_coef <- 0.5 
# dominance cofficient
dominance_coef <- 1
# length of pure drift period (number of times the population size)
number_of_times       <- 20 # see below
# Adaptation mode: "NM"=new mutation; "SV"=standing variation
selection_mode <- "SV"
initial_frequency <- 0.8
fixed_advantageous_allele <- "derived"
# number of generations between samples
selection_period_duration <- 25
# sample size
sample_size      <- c(50,50) #c(40,60)   # number of diploid individuals sampled
sample_size_loci <- 10000    # number of loci sampled for demographic inference)
#threshold for mimimum allele frequency
MAF_threshold   <- 0.05
# replicate ID
#replic  <- 0
number_of_replicates <- 2
drift_period_between_replicates <- 20
# if TRUE: estimate Ne from simulations only, do not perform neutrality tests
Ne_only <- F
# Do not output progress messages
quiet   <- F
# Do not output whole population
no_whole_pop_out <- T
# Remove simulation files
remove_files     <- F
# Number of threads for drifttest
num_of_threads <- 16

palette(cbPalette2)

# VALUES FROM COMMAND LINE CALL ARGUMENTS
# gets parameter and setting values from command line (using package 'batch')
parseCommandArgs()

# set compound parameter values that may have change through command line call
#------------------------------------------------------------------------------

# set random seed
set.seed(seed4random)

if (length(sample_size)==1) sample_size <- c(sample_size,sample_size)
if (length(sample_size)>2)  sample_size <- sample_size[1:2]

# length of pure drift period
if (file.exists(paste0("results/",simID,"/populationATequilibrium"))){
  drift_period_duration <- drift_period_between_replicates
}else{
  drift_period_duration <- number_of_times*N
}

if (selection_mode=="NM"){
  advantageous_allele <- "derived"
}else if (selection_mode!="SV"){
  selection_mode<-"SV"
  warning("Adaptation mode undefined, using standing variation") 
}

# FILE NAMING

# setwd("../")
# directory to file output
working_dir <- paste0("results/",simID) 
system( paste("mkdir",working_dir) )
replicate_dir <- paste0("results/",simID,"/replicates") 
system( paste("mkdir",replicate_dir) )
# name for log file for each replicate
log_file          <- paste0(working_dir,"/",simID,"_log.txt")
F_stats_file      <- paste0(working_dir,"/",simID,"_F_stats.RData") 
results_file      <- paste0(working_dir,"/",simID,"_results.RData") 
parameter_file    <- paste0(working_dir,"/",simID,"_params.RData") 


#log_file          <- paste0(working_dir,"/",simID,"_",replic,"_log.txt")
#whole_pop_file    <- paste0(working_dir,"/",simID,"_",replic,"_whole_pop.RData") 
#data_file         <- paste0(working_dir,"/",simID,"_",replic,"_data.RData") 
#F_stats_file      <- paste0(working_dir,"/",simID,"_",replic,"_F_stats.RData") 
#results_file      <- paste0(working_dir,"/",simID,"_",replic,"_results.RData") 
#parameter_file    <- paste0(working_dir,"/",simID,"_",replic,"_params.RData") 
#traject_file      <- paste0(working_dir,"/",simID,"_",replic,"_trajectory.RData") 
#drifttest_infile  <- paste0(working_dir,"/",simID,"_",replic,"_drifttest.dat") 
#drifttest_outfile <- paste0(working_dir,"/",simID,"_",replic,"_drifttest.out") 

