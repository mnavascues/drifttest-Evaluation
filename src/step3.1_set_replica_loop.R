
if (Ne_only){
  results <- data.frame(Fc_hat = array(NA,number_of_replicates),
                        NeHatFc = array(NA,number_of_replicates),
                        Fst_hat = array(NA,number_of_replicates),
                        NeHatFst  = array(NA,number_of_replicates),
                        Fis_hat = array(NA,number_of_replicates))

}else{
  pdf(file=paste0(working_dir,"/",simID,"_Manhattan_Fst.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_Manhattan_z2score.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_Manhattan_Pvalue.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_Manhattan_Qvalue.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_QQ_z2score.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_QQ_z2score_onlyneutral.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_QQ_drifttest.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_QQ_drifttest_onlyneutral.pdf"))
  pdf(file=paste0(working_dir,"/",simID,"_Trajectories.pdf"))
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  
  results_per_locus <- data.frame(Fst                   = array(NA, sample_size_loci*number_of_replicates),
                                  p_value               = array(NA, sample_size_loci*number_of_replicates),
                                  q_value               = array(NA, sample_size_loci*number_of_replicates),
                                  distance              = array(NA, sample_size_loci*number_of_replicates),
                                  neutral               = array(FALSE, sample_size_loci*number_of_replicates),
                                  outliers_topSNPs      = array(FALSE, sample_size_loci*number_of_replicates))
  
  
  
  first_trajectory<-T
  
  results <- data.frame(Fst_hat = array(NA,number_of_replicates),
                        Fis_hat = array(NA,number_of_replicates),
                        Ne_hat  = array(NA,number_of_replicates))
  
}

replic <- 1

