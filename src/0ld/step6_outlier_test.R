# Fst outlier neutrality test
if(Ne_only){
  num_of_sim_test <- 1
  results <- FST_outlier_test(SNP_data$genotype_data,
                             MAF_threshold=MAF_threshold,
                             delta_T=selection_period_duration,
                             num_of_sim_test=num_of_sim_test)
  Fst <- results$results_total$F_ST
}else{
  results <- FST_outlier_test(SNP_data$genotype_data,
                             MAF_threshold=MAF_threshold,
                             delta_T=selection_period_duration,
                             num_of_sim_test=num_of_sim_test)
  save(results,file=results_file)
  
  Fst <- results$results_total$F_ST
  p_value_threshold <- 5/num_of_sim_test
  any_p_value_2_repeat   <- any(results$results_by_locus$p_value < p_value_threshold, na.rm = T)
  
  while (any_p_value_2_repeat & num_of_sim_test<1e3){
    
    
    which_locus_to_repeat <- which(results$results_by_locus$p_value<p_value_threshold)
    
    num_of_sim_test <- num_of_sim_test * 10
    print(paste("Estimating p-value from",num_of_sim_test,"simulations"))
    
    temp_res <- FST_outlier_test(SNP_data$genotype_data[,c(1,2,which_locus_to_repeat+2)],
                                MAF_threshold = MAF_threshold,
                                delta_T = selection_period_duration,
                                num_of_sim_test = num_of_sim_test,
                                Fst=Fst)
    gc()
    
    results$results_by_locus$p_value[which_locus_to_repeat] <- temp_res$results_by_locus$p_value
    
    p_value_threshold <- 5/num_of_sim_test
    any_p_value_2_repeat   <- any(results$results_by_locus$p_value < p_value_threshold, na.rm = T)
    
    remove(temp_res)
    gc()
    
    save(results,file=results_file)
  }
  
}

