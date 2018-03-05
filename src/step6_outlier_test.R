if (Ne_only){
  Fstats <- compute_Fstats(SNP_data$genotype_data)

  results$Fc_hat[replic]   <- Fstats$F_C
  results$NeHatFc[replic]  <- EstimateNe.F_C(Fstats$F_C,selection_period_duration,sample_size[1],sample_size[2])
  results$Fst_hat[replic]  <- Fstats$F_ST
  results$NeHatFst[replic] <- EstimateNe.F_ST(Fstats$F_ST,selection_period_duration)
  results$Fis_hat[replic]  <- Fstats$F_IS
  
  replic <- replic+1

}else if (advantageous_allele_not_lost){
  system(paste("./bin/drifttest",
               "-tau", selection_period_duration,
               "-seed", round(runif(1,0,10000)),
               "-threads",num_of_threads,
               "-maf",MAF_threshold,
               "-infile", drifttest_infile,
               "-outfile",drifttest_outfile) )
  
  drifttest_locus_by_locus <- paste0(drifttest_outfile,"locus_by_locus")
  drifttest_multilocus     <- paste0(drifttest_outfile,"multilocus")

  res_driftest <- read.table(drifttest_locus_by_locus,header=T)
  res_driftest_multilocus <- read.table(drifttest_multilocus,header=T)

  results$Fst_hat[replic] <- res_driftest_multilocus$Fst 
  results$Fis_hat[replic] <- res_driftest_multilocus$Fis 
  results$Ne_hat[replic]  <- res_driftest_multilocus$Ne 
  
  success <- try(qvalue( res_driftest$p_value),silent = TRUE)
  if(class(success)=="try-error"){
    cat(paste0("q-values cannot be calculated for replicate ",replic,"\n"))
    q_values <- array(NA,length(nrow(res_driftest)))
    
  }else{
    q_values <- qvalue( res_driftest$p_value)$qvalues
  }
  
  
  
  fst                               <- array(NA,nrow(res_driftest))
  fst[!is.na(res_driftest$p_value)] <- res_driftest$F_ST[!is.na(res_driftest$p_value)]

  
  # define sets of loci to calculate FPR and selection footprint
  #----------------------------------------------------------------
  selected_site_position <- m2$x
  chr1 <- which(SNP_data$SNP_table$x<=genome_length/2)
  chr2 <- which(SNP_data$SNP_table$x>genome_length/2)
  if (selected_site_position <= genome_length/2){
    selected_chromosome <- 1
    neutral_chromosome  <- 2
  }else{
    selected_chromosome <- 2
    neutral_chromosome  <- 1
  }
  loci4false_positive      <- get(paste0("chr",neutral_chromosome)) 
  loci4selection_footprint <- get(paste0("chr",selected_chromosome))
  neutral_loci                      <- array(FALSE,sample_size_loci)
  neutral_loci[loci4false_positive] <- TRUE
  #selected_loci                     <- array(FALSE,sample_size_loci)
  #selected_loci[intersect(which(abs(selected_site_position-SNP_data$SNP_table$x)<=1e7),
  #                        loci4selection_footprint)] <- TRUE
  
  #----------------------------------------------------------------
  
  
  
  # Manhattan plot (Fst)
  #############################        
  plot(SNP_data$SNP_table$x[chr2],
       fst[chr2],
       xlim=c(0,genome_length),
       ylim=c(0,1),
       xlab="",
       main=paste("Simulation",simID,"replicate",replic),
       ylab=expression(italic("F")["ST"]),
       axes=F,
       pch=16,
       #cex=0.5,
       col="dark grey")
  Axis(side=2)
  axis(side=1,at=c(genome_length/4,genome_length*3/4),labels=c("chr 1","chr 2"))
  box()
  points(SNP_data$SNP_table$x[chr1],
         fst[chr1],
         pch=16,
         #cex=0.5,
         col="black")
  abline(v=selected_site_position,col="red")

  
    
  fst[fst<=0]<-NA
  z_scores <- sqrt(fst*(sum(sample_size)-2)/(1-fst))
  lambda1 <- (sum(sample_size)-2)*res_driftest_multilocus$Fst/(1-res_driftest_multilocus$Fst)
  #lambda2 <- median((z_scores)^2,na.rm=T)/qchisq(1/2,df=1)
  z2_p_values <- pchisq(z_scores^2/lambda1, df = 1, lower = FALSE)

  # Manhattan plot (z2-score)
  #############################        
  dev.set(which = dev.next())
  plot(SNP_data$SNP_table$x[chr2],
       (z_scores^2)[chr2],
       xlim=c(0,genome_length),
       ylim=c(0,max((z_scores[!is.infinite(z_scores)]^2),na.rm=T)),
       xlab="",
       main=paste("Simulation",simID,"replicate",replic),
       ylab=expression(italic("z")^2*"-score"),
       axes=F,
       pch=16,
       #cex=0.5,
       col="dark grey")
  Axis(side=2)
  axis(side=1,at=c(genome_length/4,genome_length*3/4),labels=c("chr 1","chr 2"))
  box()
  points(SNP_data$SNP_table$x[chr1],
         (z_scores^2)[chr1],
         pch=16,
         #cex=0.5,
         col="black")
  abline(v=selected_site_position,col="red")

  # Manhattan plot (p-value)
  #############################        
  dev.set(which = dev.next())
  plot(SNP_data$SNP_table$x[chr2],
       -log10(res_driftest$p_value[chr2]),
       xlim=c(0,genome_length),
       ylim=c(0,7),
       xlab="",
       main=paste("Simulation",simID,"replicate",replic),
       ylab=expression(-log[10]*"("*italic("p")*"-value)"),
       axes=F,
       pch=16,
       #cex=0.5,
       col="dark grey")
  Axis(side=2)
  axis(side=1,at=c(genome_length/4,genome_length*3/4),labels=c("chr 1","chr 2"))
  box()
  points(SNP_data$SNP_table$x[chr1],
         -log10(res_driftest$p_value[chr1]),
         pch=16,
         #cex=0.5,
         col="black")
  abline(v=selected_site_position,col="red")
  
  # Manhattan plot (q-value)
  #############################        
  dev.set(which = dev.next())
  plot(SNP_data$SNP_table$x[chr2],
       -log10(q_values[chr2]),
       xlim=c(0,genome_length),
       ylim=c(0,5),
       xlab="",
       main=paste("Simulation",simID,"replicate",replic),
       ylab=expression(-log[10]*"("*italic("q")*"-value)"),
       axes=F,
       pch=16,
       #cex=0.5,
       col="dark grey")
  Axis(side=2)
  axis(side=1,at=c(genome_length/4,genome_length*3/4),labels=c("chr 1","chr 2"))
  box()
  points(SNP_data$SNP_table$x[chr1],
         -log10(q_values[chr1]),
         pch=16,
         #cex=0.5,
         col="black")
  abline(v=selected_site_position,col="red")
  
  # QQ PLOT z2-scores
  ###########################        
  dev.set(which = dev.next())
  if (length(which(z2_p_values==0))>0){
    
    z2_p_values_4_plot <- z2_p_values[!is.na(z2_p_values)]
    uniform.quantiles <- qunif((1:length(z2_p_values_4_plot))/(length(z2_p_values_4_plot) +1))
    
    plot(sort(uniform.quantiles),
         sort(z2_p_values_4_plot),
         type = "l",
         xlim = c(0,1),
         ylim = c(0,1),
         xlab = expression(Expected ~ ~(italic(p))),
         ylab = expression(Observed ~ ~(italic(p))),
         main = paste("Simulation",simID,"replicate",replic,"(note not at log scale)"))
    abline(0,1,col="red")
  }else{
    qqman::qq(z2_p_values,
              type="l",
              xlim=c(0,ceiling(max(-log10(z2_p_values),na.rm=T))),
              ylim=c(0,ceiling(max(-log10(z2_p_values),na.rm=T))),
              main=paste("Simulation",simID,"replicate",replic))
  }

  # QQ PLOT z2-scores ONLY NEUTRAL
  ##################################       
  dev.set(which = dev.next())
  if (length(which(z2_p_values==0))>0){
    
    z2_p_values_4_plot <- z2_p_values[intersect(loci4false_positive, which(!is.na(z2_p_values)))]
    uniform.quantiles <- qunif((1:length(z2_p_values_4_plot))/(length(z2_p_values_4_plot) +1))
    
    plot(sort(uniform.quantiles),
         sort(z2_p_values_4_plot),
         type = "l",
         xlim = c(0,1),
         ylim = c(0,1),
         xlab = expression(Expected ~ ~(italic(p))),
         ylab = expression(Observed ~ ~(italic(p))),
         main = paste("Simulation",simID,"replicate",replic,"(note not at log scale)"))
    abline(0,1,col="red")
  }else{
    qqman::qq(z2_p_values[loci4false_positive],
              type="l",
              xlim=c(0,ceiling(max(-log10(z2_p_values),na.rm=T))),
              ylim=c(0,ceiling(max(-log10(z2_p_values),na.rm=T))),
              main=paste("Simulation",simID,"replicate",replic))
  }
  
  # QQ PLOT drifttest
  ###########################        
  dev.set(which = dev.next())
  qqman::qq(res_driftest$p_value,
            type="l",
            xlim=c(0,ceiling(max(-log10(res_driftest$p_value),na.rm=T))),
            ylim=c(0,ceiling(max(-log10(res_driftest$p_value),na.rm=T))),
            main=paste("Simulation",simID,"replicate",replic))
  
  # QQ PLOT drifttest ONLY NEUTRAL
  ###################################        
  dev.set(which = dev.next())
  qqman::qq(res_driftest$p_value[loci4false_positive],
            type="l",
            xlim=c(0,ceiling(max(-log10(res_driftest$p_value),na.rm=T))),
            ylim=c(0,ceiling(max(-log10(res_driftest$p_value),na.rm=T))),
            main=paste("Simulation",simID,"replicate",replic))
  
  
  # Trajectory of derived allele (selected locus)
  ################################################        
  dev.set(which = dev.next())
  
  if (length(trajectory)>0){
    if (advantageous_allele=="derived"){
      trajectory_advantageus <- trajectory/N/2
    }else if (advantageous_allele=="ancestral"){
      trajectory_advantageus <- (2*N - trajectory)/N/2
    }
    if(first_trajectory){
      plot(trajectory_advantageus,
           xlab="generation",
           ylab="derived allele frequency",
           ylim=c(0,1),
           type="l",
           col=1)
      first_trajectory<-F
    }else{
      lines(seq_along(trajectory_advantageus),
            trajectory_advantageus,
            col=replic)        
    }
  }
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())

  results_per_locus$Fst[sample_size_loci*(replic-1)+(1:sample_size_loci)]     <- res_driftest$F_ST
  results_per_locus$p_value[sample_size_loci*(replic-1)+(1:sample_size_loci)] <- res_driftest$p_value
  results_per_locus$q_value[sample_size_loci*(replic-1)+(1:sample_size_loci)] <- q_values

  results_per_locus$distance[sample_size_loci*(replic-1)+(1:sample_size_loci)] <- abs(selected_site_position-SNP_data$SNP_table$x)
  results_per_locus$neutral[sample_size_loci*(replic-1)+(1:sample_size_loci)]  <- neutral_loci
  #results_per_locus$selected[sample_size_loci*(replic-1)+(1:sample_size_loci)] <- selected_loci
  
  results_per_locus$outliers_topSNPs[sample_size_loci*(replic-1)+(1:sample_size_loci)]      <- res_driftest$p_value <= sort(res_driftest$p_value)[50]
  #results_per_locus$outliers_p_value0.001[sample_size_loci*(replic-1)+(1:sample_size_loci)] <- res_driftest$p_value < 1e-3
  #results_per_locus$outliers_q_value0.05[sample_size_loci*(replic-1)+(1:sample_size_loci)]  <- q_values < 0.05

  results_per_replica=results_per_locus[sample_size_loci*(replic-1)+(1:sample_size_loci),]
  
  #save data from replicate
  save(SNP_data,
       trajectory,
       results_per_replica,
       m2,
       z_scores,
       lambda1,
       #lambda2,
       z2_p_values,
       file=data_file)
  
  
  
  replic <- replic+1
}



















