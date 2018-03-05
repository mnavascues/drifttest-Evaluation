if (Ne_only){
  
  # Global estimates of Fst, Fis and Ne
  #-----------------------------------------
  pdf(file=paste0("results/",simID,"/Fst_hat.pdf"))
  boxplot(results$Fst_hat,
          ylim=c(0,0.2),
          main=paste("Simulation",simID,""),
          ylab=expression(italic(F)[ST]))
  abline(h=E_Fst,col="red")
  dev.off()
  pdf(file=paste0("results/",simID,"/Fis_hat.pdf"))
  boxplot(results$Fis_hat,
          ylim=c(0,1),
          main=paste("Simulation",simID,""),
          ylab=expression(italic(F)[IS]))
  abline(h=E_Fis,col="red")
  dev.off()
  pdf(file=paste0("results/",simID,"/Ne_hat.pdf"))
  boxplot(results$NeHatFst,
          ylim=c(0,2000),
          main=paste("Simulation",simID,""),
          ylab=expression(italic(N)[e]))
  abline(h=E_Ne,col="red")
  dev.off()
  
  save(results,file=results_file)
  
}else{
    
  # close previous graphic devices (pdf files for Manhattan and QQ plots)
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.cur())
  
  # Global estimates of Fst, Fis and Ne
  #-----------------------------------------
  pdf(file=paste0("results/",simID,"/Fst_hat.pdf"))
  boxplot(results$Fst_hat,
          ylim=c(0,0.2),
          main=paste("Simulation",simID,""),
          ylab=expression(italic(F)[ST]))
  abline(h=E_Fst,col="red")
  dev.off()
  pdf(file=paste0("results/",simID,"/Fis_hat.pdf"))
  boxplot(results$Fis_hat,
          ylim=c(0,1),
          main=paste("Simulation",simID,""),
          ylab=expression(italic(F)[IS]))
  abline(h=E_Fis,col="red")
  dev.off()
  pdf(file=paste0("results/",simID,"/Ne_hat.pdf"))
  boxplot(results$Ne_hat,
          ylim=c(0,2000),
          main=paste("Simulation",simID,""),
          ylab=expression(italic(N)[e]))
  abline(h=E_Ne,col="red")
  dev.off()
  
  
  results_per_locus <- cbind(results_per_locus,
                             centimorgan=bp2centimorgan(results_per_locus$distance,r))
  
  
  
  
  positives <- length(intersect(which(results_per_locus$outliers_topSNPs),
                                which(results_per_locus$neutral)))
  total <- length(intersect(which(!is.na(results_per_locus$outliers_topSNPs)),
                            which(results_per_locus$neutral)))
  FPR_top50Pvalue          <- positives/total
  
  positives <- length(intersect(which(results_per_locus$outliers_topSNPs),
                                which(results_per_locus$distance==0)))
  total <- length(intersect(which(!is.na(results_per_locus$outliers_topSNPs)),
                            which(results_per_locus$distance==0)))
  POWER_top50Pvalue          <- positives/total
  
  
  
  
  positives <- length(intersect(which(results_per_locus$p_value <= 1e-3),
                                which(results_per_locus$neutral)))
  total <- length(intersect(which(!is.na(results_per_locus$p_value)),
                            which(results_per_locus$neutral)))
  FPR_thresholdPvalue0.001 <- positives/total
  
  positives <- length(intersect(which(results_per_locus$p_value <= 1e-3),
                                which(results_per_locus$distance==0)))
  total <- length(intersect(which(!is.na(results_per_locus$p_value)),
                            which(results_per_locus$distance==0)))
  POWER_thresholdPvalue0.001 <- positives/total
  
  
  
  positives <- length(intersect(which(results_per_locus$q_value <= 0.05),
                                which(results_per_locus$neutral)))
  total <- length(intersect(which(!is.na(results_per_locus$q_value)),
                            which(results_per_locus$neutral)))
  FPR_thresholdQvalue0.05  <- positives/total
  
  positives <- length(intersect(which(results_per_locus$q_value <= 0.05),
                                which(results_per_locus$distance==0)))
  total <- length(intersect(which(!is.na(results_per_locus$q_value)),
                            which(results_per_locus$distance==0)))
  POWER_thresholdQvalue0.05  <- positives/total
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ref_points <- 10^seq(from=-2,to=1.7,by=0.01)
  ds<-3
  
  SelectionFootprint_top50Pvalue          <- array(NA,dim=length(ref_points))
  SelectionFootprint_thresholdPvalue0.001 <- array(NA,dim=length(ref_points))
  SelectionFootprint_thresholdQvalue0.05  <- array(NA,dim=length(ref_points))
  for (bin in seq_along(ref_points)){
    lower_limit <- max(0,ref_points[bin]-ds)
    upper_limit <- ref_points[bin]+ds
    loci_in_bin <- intersect(which(results_per_locus$centimorgan>lower_limit),
                             which(results_per_locus$centimorgan<upper_limit))
    
    loci_in_bin <- intersect(loci_in_bin,
                             which(!results_per_locus$neutral))
  
    loci_in_bin <- intersect(loci_in_bin,
                             which(!is.na(results_per_locus$outliers_topSNPs)))
    
    #if(!quiet) print(length(loci_in_bin))
  
    dist_cM  <- results_per_locus$centimorgan[loci_in_bin]
    dist_w   <- sqrt( ( dist_cM - ref_points[bin] )^2 )
    weight   <- (35/32) * (1 - (dist_w/ds)^2)^3
  
    positives <- vector("numeric",length(loci_in_bin))
    positives[results_per_locus$outliers_topSNPs[loci_in_bin]] <- 1
    #res <- glm(positives~dist_cM,family=binomial,weights=weight)
    SelectionFootprint_top50Pvalue[bin] <- sum(positives*weight)/sum(weight)
  
    positives <- vector("numeric",length(loci_in_bin))
    positives[which(results_per_locus$p_value[loci_in_bin] <= 1e-3)] <- 1
    #res <- glm(positives~dist_cM,family=binomial,weights=weight)
    SelectionFootprint_thresholdPvalue0.001[bin] <- sum(positives*weight)/sum(weight)
    
    positives <- vector("numeric",length(loci_in_bin))
    positives[which(results_per_locus$q_value[loci_in_bin] <= 0.05)] <- 1
    #res <- glm(positives~dist_cM,family=binomial,weights=weight)
    SelectionFootprint_thresholdQvalue0.05[bin] <- sum(positives*weight)/sum(weight)
  }
  pdf(file=paste0("results/",simID,"/SelectionFootprint_top50pvalue.pdf"))
  plot(  ref_points,
         SelectionFootprint_top50Pvalue,
         ylim=c(0,0.5),
         xlim=c(0,40),
         type="l",
         #log="x",
         main=paste("Simulation",simID),
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)") 
  abline(h=FPR_top50Pvalue,col="red")
  abline(h=POWER_top50Pvalue,col="blue")
  dev.off()
  
  
  pdf(file=paste0("results/",simID,"/SelectionFootprint_thresholdPvalue0.001.pdf"))
  plot(  ref_points,
         SelectionFootprint_thresholdPvalue0.001,
         ylim=c(0,0.5),
         xlim=c(0,40),
         type="l",
         main=paste("Simulation",simID),
         ylab="Proportion of positives",
         #log="x",
         xlab="distance to selected locus (cM)") 
  abline(h=FPR_thresholdPvalue0.001,col="red")
  abline(h=POWER_thresholdPvalue0.001,col="blue")
  dev.off()
  
  pdf(file=paste0("results/",simID,"/SelectionFootprint_thresholdQvalue0.05.pdf"))
  plot(  ref_points,
         SelectionFootprint_thresholdQvalue0.05,
         ylim=c(0,0.5),
         xlim=c(0,40),
         type="l",
         #log="x",
         main=paste("Simulation",simID),
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)") 
  abline(h=FPR_thresholdQvalue0.05,col="red")
  abline(h=POWER_thresholdQvalue0.05,col="blue")
  dev.off()
  
  
  
  SelectionFootprint_top50Pvalue <- cbind(distance      = ref_points,
                                          positive_rate = SelectionFootprint_top50Pvalue)
  
  SelectionFootprint_thresholdPvalue0.001 <- cbind(distance      = ref_points,
                                                   positive_rate = SelectionFootprint_thresholdPvalue0.001)
  
  SelectionFootprint_thresholdqvalue0.05 <- cbind(distance      = ref_points,
                                                   positive_rate = SelectionFootprint_thresholdPvalue0.001)
  
  
  threshold4ROC <- seq(0,1,0.001)
  FPR <- TPR <- array(NA,length(threshold4ROC))
  nuetral_loci <- intersect(which(results_per_locus$neutral),
                            which(!is.na(results_per_locus$p_value)))
  selected_loci <- intersect(intersect(which(!results_per_locus$neutral),
                                       which(results_per_locus$centimorgan<=1)),
                             which(!is.na(results_per_locus$p_value)))
  for (T.hold in seq_along(threshold4ROC)){
    FPR[T.hold] <- sum(results_per_locus$p_value[nuetral_loci] <= threshold4ROC[T.hold])/length(nuetral_loci)
    TPR[T.hold] <- sum(results_per_locus$p_value[selected_loci] <= threshold4ROC[T.hold])/length(selected_loci)
  }
  pdf(file=paste0("results/",simID,"/ROC_curve.pdf"))
  plot(FPR,TPR,
       main=paste("Simulation",simID),
       ylab="True positive rate",
       xlab="False positive rate",
       type="l") 
  
  abline(0,1,col="red")
  dev.off()
  ROC<-data.frame(FPR,TPR)
  
  
  
  
  
  save(SelectionFootprint_top50Pvalue,
       SelectionFootprint_thresholdPvalue0.001,
       SelectionFootprint_thresholdqvalue0.05,
       FPR_top50Pvalue,
       FPR_thresholdPvalue0.001,
       FPR_thresholdQvalue0.05,
       POWER_top50Pvalue,
       POWER_thresholdPvalue0.001,
       POWER_thresholdQvalue0.05,
       ROC,
       results,
       results_per_locus,
       file=results_file)
}
