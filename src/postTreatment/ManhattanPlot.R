#######################################################
#
# ManhattanPlot & QQplot
#
#######################################################

source("src/fun/ColorBlindPalette.R")
palette(cbPalette2)

# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
library(qvalue)

project <- "NM"
number_of_scenarios <- 50
number_of_replicates <- 1000
number_of_loci <- 10000
num_of_rep4test <- 30
library(qqman)

# read simulations description table
sim_table <- read.table(paste0("results/",project,"/simparams.txt"),header=T)

scen <- 41
#for (scen in 1){ #seq_len(nrow(sim_table))){
  
  population_size           <- sim_table$population_size[scen]
  selection_mode            <- sim_table$selection_mode[scen]
  sel_coef                  <- sim_table$sel_coef[scen]
  sigma                     <- sim_table$sigma[scen]
  selection_period_duration <- sim_table$selection_period_duration[scen]
  
  first_locus_with_trajectory <- T
  counter_allele_not_lost <- 0
  replic <- 0
  
  pdf(file=paste0("results/",project,"/",project,scen,"ManhattanPlotPvalue.pdf"))
  pdf(file=paste0("results/",project,"/",project,scen,"ManhattanPlotQvalue.pdf"))
  pdf(file=paste0("results/",project,"/",project,scen,"QQPlot_drifttest.pdf"))
  pdf(file=paste0("results/",project,"/",project,scen,"TrajectoryPlot.pdf"))
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  dev.set(which = dev.prev())
  neutral_loci<-hitchhiked_loci<-NULL
  
  while (counter_allele_not_lost < num_of_rep4test & replic < number_of_replicates) {
    replic <- replic+1
    # replic <- 14
    print(paste("scenario",scen,"; replicate", replic))
    
    RData_results    <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_results.RData")
    RData_data       <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_data.RData")
    RData_params     <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_params.RData")
    RData_trajectory <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_trajectory.RData")
    
    drifttest_locus_by_locus <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_drifttest_locus_by_locus")
    drifttest_multilocus <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_drifttest_multilocus")

    
    #load(RData_params)
    load(RData_results)
    load(RData_trajectory)
    if (advantageous_allele_not_lost) counter_allele_not_lost <- counter_allele_not_lost+1
    
    if (advantageous_allele_not_lost & file.exists(drifttest_locus_by_locus)){
      
      
      load(RData_data)
      res_driftest <- read.table(drifttest_locus_by_locus,header=T)
      
      q_values <- qvalue( res_driftest$p_value)$qvalues
      
      

      selected_site_position <- m2$x
      
      chr1 <- which(SNP_data$SNP_table$x<=5e+8/2)
      chr2 <- which(SNP_data$SNP_table$x>5e+8/2)
      
      
      # Manhattan plot (p-value)
      #############################        
      plot(SNP_data$SNP_table$x[chr2],
           -log10(res_driftest$p_value[chr2]),
           xlim=c(0,5e+8),
           ylim=c(0,10),
           xlab="",
           main=paste0("Simulation ",project,scen," replicate ",replic),
           ylab=expression(-log[10]*"(p-value)"),
           axes=F,
           pch=16,
           #cex=0.5,
           col="dark grey")
      Axis(side=2)
      axis(side=1,at=c(5e+8/4,5e+8*3/4),labels=c("chr 1","chr 2"))
      box()
      
      points(SNP_data$SNP_table$x[chr1],
             -log10(res_driftest$p_value[chr1]),
             pch=16,
             #cex=0.5,
             col="black")
      abline(v=selected_site_position,col="red")
      #abline(v=selected_site_position + 10 / (1e-8 * 100),col="red")
      #abline(v=selected_site_position - 10 / (1e-8 * 100),col="red")

      
      # Manhattan plot (q-value)
      #############################        
      dev.set(which = dev.next())
      plot(SNP_data$SNP_table$x[chr2],
           -log10(q_values[chr2]),
           xlim=c(0,5e+8),
           ylim=c(0,5),
           xlab="",
           main=paste0("Simulation ",project,scen," replicate ",replic),
           ylab=expression(-log[10]*"(q-value)"),
           axes=F,
           pch=16,
           #cex=0.5,
           col="dark grey")
      Axis(side=2)
      axis(side=1,at=c(5e+8/4,5e+8*3/4),labels=c("chr 1","chr 2"))
      box()
      
      points(SNP_data$SNP_table$x[chr1],
             -log10(q_values[chr1]),
             pch=16,
             #cex=0.5,
             col="black")
      abline(v=selected_site_position,col="red")
      
      
      # QQ PLOT
      ###########################        
      dev.set(which = dev.next())
      qqman::qq(res_driftest$p_value,
                type="l",
                xlim=c(0,6),
                ylim=c(0,6),
                main=paste0("Simulation ",project,scen," replicate ",replic))
      
      
      
      
            
      # Trajectory of derived allele (selected locus)
      ################################################        
      dev.set(which = dev.next())
      if(counter_allele_not_lost==1){
        plot(trajectory/population_size/2,
             xlab="generation",
             ylab="derived allele frequency",
             ylim=c(0,1),
             type="l",
             col=1)        
      }else{
        lines(seq_along(trajectory),
              trajectory/population_size/2,
              col=counter_allele_not_lost)        
      }
      dev.set(which = dev.prev())
      dev.set(which = dev.prev())
      dev.set(which = dev.prev())
      
      
      
      if (selected_site_position <= 5e+8/2){
        selected_chromosome <- 1
        neutral_chromosome  <- 2
      }else{
        selected_chromosome <- 2
        neutral_chromosome  <- 1
      }
      loci4false_positive      <- get(paste0("chr",neutral_chromosome)) 
      loci4selection_footprint <- get(paste0("chr",selected_chromosome)) 

      
      outliers_topSNPs      <- res_driftest$p_value <= sort(res_driftest$p_value)[50]
      outliers_p_value0.001 <- res_driftest$p_value < 1e-3

            
      neutral_loci <- rbind(neutral_loci,
                            cbind(res_driftest[loci4false_positive,2:3],
                                  outliers_topSNPs=outliers_topSNPs[loci4false_positive],
                                  outliers_p_value0.001=outliers_p_value0.001[loci4false_positive]))
      
      hitchhiked_loci <- rbind(hitchhiked_loci,
                               cbind(res_driftest[loci4selection_footprint,2:3],
                                     distance=selected_site_position-SNP_data$SNP_table$x[loci4selection_footprint],
                                     outliers_topSNPs=outliers_topSNPs[loci4selection_footprint],
                                     outliers_p_value0.001=outliers_p_value0.001[loci4selection_footprint]))
      
      
      
      
      
      
    }
    
    

  }
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.next())
  dev.off(which = dev.cur())
  
  
  neutral_loci <- neutral_loci[!is.na(neutral_loci$outliers_topSNPs),]
  hitchhiked_loci <- hitchhiked_loci[!is.na(hitchhiked_loci$outliers_topSNPs),]
  
  FPR_top50Pvalue          <- sum(neutral_loci$outliers_topSNPs)/nrow(neutral_loci)
  FPR_thresholdPvalue0.001 <- sum(neutral_loci$outliers_p_value0.001)/nrow(neutral_loci)

  
  bin_size       <- 500000
  step_size      <- 250000
  ref_points <- seq(from=0,to=4e7,by=step_size)
  
  SelectionFootprint_top50Pvalue          <- array(NA,dim=length(ref_points))
  SelectionFootprint_thresholdPvalue0.001 <- array(NA,dim=length(ref_points))
  for (bin in seq_along(ref_points)){
    lower_limit <- ref_points[bin]
    upper_limit <- ref_points[bin]+bin_size
    loci_in_bin <- intersect(which(abs(hitchhiked_loci[,"distance"])>=lower_limit),
                             which(abs(hitchhiked_loci[,"distance"])<upper_limit))
    #print(length(loci_in_bin))
    SelectionFootprint_top50Pvalue[bin]          <- sum(hitchhiked_loci[loci_in_bin,"outliers_topSNPs"])/length(loci_in_bin)
    SelectionFootprint_thresholdPvalue0.001[bin] <- sum(hitchhiked_loci[loci_in_bin,"outliers_p_value0.001"])/length(loci_in_bin)
  }
  pdf(file=paste0("results/",project,"/",project,scen,"SelectionFootprint_top50pvalue.pdf"))
  plot(  (ref_points+bin_size/2)/1000000,
         SelectionFootprint_top50Pvalue,ylim=c(0,0.5),
         type="l",
         main=paste0("Simulation ",project,scen),
         ylab="Proportion of positives",
         xlab="distance to selected locus (Mbp)") 
  abline(h=FPR_topSNPs,col="red")
  dev.off()

  
  pdf(file=paste0("results/",project,"/",project,scen,"SelectionFootprint_thresholdPvalue0.001.pdf"))
  plot(  (ref_points+bin_size/2)/1000000,
         SelectionFootprint_thresholdPvalue0.001,ylim=c(0,0.5),
         type="l",
         main=paste0("Simulation ",project,scen),
         ylab="Proportion of positives",
         xlab="distance to selected locus (Mbp)") 
  abline(h=FPR_topSNPs,col="red")
  dev.off()
  
  SelectionFootprint_top50Pvalue <- cbind(distance      = (ref_points+bin_size/2)/1000000,
                                          positive_rate = SelectionFootprint_top50Pvalue)
  
  SelectionFootprint_thresholdPvalue0.001 <- cbind(distance      = (ref_points+bin_size/2)/1000000,
                                                   positive_rate = SelectionFootprint_thresholdPvalue0.001)
  
  
  save(neutral_loci,
       hitchhiked_loci,
       SelectionFootprint_top50Pvalue,SelectionFootprint_thresholdPvalue0.001,
       FPR_top50Pvalue,FPR_thresholdPvalue0.001,
       file=paste0("results/",project,"/",project,scen,"results_drifttest.RData"))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#}








