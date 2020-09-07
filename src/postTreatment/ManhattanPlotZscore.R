#######################################################
#
# ManhattanPlot & QQplot
#
#######################################################

source("src/fun/ColorBlindPalette.R")
palette(cbPalette2)

# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
#library(qvalue)

project <- "SV"
number_of_scenarios <- 50
number_of_replicates <- 1000
number_of_loci <- 10000
num_of_rep4test <- 200

# read simulations description table
sim_table <- read.table(paste0("results/",project,"/simparams.txt"),header=T)

scen <- 31
for (scen in 31:50){#seq_len(nrow(sim_table))){

    population_size           <- sim_table$population_size[scen]
    selection_mode            <- sim_table$selection_mode[scen]
    sel_coef                  <- sim_table$sel_coef[scen]
    sigma                     <- sim_table$sigma[scen]
    selection_period_duration <- sim_table$selection_period_duration[scen]
    
    first_locus_with_trajectory <- T
    counter_allele_not_lost <- 0
    replic <- 0
    
    pdf(file=paste0("results/",project,"/",project,scen,"ManhattanPlot_Fst.pdf"))
    pdf(file=paste0("results/",project,"/",project,scen,"ManhattanPlot_z2score.pdf"))
    pdf(file=paste0("results/",project,"/",project,scen,"QQPlot_z2score.pdf"))
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
      
      drifttest_locus_by_locus <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_FSTonly_locus_by_locus")
      drifttest_multilocus <- paste0("results/",project,"/",project,scen,"/",project,scen,"_",replic,"_FSTonly_multilocus")
      
      
      #load(RData_params)
      load(RData_results)
      load(RData_trajectory)
      if (advantageous_allele_not_lost) counter_allele_not_lost <- counter_allele_not_lost+1
      
      if (advantageous_allele_not_lost & file.exists(drifttest_locus_by_locus)){
        
        
        load(RData_data)
        res_driftest <- read.table(drifttest_locus_by_locus,header=T)
        res_driftest_multilocus <- read.table(drifttest_multilocus,header=T)
        
        fst <- array(NA,nrow(res_driftest))
        fst[is.nan(res_driftest$p_value)] <- res_driftest$F_ST[is.nan(res_driftest$p_value)]
        
        
        
        selected_site_position <- m2$x
        
        chr1 <- which(SNP_data$SNP_table$x<=5e+8/2)
        chr2 <- which(SNP_data$SNP_table$x>5e+8/2)
        
        
        # Manhattan plot (Fst)
        #############################        
        plot(SNP_data$SNP_table$x[chr2],
             fst[chr2],
             xlim=c(0,5e+8),
             ylim=c(0,1),
             xlab="",
             main=paste0("Simulation ",project,scen," replicate ",replic),
             ylab=expression(italic("F")["ST"]),
             axes=F,
             pch=16,
             #cex=0.5,
             col="dark grey")
        Axis(side=2)
        axis(side=1,at=c(5e+8/4,5e+8*3/4),labels=c("chr 1","chr 2"))
        box()
        
        points(SNP_data$SNP_table$x[chr1],
               fst[chr1],
               pch=16,
               #cex=0.5,
               col="black")
        abline(v=selected_site_position,col="red")
        #abline(v=selected_site_position + 10 / (1e-8 * 100),col="red")
        #abline(v=selected_site_position - 10 / (1e-8 * 100),col="red")
        
        fst[fst<=0]<-NA
        sample_size <- 100
        z_scores <- sqrt(fst*(sample_size-2)/(1-fst))
        lambda1 <- (sample_size-2)*res_driftest_multilocus$Fst/(1-res_driftest_multilocus$Fst)
        #lambda2 <- median((z_scores)^2,na.rm=T)/qchisq(1/2,df=1)
        p_values = pchisq(z_scores^2/lambda1, df = 1, lower = FALSE)
        #q_values <- qvalue(p_values)$qvalues
        
        
        
        # Manhattan plot (z2-score)
        #############################        
        dev.set(which = dev.next())
        plot(SNP_data$SNP_table$x[chr2],
             (z_scores^2)[chr2],
             xlim=c(0,5e+8),
             ylim=c(0,max((z_scores[!is.infinite(z_scores)]^2),na.rm=T)),
             xlab="",
             main=paste0("Simulation ",project,scen," replicate ",replic),
             ylab=expression(italic("z")*"-score"),
             axes=F,
             pch=16,
             #cex=0.5,
             col="dark grey")
        Axis(side=2)
        axis(side=1,at=c(5e+8/4,5e+8*3/4),labels=c("chr 1","chr 2"))
        box()
        
        points(SNP_data$SNP_table$x[chr1],
               (z_scores^2)[chr1],
               pch=16,
               #cex=0.5,
               col="black")
        abline(v=selected_site_position,col="red")
        

        # QQ PLOT
        ###########################        
        dev.set(which = dev.next())
        qqman::qq(p_values,
                  type="l",
                  xlim=c(0,6),
                  ylim=c(0,6),
                  main=paste0("Simulation ",project,scen," replicate ",replic))

        
        # Trajectory of derived allele (selected locus)
        ################################################        
        dev.set(which = dev.next())
        if (length(trajectory)!=0){
          if(first_locus_with_trajectory){
            plot(trajectory/population_size/2,
                 xlab="generation",
                 ylab="derived allele frequency",
                 main=paste0("Simulation ",project,scen),
                 ylim=c(0,1),
                 type="l",
                 col=1)        
            first_locus_with_trajectory <- F
          }else{
            lines(seq_along(trajectory),
                  trajectory/population_size/2,
                  col=counter_allele_not_lost)        
          }
          
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
        
        
        outliers_topSNPs      <- p_values <= sort(p_values)[50]
        #outliers_p_value0.001 <- p_values < 1e-3
        #outliers_q_value0.01   <- q_values<0.01
        
        
        neutral_loci <- rbind(neutral_loci,
                              cbind(F_ST = res_driftest[loci4false_positive,2],
                                    p_value = p_values[loci4false_positive],
                                    outliers_topSNPs=outliers_topSNPs[loci4false_positive]))
                                    #outliers_p_value0.001=outliers_p_value0.001[loci4false_positive]))
                                    #outliers_q_value0.01=outliers_q_value0.01[loci4false_positive]))
        
        hitchhiked_loci <- rbind(hitchhiked_loci,
                                 cbind(F_ST = res_driftest[loci4selection_footprint,2],
                                       p_value = p_values[loci4selection_footprint],
                                       distance=selected_site_position-SNP_data$SNP_table$x[loci4selection_footprint],
                                       outliers_topSNPs=outliers_topSNPs[loci4selection_footprint]))
                                       #outliers_p_value0.001=outliers_p_value0.001[loci4selection_footprint]))
                                       #outliers_q_value0.01=outliers_q_value0.01[loci4selection_footprint]))
        
        
        
        
        
        
      }
      
      
      
    }
    dev.off(which = dev.next())
    dev.off(which = dev.next())
    dev.off(which = dev.next())
    dev.off(which = dev.cur())
    
    
    neutral_loci <- neutral_loci[!is.na(neutral_loci[,"outliers_topSNPs"]),]
    hitchhiked_loci <- hitchhiked_loci[!is.na(hitchhiked_loci[,"outliers_topSNPs"]),]
    
    FPR_top50Fst      <- sum(neutral_loci[,"outliers_topSNPs"])/nrow(neutral_loci)
    
    bin_size       <- 500000
    step_size      <- 250000
    SelectionFootprint_top50Fst      <- array(NA,dim=number_of_bins)
    ref_points <- seq(from=0,to=4e7,by=step_size)
    for (bin in seq_along(ref_points)){
      lower_limit <- ref_points[bin]
      upper_limit <- ref_points[bin]+bin_size
      loci_in_bin <- intersect(which(abs(hitchhiked_loci[,"distance"])>=lower_limit),
                               which(abs(hitchhiked_loci[,"distance"])<upper_limit))
      #print(length(loci_in_bin))
      SelectionFootprint_top50Fst[bin]      <- sum(hitchhiked_loci[loci_in_bin,"outliers_topSNPs"])/length(loci_in_bin)
    }
    pdf(file=paste0("results/",project,"/",project,scen,"SelectionFootprint_top50Fst.pdf"))
    plot(  (ref_points+bin_size/2)/1000000,
           SelectionFootprint_top50Fst,ylim=c(0,0.5),
           type="l",
           main=paste0("Simulation ",project,scen),
           ylab="Proportion of positives",
           xlab="distance to selected locus (Mbp)") 
    abline(h=FPR_topSNPs,col="red")
    dev.off()
    
    SelectionFootprint_top50Fst <- cbind(distance      = (ref_points+bin_size/2)/1000000,
                                         positive_rate = SelectionFootprint_top50Fst)
    
    
    save(neutral_loci,
         hitchhiked_loci,
         SelectionFootprint_top50Fst,
         FPR_top50Fst,
         file=paste0("results/",project,"/",project,scen,"results_topFst.RData"))
    
    
    



}








