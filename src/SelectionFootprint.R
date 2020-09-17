
source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)

num_of_bootstrap <-1000

load(file="results/simtable.RData")

pdf(file="results/SelectionFootprint.pdf",width=7,height=5.5)
layout(matrix(c(1,2,3,4), 2,2,byrow = TRUE))


palette(cbPalette1)
#pdf(file="results/SelectionFootprintNM.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- intersect(which(sim_table$selection_mode=="NM"),
                       which(sim_table$sel_coef==0.5))
scenarios <- intersect(scenarios,
                       which(sim_table$selection_period_duration==25))

#scenarios <- c(1,2,3,4,6,7,8,9)
sigma_values <-  c(0.00,0.50,0.80,0.90,0.95,1.00) #sort(sim_table$sigma[scenarios])
  

for (i in seq_along(sigma_values)){
  
  sigma <- sigma_values[i]
  simID <- sim_table$simID[intersect(which(sim_table$sigma==sigma),scenarios)]
  
  SelectionFootprint <- readRDS(file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  #load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  #ref_points <- 10^seq(from=-2,to=1.7,by=0.1)
  #ds <- ref_points*0.8
  
  #SelectionFootprint <- array(NA,dim=length(ref_points))
  #lower_SelectionFootprint <- array(NA,dim=length(ref_points))
  #upper_SelectionFootprint <- array(NA,dim=length(ref_points))
  #for (bin in seq_along(ref_points)){
  #  lower_limit <- max(0,ref_points[bin]-ds[bin])
  #  upper_limit <- ref_points[bin]+ds[bin]
  #  loci_in_bin <- intersect(which(results_per_locus$centimorgan>lower_limit),
  #                           which(results_per_locus$centimorgan<upper_limit))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!results_per_locus$neutral))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!is.na(results_per_locus$outliers_topSNPs)))
  #  
  #  if(length(loci_in_bin)>500){
  #    dist_cM  <- results_per_locus$centimorgan[loci_in_bin]
  #    dist_w   <- sqrt( ( dist_cM - ref_points[bin] )^2 )
  #    weight   <- (35/32) * (1 - (dist_w/ds[bin])^2)^3
  #    
  #    positives <- vector("numeric",length(loci_in_bin))
  #    positives[which(results_per_locus$p_value[loci_in_bin] <= 1e-3)] <- 1
  #    SelectionFootprint[bin] <- sum(positives*weight)/sum(weight)
  #    
  #    bootstrap_values <- array(NA,num_of_bootstrap)
  #    for (b in seq_len(num_of_bootstrap)){
  #      boot_sample <- sample(length(loci_in_bin),replace=T)
  #      bootstrap_values[b] <- sum(positives[boot_sample]*weight[boot_sample])/sum(weight[boot_sample])
  #    }
  #    lower_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.025)
  #    upper_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.975)
  #  }
  #}
  #SelectionFootprint <- cbind(distance      = ref_points,
  #                            positive_rate = SelectionFootprint,
  #                            upper95       = upper_SelectionFootprint,
  #                            lower95       = lower_SelectionFootprint)
  #saveRDS(SelectionFootprint,
  #        file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))

  if (i==1){
    plot(SelectionFootprint[,"distance"],
         SelectionFootprint[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,0.24),
         xlim=c(1,50),
         log="x",
         type="l",
         ylab="Proportion of positives",
         xlab="",#"distance to selected locus (cM)"
         )
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }else{
    lines(SelectionFootprint[,"distance"],
          SelectionFootprint[,"positive_rate"],
          col=i)
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }
}

legend(x=10,y=0.24,
       legend=c(expression(sigma *"="* 0.00),
                expression(sigma *"="* 0.50),
                expression(sigma *"="* 0.80),
                expression(sigma *"="* 0.90),
                expression(sigma *"="* 0.95),
                expression(sigma *"="* 1.00)),
       lty=1,
       col=1:6,
       cex=0.75,
       bty="n")
text(x=40,y=0.23,label="a",cex=1.5)
#dev.off()






palette(cbPalette1)


#pdf(file="results/SelectionFootprintSV.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- intersect(which(sim_table$selection_mode=="SV"),
                       which(sim_table$sel_coef==0.5))
scenarios <- intersect(scenarios,
                       which(sim_table$selection_period_duration==25))
scenarios <- intersect(scenarios,
                       which(sim_table$initial_frequency==-1))

#scenarios <- c(1,6,3,9)
sigma_values <-  c(0.00,0.50,0.80,0.90,0.95,1.00) #sort(sim_table$sigma[scenarios])

for (i in seq_along(sigma_values)){
  
  sigma <- sigma_values[i]
  simID <- sim_table$simID[intersect(which(sim_table$sigma==sigma),scenarios)]
  
  SelectionFootprint <- readRDS(file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  
  #load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  #ref_points <- 10^seq(from=-2,to=1.7,by=0.1)
  #ds <- ref_points*0.8
  
  #SelectionFootprint <- array(NA,dim=length(ref_points))
  #lower_SelectionFootprint <- array(NA,dim=length(ref_points))
  #upper_SelectionFootprint <- array(NA,dim=length(ref_points))
  #for (bin in seq_along(ref_points)){
  #  lower_limit <- max(0,ref_points[bin]-ds[bin])
  #  upper_limit <- ref_points[bin]+ds[bin]
  #  loci_in_bin <- intersect(which(results_per_locus$centimorgan>lower_limit),
  #                           which(results_per_locus$centimorgan<upper_limit))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!results_per_locus$neutral))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!is.na(results_per_locus$outliers_topSNPs)))
  #  
  #  if(length(loci_in_bin)>500){
  #    dist_cM  <- results_per_locus$centimorgan[loci_in_bin]
  #    dist_w   <- sqrt( ( dist_cM - ref_points[bin] )^2 )
  #    weight   <- (35/32) * (1 - (dist_w/ds[bin])^2)^3
  #    
  #    positives <- vector("numeric",length(loci_in_bin))
  #    positives[which(results_per_locus$p_value[loci_in_bin] <= 1e-3)] <- 1
  #    SelectionFootprint[bin] <- sum(positives*weight)/sum(weight)
  #    
  #    bootstrap_values <- array(NA,num_of_bootstrap)
  #    for (b in seq_len(num_of_bootstrap)){
  #      boot_sample <- sample(length(loci_in_bin),replace=T)
  #      bootstrap_values[b] <- sum(positives[boot_sample]*weight[boot_sample])/sum(weight[boot_sample])
  #    }
  #    lower_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.025)
  #    upper_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.975)
  #  }
  #}
  #SelectionFootprint <- cbind(distance      = ref_points,
  #                            positive_rate = SelectionFootprint,
  #                            upper95       = upper_SelectionFootprint,
  #                            lower95       = lower_SelectionFootprint)
  #saveRDS(SelectionFootprint,
  #        file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  if (i==1){
    plot(SelectionFootprint[,"distance"],
         SelectionFootprint[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,0.24),
         xlim=c(1,50),
         log="x",
         type="l",
         ylab="",#"Proportion of positives",
         xlab="",#"distance to selected locus (cM)"
         )
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }else{
    lines(SelectionFootprint[,"distance"],
          SelectionFootprint[,"positive_rate"],
          col=i)
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }
  
}

legend(x=10,y=0.24,
       legend=c(expression(sigma *"="* 0.00),
                expression(sigma *"="* 0.50),
                expression(sigma *"="* 0.80),
                expression(sigma *"="* 0.90),
                expression(sigma *"="* 0.95),
                expression(sigma *"="* 1)),
       lty=1,
       col=1:6,
       cex=0.75,
       bty="n")
text(x=40,y=0.23,label="b",cex=1.5)

#dev.off()















  palette(cbPalette1)
#pdf(file="results/SelectionFootprintSelCoef.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- c(49,50,51,52)#,1)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]

  SelectionFootprint <- readRDS(file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  #load(file=paste0("results/",simID,"/",simID,"_results.RData"))

  #ref_points <- 10^seq(from=-2,to=1.7,by=0.1)
  #ds <- ref_points*0.8
  #
  #SelectionFootprint <- array(NA,dim=length(ref_points))
  #lower_SelectionFootprint <- array(NA,dim=length(ref_points))
  #upper_SelectionFootprint <- array(NA,dim=length(ref_points))
  #for (bin in seq_along(ref_points)){
  #  lower_limit <- max(0,ref_points[bin]-ds[bin])
  #  upper_limit <- ref_points[bin]+ds[bin]
  #  loci_in_bin <- intersect(which(results_per_locus$centimorgan>lower_limit),
  #                           which(results_per_locus$centimorgan<upper_limit))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!results_per_locus$neutral))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!is.na(results_per_locus$outliers_topSNPs)))
  #  
  #  if(length(loci_in_bin)>500){
  #    dist_cM  <- results_per_locus$centimorgan[loci_in_bin]
  #    dist_w   <- sqrt( ( dist_cM - ref_points[bin] )^2 )
  #    weight   <- (35/32) * (1 - (dist_w/ds[bin])^2)^3
  #    
  #    positives <- vector("numeric",length(loci_in_bin))
  #    positives[which(results_per_locus$p_value[loci_in_bin] <= 1e-3)] <- 1
  #    SelectionFootprint[bin] <- sum(positives*weight)/sum(weight)
  #    
  #    bootstrap_values <- array(NA,num_of_bootstrap)
  #    for (b in seq_len(num_of_bootstrap)){
  #      boot_sample <- sample(length(loci_in_bin),replace=T)
  #      bootstrap_values[b] <- sum(positives[boot_sample]*weight[boot_sample])/sum(weight[boot_sample])
  #    }
  #    lower_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.025)
  #    upper_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.975)
  #  }
  #}
  #SelectionFootprint <- cbind(distance      = ref_points,
  #                            positive_rate = SelectionFootprint,
  #                            upper95       = upper_SelectionFootprint,
  #                            lower95       = lower_SelectionFootprint)
  #saveRDS(SelectionFootprint,
  #        file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  if (i==1){
    plot(SelectionFootprint[,"distance"],
         SelectionFootprint[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,0.24),
         xlim=c(1,50),
         log="x",
         type="l",
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)"
         )
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }else{
    lines(SelectionFootprint[,"distance"],
          SelectionFootprint[,"positive_rate"],
          col=i)
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }
}

legend(x=10,y=0.24,
       legend=c(expression("s="* 0.1),
                expression("s="* 0.2),
                expression("s="* 0.3),
                expression("s="* 0.4)),
       lty=1,
       col=1:4,
       cex=0.75,
       bty="n")
text(x=40,y=0.23,label="c",cex=1.5)

#dev.off()


































#pdf(file="results/SelectionFootprintTime5.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- c(58,57,55,56,60)# c(58,57,1,55,56,60)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  
  SelectionFootprint <- readRDS(file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  #load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  #
  #ref_points <- 10^seq(from=-2,to=1.7,by=0.1)
  #ds <- ref_points*0.8
  #
  #SelectionFootprint <- array(NA,dim=length(ref_points))
  #lower_SelectionFootprint <- array(NA,dim=length(ref_points))
  #upper_SelectionFootprint <- array(NA,dim=length(ref_points))
  #for (bin in seq_along(ref_points)){
  #  lower_limit <- max(0,ref_points[bin]-ds[bin])
  #  upper_limit <- ref_points[bin]+ds[bin]
  #  loci_in_bin <- intersect(which(results_per_locus$centimorgan>lower_limit),
  #                           which(results_per_locus$centimorgan<upper_limit))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!results_per_locus$neutral))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!is.na(results_per_locus$outliers_topSNPs)))
  #  
  #  if(length(loci_in_bin)>500){
  #    dist_cM  <- results_per_locus$centimorgan[loci_in_bin]
  #    dist_w   <- sqrt( ( dist_cM - ref_points[bin] )^2 )
  #    weight   <- (35/32) * (1 - (dist_w/ds[bin])^2)^3
  #    
  #    positives <- vector("numeric",length(loci_in_bin))
  #    positives[which(results_per_locus$p_value[loci_in_bin] <= 1e-3)] <- 1
  #    SelectionFootprint[bin] <- sum(positives*weight)/sum(weight)
  #    
  #    bootstrap_values <- array(NA,num_of_bootstrap)
  #    for (b in seq_len(num_of_bootstrap)){
  #      boot_sample <- sample(length(loci_in_bin),replace=T)
  #      bootstrap_values[b] <- sum(positives[boot_sample]*weight[boot_sample])/sum(weight[boot_sample])
  #    }
  #    lower_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.025)
  #    upper_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.975)
  #  }
  #}
  #SelectionFootprint <- cbind(distance      = ref_points,
  #                            positive_rate = SelectionFootprint,
  #                            upper95       = upper_SelectionFootprint,
  #                            lower95       = lower_SelectionFootprint)
  #saveRDS(SelectionFootprint,
  #        file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  if (i==1){
    plot(SelectionFootprint[,"distance"],
         SelectionFootprint[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,0.24),
         xlim=c(1,50),
         log="x",
         type="l",
         ylab="",#"Proportion of positives",
         xlab="distance to selected locus (cM)")
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }else{
    lines(SelectionFootprint[,"distance"],
          SelectionFootprint[,"positive_rate"],
          col=i)
    polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
            y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
            col=makeTransparent(i,alpha=0.1),
            border=NA)
  }
}

legend(x=10,y=0.24,
       legend=c(expression(tau *"="* 5),
                expression(tau *"="* 10),
                #expression(tau *"="* 25),
                expression(tau *"="* 50),
                expression(tau *"="* 100),
                expression(tau *"="* 200)),
       lty=1,
       col=1:5,
       cex=0.75,
       bty="n")
text(x=40,y=0.23,label="d",cex=1.5)

dev.off()







#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################




pdf(file="results/SelectionFootprintTime1.pdf",width=4,height=4)
par(mar=c(4,4,0.2,0.2)+0.1)

scenarios <- c(49,53,54,59)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  
  load(file=paste0("results/",simID,"/",simID,"_results.RData"))
  
  if (i==1){
    plot(SelectionFootprint_thresholdPvalue0.001[,"distance"],
         SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
         cex.lab=0.75,
         cex.axis=0.6,
         col=i,
         ylim=c(0,1),
         #xlim=c(0,40),
         log="x",
         type="l",
         ylab="Proportion of positives",
         xlab="distance to selected locus (cM)")
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }else{
    lines(SelectionFootprint_thresholdPvalue0.001[,"distance"],
          SelectionFootprint_thresholdPvalue0.001[,"positive_rate"],
          #lty=i,
          col=i)
    points(60,FPR_thresholdPvalue0.001,col=i)
    points(0.009,POWER_thresholdPvalue0.001,col=i)
    
  }
  
}

legend(x=1,y=1,
       legend=c(expression(tau *"="* 25),
                expression(tau *"="* 50),
                expression(tau *"="* 100),
                expression(tau *"="* 200)),
       lty=1,
       col=1:4,
       cex=0.75,
       bty="n")

dev.off()













palette(cbPalette1)

pdf(file="results/SelectionFootprintSVsuppl.pdf",width=4,height=3)
par(mar=c(4,4,0.2,0.2)+0.1)


#plot new mutation (scenario007)
#load(file=paste0("results/scenario007/scenario007_results.RData"))
SelectionFootprint <- readRDS(file=paste0("results/scenario007/scenario007_SelectionFootprint.RDS"))

plot(SelectionFootprint[,"distance"],
     SelectionFootprint[,"positive_rate"],
     cex.lab=0.75,
     cex.axis=0.6,
     col=1,
     ylim=c(0,0.24),
     xlim=c(1,50),
     lty=2,
     log="x",
     type="l",
     ylab="Proportion of positives",
     xlab="distance to selected locus (cM)")
polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
        y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
        col=makeTransparent(1,alpha=0.1),
        border=NA)

scenarios <- c(31,35,39)
for (i in seq_along(scenarios)){
  
  simID <- sim_table$simID[scenarios[i]]
  SelectionFootprint <- readRDS(file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  #load(file=paste0("results/",simID,"/",simID,"_results.RData"))
 
  #ref_points <- 10^seq(from=-2,to=1.7,by=0.1)
  #ds <- ref_points*0.8
  
  #SelectionFootprint <- array(NA,dim=length(ref_points))
  #lower_SelectionFootprint <- array(NA,dim=length(ref_points))
  #upper_SelectionFootprint <- array(NA,dim=length(ref_points))
  #for (bin in seq_along(ref_points)){
  #  lower_limit <- max(0,ref_points[bin]-ds[bin])
  #  upper_limit <- ref_points[bin]+ds[bin]
  #  loci_in_bin <- intersect(which(results_per_locus$centimorgan>lower_limit),
  #                           which(results_per_locus$centimorgan<upper_limit))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!results_per_locus$neutral))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!is.na(results_per_locus$outliers_topSNPs)))
  #  
  #  if(length(loci_in_bin)>500){
  #    dist_cM  <- results_per_locus$centimorgan[loci_in_bin]
  #    dist_w   <- sqrt( ( dist_cM - ref_points[bin] )^2 )
  #    weight   <- (35/32) * (1 - (dist_w/ds[bin])^2)^3
  #    
  #    positives <- vector("numeric",length(loci_in_bin))
  #    positives[which(results_per_locus$p_value[loci_in_bin] <= 1e-3)] <- 1
  #    SelectionFootprint[bin] <- sum(positives*weight)/sum(weight)
  #    
  #    bootstrap_values <- array(NA,num_of_bootstrap)
  #    for (b in seq_len(num_of_bootstrap)){
  #      boot_sample <- sample(length(loci_in_bin),replace=T)
  #      bootstrap_values[b] <- sum(positives[boot_sample]*weight[boot_sample])/sum(weight[boot_sample])
  #    }
  #    lower_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.025)
  #    upper_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.975)
  #  }
  #}
  #SelectionFootprint <- cbind(distance      = ref_points,
  #                            positive_rate = SelectionFootprint,
  #                            upper95       = upper_SelectionFootprint,
  #                            lower95       = lower_SelectionFootprint)
  #saveRDS(SelectionFootprint,
  #        file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  lines(SelectionFootprint[,"distance"],
          SelectionFootprint[,"positive_rate"],
          col=i+1)
  polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
          y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
          col=makeTransparent(i+1,alpha=0.1),
          border=NA)
}
scenarios <- c(40,44,48)
for (i in seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]

  SelectionFootprint <- readRDS(file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))

  #load(file=paste0("results/",simID,"/",simID,"_results.RData"))

  #ref_points <- 10^seq(from=-2,to=1.7,by=0.1)
  #ds <- ref_points*0.8
  
  #SelectionFootprint <- array(NA,dim=length(ref_points))
  #lower_SelectionFootprint <- array(NA,dim=length(ref_points))
  #upper_SelectionFootprint <- array(NA,dim=length(ref_points))
  #for (bin in seq_along(ref_points)){
  #  lower_limit <- max(0,ref_points[bin]-ds[bin])
  #  upper_limit <- ref_points[bin]+ds[bin]
  #  loci_in_bin <- intersect(which(results_per_locus$centimorgan>lower_limit),
  #                           which(results_per_locus$centimorgan<upper_limit))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!results_per_locus$neutral))
  #  
  #  loci_in_bin <- intersect(loci_in_bin,
  #                           which(!is.na(results_per_locus$outliers_topSNPs)))
  #  
  #  if(length(loci_in_bin)>500){
  #    dist_cM  <- results_per_locus$centimorgan[loci_in_bin]
  #    dist_w   <- sqrt( ( dist_cM - ref_points[bin] )^2 )
  #    weight   <- (35/32) * (1 - (dist_w/ds[bin])^2)^3
  #    
  #    positives <- vector("numeric",length(loci_in_bin))
  #    positives[which(results_per_locus$p_value[loci_in_bin] <= 1e-3)] <- 1
  #    SelectionFootprint[bin] <- sum(positives*weight)/sum(weight)
  #    
  #    bootstrap_values <- array(NA,num_of_bootstrap)
  #    for (b in seq_len(num_of_bootstrap)){
  #      boot_sample <- sample(length(loci_in_bin),replace=T)
  #      bootstrap_values[b] <- sum(positives[boot_sample]*weight[boot_sample])/sum(weight[boot_sample])
  #    }
  #    lower_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.025)
  #    upper_SelectionFootprint[bin] <- quantile(bootstrap_values,probs=0.975)
  #  }
  #}
  #SelectionFootprint <- cbind(distance      = ref_points,
  #                            positive_rate = SelectionFootprint,
  #                            upper95       = upper_SelectionFootprint,
  #                            lower95       = lower_SelectionFootprint)
  #saveRDS(SelectionFootprint,
  #        file=paste0("results/",simID,"/",simID,"_SelectionFootprint.RDS"))
  
  lines(SelectionFootprint[,"distance"],
        SelectionFootprint[,"positive_rate"],
        lty=2,
        col=i+1)
  polygon(x=c(SelectionFootprint[,"distance"],rev(SelectionFootprint[,"distance"])),
          y=c(SelectionFootprint[,"lower95"],rev(SelectionFootprint[,"upper95"])),
          col=makeTransparent(i+1,alpha=0.1),
          border=NA)
}





legend(x=2,y=0.24,
       title="Advantageuos allele:",
       legend=c(expression("new mutation (derived, "*pi[0]==1/2*italic(N)*")"),
                expression("ancestral, "*pi[0]==0.1),
                expression("derived, "*pi[0]==0.1),
                expression("ancestral, "*pi[0]==0.5),
                expression("derived, "*pi[0]==0.5),
                expression("ancestral, "*pi[0]==0.9),
                expression("derived, "*pi[0]==0.9)),
       lty=c(2,1),
       #pch=c(1,4,1,4,1,4,1,4), 
       col=c(1,2,2,3,3,4,4),
       cex=0.75,
       bty="n")

dev.off()









