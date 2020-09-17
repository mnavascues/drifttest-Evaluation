load(file="results/simtable.RData")
    
source("src/fun/F_stats_tools.R")
      
source("src/fun/ColorBlindPalette.R")
palette(cbPalette1)
    
          
number_of_replicates<-100
scenarios_ancestral <- c(31,35,39)
scenarios_derived   <- c(40,44,48)
      
scenarios <- c(scenarios_ancestral,scenarios_derived)


pdf(file="results/StandingVariationDiversity.pdf",width=4,height=6)
layout(matrix(c(1,2), 2,1,byrow = TRUE))
par(mar=c(4,4,0.2,0.2)+0.1)
seq_distances_cM <- readRDS(file=paste0("results/distancesCM.RDS"))
plot(seq_distances_cM,
     rep(1,length(seq_distances_cM)),
     log="x",
     type="l",
     ylim=c(0,0.8),
     lty=2,
     col="white",
     ylab="Differentiation",
     xlab="",
     cex.axis=0.7)
text(x=50,y=0.78,label="a",cex=1.5)


for (i in 1:3){#seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]
  
  mean_FST       <- readRDS(paste0("results/",simID,"/mean_FST.RDS"))
  #sd_FST <- readRDS(file=paste0("results/",simID,"/sd_FST.RDS"))
  
  #polygon(x=c(seq_distances_cM,rev(seq_distances_cM)),
  #        y=c(mean_FST+sd_FST,rev(mean_FST-sd_FST)),
  #        col=makeTransparent(i+1,alpha=0.1),
  #        border=NA)  
  lines(seq_distances_cM,
        mean_FST,
        col=i+1)
  
}



for (i in 4:6){#seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]
  
  mean_FST       <- readRDS(paste0("results/",simID,"/mean_FST.RDS"))
  #sd_FST <- readRDS(file=paste0("results/",simID,"/sd_FST.RDS"))
  
  #polygon(x=c(seq_distances_cM,rev(seq_distances_cM)),
  #        y=c(mean_FST+sd_FST,rev(mean_FST-sd_FST)),
  #        col=makeTransparent(i-2,alpha=0.1),
  #        border=NA)  
  lines(seq_distances_cM,
        mean_FST,
        col=i-2,
        lty=2)
  
}  






#distances <-readRDS(file=paste0("results/distances.RDS"))
plot(seq_distances_cM,
     rep(1,length(seq_distances_cM)),
     log="x",
     type="l",
     ylim=c(0,3e-6),
     lty=2,
     col="white",
     ylab="Diversity",
     xlab="Distance (cM)",
     cex.axis=0.7)
text(x=50,y=2.9e-6,label="b",cex=1.5)


for (i in 1:3){#seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]
  
  mean_H <- readRDS(paste0("results/",simID,"/mean_H.RDS"))
  #sd_H   <- readRDS(file=paste0("results/",simID,"/sd_H.RDS"))
  
  #polygon(x=c(seq_distances_cM,rev(seq_distances_cM)),
  #        y=c(mean_H+sd_H,rev(mean_H-sd_H)),
  #        col=makeTransparent(i+1,alpha=0.1),
  #        border=NA)  
  lines(seq_distances_cM,
        mean_H,
        col=i+1)
  
}



for (i in 4:6){#seq_along(scenarios)){
  simID <- sim_table$simID[scenarios[i]]
  
  mean_H <- readRDS(paste0("results/",simID,"/mean_H.RDS"))
  #sd_H   <- readRDS(file=paste0("results/",simID,"/sd_H.RDS"))
  
  #polygon(x=c(seq_distances_cM,rev(seq_distances_cM)),
  #        y=c(mean_H+sd_H,rev(mean_H-sd_H)),
  #        col=makeTransparent(i-2,alpha=0.1),
  #        border=NA)  
  lines(seq_distances_cM,
        mean_H,
        col=i-2,
        lty=2)
  
}  

legend(x=0.03,y=3e-6,
       title="Advantageuos allele:",
       legend=c(expression("ancestral, "*pi[0]==0.1),
                expression("derived, "*pi[0]==0.1),
                expression("ancestral, "*pi[0]==0.5),
                expression("derived, "*pi[0]==0.5),
                expression("ancestral, "*pi[0]==0.9),
                expression("derived, "*pi[0]==0.9)),
       lty=c(1,2),
       #pch=c(1,4,1,4,1,4,1,4), 
       col=c(2,2,3,3,4,4),
       cex=0.5,
       bty="n")

dev.off()




######### DO NOT RUN FROM HERE










      
genome_length <- 5e8
midpoint <- genome_length/2
source("src/fun/manipulate_output.R")
      
window_size <- 100000
seq_distances_bp <- 3*10^seq(4,7,0.2) 
seq_distances_cM <- bp2centimorgan(seq_distances_bp,1e-8)
saveRDS(seq_distances_bp,file=paste0("results/distancesBP.RDS"))
saveRDS(seq_distances_cM,file=paste0("results/distancesCM.RDS"))
      
      

for (i in seq_along(scenarios)){
  ancestral_advantageus <- FALSE
  if (any(scenarios[i]==scenarios_ancestral)) ancestral_advantageus <- TRUE
        
  simID <- sim_table$simID[scenarios[i]]
      
  FST_per_window <- matrix(data=NA,nrow=number_of_replicates,ncol=length(seq_distances_bp))
  H_per_window <- matrix(data=NA,nrow=number_of_replicates,ncol=length(seq_distances_bp))
  for (replicate in seq_len(number_of_replicates)){
    load(file=paste0("results/",simID,"/replicates/",simID,"_",replicate,"_whole_pop.RData"))
        
    if (m2$x<midpoint) {selected_chr <- 1} else {selected_chr <- 2}
      
    if (selected_chr==1){
      kept_loci <- whole_pop_data$SNP_table$x<midpoint
    }else if (selected_chr==2){
      kept_loci <- whole_pop_data$SNP_table$x>=midpoint
    }
          
    hap <- whole_pop_data$sampled_haplotypes_1[,kept_loci]
    SNP_table <- whole_pop_data$SNP_table[kept_loci,]
    distance_bp <- abs(SNP_table$x-m2$x)
    distance_cM <- bp2centimorgan(distance_bp,1e-8)
    SNP_table <- cbind(SNP_table,distance_bp,distance_cM)
          
          
    #dim(hap)
    position_under_selection <- which(SNP_table$type=="m2")
    if (ancestral_advantageus){
      haplotypes_to_keep <- which(hap[,position_under_selection]==0)
    }else{
      haplotypes_to_keep <- which(hap[,position_under_selection]==1)
    }
    selected_hap   <- hap[haplotypes_to_keep,]
    n<-nrow(selected_hap)
    unselected_hap <- hap[-haplotypes_to_keep,]
    
    selected_hap_counts <- apply(selected_hap,2,sum)
    selected_hap_counts <- cbind(derived=selected_hap_counts,
                                 ancestral=nrow(selected_hap)-selected_hap_counts)
    unselected_hap_counts <- apply(unselected_hap,2,sum)
    unselected_hap_counts <- cbind(derived=unselected_hap_counts,
                                 ancestral=nrow(unselected_hap)-unselected_hap_counts)
    
    counts<-cbind(selected_hap_counts,unselected_hap_counts)

    for (d in seq_along(seq_distances_bp)){
      loci_on_haplotype <- which(SNP_table$distance_bp>seq_distances_bp[d] &
                                 SNP_table$distance_bp<seq_distances_bp[d]+window_size)
      set_higher <- intersect(loci_on_haplotype,
                              which(SNP_table$x>SNP_table$x[position_under_selection]))
      set_lower <- intersect(loci_on_haplotype,
                             which(SNP_table$x<SNP_table$x[position_under_selection]))
      if (length(set_higher)!=0){
        loci_on_haplotype <- set_higher
      }else if(length(set_lower)!=0){
        loci_on_haplotype <- set_lower
      }
      
      if (length(loci_on_haplotype)>0){
        if (length(loci_on_haplotype)==1){
          FST_per_window[replicate,d] <- single_locus_compute_FST_from_counts(matrix(counts[loci_on_haplotype,],ncol=4))
        }else{
          FST_per_window[replicate,d] <- compute_FST_from_counts(counts[loci_on_haplotype,])
        }
        freqs<-selected_hap_counts[loci_on_haplotype,"derived"]/n
        freqs <- freqs[freqs!=0]
        H_per_window[replicate,d] <- mean((1-(freqs^2+(1-freqs)^2))*n/(n-1))/window_size
        
      }else{
        FST_per_window[replicate,d] <- 0
        H_per_window[replicate,d]   <- 0 
      }
    }
    
  }
  sd_FST   <- apply(FST_per_window, 2, sd,   na.rm=T)
  mean_FST <- apply(FST_per_window, 2, mean, na.rm=T) 
  sd_H     <- apply(H_per_window,   2, sd,   na.rm=T)
  mean_H   <- apply(H_per_window,   2, mean, na.rm=T) 
  
  saveRDS(mean_FST,file=paste0("results/",simID,"/mean_FST.RDS"))
  saveRDS(sd_FST,file=paste0("results/",simID,"/sd_FST.RDS"))
  saveRDS(mean_H,file=paste0("results/",simID,"/mean_H.RDS"))
  saveRDS(sd_H,file=paste0("results/",simID,"/sd_H.RDS"))

}

