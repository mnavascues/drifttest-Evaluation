#######################################################
#
# Calculating footprint of selection along chromosome
#
#######################################################


# read simulations description table
sim_table <- read.table("Simulations/Fc.simparams",header=T)

#pdf(file="Figures.pdf",paper="a4r")

number_of_replicates <- 100
threshold_power <- 0.05

power_only <- F 

power <- array(NA,nrow(sim_table))
selection_footprint_list <- list()

for (simulation in 1:nrow(sim_table) ){
  
  sim <- sim_table$simID[simulation]

  message(paste("Calculating footprint of selection for scenario",sim,"\n"))
  
  
  q_values_m2    <- array(NA,number_of_replicates)
  distance_cM    <- numeric()
  q_values_m1    <- numeric()
  
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      q_values_m2[replic] <- SNP_list$FC_q_value[which(SNP_list$type=="m2")]

      loci_to_sample <- intersect(which(SNP_list$type!="m2"),intersect(which(!is.na(SNP_list$distance_cM)),which(!is.na(SNP_list$FC_q_value))))
      distance_cM <- c(distance_cM,SNP_list$distance_cM[loci_to_sample])
      q_values_m1  <- c(q_values_m1,SNP_list$FC_q_value[loci_to_sample])
    }
  }
  q_values_m2 <- q_values_m2[which(!is.na(q_values_m2))]
  power[simulation] <- length(which(q_values_m2<threshold_power))/length(q_values_m2)

  distance_cM_ordered <- numeric()
  for (i in c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
              0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
              1,2,3,4,5,6,7,8,9,10,
              11,12,13,14,15,16,17,18,19,
              20,22,24,26,28,
              30,32,34,36,38,
              40,42,44,46,48,
              50,60,70,80,90,
              100,150,200)){
    distance_cM_ordered <- c(distance_cM_ordered , max(distance_cM[intersect( which(distance_cM>(i-1)),  which(distance_cM<=i) )]) )
  }
    
  distance_cM_ordered <- c(min(distance_cM), distance_cM_ordered, max(distance_cM))
  selection_footprint <- array(NA,length(distance_cM_ordered) )
  for (i in 1:length(distance_cM_ordered) ){
    ds     <- 1 # distance_cM_ordered[i]*0.5
    min_d  <- distance_cM_ordered[i]-ds
    max_d  <- distance_cM_ordered[i]+ds
    loci_in_distance_range <- intersect( which(distance_cM>min_d),
                                         which(distance_cM<max_d) )
    positives <- vector("numeric",length(loci_in_distance_range))
    positives[which(q_values_m1[loci_in_distance_range]<0.05)] <- 1
    d_values <- distance_cM[loci_in_distance_range]
    dist   <- sqrt( ( d_values - distance_cM_ordered[i] )^2 )
    weight <- 1 - (dist/ds)^2
    selection_footprint[i] <- glm(positives~d_values,family=binomial,weights=weight)$fitted.values[which(d_values==distance_cM_ordered[i])]
  }
  distance_cM_ordered <- c(0,distance_cM_ordered)
  selection_footprint <- c(power[simulation],selection_footprint)
  
  selection_footprint_list[[length(selection_footprint_list)+1]] <- cbind(distance_cM_ordered,selection_footprint)
  
  
  #save(sim_table,selection_footprint_list,file="Results/SelectionFootprint.RData")
  
}

sim_table <- cbind(sim_table,power)

save(sim_table,selection_footprint_list,file="Results/SelectionFootprint.RData")

















## Plot power


#load(file="Results/SelectionFootprint.RData")

#color blind palette
require(graphics)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette(cbPalette)

selection_mode <- "NM"
selection_period_duration <- 25
simulations2plot_total <- intersect(which(sim_table$selection_period_duration==selection_period_duration),
                                    which(sim_table$selection_mode==selection_mode))
sel_coef               <- unique(sim_table$sel_coef[simulations2plot_total])

#pdf_file_name <- paste0("Results/Power_dT",selection_period_duration,"_",selection_mode,".pdf")
#pdf(file=pdf_file_name)
#jpg_file_name <- paste0("Results/Power_dT",selection_period_duration,"_",selection_mode,".jpg")
#jpeg(file=jpg_file_name)



simulations2plot <- intersect(simulations2plot_total,which(sim_table$sel_coef==sel_coef[1]))
sigma_values     <- sim_table$sigma[simulations2plot]
power            <- sim_table$power[simulations2plot]


plot(sigma_values,power,
     type="l",
     col=1,
     ylim=c(0,1),
     ylab = expression("Power, "*italic(W)),
     xlab = expression("Selfing rate, "*sigma),
     main=substitute( italic(N)*"="*500*", "*italic(t)*"="*dt  ,list(dt=selection_period_duration)) ,
     lwd=3)
for (i in 2:length(sel_coef)){
  simulations2plot <- intersect(simulations2plot_total,which(sim_table$sel_coef==sel_coef[i]))
  sigma_values     <- sim_table$sigma[simulations2plot]
  power            <- sim_table$power[simulations2plot]
  lines(sigma_values,power,col=i,lwd=3)
}
if(selection_mode=="NM") texto<-"New mutation"
if(selection_mode=="SV") texto<-"Standing variation"
text(0.5,0.5,adj=0.5,texto,cex=3,col="grey")
#dev.off()



plot(sigma_values,power_SV,
     type="l",
     col=1,
     ylim=c(0,1),
     ylab = "Proportion of positives (q-value<0.05)",
     xlab = expression("Selfing rate, "*sigma),
     main = "Adaptation from standing variation",
     lwd=3)
lines(sigma_values,FDR_SV,lwd=3,col="grey")


plot(sigma_values,power_NM,
     type="l",
     col=1,
     ylim=c(0,1),
     ylab = "Proportion of positives (q-value<0.05)",
     xlab = expression("Selfing rate, "*sigma),
     main = "Adaptation from new mutation",
     lwd=3)
lines(sigma_values,FDR_NM,lwd=3,col="grey")






## Plot Selection footprint in funtion of selection coefficient (slection strength)
require(graphics)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette(cbPalette)

selection_period_duration <- 25
sigma_value               <- 1
simulations2plot          <- intersect(which(sim_table$sigma==sigma_value),which(sim_table$selection_period_duration==selection_period_duration))
sel_coef                  <- sim_table$sel_coef[simulations2plot]

#pdf_file_name <- paste0("Results/SelectionFootprint_sigma",sigma_value,"_dT",selection_period_duration,".pdf")
#pdf(file=pdf_file_name)
#jpg_file_name <- paste0("Results/SelectionFootprint_sigma",sigma_value,"_dT",selection_period_duration,".jpg")
#jpeg(file=jpg_file_name)



single_simulation2plot <- intersect(simulations2plot,which(sim_table$sel_coef==sel_coef[1]))
distancia        <- selection_footprint_list[[single_simulation2plot]][,1]
huella_selectiva <- selection_footprint_list[[single_simulation2plot]][,2]
plot(distancia,huella_selectiva,
     type="l",
     col=1,
     ylim=c(0,1),
     xlim=c(0,20),
     ylab = "Selection footprint",
     xlab = "Distance (cM)",
     main=substitute( expression(italic(N)*"="*500*", "*sigma*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sigma_value,dt=selection_period_duration)) ,
     lwd=3)
for (i in 2:length(sel_coef)){
  single_simulation2plot <- intersect(simulations2plot,which(sim_table$sel_coef==sel_coef[i]))
  distancia        <- selection_footprint_list[[single_simulation2plot]][,1]
  huella_selectiva <- selection_footprint_list[[single_simulation2plot]][,2]
  lines(distancia,huella_selectiva,col=i,lwd=3)
}
#dev.off()





























## Plot Selection footprint in funtion of selection coefficient (time)
require(graphics)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette(cbPalette)

sel_coef                  <- 0.5
sigma_value               <- 0.9
simulations2plot          <- intersect(which(sim_table$sigma==sigma_value),which(sim_table$sel_coef==sel_coef))
selection_period_duration <- sim_table$selection_period_duration[simulations2plot]

pdf_file_name <- paste0("SelectionFootprint_sigma",sigma_value,"_s",sel_coef,".pdf")
pdf(file=pdf_file_name)



single_simulation2plot <- intersect(simulations2plot,which(sim_table$selection_period_duration==selection_period_duration[1]))
distancia        <- selection_footprint_list[[single_simulation2plot]][,1]
huella_selectiva <- selection_footprint_list[[single_simulation2plot]][,2]
plot(distancia,huella_selectiva,
     type="l",
     col=1,
     ylim=c(0,1),
     xlim=c(0,20),
     ylab = "Selection footprint",
     xlab = "Distance (cM)",
     main=substitute( expression(italic(N)*"="*500*", "*sigma*"="*ss*", "*italic(s)*"="*select ) ,list(ss=sigma_value,select=sel_coef)) ,
     lwd=3)
for (i in 2:length(selection_period_duration)){
  single_simulation2plot <- intersect(simulations2plot,which(sim_table$selection_period_duration==selection_period_duration[i]))
  distancia        <- selection_footprint_list[[single_simulation2plot]][,1]
  huella_selectiva <- selection_footprint_list[[single_simulation2plot]][,2]
  lines(distancia,huella_selectiva,col=i,lwd=3)
}
dev.off()








## Plot Selection footprint in funtion of selection coefficient (selfing)
require(graphics)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette(cbPalette)

selection_mode <- "SV"
sel_coef                  <- 0.5
selection_period_duration <- 25
simulations2plot          <- intersect(which(sim_table$selection_period_duration==selection_period_duration),
                                       which(sim_table$sel_coef==sel_coef))
simulations2plot          <- intersect(simulations2plot,
                                       which(sim_table$selection_mode==selection_mode))
sigma_values               <- sim_table$sigma[simulations2plot]

#pdf_file_name <- paste0("Results/SelectionFootprint_s",sel_coef,"_dT",selection_period_duration,"_",selection_mode,".pdf")
#pdf(file=pdf_file_name)



single_simulation2plot <- intersect(simulations2plot,which(sim_table$sigma==sigma_values[1]))
distancia        <- selection_footprint_list[[single_simulation2plot]][,1]
huella_selectiva <- selection_footprint_list[[single_simulation2plot]][,2]
plot(distancia[-1],huella_selectiva[-1],
     type="l",
     col=1,
     ylim=c(0,0.15),
     xlim=c(0,30),
     ylab = "Probability of positive test (with q-value<0.05)",
     xlab = "Distance (cM) from locus under selection",
     #main=substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) ,
     lwd=3)
for (i in 2:7){
  single_simulation2plot <- intersect(simulations2plot,which(sim_table$sigma==sigma_values[i]))
  distancia        <- selection_footprint_list[[single_simulation2plot]][,1]
  huella_selectiva <- selection_footprint_list[[single_simulation2plot]][,2]
  lines(distancia[-1],huella_selectiva[-1],col=i,lwd=3)
}
#dev.off()
legend(10,0.13,
       legend=c(expression(sigma*"=0.00"),
                expression(sigma*"=0.50"),
                expression(sigma*"=0.75"),
                expression(sigma*"=0.90"),
                expression(sigma*"=0.95"),
                expression(sigma*"=0.99"),
                expression(sigma*"=1.00")),
       col=1:7,
       lwd=3,bty="n")


