#######################################################
#
# ManhattanPlot
#
#######################################################


# read simulations description table
sim_table <- read.table("Simulations/Fc.simparams",header=T)


number_of_replicates <- 100
population_size <- 500
selection_mode <- "NM"
sel_coef <- 0.5
sigma <- 0.0
selection_period_duration <- 25

simulation2plot <- intersect(which(sim_table$sel_coef==sel_coef),
                              which(sim_table$selection_period_duration==selection_period_duration))
simulation2plot <- intersect(simulation2plot,
                              which(sim_table$selection_mode==selection_mode))
simulation2plot <- intersect(simulation2plot,
                             which(sim_table$sigma==sigma))

sim <- sim_table$simID[simulation2plot]
 

for (replic in 1:number_of_replicates){
  success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
  
  if(class(success)=="try-error"){
    cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
  }else{

    selected_site <- which(SNP_list$type=="m2")
    selected_site_position <- SNP_list$x[selected_site]
    
    selected_region <- which(SNP_list$distance_cM<=10)
    
    
    
    plot(SNP_list$x[which(SNP_list$x>5e+8/2)],
         -log10(SNP_list$FC_q_value)[which(SNP_list$x>5e+8/2)],
         xlim=c(0,5e+8),
         ylim=c(0,2),
         xlab="",
         main=replic,
         ylab=expression(-log[10]*"(q-value)"),
         axes=F,
         pch=16,
         #cex=0.5,
         col="dark grey")
    Axis(side=2)
    axis(side=1,at=c(5e+8/4,5e+8*3/4),labels=c("chr 1","chr 2"))
    box()
    points(SNP_list$x[which(SNP_list$x<=5e+8/2)],
           -log10(SNP_list$FC_q_value)[which(SNP_list$x<=5e+8/2)],
           pch=16,
           #cex=0.5,
           col="black")
    abline(h=-log10(0.05),col="blue")
    abline(v=selected_site_position + 10 / (1e-8 * 100),col="red")
    abline(v=selected_site_position - 10 / (1e-8 * 100),col="red")
    
    #points(SNP_list$x[selected_region],
    #       -log10(SNP_list$FC_q_value)[selected_region],
    #       pch=16,
    #       #cex=0.5,
    #       col="red")
    
    if (SNP_list$FC_q_value[selected_site]==0){
      points(SNP_list$x[selected_site],
             2,
             pch=18,
             #cex=0.6,
             col="red")
    }else{
      points(SNP_list$x[selected_site],
             -log10(SNP_list$FC_q_value[selected_site]),
             pch=18,
             #cex=0.6,
             col="red")
    }
    
    
    
    
  }
}

replic<- 93






for (simulation in 1:nrow(sim_table) ){
  sim <- sim_table$simID[simulation]
  message(paste("Calculating footprint of selection for scenario",sim,"\n"))
  

  expected_Ne[simulation] <- (2-sim_table$sigma[simulation])*N
}

sim_table <- cbind(sim_table,expected_Ne)

save(sim_table,Ne_estimates,file="Results/EffectivePopSize.RData")















#load(file="Results/EffectivePopSize.RData")
source("percentile.boxplot.R")


# Make figure in function of selfing rate

selection_mode <- "SV"
sel_coef                  <- 0.0
selection_period_duration <- 25
simulations2plot <- intersect(which(sim_table$sel_coef==sel_coef),
                              which(sim_table$selection_period_duration==selection_period_duration))
simulations2plot <- intersect(simulations2plot,
                              which(sim_table$selection_mode==selection_mode))
sigma_values              <- sim_table$sigma[simulations2plot]

#pdf_file_name <- paste0("Results/EffectivePopSize_s",sel_coef,"_dT",selection_period_duration,"_",selection_mode,".pdf")
#pdf(file=pdf_file_name)
#jpg_file_name <- paste0("Results/EffectivePopSize_s",sel_coef,"_dT",selection_period_duration,"_",selection_mode,".jpg")
#jpeg(file=jpg_file_name)


boxplot(  Ne_estimates[,simulations2plot[1]],
          at=sigma_values[1],
          outline=FALSE,
          boxwex=0.02,
          names=sigma_values[1],
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = expression("selfing rate, "*sigma),
          ylim = c(0,1500),
          xlim = c(min(sigma_values),max(sigma_values)),
          main = "Neutral scenario")#substitute( expression(italic(N)*"="*500*", "*italic(s)*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sel_coef,dt=selection_period_duration)) )

for (i in 2:length(sigma_values)){
  boxplot(  Ne_estimates[,simulations2plot[i]],
            at=sigma_values[i],
            outline=FALSE,
            axes = FALSE,
            boxwex=0.02,
            names=sigma_values[i],
            add=T)
}
lines(sigma_values,sim_table$expected_Ne[simulations2plot],col="red")
axis(side=1,at=sigma_values)
if(selection_mode=="NM") texto<-"New mutation"
if(selection_mode=="SV") texto<-"Standing variation"
text(0.5,1200,adj=0.5,texto,cex=3,col="grey")
#dev.off()




















# Make figure in function of selection coefficient

selection_mode <- "NM"
sigma_value               <- 0.0
selection_period_duration <- 25
simulations2plot          <- intersect(which(sim_table$sigma==sigma_value),
                                       which(sim_table$selection_period_duration==selection_period_duration))
simulations2plot          <- intersect(simulations2plot,
                                       which(sim_table$selection_mode==selection_mode))
sel_coef                  <- sim_table$sel_coef[simulations2plot]

#pdf_file_name <- paste0("Results/EffectivePopSize_sigma",sigma_value,"_dT",selection_period_duration,"_",selection_mode,".pdf")
#pdf(file=pdf_file_name)
jpg_file_name <- paste0("Results/EffectivePopSize_sigma",sigma_value,"_dT",selection_period_duration,"_",selection_mode,".jpg")
jpeg(file=jpg_file_name)


boxplot(  Ne_estimates[,simulations2plot[1]],
          at=1,
          outline=FALSE,
          names=sel_coef[1],
          ylab = expression("effective population size, "*2*italic(Ne)),
          xlab = expression("selection coefficient, "*italic(s)),
          ylim = c(0,1500),
          xlim = c(0.5,length(sel_coef)+0.5),
          main = substitute( expression(italic(N)*"="*500*", "*sigma*"="*ss*", "*italic(t)*"="*dt ) ,list(ss=sigma_value,dt=selection_period_duration)) )

for (i in 2:length(sel_coef)){
  boxplot(  Ne_estimates[,simulations2plot[i]],
            at=i,
            outline=FALSE,
            axes = FALSE,
            names=sel_coef[i],
            add=T)
}
axis(side=1,at=1:length(sel_coef),labels=sel_coef)
abline(h=sim_table$expected_Ne[simulations2plot][1],col="red")
if(selection_mode=="NM") texto<-"New mutation"
if(selection_mode=="SV") texto<-"Standing variation"
text(2.5,1200,adj=0.5,texto,cex=3,col="grey")
dev.off()





