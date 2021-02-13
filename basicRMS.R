#basicRMSD of multi state PDB. Maybe for a GROMACS/NAMD run converted to pdb or NMR ensemble.
#To do.(1) Extend colouring options, styles in function for plot 
#      (2) Improve speed. Make pdb_single load first frame instead of whole file  
#      (3) Extend flexibility 
#      (4) Hard-coded needs fixing
#NOTE: First frame is static

#smoothNg = function(directory = "/home/dremyes/WaitingRoom/trjs", 
#                    trajFile = "Ng10nsIII.pdb",
#                    outFile  = "Ng10nsIIIMulti.png",
#                    plotTitle = "4k3b Ng 10ns Run III", 
#                    simulationLength = 10){
  directory = "/home/dremyes/WaitingRoom/trjs/cjbamd"
  
  #####################ITASS##############################  
  statFile = "BamDITass10ns.pdb"
  trajFile = "BamDITass10ns.pdb"
  outFile  = "BamDITass10nsRMSD.png"
  plotTitle = "BamD_I-TASSER, 10ns Run span = .01, (0ns Reference Frame)" 
  simulationLength = 10
  
  
  
  
  
  library(bio3d)
  library(ggplot2)
  setwd <- setwd(directory)
  pdb_single <-read.pdb(statFile)
  #All frames 
  pdb_multi <- read.pdb(trajFile, multi = TRUE)


  ITASSER.fixed.ca.inds <- atom.select(pdb_single, elety="CA")
  ITASSER.mobile.ca.inds <- atom.select(pdb_multi, elety="CA")
  ITASSER.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = ITASSER.fixed.ca.inds$xyz, mobile.inds = ITASSER.mobile.ca.inds$xyz)
  ITASSER.rmsd    <- rmsd(ITASSER.xyz[1,ITASSER.fixed.ca.inds$xyz], ITASSER.xyz[,ITASSER.mobile.ca.inds$xyz])
  ITASSER.rmsData <- data.frame(ITASSER.rmsd)
  ITASSER.rmsData[]
  ########Harmonic Freq' Analysis#########
         write.csv(traj, "happy.csv") 
  traj <- read.csv("happy.csv")
  ########trend###########################
         plot(traj)
  trend <- lm(RMSD ~ (Frame), data = traj)
       abline(trend, col="red")
  ########Detrended########################
  detrended.trajectory <- trend$residuals
  plot(detrended.trajectory, type="l", main="Detrended time series")  

  ########HHARMONICS######################
  f.data <- GeneCycle::periodogram(detrended.trajectory)
  harmonics <- 1:20 
  plot(f.data$freq[harmonics]*length(detrended.trajectory), 
       f.data$spec[harmonics]/sum(f.data$spec), 
       xlab="Harmonics (Hz)", ylab="Amplitute Density", type="h")
  
  
       
  #All the RMSD vector lists have been created, assign a time POTRA domain id to the rmsd and rbind them into a long list. Alt method in column df below.
  #Add the time variable as well
  numberFrames = length(I-TASSER.rmsd)
  timeFactor = simulationLength/numberFrames
  print(timeFactor)
  
  #Add time column
  #rmsData$time <- (1:length(I-TASSER.rmsd))*(timeFactor)
  I.TASSER.rmsData$time <- (1:length(I-TASSER.rmsd))*(timeFactor)
  #Generate lists for the names
  protein1 <- rep("I-TASSER", length(I-TASSER.rmsd))

  #Tack them onto each dataframe
  I.TASSER.rmsData$domain <- as.factor(protein1)
  
  #Make column names all the same 'rmsd'

  colnames(I-TASSER.rmsData) <- c("rmsd", "time", "protein")
  I.TASSER.rmsData[]


  
  
  
  plotIt <- ggplot(data = I-TASSER.rmsData, aes(x=time, y=rmsd)) + 
    stat_smooth(span = 0.01) +
    geom_point(aes(colour = rmsd)) +
    scale_colour_gradient2(low = "blue", mid = "green", high ="red", midpoint = 1.8) +
    facet_grid(protein ~ .) +
    ggtitle(plotTitle) +
    xlab("Time (ns)") + ylab("RMSD (Angstrom)")
  
  ggplot_build(plotIt)
  ggsave(outFile)
  
  print(plotIt)

#####################PHYRE2##############################  
  statFile = "BamDPhyre10ns.pdb"
  trajFile = "BamDPhyre10ns.pdb"
  outFile  = "BamDPhyre10nsRMSD.png"
  plotTitle = "BamD_Phyre2, 10ns Run span = .01, (0ns Reference Frame)" 
  simulationLength = 10
  
  library(bio3d)
  library(ggplot2)
  setwd <- setwd(directory)
  pdb_single <-read.pdb(statFile)
  #All frames 
  pdb_multi <- read.pdb(trajFile, multi = TRUE)
  
  
  PHYRE2.fixed.ca.inds <- atom.select(pdb_single, elety="CA")
  PHYRE2.mobile.ca.inds <- atom.select(pdb_multi, elety="CA")
  PHYRE2.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = PHYRE2.fixed.ca.inds$xyz, mobile.inds = PHYRE2.mobile.ca.inds$xyz)
  PHYRE2.rmsd    <- rmsd(PHYRE2.xyz[1,PHYRE2.fixed.ca.inds$xyz], PHYRE2.xyz[,PHYRE2.mobile.ca.inds$xyz])
  PHYRE2.rmsData <- data.frame(PHYRE2.rmsd)
  PHYRE2.rmsData[]
  
  #All the RMSD vector lists have been created, assign a time POTRA domain id to the rmsd and rbind them into a long list. Alt method in column df below.
  #Add the time variable as well
  numberFrames = length(PHYRE2.rmsd)
  timeFactor = simulationLength/numberFrames
  print(timeFactor)
  
  #Add time column
  #rmsData$time <- (1:length(PHYRE2.rmsd))*(timeFactor)
  PHYRE2.rmsData$time <- (1:length(PHYRE2.rmsd))*(timeFactor)
  
  #Generate lists for the names
  protein1 <- rep("PHYRE2", length(PHYRE2.rmsd))
  
  #Tack them onto each dataframe
  PHYRE2.rmsData$domain <- as.factor(protein1)
  
  #Make column names all the same 'rmsd'
  
  colnames(PHYRE2.rmsData) <- c("rmsd", "time", "protein")
  PHYRE2.rmsData[]
  
  
  
  plotIt <- ggplot(data = PHYRE2.rmsData, aes(x=time, y=rmsd)) + 
    stat_smooth(span = 0.01) +
    geom_point(aes(colour = rmsd)) +
    scale_colour_gradient2(low = "blue", mid = "green", high ="red", midpoint = 1.8) +
    facet_grid(protein ~ .) +
    ggtitle(plotTitle) +
    xlab("Time (ns)") + ylab("RMSD (Angstrom)")
  
  ggplot_build(plotIt)
  ggsave(outFile)
  
  print(plotIt)
  
  
  #####################MODELLER##############################  
  statFile = "BamDMOD10ns.pdb"
  trajFile = "BamDMOD10ns.pdb"
  outFile  = "BamDMOD10nsRMSD.png"
  plotTitle = "BamD_MODELLER, 10ns Run span = .01, (0ns Reference Frame)" 
  simulationLength = 10
  
  library(bio3d)
  library(ggplot2)
  setwd <- setwd(directory)
  pdb_single <-read.pdb(statFile)
  #All frames 
  pdb_multi <- read.pdb(trajFile, multi = TRUE)
  
  
  MODELLER.fixed.ca.inds <- atom.select(pdb_single, elety="CA")
  MODELLER.mobile.ca.inds <- atom.select(pdb_multi, elety="CA")
  MODELLER.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = MODELLER.fixed.ca.inds$xyz, mobile.inds = MODELLER.mobile.ca.inds$xyz)
  MODELLER.rmsd    <- rmsd(MODELLER.xyz[1,MODELLER.fixed.ca.inds$xyz], MODELLER.xyz[,MODELLER.mobile.ca.inds$xyz])
  MODELLER.rmsData <- data.frame(MODELLER.rmsd)
  MODELLER.rmsData[]
  
  #All the RMSD vector lists have been created, assign a time POTRA domain id to the rmsd and rbind them into a long list. Alt method in column df below.
  #Add the time variable as well
  numberFrames = length(MODELLER.rmsd)
  timeFactor = simulationLength/numberFrames
  print(timeFactor)
  
  #Add time column
  #rmsData$time <- (1:length(MODELLER.rmsd))*(timeFactor)
  MODELLER.rmsData$time <- (1:length(MODELLER.rmsd))*(timeFactor)
  
  #Generate lists for the names
  protein1 <- rep("MODELLER", length(MODELLER.rmsd))
  
  #Tack them onto each dataframe
  MODELLER.rmsData$domain <- as.factor(protein1)
  
  #Make column names all the same 'rmsd'
  
  colnames(MODELLER.rmsDataMOD) <- c("rmsd", "time", "protein")
  MODELLER.rmsDataMOD[]
  
  
  
  plotIt <- ggplot(data = MODELLER.rmsData, aes(x=time, y=rmsd)) + 
    stat_smooth(span = 0.01) +
    geom_point(aes(colour = rmsd)) +
    scale_colour_gradient2(low = "blue", mid = "green", high ="red", midpoint = 1.8) +
    facet_grid(protein ~ .) +
    ggtitle(plotTitle) +
    xlab("Time (ns)") + ylab("RMSD (Angstrom)")
  
  ggplot_build(plotIt)
  ggsave(outFile)
  
  print(plotIt)
  
  #####################PHYRE2##############################  
  #written by MJHardcastle, Uni of Nottingam 2016 Copyright BonevWilliamsLab
  #For help running any of this code please contact me at
  #stxmjh@nottingham.ac.uk  
  
