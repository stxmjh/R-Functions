#Megaplot. Calculates RMSD for individual elements of the protein.
#N.gono VERSION
library(bio3d)
library(ggplot2)
#First frame

smoothCj = function(
  directory = "/home/dremyes/WaitingRoom/trjs/namdsucks"
                    trajFile = "catVMD.pdb"
                    outFile  = "CjBamA111nsIMulti.png"
                    plotTitle = "I-TASSER Cj 111ns Run I"
                    simulationLength = 111
  library(bio3d)
  library(ggplot2)
  setwd <- setwd(directory)
  pdb_single <-read.pdb(trajFile)
  #All frames 
  pdb_multi <- read.pdb(trajFile, multi = TRUE)
  
  #Set POTRA Domain Elements. NOTE: Must be selected manually until rule is solved.
  POTRA1 <- atom.select(pdb_single, resno =  45:113)
  POTRA2 <- atom.select(pdb_single, resno = 114:193)
  POTRA3 <- atom.select(pdb_single, resno = 194:286)
  POTRA4 <- atom.select(pdb_single, resno = 287:370)
  POTRA5  <- atom.select(pdb_single, resno = 371:445)
  BARREL  <- atom.select(pdb_single, resno = 446:766)
  
  
  #Alternative test 1)define set 
  
  #P1.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  45:113)
  #P1.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  45:113)
  #P1.rmsd <- rmsd(pdb_single, pdb_multi, a.inds = P1.fixed.ca.inds, b.inds = P1.mobile.ca.inds, fit =TRUE)
  atoms.fixed.ca.inds  <- atom.select(pdb_single, elety="CA")
  atoms.mobile.ca.inds <- atom.select(pdb_multi, elety="CA")
  atoms.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = atoms.fixed.ca.inds$xyz, mobile.inds = atoms.mobile.ca.inds$xyz)
  atoms.rmsd    <- rmsd(atoms.xyz[1,atoms.fixed.ca.inds$xyz], atoms.xyz[,atoms.mobile.ca.inds$xyz])
  atoms.rmsData <- data.frame(atoms.rmsd)
  traj          <- atoms.rmsData
  ########Harmonic Freq' Analysis########################
  write.csv(traj, "rmsd.csv") 
  traj <- read.csv("rmsd.csv")
  ########trend##########################################
  pdf('trend.pdf')
  plot(traj)
  trend <- lm(RMSD ~ (Frame), data = traj)
  abline(trend, col="red")
  dev.off()
  ########Detrended########################
  detrended.trajectory <- trend$residuals
  pdf('detrend.pdf')
  plot(detrended.trajectory, type="l", main="Detrended time series")  
  dev.off()
  ########HHARMONICS######################
  f.data <- GeneCycle::periodogram(detrended.trajectory)
  harmonics <- 1:20 
  pdf('harmonics.pdf')
  plot(f.data$freq[harmonics]*length(detrended.trajectory), 
       f.data$spec[harmonics]/sum(f.data$spec), 
       xlab="Harmonics (Hz)", ylab="Amplitute Density", type="h")
  dev.off()
 write.csv(harmonics, "harmoData111.csv")
  
  
  #ca.inds <- atom.select(pdb_single, elety="CA", resno=45:113)
  #ca.inds    <- atom.select(pdb_single, elety ="CA")
  
  P1.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  45:113)
  P1.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  45:113)
  P1.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = P1.fixed.ca.inds$xyz, mobile.inds = P1.mobile.ca.inds$xyz)
  P1.rmsd    <- rmsd(P1.xyz[1,P1.fixed.ca.inds$xyz], P1.xyz[,P1.mobile.ca.inds$xyz])
  P1.rmsData <- data.frame(P1.rmsd)
  P1.rmsData[]
  
  P2.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  114:193)
  P2.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  114:193)
  P2.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = P2.fixed.ca.inds$xyz, mobile.inds = P2.mobile.ca.inds$xyz)
  P2.rmsd    <- rmsd(P2.xyz[1,P2.fixed.ca.inds$xyz], P2.xyz[,P2.mobile.ca.inds$xyz])
  P2.rmsData <- data.frame(P2.rmsd)
  P2.rmsData[]
  
  P3.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  194:286)
  P3.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  194:286)
  P3.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = P3.fixed.ca.inds$xyz, mobile.inds = P3.mobile.ca.inds$xyz)
  P3.rmsd    <- rmsd(P3.xyz[1,P3.fixed.ca.inds$xyz], P3.xyz[,P3.mobile.ca.inds$xyz])
  P3.rmsData <- data.frame(P3.rmsd)
  P3.rmsData[]
  
  P4.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  287:370)
  P4.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  287:370)
  P4.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = P4.fixed.ca.inds$xyz, mobile.inds = P4.mobile.ca.inds$xyz)
  P4.rmsd    <- rmsd(P4.xyz[1,P4.fixed.ca.inds$xyz], P4.xyz[,P4.mobile.ca.inds$xyz])
  P4.rmsData <- data.frame(P4.rmsd)
  P4.rmsData[]
  
  P5.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  371:445)
  P5.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  371:445)
  P5.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = P5.fixed.ca.inds$xyz, mobile.inds = P5.mobile.ca.inds$xyz)
  P5.rmsd    <- rmsd(P5.xyz[1,P5.fixed.ca.inds$xyz], P5.xyz[,P5.mobile.ca.inds$xyz])
  P5.rmsData <- data.frame(P5.rmsd)
  P5.rmsData[]
  
  B1.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  446:766)
  B1.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  446:766)
  B1.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = B1.fixed.ca.inds$xyz, mobile.inds = B1.mobile.ca.inds$xyz)
  B1.rmsd    <- rmsd(B1.xyz[1,B1.fixed.ca.inds$xyz], B1.xyz[,B1.mobile.ca.inds$xyz])
  B1.rmsData <- data.frame(B1.rmsd)
  B1.rmsData[]
  
  #All the RMSD vector lists have been created, assign a time POTRA domain id to the rmsd and rbind them into a long list. Alt method in column df below.
  #Add the time variable as well
  numberFrames = length(P1.rmsd)
  timeFactor = simulationLength/numberFrames
  print(timeFactor)
  
  #Add time column
  #rmsData$time <- (1:length(P1.rmsd))*(timeFactor)
  P1.rmsData$time <- (1:length(P1.rmsd))*(timeFactor)
  P2.rmsData$time <- (1:length(P2.rmsd))*(timeFactor)
  P3.rmsData$time <- (1:length(P3.rmsd))*(timeFactor)
  P4.rmsData$time <- (1:length(P4.rmsd))*(timeFactor)
  P5.rmsData$time <- (1:length(P5.rmsd))*(timeFactor)
  B1.rmsData$time <- (1:length(B1.rmsd))*(timeFactor)
  P1.rmsData[]
  #Generate lists for the names
  potra1 <- rep("POTRA1", 1111)
  potra2 <- rep("POTRA2", 1111)
  potra3 <- rep("POTRA3", 1111)
  potra4 <- rep("POTRA4", 1111)
  potra5 <- rep("POTRA5", 1111)
  barrel <- rep("BARREL", 1111)
  
  #Tack them onto each dataframe
  P1.rmsData$domain <- as.factor(potra1)
  P2.rmsData$domain <- as.factor(potra2)
  P3.rmsData$domain <- as.factor(potra3)
  P4.rmsData$domain <- as.factor(potra4)
  P5.rmsData$domain <- as.factor(potra5)
  B1.rmsData$domain <- as.factor(barrel)
  
  
  
  
  #Make column names all the same 'rmsd'
  colnames(P1.rmsData) <- c("rmsd", "time")
  colnames(P2.rmsData) <- c("rmsd", "time")
  colnames(P3.rmsData) <- c("rmsd", "time")
  colnames(P4.rmsData) <- c("rmsd", "time")
  colnames(P5.rmsData) <- c("rmsd", "time")
  colnames(B1.rmsData) <- c("rmsd", "time")
  
  #rbind them all into one long list with
  bindData <- rbind(P1.rmsData, P2.rmsData)
  bindData <- rbind(bindData, P3.rmsData)
  bindData <- rbind(bindData, P4.rmsData)
  bindData <- rbind(bindData, P5.rmsData)
  bindData <- rbind(bindData, B1.rmsData)
  bindData[]
  #bindData$domain_f = factor(bindData$domain, levels=c('POTRA1','POTRA2','POTRA3','POTRA4','POTRA5','BARREL'))
  
  colnames(bindData) <- c("rmsd", "time", "domain")
  bindData[]
  bindData$domain <- factor(bindData$domain, levels = c("BARREL", "POTRA5", "POTRA4", "POTRA3", "POTRA2", "POTRA1"))
  bindData[]
  
  #p = qplot(time, rmsd, 
  #          data=bindData, 
  #          geom=c("point", "smooth")) +
  #  facet_grid(domain ~ .) +
  #  labs(title = plotTitle, x = "Time (ns)", y = "RMSD (Angstrom)")
  
  
  plotIt <- ggplot(data = bindData, aes(x=time, y=rmsd)) + 
    stat_smooth(span = 0.1) +
    geom_point(aes(colour = rmsd)) +
    scale_colour_gradient2(low = "blue", mid = "green", high ="red", midpoint = 1.8) +
    facet_grid(domain ~ .) +
    ggtitle(plotTitle) +
    xlab("Time (ns)") + ylab("RMSD (Angstrom)")
  
  ggplot_build(plotIt)
  ggsave(outFile)
  
  print(plotIt)
  bindData[]
}