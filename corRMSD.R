
#Megaplot. Calculates RMSD for individual elements of the protein.
#N.gono VERSION
library(bio3d)
library(ggplot2)
library(Hmisc)
#First frame

directory = "/home/dremyes/WaitingRoom/trjs"
trajFile = "Ng10nsIII.pdb"
outFile = "SomethingAwesome.png"
plotTitle = "SomethingAwesome"
simulationLength = 10 


setwd <- setwd(directory)
pdb_single <-read.pdb(trajFile)
  #All frames 
pdb_multi <- read.pdb(trajFile, multi = TRUE)
POTRA1 <- atom.select(pdb_single, resno =  1:70)


Full.fixed.ca.inds <- atom.select(pdb_single, elety="CA")
Full.mobile.ca.inds <- atom.select(pdb_multi, elety="CA")
Full.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = Full.fixed.ca.inds$xyz, mobile.inds = Full.mobile.ca.inds$xyz)

P1.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno = 1:70)
P1.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno = 1:70)
P1.rmsd    <- rmsd(P1.xyz[1,P1.fixed.ca.inds$xyz], P1.xyz[,P1.mobile.ca.inds$xyz])
P1.rmsData <- data.frame(P1.rmsd)
P1.rmsData[]
numberFrames = length(P1.rmsd)
timeFactor = simulationLength/numberFrames
print(timeFactor)

P1.rmsData$time <- (1:length(P1.rmsd))*(timeFactor)
potra1 <- rep("POTRA1", 101)
P1.rmsData$domain <- potra1
colnames(P1.rmsData) <- c("rmsd", "time")

P1.rmsData[]
################################################################################################
#P1 Fit Function Done#
################################################################################################
POTRA3 <- atom.select(pdb_single, resno = 155:244)
P3.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  155:244)
P3.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  155:244)
P3.xyz     <- fit.xyz(fixed = pdb_single$xyz, mobile=pdb_multi, fixed.inds = P3.fixed.ca.inds$xyz, mobile.inds = P3.mobile.ca.inds$xyz)
P3.rmsd    <- rmsd(P3.xyz[1,P3.fixed.ca.inds$xyz], P3.xyz[,P3.mobile.ca.inds$xyz])
P3.rmsData <- data.frame(P3.rmsd)
P3.rmsData[]
P3.rmsData$time <- (1:length(P3.rmsd))*(timeFactor)
potra3 <- rep("POTRA3", 101)
P3.rmsData$domain <- potra3
colnames(P3.rmsData) <- c("rmsd", "time")

corData <- data.frame()
corData[]
corData <- rbind(P1.rmsData, P3.rmsData) 
colnames(corData) <- c("rmsd", "time", "domain")
corData[]

cor(corData, use="complete.obs", method="spearman") 

p = qplot(time, rmsd, 
          data=corData, 
          geom=c("point", "smooth")) +
  facet_grid(domain ~ .) +
  labs(title = plotTitle, x = "Time (ns)", y = "RMSD (Angstrom)")


print(p)
ggsave(outFile)
bindData[]