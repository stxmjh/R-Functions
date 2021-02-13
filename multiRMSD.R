#Megaplot. Calculates RMSD for individual elements of the protein.
library(bio3d)
library(ggplot2)
#First frame

multiRMSD <- function(directory = "/home/dremyes/WaitingRoom", trajFile = "step7_production.II.protein.pdb", outFile, plotTitle, simulationLength = 10){ 
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
POTRA5 <- atom.select(pdb_single, resno = 371:445)
BARREL <- atom.select(pdb_single, resno = 446:766)


#Alternative test 1)define set 

#P1.fixed.ca.inds <- atom.select(pdb_single, elety="CA", resno =  45:113)
#P1.mobile.ca.inds <- atom.select(pdb_multi, elety="CA", resno =  45:113)
#P1.rmsd <- rmsd(pdb_single, pdb_multi, a.inds = P1.fixed.ca.inds, b.inds = P1.mobile.ca.inds, fit =TRUE)





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
P3.rmsd    <- rmsd(P3.xyz[1,P3.fixed.ca.inds$xyz], P2.xyz[,P3.mobile.ca.inds$xyz])
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

#Generate lists for the names
potra1 <- rep("POTRA1", 101)
potra2 <- rep("POTRA2", 101)
potra3 <- rep("POTRA3", 101)
potra4 <- rep("POTRA4", 101)
potra5 <- rep("POTRA5", 101)
barrel <- rep("BARREL", 101)

#Tack them onto each dataframe
P1.rmsData$domain <- potra1
P2.rmsData$domain <- potra2
P3.rmsData$domain <- potra3
P4.rmsData$domain <- potra4
P5.rmsData$domain <- potra5
B1.rmsData$domain <- barrel


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

colnames(bindData) <- c("rmsd", "time", "domain")
bindData[]

p = qplot(time, rmsd, 
          data=bindData, 
          geom=c("point", "smooth")) +
          facet_grid(domain ~ .) +
          labs(title = plotTitle, x = "Time (ns)", y = "RMSD (Angstrom)")

  
print(p)
ggsave(outFile)
bindData[]
}









#numberFrames = length(P1.rmsd)
#timeFactor = 10/numberFrames
#print(timeFactor)  
#Add a time column
#rmsData$time <- (1:length(P1.rmsd))*(timeFactor)
#rmsData[]


#rmsData <- data.frame(P1.rmsData)
#rmsData$rmsd <- P2.rmsd 
#rmsData$rmsd <- P3.rmsd
#rmsData$rmsd <- P4.rmsd
#rmsData$rmsd <- P5.rmsd
#rmsData$rmsd <- B1.rmsd
#rmsData[]



#p = qplot(Concentration, Percent.of.control, 
#          data=rmsData, 
#          geom=c("point", "smooth"), colour=Response.type) +
#  scale_x_log10() +
#  facet_grid(Compound ~ Cell.line) +
#  coord_cartesian(ylim=c(-10, 110))
#print(p)##


#Define the time factor in nanoseconds
#/test simulationLength = 10
#numberFrames = length(rmsData)
#timeFactor = simulationLength/numberFrames
#timeFactor = 10/numberFrames
#print(timeFactor)  
#Add a time column
#rmsData$time <- (1:length(rmsd))*(timeFactor)
#plotIt <- ggplot(data = rmsData, aes(x=time, y=rmsd)) + 
#  geom_point(aes(colour = rmsd)) +
#  scale_colour_gradient2(low = "blue", mid = "blue", high ="red", midpoint = 4) +
#  geom_smooth(colour = "black", size = 1) +
#  ggtitle("Root mean square deviation variation") +
#  xlab("Time (ns)") + ylab("RMSD (Angstrom)")#

#ggplot_build(plotIt)#



