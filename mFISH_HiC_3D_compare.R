"mFISH_HiC_3D_compare" <- function(mFISH.3D.data,mFISH.genome.data,HiC.3D.data,symmetric = F){

## mFISH.3D.data: multiplex FISH 3D coordinates and replicates cf Table S4 Wang et al., 2016.
## mFISH.genome.data: multiplex FISH genome coordinates cf Table S1 Wang et al., 2016.

## HiC.3D.data: genomic and 3D coordinates from reconstruction based on HiC or other conformational assay.  Format: start,end,x,y,z where start,end give genome coords
## all data are chromosome specific
## symmetric argument to procrustes alignment function : set to F so mean mFISH structure serves as common referent

library(vegan)

names(mFISH.3D.data) <- c("Serial.Num","TAD.id","x","y","z")

nSerNum <- length(unique(mFISH.3D.data$Serial.Num)) ## number of replicates
nTADs <- length(unique(mFISH.3D.data$TAD.id)) ## number of FISH probes (located at TADs)

## center each mFISH replicate

ctr.mF <- numeric()

for(i in 1:nSerNum){
	mF.i <- mFISH.3D.data[mFISH.3D.data$Serial.Num == i,3:5]
	ctr.mF.i <- scale(mF.i,center=T,scale=F)
	ctr.mF <- rbind(ctr.mF,ctr.mF.i)
	}

ctr.mF <- cbind(mFISH.3D.data[,1:2],ctr.mF)
names(ctr.mF) <- c("Serial.Num","TAD.id","x","y","z")

## obtain mean (over replicates) 3D structure
## depends on multiplex FISH replicates having and/or being aligned to a common reference frame 

mean.mF.coords <- aggregate(ctr.mF,by=list(ctr.mF$TAD.id),FUN=mean,na.rm=T)

## impute individual replicate/TAD missing 3D coords using mean structure as template

imp.mF.coords <- ctr.mF

for(i in 1:nSerNum){
	for(j in 1:nTADs)
		if(is.na(ctr.mF[(i-1)*nTADs + j,3])) imp.mF.coords[(i-1)*nTADs + j,3:5] <- mean.mF.coords[j,4:6]
}

## compute procrustes sum-of-squared deviations of replicates from mean multiplex FISH structure

procrustes.mean.mF.ss <- rep(NA,nSerNum)
procrustes.mean.mF.scale <- rep(NA,nSerNum)

for(i in 1:nSerNum){
	procrust.mean.i <- procrustes(mean.mF.coords[,4:6],imp.mF.coords[((i-1)*nTADs + 1):	(i*nTADs),3:5],symmetric=symmetric)
	procrustes.mean.mF.ss[i] <- procrust.mean.i$ss
	procrustes.mean.mF.scale[i] <- procrust.mean.i$scale
	}

names(HiC.3D.data) <- c("start","end","x","y","z")		
names(mFISH.genome.data)[4:5] <- c("Start.genomic.coordinate","End.genomic.coordinate")

## bin HiC data to align with (generally lower resolution) multiplex FISH data

HiC.genome.matched.coords <- matrix(NA,nrow=nrow(mFISH.genome.data),ncol=3)		

for(i in 1:nrow(mFISH.genome.data)){
	mat.HiC <- HiC.3D.data[(HiC.3D.data$start >= 	mFISH.genome.data$Start.genomic.coordinate[i]) & (HiC.3D.data$end <= mFISH.genome.data$End.genomic.coordinate[i]), 3:5 ]
	HiC.genome.matched.coords[i,] <- apply(mat.HiC,2,mean)
	}

## compute procrustes sum-of-squared deviations of HiC reconstruction from mean multiplex FISH structure

procrustes.mean.mF.HiC <- procrustes(mean.mF.coords[,4:6],HiC.genome.matched.coords,symmetric=symmetric)
p.ss <- procrustes.mean.mF.HiC$ss

## histogram of sum-of-squared deviations of multiplex FISH replicates from mean multiplex FISH structure with sum-of-squared deviations of HiC reconstruction from mean multiplex FISH structure indicated

pdf(file=paste("histogram-comp-mF-HiC","pdf",sep="."))
hist(na.omit(procrustes.mean.mF.ss),col=4,main="",xlab="mFISH Procrustes ss",xlim=c(min(p.ss,procrustes.mean.mF.ss)-0.01,max(p.ss,procrustes.mean.mF.ss)+0.01))
abline(v=p.ss,col=2)

dev.off()
	
reslt <- list(mean.mF.coords = mean.mF.coords, imp.mF.coords = imp.mF.coords, HiC.genome.matched.coords = HiC.genome.matched.coords, procrustes.mean.mF.ss = procrustes.mean.mF.ss, procrustes.mean.mF.scale = procrustes.mean.mF.scale, procrustes.mean.mF.HiC = procrustes.mean.mF.HiC, nSerNum = nSerNum, nTADs = nTADs)

reslt
}

