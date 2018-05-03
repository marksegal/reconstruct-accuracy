"GAM_perm_slice_principal_comp" <- function(GAM.3D.data,GAM.segment.data,NP.num=100,nperm=1000,scale=T){

## GAM.3D.data : 3D coordinates of reconstruction based on GAM proximities (normalized linkage disequilibrium) for a given chromosome, as obtained for example from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64881.  Format is start,end,x,y,z where start,end specify genome coordinates corresponding to the 3D reconstruction given by x,y,z.
## GAM.segment.data : binary indicators of whether a given locus (defined by start,end genomic coordinates corresponding to those in GAM.3D.data and hence for the same chromosome) was detected in each of the NP_1,..., NP_n nuclear profiles. Format is start,end, NP_1,...,NP_n as obtained, for example, from chromosome subsetted segmentation files at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64881
## NP.num: minimum number of within-profile detections in order for the profile to be used in reconstruction accuracy assessment
## nperm: number of permutations of detection indicators
## scale: perform principal components on scaled (unit variance) data

reslt <- matrix(NA,nrow=nperm,ncol=18)

## obtain profiles with sufficient detections

slice.cts <- apply(GAM.segment.data[,-(1:2)],2,sum)
NP.set <- which(slice.cts >= sort(slice.cts,decreasing=T)[NP.num])

## determine PC1+PC2, PC1, PC2 explained variances under permutation

for(j in 1: nperm){
	seg.mat <- apply(GAM.segment.data[,-(1:2)],2,sample)
		
	slice.pct <- first.pct <- second.pct <- numeric()
	
	for(i in 1:NP.num){	
		data.i <- GAM.3D.data[which(seg.mat[,NP.set[i]] == 1),]
		prcomp.i <- prcomp(data.i[,3:5],scale=scale)	
		slice.pct.var <- (prcomp.i$sdev[1]^2 + prcomp.i$sdev[2]^2) / sum(prcomp.i$sdev^2)
		slice.pct <- c(slice.pct.var,slice.pct)
		first.pct.var <- prcomp.i$sdev[1]^2  / sum(prcomp.i$sdev^2)
		first.pct <- c(first.pct.var,first.pct)
        second.pct.var <- prcomp.i$sdev[2]^2  / sum(prcomp.i$sdev^2)
        second.pct <- c(second.pct.var,second.pct)
		}
	
	reslt[j,] <- c(summary(first.pct),summary(slice.pct),summary(second.pct))
	}
				
## determine PC1+PC2, PC1, PC2 explained variances from original data
				
slice.pct <- first.pct <- second.pct <- numeric()
		
for(i in 1:NP.num){	
	data.i <- GAM.3D.data[which(seg.mat[,NP.set[i]] == 1),]
	prcomp.i <- prcomp(data.i[,3:5],scale=scale)	
	slice.pct.var <- (prcomp.i$sdev[1]^2 + prcomp.i$sdev[2]^2) / sum(prcomp.i$sdev^2)
	slice.pct <- c(slice.pct.var,slice.pct)
	first.pct.var <- prcomp.i$sdev[1]^2  / sum(prcomp.i$sdev^2)
	first.pct <- c(first.pct.var,first.pct)
    second.pct.var <- prcomp.i$sdev[2]^2  / sum(prcomp.i$sdev^2)
    second.pct <- c(second.pct.var,second.pct)
	}
		
reslt <- rbind(c(summary(first.pct),summary(slice.pct),summary(second.pct)),reslt)

## obtain summaries of explained variances over permutations 

mean.means <- apply(reslt[-1,c(4,10,16)],2,mean)
median.medians <- apply(reslt[-1,c(3,9,15)],2,median)

sd.means <- apply(reslt[-1,c(4,10,16)],2,sd)
mad.medians <- apply(reslt[-1,c(3,9,15)],2,mad)

## combine with explained variance summaries of actual reconstruction
var.summ <- as.data.frame(rbind(reslt[1,c(4,10,16)],mean.means,sd.means,reslt[1,c(3,9,15)],median.medians,mad.medians))

names(var.summ) <- c("PC1","PC1+PC2","PC2")
row.names(var.summ)[1] <- "original.mean" 
row.names(var.summ)[4] <- "original.median"

var.summ
}

