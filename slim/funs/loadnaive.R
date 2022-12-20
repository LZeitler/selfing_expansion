require(data.table)
args = commandArgs(trailingOnly=TRUE) # 1=scoeffs, 2=gmatrix.temp, 3=output (gmatrix)

s<-fread(args[1],header=F,data.table=F)
m<-fread(args[2],header=F,data.table=F,sep="\t")
nind <- length(names(m))

## naive load, recessive and additive models, only deleterious mutations
mdel=m[s<0,]

## additive
madd=matrix(NA, nrow = nrow(mdel), ncol = ncol(mdel))
madd[mdel=="0|0"] <- 0
madd[mdel=="1|1"] <- 2
madd[mdel=="0|1"] <- 1
madd[mdel=="1|0"] <- 1

## recessive
mrec=matrix(NA, nrow = nrow(mdel), ncol = ncol(mdel))
mrec[mdel=="0|0"] <- 0
mrec[mdel=="1|1"] <- 1
mrec[mdel=="0|1"] <- 0
mrec[mdel=="1|0"] <- 0

d <- data.frame(nind,
                meanAddLoad=mean(colSums(madd)),
                meanRecLoad=mean(colSums(mrec)),
                loci=nrow(mdel),
                polyloci=nrow(mrec[rowSums(madd)!=nind*2,]))

fwrite(d,args[3],sep='\t',col.names=F)
