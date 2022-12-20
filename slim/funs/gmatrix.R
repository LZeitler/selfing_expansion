require(data.table)
args = commandArgs(trailingOnly=TRUE) # 1=dcoeffs, 2=gmatrix.temp, 3=output (gmatrix)

d<-fread(args[1],header=F,data.table=F)
m<-fread(args[2],header=F,data.table=F,sep="\t")

md <- replicate(ncol(m),d[,1])
md[m=="0|0"] <- 0
md[m=="1|1"] <- 2

fwrite(md,args[3],sep='\t',col.names=F)
