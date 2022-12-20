require(data.table)
args = commandArgs(trailingOnly=TRUE) # 1=scoeffs, 2=gmatrix.temp, 3=output

s<-fread(args[1],header=F,data.table=F)
m<-fread(args[2],header=F,data.table=F,sep="\t")

## freq calculation
ind <- ncol(m)

cnt <- sapply(1:nrow(m), function(r) {
    sum(as.numeric(c(substr(m[r,],1,1),substr(m[r,],3,3))))
})

c=data.frame(sel=unname(s),count=cnt,freq=cnt/(2*ind),populationSize=ind)
fwrite(c,args[3],sep='\t')
