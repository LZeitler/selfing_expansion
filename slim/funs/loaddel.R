require(data.table)
args = commandArgs(trailingOnly=TRUE) # 1=scoeffs, 2=gmatrix, 3=output

s<-fread(args[1],header=F,data.table=F)
m<-fread(args[2],data.table=F)

## only deleterious sites go into
## load calculation, using S coefficients and matrix of gt*dom
if (all(dim(m)>1)) {
    ms=m[s<0,]*s[s<0,1]
    v=apply(ms,2,function(x) prod(1+x))
} else {
    v=data.frame(f='NA')
}
fwrite(data.frame(v),args[3],col.names=F,sep='\t')
