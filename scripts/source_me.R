require(data.table)
require(tidyr)
require(dplyr)
require(ggplot2)
require(cowplot)
require(grid)
require(gridExtra)
require(viridis)
require(colorspace)
require(reshape2)
require(ggrepel)
require(RColorBrewer)
options(max.print=999)

## install.packages(c("data.table","tidyr","dplyr","ggplot2","cowplot","grid","gridExtra","viridis","colorspace","reshape2","ggrepel") 

## cowplot default
theme_set(theme_cowplot())

## a function to save plots quickly in pdf and png
saveplot <- function(x, name, width=10, height=6){
    ggsave(paste0(name, ".pdf"),
           x, width = width, height = height, limitsize = F)
    ggsave(paste0(name, ".png"),
           x, width = width, height = height, limitsize = F,
	   type="cairo",
           bg = "white")
}

## a function to get complements
complement <- function(x){ ## A-T and G-C are complementary
    bases=c("A","C","G","T")
    xx <- unlist(strsplit(toupper(x),NULL))
    paste(unlist(lapply(xx,function(bbb){
    if(bbb=="A") compString<-"T"
    if(bbb=="C") compString<-"G"
    if(bbb=="G") compString<-"C"
    if(bbb=="T") compString<-"A"
    if(!bbb %in% bases) compString<-"N"
    return(compString)
})),collapse="")
}

mydatestamp <- function() as.character(system('date +%Y%m%d',T))
mytimestamp <- function() as.character(system('date +%Y%m%d%H%M%S',T))

## fread and save a remote file to local without ESS being weird
freadr <- function(server,name,data.table=F,...){
    if (Sys.info()["nodename"]=="leo-thinkpad") {
        lname <- paste0("/tmp/fread_",rev(unlist(strsplit(name,"/")))[1])
        cmd <- paste0("rsync -a ",server,":",name," ",lname,"; cat ",lname)
        cat(cmd,"\n")
        fread(cmd=cmd,
              data.table=data.table,...)
    } else {
        fread(name,data.table=data.table,...)
    }   
}
