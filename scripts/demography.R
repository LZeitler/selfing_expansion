source('~/code/r/source_me.R')

setwd("../output")

load("empi/geos.RData")
load("empi/geos_regionColors.RData")

dem <- fread("empi/dadi_demography_all.txt",data.table=F)

names(dem) <-
    c("grpidx","params","theta","Ne","llopt","llmodel","vcf","modelname","file")
dem$pop <- substr(dem$vcf,60,61)
dem$numpars <- lengths(sapply(strsplit(dem$params,"\\s|\\[|\\]"), # parse number of parameters
                              function(x) as.numeric(x[!is.na(as.numeric(gsub("None",NA,x)))])))
dem <- dem %>% mutate(aic=-2*llmodel+2*numpars)

top <- 1
dembestTop <- dem %>%
    group_by(file,pop,modelname) %>%
    filter(!duplicated(aic)) %>%
    group_by(file,pop) %>% 
    filter(rank(abs(aic))%in%1:top) %>%
    mutate(rank=rank(abs(aic))
           ) %>% ungroup %>%
    select(-file,-rank,-vcf,-llopt)

## parse parameters for the best models
pars <- t(sapply(strsplit(dembestTop$params,"\\s|\\[|\\]"),
               function(x) as.numeric(x)[!is.na(as.numeric(x))]))
pars <- t(sapply(pars,"[", i = seq_len(max(lengths(pars)))))
dembestTop <- cbind(dembestTop,pars)
dembestTop$modelnameF <- factor(dembestTop$modelname,
                                levels=c("snm", "snm_inbreeding", 
                                "two_epoch", "two_epoch_inbreeding",
                                "bottlegrowth", "three_epoch"))

dembestTop <- dembestTop %>%
    mutate(
        N_ec= # contemporary N
            ifelse(modelnameF%in%
                   c("two_epoch","two_epoch_inbreeding"),
                   Ne*`1`,
            ifelse(modelnameF%in%
                   c("snm","snm_inbreeding"),
                   Ne,
            ifelse(modelnameF%in%
                   c("three_epoch","bottlegrowth"),
                   Ne*`2`,NA))),
        N_bot=
            ifelse(modelnameF%in%
                   c("three_epoch","bottlegrowth"),
                   Ne*`1`,Ne),
        timeptBot=
            ifelse(modelnameF%in%
                   c("two_epoch","two_epoch_inbreeding"),
                   2*Ne*`2`,
            ifelse(modelnameF%in%
                   c("bottlegrowth"),
                   2*Ne*`3`,
            ifelse(modelnameF%in%
                   c("three_epoch"),
                   2*Ne*`3`,
            ifelse(modelnameF%in%
                   c("snm","snm_inbreeding"),
                   1,NA)))),
        timept=
            ifelse(modelnameF%in%
                   c("three_epoch"),
                   2*Ne*`4`,timeptBot),
        F_ib=
            ifelse(modelnameF%in%
                   c("snm_inbreeding","two_epoch_inbreeding"),
                   `3`,NA))

dembestTop <- inner_join(dembestTop,geos,by="pop")

cliffsA <- dembestTop %>% select(modelnameF,pop,country,F_ib=F_ib,Ne,N_ec,N_bot,timept,timeptBot,theta,region,region2.color,region2.order)

cliffsDA <- bind_rows(cliffsA %>% mutate(idx=1, N=N_ec, T=0),
                      cliffsA %>% mutate(idx=2, N=N_ec, T=timept),
                      cliffsA %>% mutate(idx=3, N=N_bot, T=timept),
                      cliffsA %>% mutate(idx=4, N=N_bot, T=timeptBot),
                      cliffsA %>% mutate(idx=5, N=Ne, T=timeptBot),                      
                      cliffsA %>% mutate(idx=6, N=Ne, T=max(cliffsA$timeptBot,na.rm=T)*1.1))

## numbers: timing of bottlenecks
cliffsA %>% filter(modelnameF=="three_epoch") %>% group_by(isedge=grepl("French|Swiss",region))  %>%
    mutate(bottleneckDuration=timeptBot-timept) %>% select(pop,bottleneckDuration,isedge) %>% 
    filter(bottleneckDuration==min(bottleneckDuration) |
           bottleneckDuration==max(bottleneckDuration))

## regression for bottlegrowth model
bottle <- filter(cliffsDA,modelnameF%in%"bottlegrowth",idx%in%c(1,3))

expon <- function(t2,n1,n2,length=10){
    a <- n1
    b <- (n2/n1)**(1/t2)
    x <- seq(0,t2,length.out=length)
    y <- a * b**x
    return(data.frame(T=x,N=y))
}

predict <- data.frame()
for (p in unique(bottle$pop)) {
    t <- bottle %>% filter(pop==p)
    d <- expon(t[2,"T"],t[1,"N"],t[2,"N"])
    p <- bind_cols(select(t[1,],-N,-T),d) %>% mutate(idx=2)
    predict <- rbind(predict,p)
}

cliffsDAF <- bind_rows(
    filter(cliffsDA,!(modelnameF=="bottlegrowth" & idx==2)),
    predict
) %>% arrange(idx)

cliffPlotA <- ggplot(filter(cliffsDAF))+
    geom_path(aes(x=T,y=N,color=region,group=pop),size=1)+
    geom_label_repel(data=filter(cliffsDAF,idx==5,pop!="St"),aes(x=T,y=N,label=pop),
                     seed=12345)+
    geom_label_repel(data=filter(cliffsDAF,idx==5,pop=="St"),aes(x=T,y=N,label=pop),
                     nudge_y=-6000, nudge_x = 2000,
                     seed=12345)+
    labs(x="generations",y=expression(N["e"]),color="Region")+
    facet_wrap(~factor(region,levels=c("Abruzzo","Apuan Alps","French Alps","Swiss Alps")))+
    theme(legend.position="none")+
    background_grid(size.major = .2)+
    scale_color_manual(
        values=regionColors[names(regionColors)%in%c("Abruzzo","Apuan Alps","French Alps","Swiss Alps")])
cliffPlotA

saveplot(cliffPlotA,"../plots/20220908_demography",10,10)




## Table. Best AIC for every fitted model
top <- 1
demographyTable <- dem %>%
    group_by(file,pop,modelname) %>%
    mutate(rank=rank(abs(aic),ties.method = "first")) %>% 
    filter(rank%in%1:top) %>%
    ungroup %>%
    select(-file,-vcf,-llopt) %>%
    arrange(pop,aic) %>%
    select(Population=pop,Model=modelname,Ne,Likelihood=llmodel,AIC=aic)

fwrite(demographyTable,"empi/20220908_demographyTable.csv",sep="\t")
