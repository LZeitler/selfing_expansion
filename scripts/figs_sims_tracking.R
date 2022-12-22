source('~/code/r/source_me.R')

setwd("../output")

load("sims/init.RData") # init files for simulations

##############
## FIGURE 2 ##
##############

## load summary stats file
f <- fread("sims/fitness-stats_summary.txt",data.table=F) %>%
    rename(gen=1,deme=2,fitnessVar=3,fitnessMean=4,fitnessMedian=5,nind=6,countDel=7,countLethal=8,
           pi=9,meanObsHetSites=10,meanObsHet=11,dfeDel=12,dfeBen=13,file=14)
f$par <- as.numeric(gsub("(.+par)(\\d+)(.+$)","\\2",f$file))
f$rep <- as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",f$file))
f$file <- NULL
f <- inner_join(f,pan,by='par')
f <- inner_join(f,par,by="par",suffix = c("",".num"))
f <- f %>% filter(gen%%5==0)

f$meanObsHetFixed <- f$meanObsHetSites/(1e7-1)   # corrected in slim script from 20220813

## define factors and calculate contrasts
contrastGen <- 100 # how many generations should be contrasted?
lagK <- contrastGen/diff(c(f$gen %>% unique %>% sort)[1:2])

stationaryDemes <- c(shiftDeme+1,
                     shiftDeme+11,
                     shiftDeme+16)
names(stationaryDemes) <- sapply(1:3,function(x) paste0("Stationary deme ",stationaryDemes[x]))
ffEveryGenMinimal <- f %>%
    group_by(gen,par,rep) %>%
    filter(nind>10)
ffEveryGenStaticDemes <- ffEveryGenMinimal %>%
    mutate(isEdge=ifelse(deme==stationaryDemes[1],names(stationaryDemes[1]),
                  ifelse(deme==stationaryDemes[2],names(stationaryDemes[2]),
                  ifelse(deme==stationaryDemes[3],names(stationaryDemes[3]),
                         NA)))) %>%
    filter(!is.na(isEdge)) %>% ungroup()
ffEveryGenEdgeDemes <- ffEveryGenMinimal %>%
    filter(deme==max(deme)) %>%
    mutate(isEdge="Edge") %>% ungroup()
ffEveryGen <- bind_rows(
    ffEveryGenStaticDemes,
    ffEveryGenEdgeDemes
) %>%
    mutate(isEdge=factor(isEdge,
                         levels=c(names(stationaryDemes),"Edge"),
                         ordered=T)) %>% 
    group_by(par,rep,isEdge) %>%
    arrange(gen) %>% 
    mutate(fitnessChangePerGen=fitnessMean/lag(fitnessMean,lagK),
           countDelChangePerGen=countDel-lag(countDel,lagK),
           diffGen=gen-lag(gen,lagK),
           selfRate.current=
               factor(ifelse(deme>shiftDeme,selfRate.num,0)),
           demeAndMating=
               paste(isEdge,selfRate.current),
           ) %>%
    group_by(par,rep,gen) %>%
    mutate(demeReached=max(deme)) %>%
    filter(diffGen==contrastGen) %>%
    ungroup


## standardize time by start (0), shift (0.5), end (1)
ffBefore <- ffEveryGen %>%
    filter(demeReached<=shiftDeme) %>%
    group_by(par,rep) %>% 
    mutate(timeRank=rank(gen)/max(rank(gen))/2)
ffAfter <- ffEveryGen %>%
    filter(demeReached>shiftDeme) %>%
    group_by(par,rep) %>% 
    mutate(timeRank=rank(gen)/max(rank(gen))/2+.5)

ffEveryGen <- bind_rows(
    ffBefore,
    ffAfter
)


colorDemeAndMating <- setNames(c(
    "#195016", # dark green: edge outcr DFE:195016
    "#79B076", # light green: edge self #B2DF8A
    "#70afd9", # light blue: interior self #A6CEE3 DFE:82AFCD
    "#1F78B4"), # dark blue: interior outcr
    unique(filter(ffEveryGen,selfRate.num%in%c(0,.95),
                  isEdge%in%c("Edge","Stationary deme 35"))$demeAndMating))

colorDemeAndMatingLines <- colorDemeAndMating
colorDemeAndMatingLines["Edge 0"] <- "#33A02C"

colorDemeAndMatingBackground <- setNames(c("white","black","black","white"),names(colorDemeAndMating))


trajectoryPlotFun <- function(data, y.variable, y.label, smoothing.span=0.2) {
    require(ggnewscale)
    p <- data %>%
        group_by(par,rep,timeRank,isEdge,demeAndMating) %>%
        group_by_at(vars(all_of(c(allfacet,"selfRate","selfRate.current"))),.add=T) %>%
        summarise(!!y.variable := mean(!!sym(y.variable),na.rm=T)) %>%
        ggplot(aes_string("timeRank",y.variable,
                          group="demeAndMating"
                          )) +
        geom_point(aes(color=demeAndMating),
                   alpha=.9,size=.7)+
        scale_color_manual(values=colorDemeAndMating)+
        new_scale_color()+
        stat_smooth(
            geom="line",
            aes(color=demeAndMating),
            alpha=0.9,method="loess",span=smoothing.span,se=F,size=1.5)+
        scale_color_manual(values=colorDemeAndMatingBackground,,na.value="black")+
        new_scale_color()+
        stat_smooth(
            geom="line",
            aes(color=demeAndMating),
            alpha=.7,method="loess",span=smoothing.span,se=F,size=1)+
        scale_color_manual(values=colorDemeAndMating)+
        labs(x="Deme of edge",y=y.label,color="Population")+
        background_grid()+
        geom_vline(xintercept = 0.5, linetype="dashed", size=1)
    return(p)
}


## heterozygosity one plot (panel A)
trajectoryPlotHet <- trajectoryPlotFun(
    filter(ffEveryGen,
           selfRate.num%in%c(0,0.95),
           mben=="mben=0.01",
           isEdge%in%c(names(stationaryDemes[2]),"Edge")),
    "meanObsHetFixed",
    "Observed heterozygosity"
)+
    theme(legend.position="none")+
    scale_y_continuous()
trajectoryPlotHet

## Lethals (panel B)
trajectoryPlotLethals <- trajectoryPlotFun(
    filter(ffEveryGen,
           selfRate.num%in%c(0,0.95),
           mben=="mben=0.01",
           isEdge%in%c(names(stationaryDemes[2]),"Edge")),
        "countLethal",
        "Count of lethal alleles"
)+
    theme(legend.position="none")+
    scale_y_log10()+
    coord_cartesian(y=c(0.05,10))
trajectoryPlotLethals



## change every generation (C, D)
plotFitnessChangeFun <- function(data,colorDemes) {
    p <- data %>%
        group_by(par,rep,timeRank,isEdge,selfRate.current,demeAndMating) %>%
        group_by_at(vars(all_of(c(allfacet,"selfRate.num"))),.add=T) %>%
        summarise(fitnessChangePerGen=mean(fitnessChangePerGen,na.rm=T)) %>% # because population can be in same
                                        # deme, even if time is passing on and multiple sampling
                                        # times exist
        ggplot() +
        geom_point(aes(timeRank,fitnessChangePerGen,
                       color=paste(demeAndMating),
                       group=paste(isEdge)),alpha=.9,size=.7)+
        scale_color_manual(values=colorDemeAndMating)+
        new_scale_color()+
        stat_smooth(
            geom='line',
            aes(timeRank,fitnessChangePerGen,group=paste(isEdge),color=paste(demeAndMating)),
            alpha=.9,method='loess',span=.2,se=F,size=1.5)+
        scale_color_manual(values=colorDemeAndMatingBackground,na.value="black")+
        new_scale_color()+
        stat_smooth(
            geom='line',
            aes(timeRank,fitnessChangePerGen,color=paste(isEdge)),
            alpha=.7,size=1,method='loess',span=.2,se=F)+
        scale_color_manual(values=colorDemes,
                           guide=guide_legend(override.aes = list(alpha=1,size=7,shape=22,linetype=0))
                           )+
        geom_vline(xintercept = 0.5, linetype="dashed", size=1)+
        geom_hline(yintercept = 1, linetype="dashed", color="grey50", size=1)+
        facet_grid(~ifelse(selfRate.num==0,"outcrossing",paste0(selfRate.num*100,'% selfing')),
                   scales="free_y",switch="x")+
        labs(y=as.expression(bquote(frac(fitness[t],fitness[t-.(contrastGen)]))))+
        theme(legend.position="none",
              axis.title.x=element_blank()
              )+
        coord_cartesian(
            xlim=c(0.5,NA)
           ,ylim=c(0.8,1.17)
        )+
        background_grid()
    return(p)
}


plotFitnessChange0 <- plotFitnessChangeFun( # panel C
    filter(ffEveryGen,
           selfRate.num%in%c(0),
           mben=="mben=0.01",
           isEdge%in%c(names(stationaryDemes[2]),"Edge")),
    colorDemes = setNames(colorDemeAndMating[c(1,4)],c("Edge","Stationary deme 35"))
)
plotFitnessChange0

plotFitnessChange095 <- plotFitnessChangeFun( # panel D
    filter(ffEveryGen,
           selfRate.num%in%c(0.95),
           mben=="mben=0.01",
           isEdge%in%c(names(stationaryDemes[2]),"Edge")),
    colorDemes = setNames(colorDemeAndMatingLines[c(2,3)],c("Edge","Stationary deme 35"))
)
plotFitnessChange095

## numbers: how many generations does it take from deme 25 to 35 colonization (fitness drop)
filter(ffEveryGen,
       selfRate.num%in%c(0.95),
       mben=="mben=0.01",
       isEdge%in%c("Edge"),
       demeReached%in%c(25,35)) %>%
    group_by(demeReached) %>% 
    summarise(gen=mean(gen)) %>% 
    pull(gen) %>% diff

n0 <- filter(ffEveryGen,
       selfRate.num%in%c(0.95),
       mben=="mben=0.01",
       isEdge%in%c("Edge"),
       demeReached%in%c(25,35)) %>%
    group_by(rep,demeReached) %>% 
    summarise(gen=min(gen)) %>%
    group_by(rep) %>% 
    filter(length(rep)==2)  %>%
    ungroup %>% 
    mutate(d=gen-lead(gen)) %>% filter(d>0) %>% pull(d) 
min(n0) #45   
max(n0) #235  

## legend
dummyPlotChanges <- expand_grid(
    x=1:2,
    y=1:2,
    Edge=factor(c("outcrossing","95% selfing"),
                c("outcrossing","95% selfing")),
    Interior=factor(c("outcrossing","95% selfing"),
                    c("outcrossing","95% selfing"))
)

legend <- get_legend(
    ggplot(dummyPlotChanges,
           aes(
               x=x,
               y=y
           ))+
    geom_bar(aes(fill=`Interior`),stat="identity")+
    scale_fill_manual(values=unname(colorDemeAndMating[c(4,3)]),
                      guide = guide_legend(order = 2))+ #stationary: outcrossing, selfing
    new_scale_fill()+
    geom_bar(aes(fill=`Edge`),stat="identity")+
    scale_fill_manual(values=unname(colorDemeAndMating[c(1,2)]),
                      guide = guide_legend(order = 1))+ # edge: outcrossing, selfing
    theme(legend.text = element_text(size=15),
          legend.title= element_text(size=18))
)


breaks <- scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),
                             labels=c("Start","","Mating shift","","End"))

fig2 <- 
    plot_grid(
        trajectoryPlotHet+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+breaks,
        trajectoryPlotLethals+theme(axis.title.x=element_blank())+breaks,
        plot_grid(
            NULL,
            legend,
            plot_grid(
                plotFitnessChange0+theme(axis.text.x=element_blank())+breaks,
                plotFitnessChange095+breaks,
                ncol=1,
                labels=c("C","D")
            ),
            ncol=3,
            rel_widths = c(-.3,1.4,2.9)
        ),
        align="v",
        axis="b",
        labels=c("A","B"),
        ncol=1,
        rel_heights = c(1,1,2)
    )
fig2 <- ggdraw(add_sub(fig2, "Time",
                       vpadding=grid::unit(1,"lines"),
                       y=7, x=0.55, vjust=15))
fig2


saveplot(fig2,"../plots/20221202_fig2",5.99,11.7)



#### OPTIONAL AND SUPPLEMENTARY
colorDemeAndMatingAdd <- setNames(c(
    "#31a354", # dark green: edge outcr DFE:195016
    "#a1d99b", # light green: edge self #B2DF8A
    "#4292c6", # light blue: interior self #A6CEE3 DFE:82AFCD
    "#c6dbef"), # dark blue: interior outcr
    c("Edge 0.5", "Edge 1", 
      "Stationary deme 35 0.5", "Stationary deme 35 1")
    )
colorDemeAndMating <- c(colorDemeAndMatingAdd,colorDemeAndMating)

dummyPlotSup <- data.frame(color=colorDemeAndMating[sort(names(colorDemeAndMating))],
                           factor=sort(names(colorDemeAndMating))) %>%
    mutate(Population = ifelse(grepl("Stationary",factor),"Interior","Edge"))
dummyPlotSup$self <- rep(c(0,0.5,0.95,1),2)
dummyPlotSup$annotation <- paste0(dummyPlotSup$Population," ",dummyPlotSup$self*100, "%")

legendSup <- get_legend(ggplot(dummyPlotSup)+
    geom_bar(aes(self,annotation,fill=annotation),stat='identity')+
    scale_fill_manual(values=setNames(dummyPlotSup$color,dummyPlotSup$annotation))+labs(fill="Population &\nSelfing %"))



## optional supplement (het all self-rates)
y.lims <- c(3e-7,3e-4)
plotHetAll <-
    trajectoryPlotFun(
        filter(ffEveryGen,
               selfRate.num%in%c(0,0.5,0.95,1),
               mben=="mben=0.01",
               isEdge%in%c(names(stationaryDemes[2]),"Edge")),
        "meanObsHetFixed",
        "Observed Heterozygosity")+
    coord_cartesian(y=c(3e-7,3e-4))+
    theme(legend.position='none')+
    labs(x="Time")
plotHetAll <- plot_grid(plotHetAll,legendSup,ncol=2,rel_widths=c(1,0.2))
plotHetAll

saveplot(plotHetAll,"../plots/trajectoryHet_all_20221202",10,7)


plotLethalAll <-
    trajectoryPlotFunOne(
        filter(ffEveryGen,
               selfRate.num%in%c(0,0.5,0.95,1),
               mben=="mben=0.01",
               isEdge%in%c(names(stationaryDemes[2]),"Edge")),
        "countLethal",
        "Count of lethal alleles")+
    theme(legend.position='none')+
    labs(x="Time")+
    scale_y_log10()
plotLethalAll <- plot_grid(plotLethalAll,legendSup,ncol=2,rel_widths=c(1,0.2))
plotLethalAll

saveplot(plotLethalAll,"../plots/trajectoryLethal_all_20221202",10,7)


#### numbers

## ObsHet, outcrossing 24, edge
filter(ffEveryGen,deme==24,mben.num==0.01,isEdge=="Edge") %>%
    pull(meanObsHetFixed) %>% mean %>% format(scientific = T)

## ObsHet, selfing 35, edge
filter(ffEveryGen,deme==35,mben.num==0.01,isEdge=="Edge",selfRate.num==0.95) %>%
    pull(meanObsHetFixed) %>% mean %>% format(scientific = T)

## ObsHet, outcrossing 35, edge, 
filter(ffEveryGen,deme==35,mben.num==0.01,isEdge=="Edge",selfRate.num==0) %>%
    pull(meanObsHetFixed) %>% mean %>% format(scientific = T)

## ObsHet, selfing 35, end of sim (interior)
filter(ffEveryGen,deme==35,mben.num==0.01,selfRate.num==0.95) %>%
    group_by(rep) %>%
    filter(gen==max(gen)) %>% 
    pull(meanObsHetFixed) %>% mean %>% format(scientific = T)

## ObsHet, outcrossing 35, end of sim (interior)
filter(ffEveryGen,deme==35,mben.num==0.01,selfRate.num==0) %>%
    group_by(rep) %>%
    filter(gen==max(gen)) %>% 
    pull(meanObsHetFixed) %>% mean %>% format(scientific = T)

## Count of Lethals, mean before shift, at the end
n1 <- filter(ffEveryGen,deme==24,mben.num==0.01) %>%
    group_by(rep) %>%
    ## filter(gen==max(gen)) %>%   # at the end
    filter(isEdge=="Edge") %>% # when it is a edge
    group_by(isEdge,demeReached) %>%
    summarise(l=mean(countLethal),n()) %>% data.frame
n1

## Count of Lethals, Edge at the end
n4 <- filter(ffEveryGen,deme==50,mben.num==0.01,isEdge=="Edge") %>%
    group_by(isEdge,selfRate.current) %>%
    summarise(l=mean(countLethal)) %>% data.frame
n4

## percent change in lethals
n4$percL <- (n4$l-n1$l)/n1$l*100
n4$startL <- n1$l
n4

## Lethals: how much time required for purging in selfRate=1
n5 <- filter(ffEveryGen,demeReached%in%c(25,28),
             mben.num==0.01,selfRate.num==1,isEdge=="Edge")%>%
    group_by(rep,demeReached) %>%
    summarise(gen=min(gen)) %>%
    group_by(rep) %>%
    mutate(d=lead(gen)-gen) %>%
    ungroup %>%
    summarise(min(d,na.rm=T),max(d,na.rm=T))
n5

