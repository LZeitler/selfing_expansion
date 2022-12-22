source('~/code/r/source_me.R')

setwd("../output")

load("sims/init.RData") # init files for simulations

####
#### FIGURE 1
#### UPDATED 20221208
####

## load sims
f <- fread("sims/fitness-stats_summary.txt",data.table=F) %>%
    rename(gen=1,deme=2,fitnessVar=3,fitnessMean=4,fitnessMedian=5,nind=6,countDel=7,countLethal=8,
           pi=9,meanObsHetSites=10,meanObsHet=11,dfeDel=12,dfeBen=13,file=14)

f$par <- as.numeric(gsub("(.+par)(\\d+)(.+$)","\\2",f$file))
f$rep <- as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",f$file))
f$file <- NULL
f <- inner_join(f,pan,by='par')
f <- inner_join(f,par,by="par",suffix = c("",".num"))

simData <- f %>%
    group_by(par,rep) %>%
    filter(mben.num==0.01) %>% # with beneficials
    filter(gen==max(gen)) %>% 
    filter(deme==1 |
           deme==max(deme)) %>%
    mutate(isEdge=ifelse(deme==1,"Core",
                  ifelse(deme==max(deme),"Edge",NA))) %>%
    mutate(selfRate=ifelse(isEdge=="Core","selfRate=0",paste(selfRate))) %>% 
    ungroup

colors <- grey.colors(n=4,start=.1,end=.9)


## sims speed
fmax <- f %>%
    group_by(par,rep) %>%
    group_by_at(vars(all_of(c(allfacet,"selfRate"))),.add=T) %>%
    summarise(maxGen=max(gen)) %>%
    mutate(selfRate=factor(selfRate,
                           levels=c("selfRate=0","selfRate=0.5","selfRate=0.95","selfRate=1"),
                           labels=c("outcrossing","50%","95%","100%"),
                           ordered=T)) %>% 
    filter(mben=="mben=0.01")

maxDeme <- f %>% summarise(max(deme)) %>% unlist %>% unname

## numbers
fmax %>% group_by(selfRate) %>% summarise(mean(maxGen),sd(maxGen))


## panel C max Time
maxTimePlot <-fmax %>%
    ggplot()+
    geom_boxplot(aes(selfRate,maxGen,fill=selfRate))+
    theme(axis.title.x=element_blank())+
    labs(y=paste0("Generations to deme ",maxDeme))+
    background_grid(major="y")+
    scale_fill_manual(values=c(colors[1:4]))
maxTimePlot



## sims pi
## numbers
simData %>%
    group_by(selfRate,isEdge) %>%
    summarise(format(mean(pi)))

diversityPlotSims <- simData %>%
    filter(selfRate.num%in%c(0,.5,.95,1)) %>%
    mutate(selfRate=factor(selfRate,
                           levels=c("selfRate=0","selfRate=0.5","selfRate=0.95","selfRate=1"),
                           labels=c("outcr.","50%","95%","100%"),
                           ordered=T)) %>% 
    ggplot()+
    geom_boxplot(aes(x=selfRate,y=pi,fill=selfRate))+
    scale_fill_manual(values=colors[1:4])+
    background_grid(major="y")+
    coord_cartesian(ylim=c(0,10e-5),
                    clip="off")+
    labs(y="Nucleotide diversity")+
    facet_grid(~isEdge,switch="x",scales = "free", space = "free")+
    theme(axis.title.x=element_blank())
diversityPlotSims


## sims fitness
fitMaxCore <- simData %>% # for rel fit calc
    filter(isEdge=="Core") %>% 
    group_by(par,rep) %>%
    group_by_at(vars(all_of(c(allfacet,"selfRate"))),.add=T) %>%
    filter(fitnessMean==max(fitnessMean)) %>%
    select(fitnessMaxCore=fitnessMean) %>%
    ungroup %>% select(-selfRate)

simData <- simData %>%
    filter(selfRate.num%in%c(0,.5,.95,1)) %>%
    inner_join(fitMaxCore) %>% 
    mutate(fitnessMean=fitnessMean/fitnessMaxCore) %>% # rel fit
    mutate(selfRate=factor(selfRate,
                           levels=c("selfRate=0","selfRate=0.5","selfRate=0.95","selfRate=1"),
                           labels=c("outcr.","50%","95%","100%"),
                           ordered=T))
simFitnessPlot <- simData %>% 
    ggplot(aes(selfRate,fitnessMean,fill=selfRate))+
    geom_boxplot()+
    scale_fill_manual(values=colors[1:4])+
    facet_grid(~isEdge,switch="x",scales = "free", space = "free")+
    background_grid(major="y")+
    labs(y="Relative fitness")+
    theme(axis.title.x=element_blank())
simFitnessPlot

## saveplot(simFitnessPlot+scale_fill_viridis_d(alpha=.8),
##          "../plots/loadBoxRecessivePres_20221212",
##          width=8,height=4.5)


## numbers
numbers <- simData %>%
    group_by(isEdge,selfRate) %>%
    summarise(mf=mean(fitnessMean),sd=sd(fitnessMean))

numbers$change <- (numbers$mf-numbers$mf[1])/numbers$mf[1]*100
numbers


## overview sims setup
require(magick)
modelPlot <- image_read_pdf("../plots/expansion_sim_comic1_bw_fixLeg_big.pdf")
modelPlot <- ggdraw()+draw_image(modelPlot)


## assemble main figure
## run fig4_empirical.R before running this
fig1 <- plot_grid(
    plot_grid(
        NULL,
        modelPlot,
        maxTimePlot+theme(legend.position="none"),
        ## loadVsFitnessPlotRec,
        nrow=1,
        labels=c("A","","B"),
        scale=1.05,
        rel_widths=c(.1,1.5,1)
    ),
    NULL,
    plot_grid(
        ## maxTimePlot+theme(legend.position="none"),
        diversityPlotSims+theme(legend.position="none"),
        simFitnessPlot+theme(legend.position="none"),
        nrow=1,
        align="vh",
        labels=c("C","D"),
        scale=1
    ),
    axis="t",
    nrow=3,
    rel_heights=c(1,.1,.8)
)
fig1


saveplot(fig1,"../plots/20221208_fig1",10,6.5)
