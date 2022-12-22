source('~/code/r/source_me.R')

setwd("../output")

load("empi/geos.RData") # init data for empirical data
load("empi/geos_regionColors.RData")
load("sims/init.RData") # init files for simulations

####
#### FIGURE 4
####

## empirical pi
## load empirical results

pixyData <- fread("empi/20220317_pixy_pi.txt",data.table=F)

pixyData <- left_join(pixyData,geos,by="pop") %>%
    mutate(selfRate=ifelse(country%in%c("Italy","Spain","Greece"),
                           "selfRate=0",
                           "selfRate=0.95"),
           isEdge=ifelse(country%in%c("Italy","Greece"),
                         "Core",
                         "Edge"))

## numbers
pixyData %>%
    group_by(region) %>%
    summarise(mean(avg_pi))

diversityPlotEmpi <- pixyData %>%
    filter(region%in%c("Apuan Alps",'Abruzzo','Swiss Alps','French Alps')) %>%
    ggplot()+
    geom_violin(aes(x=region,y=avg_pi,fill=region),alpha=.8)+
    stat_summary(aes(x=region,y=avg_pi,group=region),
                 fun="mean",geom="point",shape=18,size=5)+
    scale_fill_manual(values=c(regionColors,setNames(regionColors,paste0("\n",names(regionColors)))))+
    labs(y="Nucleotide diversity",x="Region")+
    background_grid(major="y")+
    theme(legend.position='none')
diversityPlotEmpi



## empirical inbreeding
## inbreeding coefficients Arabis data
inbreedingData <- fread("empi/inbreedingF_ngsf_20211207.txt",data.table=F,header=F)
inbreedingNames <- fread("empi/inbreedingF_names_ngsf_20211207.txt",data.table=F,header=F)
badSamples <-
    c('Br22','Br06','Cc05','St15','Am01','Br18','Br24','Pi40','Ma97',
      'Am92','Ga91','Ga95','Mv93','Ma28','Pa9','Pi9','Pi95') # samples with low quality mentioned in
                                        # ms and PCR doubles
inbreedingData <- data.frame(
    ind=inbreedingNames$V1,
    pop=substr(inbreedingNames$V1,1,2),    
    ibcoeff=inbreedingData$V1
) %>% filter(!ind%in%badSamples)

inbreedingData <- left_join(inbreedingData,geos,by="pop")

inbreedingData <- inbreedingData %>% 
    filter(!pop%in%c("CG","WT","GG")) %>% 
    mutate(pop=ifelse(pop %in% c("NU","TV","SV","RI","GE"), "Scand",
               ifelse(pop %in% c("LC","AN"), "Spain",
                      ifelse(pop == "VI", "Greece",pop))))
inbreedingData$popC <-
    factor(inbreedingData$pop,unique(inbreedingData$pop[order(inbreedingData$region2.order)]),ordered=T)


inbreedingPlotX <- inbreedingData %>%
    ggplot()+
    geom_boxplot(aes(popC,ibcoeff,fill=region),outlier.alpha=0)+
    geom_jitter(aes(popC,ibcoeff,group=region),alpha=.5)+
    scale_fill_manual(values=regionColors)+
    labs(y="Inbreeding Coefficient",fill="Region")+
    background_grid()+
    theme(
        axis.text.x=element_text(angle=90,
                                   hjust=1,
                                   vjust=.4,
                                   size=10,
                                   lineheight = .7),
        axis.title.x = element_blank(),
        legend.position='none')
inbreedingPlotX

## colorful annotation for regions
annoData1 <- inbreedingData %>% group_by(popC) %>%
    summarise(region) %>% unique %>% ungroup %>%
    mutate(r=seq_along(popC))
annoData2 <- annoData1 %>% group_by(region) %>% summarise(m=mean(r)) %>%
    mutate(text=ifelse(region%in%c("Scandinavia","Spain","Greece"),"",as.character(region)))
    
annoRasterPlot <- ggplot(annoData1) +
    geom_raster(aes(r,1,fill=region,alpha=.5))+
    scale_fill_manual(
        values=regionColors
    )+
    scale_x_continuous(expand = c(0,0))+
    geom_text(data=annoData2,
              aes(m,1,label=text))+
    theme(
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        legend.position="none"
    )

inbreedingPlot <- plot_grid(
    inbreedingPlotX,
    NULL,
    annoRasterPlot,
    ncol=1,
    align="v",
    axis="lr",
    rel_heights = c(.9,-.05,.1)
)
inbreedingPlot


## numbers
inbreedingData %>%
    group_by(pop,region) %>%
    summarise(max(ibcoeff),min(ibcoeff),mean(ibcoeff),sd=sd(ibcoeff)) %>% arrange(sd) %>% data.frame

inbreedingData %>%
    group_by(region) %>%
    summarise(max(ibcoeff),min(ibcoeff),mean(ibcoeff),sd=sd(ibcoeff)) %>% arrange(sd) %>% data.frame


## map
require(magick)
mapPlot <- image_read("../plots/map_20220620_aw_margins_noarrow.png")
mapPlot <- ggdraw()+draw_image(mapPlot)


## R'xy
load("empi/comparisons_Rxy_20221206.RData")

geosRegions <- geos %>% ungroup %>%  select(pop,region)
seventeen <- c("Am","Br","Ca","Cc","Es","Ga","Gf","Gs","Gz",
               "La","Ma","Mv","Pa","Pi","Po","Se","St")

## add annotation
comparisons <- comparisons %>%
    inner_join(geosRegions,by=c("px"="pop")) %>%
    rename(region1=region) %>% 
    inner_join(geosRegions,by=c("py"="pop")) %>%
    rename(region2=region) %>%
    mutate(Comparison=paste0(region1," x ",region2)) 

## calculate R'xy for nonsyn and lof
comparisonsDash <- comparisons %>%
    select(-lxnoty,-lynotx) %>% 
    pivot_wider(names_from = type, values_from = rxy) %>%
    mutate("R'xy_del"=sel/neu,  # calculate R'xy for del sites
           "R'xy_lof"=lof/neu) %>% # and for lof sites
    select(-neu,-sel,-lof) %>% 
    pivot_longer(cols=starts_with("R'xy"),
                 names_to="type",
                 values_to="R'xy",
                 names_pattern = "R'xy_(.*)")

## confidence intervals
myCI <- 0.95
alpha <- 1-myCI
comparisonsCI <- comparisonsDash %>%
    filter(run!='all') %>%
    group_by(across(c(-`R'xy`,-"run"))) %>%
    summarise(n=n(),
              jackEST=mean(`R'xy`),
              jackSE=sqrt(var(`R'xy`)/n),
              CI=qt(alpha/2,n-1,lower.tail=F)*jackSE,
              jackCI_top=c(jackEST+CI),
              jackCI_bot=c(jackEST-CI)
              )

## merge full estimates with partial CI estimates
comparisonsDash <- inner_join(
    filter(comparisonsDash,run=='all') %>% select(-run),
    comparisonsCI %>% select(-n)
)
    

## filtering by only interesting comparisons
selectCombinations <- c("Apuan Alps x Abruzzo",
                        "French Alps x Abruzzo",
                        "Swiss Alps x Abruzzo",
                        "French Alps x Apuan Alps",
                        "Swiss Alps x Apuan Alps",
                        "Swiss Alps x French Alps")

comparisonsDashFiltered <- comparisonsDash %>%
    filter(
        px%in%seventeen,
        py%in%seventeen,
        Comparison%in%selectCombinations
    ) %>% mutate(Comparison=factor(Comparison,
                                   levels=selectCombinations,
                                   labels=gsub(" x "," ×\n",selectCombinations)))

## plotting R'xy
rxyPlot <- comparisonsDashFiltered %>%
    mutate(popCombi=paste(px,"x",py),
           type=factor(type,
                       levels=c("del","lof"),
                       labels=c("deleterious","LoF"))
           ) %>% 
    ggplot()+
    geom_point(
        aes(paste(type,
            reorder(`R'xy`,`R'xy`)),
            `R'xy`,
            color=type,
            group=popCombi),
        position=position_dodge(1)
    )+
    geom_errorbar(
        aes(paste(type,
            reorder(`R'xy`,`R'xy`)),
            ymin=jackCI_bot,
            ymax=jackCI_top,
            color=type,
            group=popCombi),
        position=position_dodge(1),
        width=3)+
    facet_grid(~Comparison,switch="x",scales="free_x")+
    geom_hline(yintercept = 1, linetype="dashed", color="grey50", size=1)+
    scale_color_manual(values=c(`deleterious`="#998ec3",LoF="#f1a340"))+
    labs(y=expression(R*"'"[xy])
        ,color="Prediction"
         )+
    background_grid(major="y")+
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=10),
        legend.position=c(.35,.95),
        legend.background=element_rect(fill="white")
    )+
    guides(color=guide_legend(title.position = "left",
                              ncol=2,
                              override.aes=list(shape=15,size=8,linetype=NA)))
rxyPlot


## Assemble Figure 4.
## run fig1.R before running this
fig4 <- plot_grid(
    plot_grid(mapPlot,NULL,rel_heights=c(20,1),ncol=1,scale=.95),
    inbreedingPlot,
    diversityPlotEmpi,
    rxyPlot,
    nrow=2,
    labels=c("A","B","C","D")
)
fig4

saveplot(fig4,"../plots/20221216_fig4",13,9)





#### phyloP recessive (skipped)
#phylo <- fread("empi/phyloPLoad.inds.q80.noscore.20220425.csv",data.table=F) %>%

