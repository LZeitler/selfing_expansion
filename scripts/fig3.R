source('~/code/r/source_me.R')

setwd("../output")

load("sims/init.RData") # init files for simulations
load("empi/geos.RData")
load("empi/geos_regionColors.RData")

require(ggnewscale)
require(patchwork)

#### comprare DFE simulated and empirical

#### FIGURE 3 (sims) AND supplemental (empirical)

## DFE sims

de.realized <- fread("sims/freq-core-edge_20220518022225.txt",data.table=F) %>%
    filter(sel<0) # we are interested in deleterious sites

## define bins of selection coefficient counts
selCoefBreaks1 <- c(-10,-1,-1e-3,-1e-4,-1e-99)
de.realized$sBinX <- cut(de.realized$sel,breaks=selCoefBreaks1,include.lowest=T,right=T,ordered_result=T)
de.realized$sBin <- factor(de.realized$sBinX,
                           levels=levels(de.realized$sBinX),
                           labels=c("lethal","[-0.001,-1)","[-0.0001,-0.001)","[-1e-99,-0.0001)")
                           )
de.realized %>% group_by(sBinX) %>% summarise(unique(sBin)) # sanity check
de.realized %>% group_by(sBinX) %>% summarise(n())          # info

## DFE-like
de.realized <- de.realized %>%
    filter(pop%in%c("t3_C\nfinal edge","core")) %>%
    mutate(pop=ifelse(pop=="core","Core","Edge"),
           selfRate.num=ifelse(pop=="Core",1,selfRate.num)) %>% 
    group_by(rep,selfRate.num,pop,
             .dots=c(allfacet)) %>%
    count(sBin,.drop=F,name="n") %>% 
    mutate(Proportion=prop.table(n)) %>%
    group_by(sBin,selfRate.num,pop,
             .dots=c(allfacet)) %>%
    summarise(n_=n(),
              n50=median(n),
              Proportion50=median(Proportion),
              ProportionMin=quantile(Proportion,0.05),
              ProportionMax=quantile(Proportion,0.95)) %>%
    ungroup %>%
    mutate(
        selfRate.numF=ifelse(selfRate.num==0,"outcrossing",paste0(selfRate.num*100,'% selfing')),
        selfRate.numF=factor(selfRate.numF,levels=c(c("outcrossing","50% selfing","95% selfing","100% selfing")),ordered=T),
        grouping=factor(paste(selfRate.numF,pop),
                        levels=c("100% selfing Core", "outcrossing Edge", "50% selfing Edge",
                                 "95% selfing Edge", "100% selfing Edge"), ordered=T))
gc()

plotFunDFEsims <- function(data) {
    dodge <- position_dodge(width=0.9)
    plot <- data %>%
        ggplot()+
        geom_bar(
            aes(sBin,Proportion50,
                fill=pop,
                group=grouping),
            color="grey50",
            stat="identity",position=dodge)+
        scale_fill_manual(values=unname(regionColors[names(regionColors)%in%c("Apuan Alps","Swiss Alps")]))+
        labs(fill="Population")+
        new_scale_fill()+
        geom_bar(
            aes(sBin,Proportion50,fill=selfRate.numF,
                group=grouping),
            color="grey50",
            alpha=.5,stat="identity",position=dodge)+
        scale_fill_grey(start=0,end=0.9)+
        background_grid(major="y",minor="y")+
        labs(fill="Mating",
             x="Selection coefficient",
             y="Proportion of deleterious sites")+
        guides(fill=guide_legend(override.aes = list(alpha=1)),
               color="none")+
        geom_errorbar(
            aes(sBin,ymax=ProportionMax,ymin=ProportionMin,
                group=grouping),position=dodge, width = 0.5)
    return(plot)
}

plotDFEwithMben <- plotFunDFEsims(filter(de.realized,mben=="mben=0.01"))+
    scale_x_discrete(limits=rev(levels(de.realized$sBin)))
plotDFEwithMben

## (Alternative) small lethal inset 
plotDFEInset <- plotFunDFEsims(filter(de.realized,sBin=="lethal",mben=="mben=0.01"))+
    facet_wrap(~sBin,nrow=1,scales="free")+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.title.y=element_blank(),
          axis.title.x = element_blank(),
          plot.background = element_rect(color="black",fill=alpha("white",1)),
          ## panel.background = element_rect(fill=alpha("white",1)),
          ## panel.border = element_rect(colour = "black", fill=NA, size=.5),
          )
plotDFEInset

legendDummy <- ggplot(data.frame(x=c(0,0,0,0,0),
                                 y=-c(1,3,4,5,6),
                                 self=paste(c(1,0,.5,.95,1)),
                                 deme=c(T,F,F,F,F),
                                 text=c("outcrossing","outcrossing",
                                        "50% selfing","95% selfing", "100% selfing")))+
    geom_point(aes(x=x,y=y,color=self),shape=15,size=7)+
    scale_color_grey(start=0,end=0.9)+
    new_scale_color()+
    geom_point(aes(x=x,y=y,color=deme),alpha=.5,shape=15,size=7)+
    scale_color_manual(values=rev(unname(regionColors[names(regionColors)%in%c("Apuan Alps","Swiss Alps")])))+
    geom_text(aes(x=x+.05,y=y,label=text), hjust=0,size=5)+
    annotate("text",x=0,y=-0.2,label="Core",size=5.5,hjust=0.3)+
    annotate("text",x=0,y=-2.2,label="Edge",size=5.5,hjust=0.3)+
    xlim(-0.03,.3)+
    ylim(-6.3,0)+
    theme_void()+
    theme(legend.position="none",
          ## panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.background = element_rect(fill=alpha("white",.7),color=NA)
          )
legendDummy

saveplot( # for presentation
    plot_grid(
        plotDFEInset+theme(axis.title.y=element_text(),
                           plot.background=element_blank())+
        labs(y="Proportion of DFE", title="Lethal proportion of DFE"),
        plot_grid(NULL,legendDummy,NULL,ncol=1,rel_heights=c(1,1,1)),
        ncol=2,
        rel_widths=c(1,.3)
    ),"~/pro/selfing-sims/plots/fig3_dfe_sims_lethal_inset_20221212",10,6
)


layout <- c(
  area(t = 1, l = 0, b = 5, r = 8), #main
  area(t = 1, l = 8, b = 3, r = 9)  #inset
)

plotDFEs <- plotDFEwithMben+
    theme(legend.position="none",
          ## plot.margin = margin(5,50,5,5)
          )+
    annotation_custom(ggplotGrob(legendDummy),
                      xmin=0.5,xmax=1.45,
                      ymin=0.3,ymax=Inf)+
    plotDFEInset+
    plot_layout(design=layout)
plotDFEs

saveplot(plotDFEs,"../../paper/plots/fig3_dfe_sims_mben_ms_20221129",9,5)

#### EMPIRICAL
## fitdadi
## dadi-demography
dem <- fread("empi/dadi_demography_all.txt",data.table=F)

names(dem) <-
    c("grpidx","params","theta","Ne","llopt","llmodel","vcf","modelname","file")
dem$pop <- substr(dem$vcf,60,61)

dem$numpars <- lengths(sapply(strsplit(dem$params,"\\s|\\[|\\]"), # parse number of parameters
                              function(x) as.numeric(x[!is.na(as.numeric(gsub("None",NA,x)))])))
dem <- dem %>% mutate(aic=-2*llmodel+2*numpars)

## best X inferences AIC for each pop
top <- 100
dembestTop <- dem %>%
    group_by(file,pop,modelname) %>%
    filter(!duplicated(aic)) %>%
    group_by(file,pop) %>% 
    filter(rank(abs(aic))%in%1:top) %>%
    mutate(rank=rank(abs(aic))
           ) %>% ungroup
    
dembestTop$modelnameF <- factor(dembestTop$modelname,
                                levels=c("snm", "snm_inbreeding", 
                                "two_epoch", "two_epoch_inbreeding",
                                "bottlegrowth", "three_epoch"))



## dadi-fitdadi
dfe <- fread("empi/dadi_dfe_all.txt",data.table=F)
names(dfe) <-
    c("grpidx","shape","scale","[0,1]","[1,10]","[10,100]","[100,-99]","llmodel","Ne","modelname","LnsLs","file")

dfe$pop <- substr(dfe$file,13,14)

## choose best demography based on AIC
bestDemography <- dembestTop %>% group_by(pop) %>% filter(aic==min(aic))

## choose best dfe based on ll
bestDFEeachModel <- dfe %>%
    group_by(pop,modelname) %>% 
    filter(llmodel==max(llmodel)) %>%
    mutate(modelname=ifelse(modelname=="equilibrium","snm",modelname))

## merge by modelname
bestDemographyDFE <-
    inner_join(bestDemography,bestDFEeachModel,by=c("pop","modelname"),suffix=c("",".DFE"))

## annotate with pop infos
bestDemographyDFE <- inner_join(bestDemographyDFE,geos,by="pop")
bestDemographyDFE$popC <-
    factor(bestDemographyDFE$pop,unique(bestDemographyDFE$pop[order(bestDemographyDFE$region2.order)]),ordered=T)

## final overview
pickedDFE <- bestDemographyDFE %>% select(Ne,Ne.DFE,llmodel,modelname,pop,numpars,region) %>% unique %>% arrange(region)

## recalculate DFE from shape and rate, scale by Ne
my.h <- 0.3

dfeShape <- bestDemographyDFE %>% select(pop,popC,region,shape,scale,modelname,Ne) %>% unique %>% 
    ## select(-all_of(c("[0,1]","[1,10]","[10,100]","[100,-99]"))) %>%
    mutate(
        g1 = pgamma((1/my.h)*(0.5*1),   shape=shape, scale=scale),
        g2 = pgamma((1/my.h)*(0.5*10),  shape=shape, scale=scale),
        g3 = pgamma((1/my.h)*(0.5*100), shape=shape, scale=scale),
        g4 = pgamma((1/my.h)*(0.5*10000000000000), shape=shape, scale=scale),
        `[0,1)` = g1,          # bin 1 
        `[1,10)` = g2 - g1,    # bin 2 
        `[10,100)` = g3 - g2,     # bin 3 
        `[100,∞]` = g4 - g3    # bin 4 
    )
    
dfeShape1 <- melt(dfeShape, variable.name = "s", value.name = "ProportionX", measure.vars = c(
                   "[0,1)",
                   "[1,10)",
                   "[10,100)",
                   "[100,∞]"
                   ))


plotDFEempi1 <- ggplot(dfeShape1)+
    geom_bar(
        aes(s,ProportionX,fill=region,group=popC),stat="identity",position="dodge",color="grey50")+
    geom_text(
        aes(s,ProportionX,label=pop,group=popC), position=position_dodge(0.9),
        hjust=-0.5,
        size=3, angle=90)+
    scale_fill_manual(
        values=regionColors[names(regionColors)%in%c("Abruzzo","Apuan Alps","French Alps","Swiss Alps")]) +
    background_grid(major="y",minor="y")+
    labs(
        fill="Region",
        ## x=as.expression(bquote(italic("N")%*%"selection coefficient")),
        x=as.expression(bquote(italic("N")%*%"selection coefficient")),
        y="Proportion of deleterious sites")+
    theme(legend.position=c(0.1,0.85))
plotDFEempi1 # 20221205

saveplot(plotDFEempi1,paste0("/home/leo/projects/alpina-genomics/plots/",mytimestamp(),".dadi_dfe_binsCorrected"),10,6)

dfeShape1%>% group_by(pop,modelname) %>% summarise(sum(ProportionX)) # sanity, sums must be 1

dfeShape1



## numbers
dfeShape1 %>%
    group_by(region) %>%
    summarise(
        min(g1),
        min(g2),
        min(g3),
        min(g4),
        max(g1),
        max(g2),
        max(g3),
        max(g4)
              ) %>% data.frame

dfeShape1 %>% group_by(region,s) %>%
    summarise(mean(ProportionX))

dfeShape1 %>% group_by(s) %>%
    summarise(mean(ProportionX))

dfeShape1 %>%
    inner_join(geos) %>% 
    mutate(cc=ifelse(country=="Italy","IT","ALPS")) %>%
    group_by(cc,s) %>%
    summarise(mean(ProportionX))
