source('~/code/r/source_me.R')

setwd("../output")

load("sims/init.RData") # init files for simulations

de.realized <- fread("sims/freq-core-edge_20220518022225.txt",data.table=F) %>%
    filter(sel<0,
           mben=="mben=0.01") # we are interested in deleterious sites

## segregating vs fixed alleles
alleleData <- de.realized %>%
    filter(pop%in%c("t3_C\nfinal edge","core")) %>%
    mutate(pop=ifelse(pop=="core","Core","Edge"),
           selfRate.num=ifelse(pop=="Core",0,selfRate.num),
           Allele=ifelse(freq==1,"fixed","segregating")) %>% 
    group_by(rep,selfRate.num,pop,.dots=c(allfacet),Allele) %>%
    summarise(#sumFreq=sum(freq), # sum freq ~ additive naive load (segregating load vs fixed load)
        Count=length(Allele) # count of loci that are fixed (or segregating) per replicate
    ) %>%
    group_by(rep,selfRate.num,pop,.dots=c(allfacet)) %>%
    mutate(
        ## Proportion=prop.table(Count)
        Proportion=Count/1e7
    ) %>%
    mutate(selfRate=factor(selfRate.num,
                           levels=c(0,0.5,0.95,1),
                           labels=c("outcr.","50%","95%","100%"),
                           ordered=T))

labels <- alleleData %>%
    filter(Allele=="fixed") %>% 
    group_by(selfRate,pop) %>%
    summarise(yval=max(Proportion),
              lval=mean(Count),
              lab=paste0("bar({x}) == ", lval))
    
fixSegPlot <- ggplot(filter(alleleData,Allele=="fixed"))+
    geom_boxplot(aes(selfRate,
                     Proportion,
                     fill=selfRate
                     ))+
    ## geom_text(data=labels,
    ##           mapping=aes(selfRate,
    ##                       yval,
    ##                       label=lab
    ##                       ),
    ##           parse=T,
    ##           nudge_y=.05
    ##           )+
    facet_grid(~pop,
               switch="x",
               scales = "free", space = "free")+
    ## scale_fill_manual(values=c("#9BB453","#495F3A"))+
    scale_fill_manual(values=grey.colors(n=4,start=.1,end=.9))+
    background_grid(major="y")+
    labs(y="Proportion of sites fixed for deleterious alleles")+
    theme(axis.title.x=element_blank(),
          legend.position='none')
fixSegPlot

saveplot(fixSegPlot,"~/pro/alpina-genomics/paper/plots/fixSegBox-Sims.20221019",
         width=8,height=4.5)

## combine fixation and counts into one plot.
## run fig4_empirical.R before this
allelesMultiPlot <- plot_grid(
    NULL,
    fixSegPlot,
    loadBoxNaiveRec+theme(axis.title.x=element_blank()),
    loadBoxNaiveAdditive+theme(axis.title.x=element_blank()),
    labels=c(NA,LETTERS[1:3]),
    nrow=1,
    rel_widths=c(.05,1,1,1),
    hjust = c(0,.4,0,0)
)
allelesMultiPlot

saveplot(allelesMultiPlot,"~/pro/alpina-genomics/paper/plots/allelesMulti.20221019",
         width=16,height=4.5)
