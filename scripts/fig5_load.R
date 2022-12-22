source('~/code/r/source_me.R')

setwd("../output")

load("empi/geos.RData") # init data for empirical data
load("empi/geos_regionColors.RData")
load("sims/init.RData") # init files for simulations

####
#### FIGURE 5
#### 20221208
####


## phyloP load recessive

loadPlot <- function(data,y,y.label,facetname,hide.facet=F){
    additionalTheme <- NULL
    if (hide.facet) {
        additionalTheme <- theme(
            strip.background = element_blank(),
            strip.text.x = element_blank()
        )
    }

    ## actual data plot
    boxplot <- data %>%
        mutate(facetvar=facetname) %>% 
        ggplot()+
        geom_boxplot(
            aes_string("popC",y,fill="region",group="popC")
        )+
        facet_wrap(facetFormula(NA,"facetvar"))+
        scale_fill_manual(
            values=regionColors
        )+
        background_grid()+
        labs(
            fill="Region",
            x="Population",
            y=y.label
        )+
        theme(
            legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_text(angle=90,
                                     hjust=1,
                                     vjust=.4,
                                     size=10,
                                     lineheight = .7),
            strip.text = element_text(size=16)
        )+
        additionalTheme
    
    return(boxplot)
}


snpeff <- fread("empi/snpEffLoad.py.inds.20220427.csv",data.table=F) %>%
    mutate(pop=substr(individual,1,2)) %>%
    inner_join(geos,by="pop") %>%
    filter(!pop%in%c("CG","WT","GG")) %>%
    filter(!individual%in%badSamples) %>% 
    mutate(pop=ifelse(pop %in% c("NU","TV","SV","RI","GE"), "Scand",
               ifelse(pop %in% c("LC","AN"), "Spain",
                      ifelse(pop == "VI", "Greece",pop))))
snpeff$popC <- factor(snpeff$pop,unique(snpeff$pop[order(snpeff$region2.order)]),ordered=T)
snpeff <- snpeff %>%
    group_by(individual,popC,region) %>%
    summarise(recessiveSNPeff=sum(recessiveSNPeff),
              additiveSNPeff=sum(additiveSNPeff))


## colorful annotation for regions
annoData1 <- snpeff %>% group_by(popC) %>% summarise(region) %>% unique %>% ungroup %>%
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




loadPlotSNPeffRecessive <- loadPlot(snpeff,"recessiveSNPeff",
                                    "Recessive Genetic Load",
                                    "SNPeff", hide.facet=T)
loadPlotSNPeffRecessive

loadPlotSNPeffAdditive <- loadPlot(snpeff,"additiveSNPeff",
                                   "Additive Genetic Load",
                                   "SNPeff", hide.facet = T)
loadPlotSNPeffAdditive

snpeffBoth <- plot_grid(
    loadPlotSNPeffRecessive,
    loadPlotSNPeffAdditive,
    NULL,
    annoRasterPlot,
    nrow=4,
    labels=c("A","B",NA,NA),
    align="v",
    axis="l",
    rel_heights = c(1,1,-.05,.12)
)
snpeffBoth

saveplot(snpeffBoth,"../plots/20221208_fig5",5,8)
