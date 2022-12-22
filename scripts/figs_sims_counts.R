source('~/code/r/source_me.R')

setwd("../output")

load("empi/geos.RData") # init data for empirical data
load("empi/geos_regionColors.RData")
load("sims/init.RData") # init files for simulations


############################################################
###### SUPPLEMENTAL: Proxy Regression & Allele Counts ######
############################################################


#### Supplementary plots for proxy regression

## proxy regression for h=.3, recessive
require(ggpmisc)

naiveLoad <- bind_rows(
    fread("sims/loadnaive_t3_Z_summary.txt",data.table=F) %>% mutate(pop="Core"),
    fread("sims/loadnaive_t3_A_summary.txt",data.table=F) %>% mutate(pop="Deme before shift"),
    fread("sims/loadnaive_t3_B_summary.txt",data.table=F) %>% mutate(pop="Deme after shift"),
    fread("sims/loadnaive_t3_C_summary.txt",data.table=F) %>% mutate(pop="Edge")
) %>%
    rename(nind=1,
           additive=2,
           recessive=3,
           loci=4,
           polyloci=5,
           file=6) %>%
    mutate(par=as.numeric(gsub("(.+par)(\\d+)(.+$)","\\2",file)),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",file))) %>%
    select(-file) %>%
    inner_join(pan,by='par') %>%
    inner_join(par,by='par',suffix = c("",".num")) %>%
    mutate(selfRate=ifelse(pop%in%c("Core","Deme before shift"),"selfRate=0",paste(selfRate))) %>%
    mutate(popAndMating=paste0(pop,"\n",selfRate))


allFitness <- bind_rows(
    fread("sims/load_t3_Z_summary.txt",data.table=F) %>% mutate(pop="Core"),
    fread("sims/load_t3_A_summary.txt",data.table=F) %>% mutate(pop="Deme before shift"),
    fread("sims/load_t3_B_summary.txt",data.table=F) %>% mutate(pop="Deme after shift"),
    fread("sims/load_t3_C_summary.txt",data.table=F) %>% mutate(pop="Edge")
) %>%
    rename(fitness=1,
           file=2) %>%
    mutate(par=as.numeric(gsub("(.+par)(\\d+)(.+$)","\\2",file)),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",file))) %>%
    select(-file) %>%
    group_by(par,rep,pop) %>%
    summarise(fitness=mean(fitness)) %>% ungroup

countVsFitness <- inner_join(naiveLoad,allFitness,by=c("par","rep","pop")) %>%
    mutate(`additive model`=additive,
           `recessive model`=recessive) %>%
    pivot_longer(cols=c(`additive model`,`recessive model`),
                 names_to="Model",
                 values_to="oneOverLoad")


loadVsFitnessPlotFun <-
    function(data,model,bycolor,leg.pos=c(.75,.2),hide.facet=F,dotsize=5,reg.formula="y~x",y.lab) {
        formula <- as.formula(reg.formula)
        if (bycolor=="pop") {
            pointgeom <- geom_point(aes(color=pop),size=dotsize,alpha=.8)
            rsqAnnotation <- stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                                            after_stat(p.value.label),
                                                            sep = "*\", \"*")),
                                          formula=formula,
                                          label.x=0.02,
                                          label.y=0.98,
                                          size=5
                                          )
            smoothing <- geom_smooth(color="grey50",se=F,method="lm",formula=formula)
            scaleColor <- scale_color_manual(values=brewer.pal(n=4,"Set2"),na.value=alpha("white",0))
            additionalTheme <- labs(color="Demography")
            
        } else if (bycolor=="selfRate") {
            
            colors <- grey.colors(n=4,start=.1,end=.9)
            pointgeom <- geom_point(aes(color=as.character(selfRate)),size=5,alpha=.8)
            rsqAnnotation <- NULL
            smoothing <- NULL
            scaleColor <- scale_color_manual(values = colors)
            additionalTheme <- labs(color="Selfing Rate")

        } else stop("check bycolor")

        hide.facet.theme <- NULL
        if (hide.facet) {
            hide.facet.theme <- theme(
                strip.background = element_blank(),
                strip.text.x = element_blank()
            )
        }
        
        plot <- data %>%
            filter(mben!="mben=0") %>% # simulations with beneficials introduce more noise
            filter(Model%in%model) %>%
            mutate(selfRate=factor(selfRate,
                                   levels=c("selfRate=0","selfRate=0.5","selfRate=0.95","selfRate=1"),
                                   labels=c("0","0.5","0.95","1"))) %>% 
            mutate(pop=factor(pop,
                              levels = c(
                                  "Core",           
                                  "Deme before shift",
                                  "Deme after shift",
                                  "Edge"),
                              labels=c(
                                  "Core",           
                                  "before shift",
                                  "after shift",
                                  "Edge"
                              )
                             ,ordered=T)) %>% 
            ggplot(aes(fitness,oneOverLoad))+
            geom_point(size=dotsize+1,alpha=.8,shape=1)+
            pointgeom+
            rsqAnnotation+
            smoothing+
            facet_wrap(facetFormula(NA,"Model"),scales="free_y")+
            labs(x="Fitness",y=y.lab)+
            scaleColor +
            background_grid()+
            scale_shape_manual(values=c(16,5))+
            guides(shape = guide_legend(override.aes = list(alpha=1)),
                   color = guide_legend(override.aes = list(alpha=1)))+
            theme(legend.position = leg.pos,
                  legend.background = element_rect(fill=alpha('white',.6),
                                                   size=.3,
                                                   linetype='solid',
                                                   color='black'),
                  legend.margin = margin(r=.8,l=.8,t=.2,b=.7,unit='lines'),
                  strip.text = element_text(size=16))+
            hide.facet.theme+
            additionalTheme
        return(plot)
    }

## Recessive version for .3
loadVsFitnessPlotRec <- loadVsFitnessPlotFun(data=countVsFitness,
                                             model=c('recessive model'),
                                             bycolor='pop',
                                             leg.pos='none',
                                             hide.facet=T,
                                             y.lab="Recessive Load")
loadVsFitnessPlotRec

## Additive version for .3
loadVsFitnessPlotAdd <- loadVsFitnessPlotFun(data=countVsFitness,
                                             model=c('additive model'),
                                             bycolor='pop',
                                             leg.pos=c(.65,.8),
                                             hide.facet=T,
                                             y.lab="Additive Load")
loadVsFitnessPlotAdd


loadVsFitnessPlotBoth <- plot_grid(loadVsFitnessPlotRec,
          loadVsFitnessPlotAdd,
          nrow=1,
          labels=c("A","B"))
loadVsFitnessPlotBoth

saveplot(loadVsFitnessPlotBoth,'../plots/loadVsFitnessPlotH.3Both_20221208',12,6)



## Both versions for h=.5 sims proxy regression
load("sims/pointfive/init.RData") # init files for simulations .5 version

naiveLoad <- bind_rows(
    fread("sims/pointfive/loadnaive_t3_Z_summary.txt",data.table=F) %>% mutate(pop="Core"),
    fread("sims/pointfive/loadnaive_t3_A_summary.txt",data.table=F) %>% mutate(pop="Deme before shift"),
    fread("sims/pointfive/loadnaive_t3_B_summary.txt",data.table=F) %>% mutate(pop="Deme after shift"),
    fread("sims/pointfive/loadnaive_t3_C_summary.txt",data.table=F) %>% mutate(pop="Edge")
) %>%
    rename(nind=1,
           additive=2,
           recessive=3,
           loci=4,
           polyloci=5,
           file=6) %>%
    mutate(par=as.numeric(gsub("(.+par)(\\d+)(.+$)","\\2",file)),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",file))) %>%
    select(-file) %>%
    inner_join(pan,by='par') %>%
    inner_join(par,by='par',suffix = c("",".num")) %>%
    mutate(selfRate=ifelse(pop%in%c("Core","Deme before shift"),"selfRate=0",paste(selfRate))) %>%
    mutate(popAndMating=paste0(pop,"\n",selfRate))


allFitness <- bind_rows(
    fread("sims/pointfive/load_t3_Z_summary.txt",data.table=F) %>% mutate(pop="Core"),
    fread("sims/pointfive/load_t3_A_summary.txt",data.table=F) %>% mutate(pop="Deme before shift"),
    fread("sims/pointfive/load_t3_B_summary.txt",data.table=F) %>% mutate(pop="Deme after shift"),
    fread("sims/pointfive/load_t3_C_summary.txt",data.table=F) %>% mutate(pop="Edge")
) %>%
    rename(fitness=1,
           file=2) %>%
    mutate(par=as.numeric(gsub("(.+par)(\\d+)(.+$)","\\2",file)),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",file))) %>%
    select(-file) %>%
    group_by(par,rep,pop) %>%
    summarise(fitness=mean(fitness)) %>% ungroup

countVsFitnessAdd <- inner_join(naiveLoad,allFitness,by=c("par","rep","pop")) %>%
    mutate(`additive model`=additive,
           `recessive model`=recessive) %>%
    pivot_longer(cols=c(`additive model`,`recessive model`),
                 names_to="Model",
                 values_to="oneOverLoad")


## Recessive version for .5
loadVsFitnessPlotRecA <- loadVsFitnessPlotFun(data=countVsFitnessAdd,
                                             model=c('recessive model'),
                                             bycolor='pop',
                                             leg.pos='none',
                                             hide.facet=T,
                                             y.lab="Recessive Load")
loadVsFitnessPlotRecA

## Additive version for .5
loadVsFitnessPlotAddA <- loadVsFitnessPlotFun(data=countVsFitnessAdd,
                                             model=c('additive model'),
                                             bycolor='pop',
                                             leg.pos=c(.65,.8),
                                             hide.facet=T,
                                             y.lab="Additive Load")
loadVsFitnessPlotAddA


loadVsFitnessPlotBothA <- plot_grid(loadVsFitnessPlotRecA,
          loadVsFitnessPlotAddA,
          nrow=1,
          labels=c("A","B"))
loadVsFitnessPlotBothA

saveplot(loadVsFitnessPlotBothA,'../plots/loadVsFitnessPlotH.5Both_20221208',12,6)



#### SIMULATIONS: COUNTS OF ALLELES/LOCI
#### Supplementary for counts of deleterious alleles and loci
loadBoxFun <- function(data,model,ydata) {
    p <- data %>%
        filter(mben!="mben=0") %>% # simulations with beneficials introduce more noise
        filter(Model==model) %>% # filter by model
        mutate(selfRate=factor(selfRate,
                               levels=c("selfRate=0","selfRate=0.5","selfRate=0.95","selfRate=1"),
                               labels=c("outcr.","50%","95%","100%"),
                               ordered=T)) %>% 
        ggplot(aes_string("selfRate",ydata,fill="selfRate"))+
        geom_boxplot()+
        scale_fill_manual(values=grey.colors(n=4,start=.1,end=.9))+
        background_grid(major="y")+
        labs(x="Mating")+
        theme(legend.position="none")
    return(p)
}

## counts of deleterious alleles for simulated data. Uses same functions 
naiveAdditive <- filter(countVsFitness,pop%in%c("Edge","Core"))
loadBoxNaiveAdditive <- loadBoxFun(naiveAdditive,"additive model","additive")+
    labs(y="Mean deleterious allele count")+
    facet_grid(~pop,switch="x",scales = "free", space = "free")
loadBoxNaiveAdditive

## numbers
naiveAdditive %>% 
    filter(mben!="mben=0") %>%
    group_by(popAndMating) %>%
    summarise(mean(additive))

saveplot(loadBoxNaiveAdditive,
         "../plots/loadBoxNaiveAdditive_alleleCounts_20221015",
         width=8,height=4.5)

## counts of deleterious loci
naiveRecessive <- filter(countVsFitness,pop%in%c("Edge","Core"))
loadBoxNaiveRec <- loadBoxFun(naiveRecessive,"recessive model","recessive")+
    labs(y="Mean deleterious loci count")+
    facet_grid(~pop,switch="x",scales = "free", space = "free")
loadBoxNaiveRec

## numbers
naiveRecessive %>% 
    filter(mben!="mben=0") %>%
    group_by(popAndMating) %>%
    summarise(mean(recessive))

saveplot(loadBoxNaiveRec,
         "../plots/loadBoxNaiveRecessive_lociCounts_20221017",
         width=8,height=4.5)
