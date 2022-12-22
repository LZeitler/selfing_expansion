##  R XY compares the number of derived alleles of a particular category (LOF, missense, or
##  putatively neutral) found in one population relative to the other. If RXY = 1, both
##  populations X and Y have the same load of derived alleles, whereas if R XY > 1, then population
##  X has more derived alleles than Y and vice versa if R XY < 1
##  https://www.pnas.org/doi/pdf/10.1073/pnas.2023018118

## prepare in load_snpeff.py
setwd("/storage/homefs/lz20n032/pro/alpina-genomics/output")
.libPaths("~/R/R-4.1-manually/")
source('~/code/r/source_me.R')
load("../paper/output/empi/geos.RData")


## prepare populations, factor combinations
geosRegions <- geos %>% ungroup %>%  select(pop,region)

factors <- unique(geos$region)
factorsComb <- t(combn(factors,2))
factorsStr <- paste0(factorsComb[,1]," x ",factorsComb[,2])

seventeen <- c("Am","Br","Ca","Cc","Es","Ga","Gf","Gs","Gz","La","Ma","Mv","Pa","Pi","Po","Se","St")

## load tables
loci <- rbind(
    data.frame(fread("snpEffLoad.py.loci.20221202.csv",data.table=F),type="sel"),
    data.frame(fread("snpEffLoad.py.neutral.loci.20221202.csv",data.table=F),type="neu"),
    data.frame(fread("snpEffLoad.py.lof.loci.20221202.csv",data.table=F),type="lof")
)

## count derived alleles
counts <- loci %>%
    select(population,type,alleleFreq) %>%
    group_by(type)

## nest by populations and type, calculate Rxy
## Rxy from Do et al. 2015:
## Lx,notY: sum((cix/nix)*(1-ciy/niy))
## x, y: populations
## cix: derived count at site i and pop x
## nix: number of chromosomes at site i and pop x (c/n: derived allele freq in pop x)
## Rxy = Lx,notY/Ly,notx
allpops <- expand.grid(p1=unique(counts$population),p2=unique(counts$population))
comparisons <- data.frame()
for (r in rownames(allpops)) {
    px <- allpops[r,"p1"] # first population "x"
    py <- allpops[r,"p2"] # second population "y"
    for (t in unique(counts$type)) { # iterate through types of mutations (syn/nonsyn/etc)
        f1 <- counts %>% filter(population==px,type==t) %>%
            pull(alleleFreq) # vector of allele frequencies for population x
        f2 <- counts %>% filter(population==py,type==t) %>%
            pull(alleleFreq) # vector of allele frequencies for population y
        ## jackknife: 100 bootstraps with 100 contiguous blocks, 1 time with all data.
        len <- length(f1)
        segmentLength <- floor(len/101)
        jack <- rep(T,len)
        for (i in 0:100) {
            
            ## pick loci
            f11 <- f1[jack]
            f21 <- f2[jack]

            ## calculate L and Rxy, summarize over loci
            lxnoty <- sum(f11*(1-f21),na.rm=T)
            lynotx <- sum(f21*(1-f11),na.rm=T)
            rxy <- lxnoty/lynotx

            ## save calculations
            comparisons <- rbind(comparisons,
                                 data.frame(px,py,type=t,lxnoty,lynotx,rxy,
                                            run=ifelse(i==0,"all",i)))

            ## pick loci for jack for the next iteration
            jack <- c(rep(T,i*segmentLength),
                      rep(F,segmentLength),
                      rep(T,len-((i+1)*segmentLength)))
        }
    }
    cat(r,"\n")
}

save(comparisons, file = "comparisons_Rxy_20221206.RData")
save(loci, file = "loci_20221206.RData")



## From here in fig4_empirical
load("empi/comparisons_Rxy_20221206.RData")

## add annotation
comparisons <- comparisons %>%
    inner_join(geosRegions,by=c("px"="pop")) %>%
    rename(region1=region) %>% 
    inner_join(geosRegions,by=c("py"="pop")) %>%
    rename(region2=region) %>%
    mutate(Comparison=paste0(region1," x ",region2)) 

## R'xy
comparisonsDash <- comparisons %>%
    select(-lxnoty,-lynotx) %>% 
    pivot_wider(names_from = type, values_from = rxy) %>%
    mutate("R'xy_nonsyn"=sel/neu,  # calculate R'xy for nonsyn sites
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
    ) %>% mutate(Comparison=factor(Comparison,levels=selectCombinations))

## plotting
rxyPlot <- comparisonsDashFiltered %>%
    mutate(popCombi=paste(px,"x",py),
           type=factor(type,
                       levels=c("nonsyn","lof"),
                       labels=c("non-synonymous","LoF"))
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
    scale_color_manual(values=c(`non-synonymous`="#998ec3",LoF="#f1a340"))+
    labs(y=expression(R*"'"[xy])
        ,color="Prediction"
         )+
    background_grid(major="y")+
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=8),
        legend.position=c(.8,.9)
    )+
    guides(color=guide_legend(override.aes=list(shape=15,size=8,linetype=NA)))
rxyPlot

saveplot(rxyPlot,"../../plots/rdashxy_nns_20221207",10,6)



#############################
###### FIXATION #############
## fixed alleles per pop
load("empi/loci_20221206.RData")

nlen <- loci %>%
    group_by(population) %>% summarise(n=n())

fixed <- loci %>%
    group_by(population,type
             ) %>%
    summarise(
        fixed=sum(isDerivedFixed,na.rm=T),  # count of fixed
        fixedp=fixed/length(isDerivedFixed)# proportion of fixed over all of
                                        # this type segregating in population
    ) %>%
    inner_join(nlen) %>% 
    mutate(fixedpp=fixed/n) %>% # or by n+ns
    inner_join(geosRegions,by=c("population"="pop")) %>%
    mutate(type=factor(type,
                       levels=c("neu","sel","lof"),
                       labels=c("neu"="neutral",
                                "sel"="deleterious",
                                "lof"="LoF")
                       ))

fixedPPlot <- ggplot(filter(fixed,population%in%seventeen)) +
    geom_bar(aes(factor(population,levels=unique(geos$pop[order(geos$region2.order)])),fixedp,
                 fill=type),
             stat='identity',
             position='dodge')+
    scale_fill_manual(values=c(deleterious="#998ec3",LoF="#f1a340",
                               neutral="#542788"))+
    labs(x="Population",y="Proportion of fixed sites",fill="Prediction")+
    facet_grid(~region,switch='x',scales="free")+
    background_grid(major="y")
fixedPPlot # for supplement

saveplot(fixedPPlot,"../../plots/fixedProportion_nns_byPop_byAllSitesWithinType_20221216",9,5)

## numbers for fixation
filter(fixed,population%in%seventeen) %>%
    group_by(type,region) %>%
    summarise(mean(fixedp)*100)

filter(fixed,population%in%seventeen) %>%
    group_by(region) %>%
    summarise(mean(fixedp)*100)

filter(fixed,population%in%seventeen) %>%
    group_by(type) %>%
    summarise(mean(fixedp)*100)

filter(fixed,population%in%seventeen) %>%
    group_by(population) %>%
    summarise(m=mean(fixedp)*100) %>% arrange(-m)



#### Other plots - not in ms.
fixedPlot <- ggplot(filter(fixed,population%in%seventeen)) +
    geom_boxplot(aes(region,fixed,
                     fill=type))+
    labs(x="Region",y="Count of fixed loci",fill="Prediction")+
    background_grid(major="y")
fixedPlot

saveplot(fixedPlot,"../../plots/fixed_nns_20221205",9,5)


fixedPPlot <- ggplot(filter(fixed,population%in%seventeen)) +
    geom_boxplot(aes(region,fixedp,
                     fill=type))+
    labs(x="Region",y="Proportion of fixed loci",fill="Prediction")+
    background_grid(major="y")
fixedPPlot

saveplot(fixedPPlot,"../../plots/fixedProportion_nns_20221205",9,5)


fixedPPPlot <- ggplot(filter(fixed,population%in%seventeen)) +
    geom_bar(aes(factor(population,levels=unique(geos$pop[order(geos$region2.order)])),fixedpp,
                 fill=type),
             stat='identity',
             position='dodge')+
    labs(x="Population",y="Proportion of fixed sites",fill="Prediction")+
    facet_grid(~region,switch='x',scales="free")+
    background_grid(major="y")
fixedPPPlot

saveplot(fixedPPPlot,"../../plots/fixedProportion_nns_byPop_byAllSitesOverAllTypes_20221205",9,5)
