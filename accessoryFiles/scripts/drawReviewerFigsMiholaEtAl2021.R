library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)

#############################
### theme7point
### KB June 29 2018
### Set the default theme to 7-point font
### suitable for publication
### ARGS: none
### OUTPUTS: Sets the theme
theme7point <- function() {
  theme_base_size <- 7
  theme_set(theme_bw(base_size = theme_base_size) %+replace%
              theme(axis.text = element_text(size=theme_base_size,
                                             color='black'),
                    panel.grid = element_blank(),
                    panel.border=element_blank(),
                    axis.line=element_line(size=.2),
                    axis.ticks = element_line(size=.2)))
}

#############################
### fancy_scientific
### KB July 10 2018
### Change number to scientific notation format:
### Most useful for scale_x[y]_log10(labels=fancy_scientific)
### ARGS:
# l   number
## RETURNS: Formatted number as expression
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # remove +s
  l <- gsub("\\+", "", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

#############################
### toTPM
### KB July 03 2018
### Convert a vector of strengths into a Tags (or Fragments) per million value
### ARGS:
# x           vector of values
# noNeg       remove negative values by adding the minimum FPM to all
## RETURNS: Normalized vector in T(F)PM
toTPM <- function (x,noNeg=FALSE){
  #x <- x+abs(min(x))+1;
  v <- (x/sum(x)*1000000)
  if (noNeg & min(v) < 0){
    v <- v - min(v) + 1
  }
  return(v)
}

#############################
### plotHSstrength
### KB July 23 2018
### Make a density scatterplot of hotspot strength in two samples
### ARGS:
# sA      vector of values for x-axis
# sB      vector of values for y-axis
# nameA   Name for x-axis
# nameB   Name for y-axis
#
## RETURNS: geom
plotHSstrength <- function (sA,sB,nameA,nameB,
                            noGradient = FALSE,
                            sOK        = NULL,
                            lblAsIs    = FALSE,
                            topLeft    = FALSE,
                            xLim       = 0,
                            yLim       = 0,
                            xLimits    = NULL,
                            yLimits    = NULL,
                            stdStyle   = TRUE,
                            pointSz    = 2,
                            txtSz      = 7,
                            noLeg      = FALSE,
                            labelSize  = 8){
  
  sDataInit <- data.frame('tpmA'=toTPM(sA),'tpmB'=toTPM(sB))
  
  if (length(sOK) > 0){
    sData <- sDataInit[sOK>0,]
  }else{
    sData <- sDataInit
  }
  
  ## Get max / min coords
  xMax <- max(sData$tpmA)
  xMin <- min(sData$tpmA)
  xMinLbl <- min(sData$tpmA[sData$tpmA>0])
  
  yMax <- max(sData$tpmB)
  yMin <- min(sData$tpmB)
  
  ## Set limits if we're using them
  if (xLim){
    xMax <- xLim
  }
  
  if (yLim){
    yMax <- yLim
  }
  
  if (!is.null(xLimits)){
    xMin <- xLimits[1]
    xMax <- xLimits[2]
  }
  
  if (!is.null(yLimits)){
    yMin <- yLimits[1]
    yMax <- yLimits[2]
  }
  
  sCC <- round(cor(sData[,1:2],method='spearman')[1,2]^2,2)
  
  if (lblAsIs){
    xLblName <- nameA
    yLblName <- nameB
  }else{
    xLblName <- paste0(nameA,' ','(FPM)')
    yLblName <- paste0(nameB,' ','(FPM)')
  }
  
  okData <- sData[sData$tpmA>xMin & sData$tpmA<xMax &
                    sData$tpmB>yMin & sData$tpmB<yMax,]
  
  if (stdStyle){
    gInit <- ggplot(data=okData,aes(x=tpmA,y=tpmB)) +
      geom_point(color='grey50',size=pointSz) +
      scale_x_log10(labels=fancy_scientific) +
      scale_y_log10(labels=fancy_scientific) +
      geom_smooth(method=lm,linetype=2,colour="NA",se=F) +
      guides(alpha="none") +
      annotation_logticks(sides='lb',
                          size=.2,
                          short=unit(0.050,'cm'),
                          mid=unit(0.075,'cm'),
                          long=unit(0.100,'cm')) +
      coord_cartesian(xlim=c(xMin,max(xMax,yMax)),
                      ylim=c(xMin,max(xMax,yMax))) +
      xlab(xLblName) +
      ylab(yLblName) +
      theme(legend.position=c(1,0),
            legend.justification=c(1,0))
  }else{
    gInit <- ggplot(data=okData,aes(x=tpmA,y=tpmB)) +
      geom_point(color='grey50',size=pointSz) +
      scale_x_log10(labels=fancy_scientific) +
      scale_y_log10(labels=fancy_scientific) +
      geom_smooth(method=lm,linetype=2,colour="NA",se=F) +
      guides(alpha="none") +
      annotation_logticks(sides='lb') +
      coord_cartesian(xlim=c(xMin,max(xMax,yMax)),
                      ylim=c(xMin,max(xMax,yMax))) +
      theme_MF() +
      xlab(xLblName) +
      ylab(yLblName) +
      theme(plot.title=element_text(size=20),
            legend.position=c(1,0),
            legend.justification=c(1,0),
            legend.title=element_text(size=15))
  }
  #ggtitle(substitute(paste('Spearman ', R^"2"," = ", cc ,sep=''),list(cc = sCC)))
  
  if (noGradient){
    g <- gInit
  }else{
    g <- gInit + stat_density2d(aes(fill=..level..,
                                    alpha=..level..),
                                geom='polygon',
                                colour='NA')
    #scale_fill_continuous(low="grey30",high="red")
  }
  
  g$labels$fill <- "Hotspot density"
  
  myLbl <- substitute(paste('Spearman ', R^"2"," = ", cc ,sep=''),list(cc = sCC));
  myLbl <- paste("  R^2 == ", sCC)
  
  gg <- g + geom_text(label=myLbl,
                      x=-Inf,y=Inf,
                      check_overlap = TRUE,
                      hjust=-0.5,vjust=1,
                      size=labelSize*5/14,
                      parse=TRUE) +
    theme(axis.line = element_line(size=.2),
          legend.background=element_blank())
  
  if (noLeg){
    gg <- gg + theme(legend.position='none')
  }
  
  return(gg)
}


###############################
theme7point()
###############################

## Import hotspots table
dfHotspots <- read.table('allMergedRatHotspots.finalTable.tab',header=TRUE)
names(dfHotspots) <- c('chr','from','to','SHR(wt)','SHR(het)','SHR(ko)','WKY','BN','hotspotSHR(wt)','hotspotSHR(het)','hotspotSHR(ko)','hotspotWKY','hotspotBN')

dfHotspots$scalefactor <- 1
dfHotspots$scalefactor[dfHotspots$chr == 'chrX'] <- 2
dfHotspots$sexCS <- 'Auto'
dfHotspots$sexCS[dfHotspots$chr == 'chrX'] <- 'chrX'

## Make plot dataframe
dfMeltedHS <- reshape2::melt(dfHotspots, id.vars = c('chr','from','to','sexCS','scalefactor'),
                             measure.vars=c('hotspotSHR(wt)','hotspotSHR(het)','hotspotSHR(ko)','hotspotWKY','hotspotBN')) %>%
  mutate(variable = gsub('hotspot','',variable))

dfMeltedStr <- reshape2::melt(dfHotspots, id.vars = c('chr','from','to','sexCS','scalefactor'),
                              measure.vars=c('SHR(wt)','SHR(het)','SHR(ko)','WKY','BN')) %>%
  inner_join(dfMeltedHS,by=c('chr','from','to','sexCS','scalefactor','variable')) %>%
  rename(strain = variable ,strength = value.x, hotspot = value.y) 

## ChrX V Autosomes plot
gXHSAreStrong = ggplot(dfMeltedStr %>% filter(hotspot >0), 
       aes(y=sexCS,x=strength*scalefactor)) + 
  geom_violin(aes(fill=sexCS),lwd=.2) + 
  geom_boxplot(notch=TRUE,width=.2,outlier.size=.1,lwd=.2) + 
  facet_wrap(~strain,ncol=3) + 
  scale_x_log10() + 
  annotation_logticks(sides='b',
                      long  = unit(0.2,'cm'),
                      mid   = unit(0.1,'cm'),
                      short = unit(0.1,'cm'),
                      size=.2) + 
  scale_fill_manual(values=c('darkorange','grey50')) +
  geom_text(aes(label=paste0("  ",strain)),
            x=-Inf,y=Inf,
            hjust=0,vjust=1,
            size=8*5/14,
            check_overlap=TRUE,
            fontface='bold') + 
  ylab('') + 
  xlab('Hotspot strength (A.U.)') + 
  theme(legend.position='none',
        strip.background = element_blank(),
        strip.text=element_blank())

ggsave('Mihola_et_al_2021_strongHSonX.png',plot = gXHSAreStrong, height=4,width=5, dpi=500)
ggsave('Mihola_et_al_2021_strongHSonX.pdf',plot = gXHSAreStrong, height=4,width=5)

###########################################################################################################################################################
# Count default HS
dfHotspots$default <- dfHotspots$'hotspotSHR(ko)' * dfHotspots$'hotspotSHR(wt)' * dfHotspots$'hotspotWKY' * dfHotspots$'hotspotBN'

print (paste0(sum(dfHotspots$default)," default hotspots shared by all three strains ..."))

###########################################################################################################################################################
## Make all correlation plots
corrMat  <- matrix(data = 0,nrow=5,ncol=5)
plotLst <- list()
n <- 0
for (ni in 4:8){
  for (nj in 4:8){
    n <- n+1
    i<- names(dfHotspots)[ni]
    j<- names(dfHotspots)[nj]
    
    iHS <- dfHotspots[[i]][dfHotspots$chr != 'chrX' & dfHotspots[[paste0('hotspot',i)]]*dfHotspots[[paste0('hotspot',j)]] >0]
    jHS <- dfHotspots[[j]][dfHotspots$chr != 'chrX' & dfHotspots[[paste0('hotspot',i)]]*dfHotspots[[paste0('hotspot',j)]] >0]
    
    plotLst[[n]] <- plotHSstrength(iHS,jHS,i,j,pointSz = .1,noLeg = TRUE,labelSize = 8)
    
    corrMat[ni-3,nj-3] <- cor(iHS,jHS,method='spearman')
  }
}

gAllScatters <- ggarrange(plotlist = plotLst)

ggsave('Mihola_et_al_2021_allSSDSScatterplots.png',plot = gAllScatters, height=8,width=8, dpi=500)
ggsave('Mihola_et_al_2021_allSSDSScatterplots.pdf',plot = gAllScatters, height=8,width=8)


