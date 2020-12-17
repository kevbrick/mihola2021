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

theme7point()

z <- read.table('allCoverage.Rdata.tab',header=TRUE) %>% 
  mutate(strain = ifelse(grepl('WKY',sample),'WKY',ifelse(grepl('RIJH',sample),'BN/RIJHsd','SHR')),
         sample = gsub("_"," ",gsub("BN","BN/",sample)),
         from=from/1000000,
         to=to/1000000) %>%
  filter(!grepl('Input',sample)) 

 z$label <- paste0("'",z$sample,"'")

 z$label[z$sample == 'SHR DMC1 SSDS Prdm9KO'] <- '"SHR DMC1 SSDS "*Prdm9^"-/-"'
 z$label[z$sample == 'SHR DMC1 SSDS Prdm9KOHet'] <- '"SHR DMC1 SSDS "*Prdm9^"+/-"'
 
 z$sample <- factor(z$sample,levels=c('BN/RIJHsd H3K4m3 Liver',
                                      'BN/RIJHsd H3K4m3 Testis',
                                      'WKY H3K4m3 Testis',
                                      'SHR DMC1 SSDS Prdm9KO',
                                      'SHR DMC1 SSDS Prdm9KOHet',
                                      'SHR DMC1 SSDS wildtype',
                                      'WKY DMC1 SSDS wildtype',
                                      'BN/RIJHsd DMC1 SSDS wildtype'))

xMin <- min(z$from)
xMax <- max(z$to)

# for (s in unique(z$sample)){
#   z$coverage[z$sample == s] <- z$coverage[z$sample==s]/sum(z$coverage[z$sample==s])*100
# }

gCoverage <- ggplot(z,aes(x=(from+to)/2,y=coverage)) + 
  geom_area(aes(fill=strain),color='black',lwd=.2) + 
  facet_wrap(~sample,ncol=1,scales='free_y') + 
  geom_text(x=-Inf,y=Inf,
            size=8*5/14,
            hjust=-0.1,vjust=1,
            fontface='bold',
            check_overlap=TRUE,
            aes(label=label),parse=TRUE)+
  theme(legend.position='none',
        strip.background = element_blank(),
        strip.text=element_blank()) +
  xlab('Position on chr 17') + 
  ylab('Coverage (read depth)') + 
  scale_fill_manual(values=c('chocolate1','maroon1','dodgerblue2')) + 
  xlab('') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        panel.spacing = unit(0.75, "lines")) + 
  coord_cartesian(xlim=c(xMin,xMax),expand=FALSE) 

dfGTF <- read.table('region.gtf.forR.tab')
names(dfGTF) <- c('type','from','to','dir','hs')

dfGTF$from <- dfGTF$from/1000000
dfGTF$to <- dfGTF$to/1000000

dfT <- dfGTF %>% filter(type=='transcript') %>% mutate(ypos=1)
for (t in 2:length(dfT)){
  if (sum((dfT$to-dfT$from[3])<0)){
    dfT$ypos[t] <- dfT$ypos[t-1]-1
  }
}

dfE <- dfGTF %>% filter(type %in% c('start_codon','CDS','exon','3UTR','5UTR','transcript')) %>% mutate(ypos=1)
t<-0
for (i in 1:length(dfE$type)){
  if (dfE$type[i] == 'transcript'){
    t <- t+1
  }else{
    dfE$ypos[i] <- dfT$ypos[t]
  }
}

yTri <- min(dfT$ypos)-1.5

gGenes <- ggplot() + 
  geom_segment(data=dfT %>% filter(dir=="+"),lwd=.3,
               aes(x=from,xend=from+(to-from)*.33,y=ypos,yend=ypos),arrow = arrow(length = unit(0.05, "inches"))) +
  geom_segment(data=dfT %>% filter(dir=="+"),lwd=.3,
               aes(x=from+(to-from)*.33,xend=from+(to-from)*.66,y=ypos,yend=ypos),arrow = arrow(length = unit(0.05, "inches"))) + 
  geom_segment(data=dfT %>% filter(dir=="+"),lwd=.3,
               aes(x=from+(to-from)*.66,xend=to,y=ypos,yend=ypos)) + 
  geom_segment(data=dfT %>% filter(dir=="-"),lwd=.3,
               aes(xend=from,x=from+(to-from)*.33,y=ypos,yend=ypos)) +
  geom_segment(data=dfT %>% filter(dir=="-"),lwd=.3,
               aes(xend=from+(to-from)*.33,x=from+(to-from)*.66,y=ypos,yend=ypos),arrow = arrow(length = unit(0.05, "inches"))) + 
  geom_segment(data=dfT %>% filter(dir=="-"),lwd=.3,
               aes(xend=from+(to-from)*.66,x=to,y=ypos,yend=ypos),arrow = arrow(length = unit(0.05, "inches"))) + 
  geom_point(data=dfE %>% filter(type %in% c('start_codon')),
            aes(x=from,fill=hs>0),y=yTri,size=3,shape=24) +
  geom_rect(data=dfE %>% filter(type %in% c('exon')),
               aes(xmin=from-0.0005,
                   xmax=to+0.0005,
                   ymin=ypos+.3,
                   ymax=ypos-.3),lwd=.2)+
  geom_rect(data=dfE %>% filter(grepl('UTR',type)),
            aes(xmin=from-0.0005,
                xmax=to+0.0005,
                ymin=ypos+.5,
                ymax=ypos-.5),lwd=.2) +
  ylab('Genes') + xlab('Position on chr17 (Mb)') +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position='none') + 
  coord_cartesian(xlim=c(xMin,xMax),
                  ylim=c(yTri-1,1.5),
                  expand=FALSE) +  
  scale_fill_manual(values=c('pink','green')) 

gPlot <- ggarrange(gCoverage ,
          gGenes,ncol=1,heights=c(6,1),
          align='v')

ggsave('Mihola_et_al_2021_rn5Coverage.png',gPlot,height=7,width=7,dpi=500)
ggsave('Mihola_et_al_2021_rn5Coverage.pdf',gPlot,height=7,width=7,dpi=500)