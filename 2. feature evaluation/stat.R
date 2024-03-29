#### OED grade distribution of positive patches ####
dat1 = read.csv('./patch-C.csv', T, stringsAsFactors=F)
name = colnames(dat1)[2:17]
name = gsub('X[0-9]{1,2}', '', name)
name = paste0(1:16, '. ', name)
rm(dat1)
dat = read.csv('./label.csv', T, stringsAsFactors=F)
dat$dysplasia = gsub(' dysplasia', '', dat$dysplasia)
dat$dysplasia = factor(dat$dysplasia, levels=c('hyperplasia', 'mild',
                                               'moderate', 'severe', 'cancer'))
dat$size = 'large'
dat$size[grepl('C', dat$number)] = 'small'
temp = dat[dat$X1 == 'YY',]

library(ggplot2)
library(reshape2)
library(extrafont)
for(i in 1:16) {
  temp = dat[dat[,(4+i)] == 'YY',]
  if(nrow(temp) != 0) {
    if(length(unique(temp$size)) == 2) {
      tiff(paste0('./barplot of positive patches/', name[i], '.tiff'), res=300, height=1200, width=1800)
      print(ggplot(temp, aes(x=dysplasia, group=size)) + 
              geom_bar(aes(fill=factor(..x..)), stat='count', color='black') +
              theme_bw() + theme(panel.background=element_blank(), legend.position='none',
                                 axis.text=element_text(size=15), title=element_text(size=18),
                                 plot.title=element_text(size=20, hjust=0.5, family='SimHei'),
                                 axis.text.x=element_text(size=13)) +
              geom_text(aes(label = scales::percent(..prop.., accuracy=0.1)), 
                        stat= "count", vjust = -0.25) +
              scale_x_discrete(labels=c('1', '2', '3', '4', '5')) +
              ggtitle(name[i]) + xlab('') + ylab('Quantity') + facet_grid(~size))
      dev.off()
    } else if(unique(temp$size) == 'large') {
      tiff(paste0('./barplot of positive patches/', name[i], '.tiff'), res=300, height=1200, width=1800)
      print(ggplot(temp, aes(x=dysplasia)) + 
              geom_bar(aes(fill=factor(..x..)), stat='count', color='black') +
              theme_bw() + theme(panel.background=element_blank(), legend.position='none',
                                 axis.text=element_text(size=15), title=element_text(size=18),
                                 plot.title=element_text(size=20, hjust=0.5, family='SimHei'),
                                 axis.text.x=element_text(size=13)) +
              geom_text(aes(label = scales::percent(..count../sum(..count..), accuracy=0.1)), 
                        stat= "count", vjust = -0.25) +
              scale_x_discrete(labels=c('1', '2', '3', '4', '5')) +
              ggtitle(paste0(name[i], ' - large')) + xlab('') + ylab('Quantity'))
      dev.off()
    } else if(unique(temp$size) == 'small') {
      tiff(paste0('./barplot of positive patches/', name[i], '.tiff'), res=300, height=1200, width=1800)
      print(ggplot(temp, aes(x=dysplasia)) + 
              geom_bar(aes(fill=factor(..x..)), stat='count', color='black') +
              theme_bw() + theme(panel.background=element_blank(), legend.position='none',
                                 axis.text=element_text(size=15), title=element_text(size=18),
                                 plot.title=element_text(size=20, hjust=0.5, family='SimHei'),
                                 axis.text.x=element_text(size=13)) +
              geom_text(aes(label = scales::percent(..count../sum(..count..), accuracy=0.1)), 
                        stat= "count", vjust = -0.25) +
              scale_x_discrete(labels=c('1', '2', '3', '4', '5')) +
              ggtitle(paste0(name[i], ' - small')) + xlab('') + ylab('Quantity'))
      dev.off()
    }
  }
}
