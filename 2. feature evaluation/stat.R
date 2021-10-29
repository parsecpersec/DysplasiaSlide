#### OED grade distribution of positive patches ####
dat = read.csv('./label.csv', T, stringsAsFactors=F)
dat$dysplasia = gsub(' dysplasia', '', dat$dysplasia)
dat$dysplasia = factor(dat$dysplasia, levels=c('hyperplasia', 'mild',
                                               'moderate', 'severe', 'cancer'))
dat$size = 'large'
dat$size[grepl('C', dat$number)] = 'small'
temp = dat[dat$X1 == 'YY',]

library(ggplot2)
library(reshape2)
for(i in 1:16) {
  temp = dat[dat[,(4+i)] == 'YY',]
  if(length(unique(temp$size)) == 2) {
    tiff(paste0('./barplot of positive patches/', name[i], '.tiff'), res=300, height=1200, width=1800)
    print(ggplot(temp, aes(x=dysplasia, group=size)) + 
            geom_bar(aes(fill=factor(..x..)), stat='count', color='black') +
            theme_bw() + theme(panel.background=element_blank(), legend.position='none',
                               axis.text=element_text(size=15), title=element_text(size=18),
                               plot.title=element_text(size=20, hjust=0.5, family='SimHei'),
                               axis.text.x=element_text(size=13)) +
            geom_text(aes(label = scales::percent(2*..count../sum(..count..))), 
                      stat= "count", vjust = -0.25) +
            ggtitle(name[i]) + xlab('') + ylab('Quantity') + facet_grid(~size))
    dev.off()
  }
}
