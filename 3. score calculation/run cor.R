setwd('')
# library(ggplot2)
# library(GGally)
# library(corrgram)
library(corrplot)

#### start ####

#### slide level ####
name = 'mean'
for(f in filename) {
  if(grepl('large', f)) {
    temp = aggregate(soft1[, f], FUN=mean, by=list(path=soft1$path))
    dat = cbind(dat, temp$x)
  } else {
    temp = aggregate(soft2[, f], FUN=mean, by=list(path=soft2$path))
    dat = cbind(dat, temp$x)
  }
}

colnames(dat)[series] = filename
for(i in series) {
  dat[,i] = 100*dat[,i]    # (0-1) -> (0-100)
}

# heatmap(as.matrix(dat[,series]), scale='row')
# ggcorr(dat[,series], method=c('pairwise', 'pearson'))
# corrgram(dat[,series])
corrplot(cor(as.matrix(dat[,series]), method='pearson'), method='square', order='hclust',
         col=colorRampPalette(c('blue', 'white', 'tomato'))(200), mar=c(0, 0, 1, 0),
         tl.col='black', tl.srt=45, tl.cex=1.5, cl.cex=1.5)

#### patch level ####
# rm(list=ls())
corrplot(cor(as.matrix(large[,2:7]), method='pearson'), method='square', order='hclust',
         col=colorRampPalette(c('blue', 'white', 'tomato'))(200), mar=c(0, 0, 1, 0),
         tl.col='black', tl.srt=45, tl.cex=1.5, cl.cex=1.5)

corrplot(cor(as.matrix(small[,2:7]), method='pearson'), method='square', order='hclust',
         col=colorRampPalette(c('blue', 'white', 'tomato'))(200), mar=c(0, 0, 1, 0),
         tl.col='black', tl.srt=45, tl.cex=1.5, cl.cex=1.5)
