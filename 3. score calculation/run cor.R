setwd('')
# library(GGally)
# library(corrgram)
library(corrplot)
library(ggplot2)
library(factoextra)

#### start ####

#### slide cor ####
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

#### slide PCA ####
slide.pca = prcomp(dat[,series], center=T, scale.=T)
summary(slide.pca)
scree = data.frame(pc=1:length(slide.pca$sdev), Variances=slide.pca$sdev)
ggplot(scree, aes(x=pc, y=Variances)) + geom_line(stat='identity', size=1.2) + 
  geom_point(size=3.5, color='blue') + ggtitle('Screeplot of the Principal Components') +
  scale_x_continuous(breaks=1:12) + theme_bw() + xlab('') +
  geom_hline(yintercept=1, color='red', linetype=2, size=1) +
  theme(plot.title=element_text(size=20, hjust=0.5), axis.title=element_text(size=18),
        axis.text=element_text(size=14))

dat$dysplasia = ifelse(dat$dys_bdkq == 'Hyperplasia', 'Nondysplasia', 'Dysplasia')
fviz_pca_ind(slide.pca, geom.ind='point', pointshape=21, pointsize=2, 
             fill.ind=dat$dysplsia, col.ind='black', palette='jco', 
             addEllipses=T, label='var', col.var='black', repel=T, legend.title='Class') +
  ggtitle("2D PCA-plot from 12-feature dataset") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#### patch cor ####
# rm(list=ls())
corrplot(cor(as.matrix(large[,2:7]), method='pearson'), method='square', order='hclust',
         col=colorRampPalette(c('blue', 'white', 'tomato'))(200), mar=c(0, 0, 1, 0),
         tl.col='black', tl.srt=45, tl.cex=1.5, cl.cex=1.5)

corrplot(cor(as.matrix(small[,2:7]), method='pearson'), method='square', order='hclust',
         col=colorRampPalette(c('blue', 'white', 'tomato'))(200), mar=c(0, 0, 1, 0),
         tl.col='black', tl.srt=45, tl.cex=1.5, cl.cex=1.5)
