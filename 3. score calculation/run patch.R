setwd('')

#### generate patch file - large ####
files = list.files('./soft_patch/temp_large/')
path = unique(substr(files, 1, 9))
soft = unique(substr(files, 11, length(files)))
soft = gsub('[.]csv', '', soft)

dat = data.frame()
dat0 = data.frame()
for(p in path) {
  for(s in soft) {
    filename = paste0('./soft_patch/temp_large/', p, '_', s, '.csv')
    temp = read.csv(filename, header=F, stringsAsFactors=F, row.names='V1')
    if(s == soft[1]) {
      rownames(temp) = paste0(p, '_', rownames(temp))
      dat0 = temp
    } else {
      rownames(temp) = paste0(p, '_', rownames(temp))
      dat0 = cbind(dat0, temp)
    }
  }
  colnames(dat0) = soft
  dat = rbind(dat, dat0)
}

dat$path = substr(rownames(dat), 1, 9)
lab = read.csv('./saved.csv', T, stringsAsFactors=F)
lab = lab[lab$dys_bdkq != 'Cancer',]
dat = dat[dat$path %in% lab$path,]
dat$dys = lab$dys_bdkq[match(dat$path, lab$path)]
write.csv(dat, './soft_patch/large.csv', quote=F, row.names=T)

#### generate patch file - small ####
files = list.files('./soft_patch/temp_small/')
path = unique(substr(files, 1, 9))
soft = unique(substr(files, 11, length(files)))
soft = gsub('[.]csv', '', soft)

dat = data.frame()
dat0 = data.frame()
for(p in path) {
  for(s in soft) {
    filename = paste0('./soft_patch/temp_small/', p, '_', s, '.csv')
    temp = read.csv(filename, header=F, stringsAsFactors=F, row.names='V1')
    if(s == soft[1]) {
      rownames(temp) = paste0(p, '_', rownames(temp))
      dat0 = temp
    } else {
      rownames(temp) = paste0(p, '_', rownames(temp))
      dat0 = cbind(dat0, temp)
    }
  }
  colnames(dat0) = soft
  dat = rbind(dat, dat0)
}

dat$path = substr(rownames(dat), 1, 9)
lab = read.csv('./saved.csv', T, stringsAsFactors=F)
lab = lab[lab$dys_bdkq != 'Cancer',]
dat = dat[dat$path %in% lab$path,]
dat$dys = lab$dys_bdkq[match(dat$path, lab$path)]
write.csv(dat, './soft_patch/small.csv', quote=F, row.names=T)

#### density plot ####
library(ggplot2)
plot_dens = function(index) {
  tiff(paste0('./density_patch/', soft[index], '.tiff'), width=1800, height=1200, res=300)
  print(ggplot(dat, aes(x=dat[,index], color=dys, fill=dys)) + geom_density(size=1, alpha=0.1) +
          theme_bw() + theme(panel.background=element_blank(), title=element_text(size=18),
                             axis.text=element_text(size=15), plot.title=element_text(hjust=0.5)) +
          ggtitle(soft[index]) + xlab('Softmax') + ylab('Density'))
  dev.off()
}
for(i in 1:length(soft)) {
  plot_dens(i)
}

#### try scores ####
library(MASS)
library(ggplot2)
library(ggpubr)
eng = c('irregular epithelial stratification',
        'loss of polarity of basal cells',
        'drop-shaped rete ridges',
        'increased number of mitotic figures',
        'abnormally superficial mitotic figures',
        'premature keratinization in single cells',
        'keratin pearls within rete ridges',
        'loss of epithelial cell cohesion',
        'abnormal variation in nuclear size',
        'abnormal variation in nuclear shape',
        'abnormal variation in cell size',
        'abnormal variation in cell shape',
        'increased N:C ratio',
        'atypical mitotic figures',
        'increased number and size of nucleoli',
        'hyperchromasia')
filename = c('1_large', '2_large', '3_large', '4_large',
             '6_large', '8_small', '9_small', '10_small',
             '11_small', '12_small', '13_large', '16_small')
filename = paste0('X', filename)
series = 9:(8+length(filename))

soft1 = read.csv('./soft_patch/large.csv', T, stringsAsFactors=F)
soft2 = read.csv('./soft_patch/small.csv', T, stringsAsFactors=F)

#### 1. mean ####
# run score.R #
name = 'mean'
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
dat = dat[dat$dys_bdkq != 'Cancer', ]
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
dat$dys_bdkq = factor(dat$dys_bdkq, levels=c("Hyperplasia",
                                             "Mild dysplasia",
                                             "Moderate dysplasia",
                                             "Severe dysplasia"))

# OR assess
OR = data.frame()
for (i in 1:length(filename)) {
  fit = polr(dat$dys_bdkq ~ dat[,(i+8)], Hess=T)
  res = exp(cbind(OR=coef(fit), confint.lm(fit, level=0.95)[1,]))
  res = as.data.frame(res, stringsAsFactors=F)
  res$var = c(paste0('lower_', i), paste0('upper_', i))
  OR = rbind(OR, res)
}

OR$var = 1
OR$var[seq(1, 2*length(files), 2)] = OR$V2[seq(2, 2*length(files), 2)]
OR = OR[-seq(2, 2*length(files), 2),]
rownames(OR) = 1:length(files)
colnames(OR) = c('OR', 'lower', 'upper')
OR$include = ifelse((OR$lower-1)*(OR$upper-1)>0, T, F)
OR$var = colnames(dat)[series]

# OR sorted by size
OR2 = plyr::arrange(OR, var, decreasing=F)
OR2$var = factor(OR2$var, levels=rev(OR2$var))
# OR2[,1:3] = log(OR2[,1:3])  # if OR too large

pd = position_dodge(0.1)
ggplot(OR2, aes(x=var, y=OR, color=include)) + coord_flip() +
  ggtitle(name) + 
  geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  geom_line(position=pd, size = 0.75) + xlab('') + ylab('Odds Ratio') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) + theme_bw() +
  theme(panel.background=element_blank()) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        title = element_text(size = 18, hjust=0.5), legend.position='none',
        plot.title=element_text(hjust=0.5)) +
  scale_color_manual(values=c('grey', 'tomato'))

# calculate score
sum = rep(0, nrow(dat))
n = 0
for (i in 1:length(filename)) {
  if (T) {
    sum = sum + dat[(i+8)] *  OR$OR[(i)]
    n = n + 1
  }
}
dat$score = sum[,1] / n

fit = polr(dys_bdkq~score, Hess=T, data=dat)
summary(fit)
res = exp(cbind(OR=coef(fit), confint(fit)))
res = as.data.frame(res, stringsAsFactors=F)
res$var = rownames(res)

# Score
fit2 = polr(dys_bdkq~score+age+sex+smoke+drink, Hess=T, data=dat)
summary(fit2)
res2 = exp(cbind(OR=coef(fit2), confint(fit2)))
res2 = as.data.frame(res2, stringsAsFactors=F)
res2$var = rownames(res2)

colnames(res2) = c('OR', 'lower', 'upper', 'var')
res2$var = c('')
res2$var = factor(res2$var, levels=res2$var)
res2$mark = c(T, F, F, F, F)
pd = position_dodge(0.1)
ggplot(res2, aes(x=var, y=OR, color=mark)) + geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  ggtitle(name) + xlab('') + ylab('Odds Ratio (OR)') +
  theme_bw() + theme(panel.background=element_blank(), legend.position='none',
                     axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
                     plot.title=element_text(hjust=0.5),
                     title = element_text(size = 18, hjust=0.5)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) + 
  scale_color_manual(values=c('black', 'tomato'))

# extra
label = levels(factor(dat$dys_bdkq))
comps = combn(label, 2)
mycomp = list(comps[,1], comps[,2], comps[,3], comps[,4], comps[,5], comps[,6])
ggplot(dat, aes(x=dys_bdkq, y=score, fill=dys_bdkq)) + geom_boxplot() + 
  stat_compare_means(method='anova', label.y=85) +
  stat_compare_means(method='kruskal', label.x=2, label.y=85) + 
  stat_compare_means(comparisons=mycomp, method='wilcox') + theme_bw() +
  theme(panel.background=element_blank(), legend.position='none',
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle(name) + xlab('') + ylab('Score')

ggplot(dat, aes(x=score, fill=dys_bdkq, color=dys_bdkq)) + geom_density(size=1, alpha=0.1) + 
  theme_bw() + xlim(c(0.6 * min(dat$score), 1.25 * max(dat$score))) +
  theme(panel.background=element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle(name) + xlab('Score') + ylab('Density')

#### 2. ####
