setwd('')
library(MASS)
library(ggplot2)


#### part3 ####

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

#### input all ####
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
files = paste0('./soft/', list.files('./soft/'))
series = 9:(8+length(files))
for(i in 1:length(files)) {
  temp = read.csv(files[i], F, stringsAsFactors=F)
  dat = cbind(dat, temp$V2)
}
name = gsub('./soft/Soft_', '', files)
name = gsub('.csv', '', name)
name[grepl('^[1-9]_', name)] = paste0('0', name[grepl('^[1-9]_', name)])
colnames(dat)[series] = name
for(i in series) {
  dat[,i] = 100*dat[,i]    # softmax (0-1) -> (0-100)
}
dat = dat[dat$dys_bdkq != 'Cancer',]
dat$dys_bdkq = factor(dat$dys_bdkq, levels=c("Hyperplasia",
                                             "Mild dysplasia",
                                             "Moderate dysplasia",
                                             "Severe dysplasia"))
#### view feature score ####
library(ggplot2)
library(reshape2)
dat2 = melt(dat[,series])
dat2$size = gsub('.*_', '', dat2$variable)
dat2$size = factor(dat2$size, levels=c('small', 'large'))
lev = name[order(name)]
dat2$variable = factor(dat2$variable, levels=rev(lev))
ggplot(dat2, aes(x=variable, y=value)) + geom_boxplot(aes(fill=size), size=0.1, outlier.shape=NA) + 
  coord_flip() + theme_bw() + xlab('Scores of the Features') + ylab('') +
  ggtitle('Distribution of Feature Scores') +
  theme(panel.background=element_blank(), axis.text=element_text(size = 15, face='bold'),
        title=element_text(size = 24, hjust=0.5),
        legend.text=element_text(size=15, face='bold'))

#### OR assess ####
OR = data.frame()
for (i in 1:length(files)) {
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
  ggtitle(paste0('Ordinal Logistic Regression Analyses of OED Features')) + 
  geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  geom_line(position=pd, size = 0.75) + xlab('') + ylab('Odds Ratio') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) + theme_bw() +
  theme(panel.background=element_blank()) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        title = element_text(size = 18, hjust=0.5), legend.position='none') +
  scale_color_manual(values=c('grey', 'tomato'))
  
#### part3 fixed features ####
# included features
features = c(paste0(eng[1], ' - 1024'),
             paste0(eng[2], ' - 1024'),
             paste0(eng[3], ' - 1024'),
             paste0(eng[4], ' - 1024'),
             paste0(eng[6], ' - 1024'),
             paste0(eng[8], ' - 224'),
             paste0(eng[9], ' - 224'),
             paste0(eng[10], ' - 224'),
             paste0(eng[11], ' - 224'),
             paste0(eng[12], ' - 224'),
             paste0(eng[13], ' - 1024'),
             paste0(eng[16], ' - 224'))
filename = c('1_large', '2_large', '3_large', '4_large',
             '6_large', '8_small', '9_small', '10_small',
             '11_small', '12_small', '13_large', '16_small')
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
files = paste0('./soft/Soft_', filename, '.csv')
series = 9:(8+length(files))
for(i in 1:length(files)) {
  temp = read.csv(files[i], F, stringsAsFactors=F)
  dat = cbind(dat, temp$V2)
}
name = gsub('./soft/Soft_', '', files)
name = gsub('.csv', '', name)
name[grepl('^[1-9]_', name)] = paste0('0', name[grepl('^[1-9]_', name)])
colnames(dat)[series] = name
for(i in series) {
  dat[,i] = 100*dat[,i]    # softmax (0-1) -> (0-100)
}
dat = dat[dat$dys_bdkq != 'Cancer',]
dat$dys_bdkq = factor(dat$dys_bdkq, levels=c("Hyperplasia",
                                             "Mild dysplasia",
                                             "Moderate dysplasia",
                                             "Severe dysplasia"))

# OR assess
OR = data.frame()
for (i in 1:length(files)) {
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

# OR sorted by name and sig
OR2 = plyr::arrange(OR, include, var, decreasing=F)
OR2$var = factor(OR2$var, levels=rev(OR2$var))
OR2[,1:3] = log(OR2[,1:3])  # if OR too large

# OR sorted by size
OR2 = plyr::arrange(OR, var, decreasing=F)
OR2$var = factor(OR2$var, levels=rev(OR2$var))
OR2[,1:3] = log(OR2[,1:3])  # if OR too large

pd = position_dodge(0.1)
ggplot(OR2, aes(x=var, y=OR, color=include)) + coord_flip() +
  ggtitle(paste0('Ordinal Logistic Regression Analyses of OED Features')) + 
  geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  geom_line(position=pd, size = 0.75) + xlab('') + ylab('Odds Ratio (log scale)') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) + theme_bw() +
  theme(panel.background=element_blank()) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        title = element_text(size = 18, hjust=0.5), legend.position='none') +
  scale_color_manual(values=c('grey', 'tomato'))

#### calculate score ####
sum = rep(0, nrow(dat))
n = 0
for (i in 1:length(files)) {
  if (T) {
    sum = sum + dat[(i+8)] *  OR$OR[(i)]
    n = n + 1
  }
}
dat$score = sum[,1] / n

colnames(res2) = c('OR', 'lower', 'upper', 'var')
res2$var = c('')  # rename
res2$var = factor(res2$var, levels=res2$var)
res2$mark = c(T, F, F, F, F)
pd = position_dodge(0.1)
ggplot(res2, aes(x=var, y=OR, color=mark)) + geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  ggtitle('') + xlab('') + ylab('Odds Ratio (OR)') +
  theme_bw() + theme(panel.background=element_blank(), legend.position='none',
                     axis.text.x = element_text(size = 15),
                     axis.text.y = element_text(size = 15),
                     title = element_text(size = 18, hjust=0.5)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) + 
  scale_color_manual(values=c('black', 'tomato'))

#### extra ####
library(ggpubr)
label = levels(factor(dat$dys_bdkq))
comps = combn(label, 2)
mycomp = list(comps[,1], comps[,2], comps[,3], comps[,4], comps[,5], comps[,6])
ggplot(dat, aes(x=dys_bdkq, y=score)) + geom_boxplot() + stat_compare_means(method='anova', label.y=85) +
  stat_compare_means(method='kruskal', label.x=2, label.y=85) + 
  stat_compare_means(comparisons=mycomp, method='wilcox')

#### loop score loop boxplot ####
for(i in series) {
  fit0 = orm(dys_bdkq~dat[,i]+age+sex+smoke+drink, data=dat)
  res0 = exp(cbind(OR=coef(fit0), confint(fit0)))
  res0 = as.data.frame(res0, stringsAsFactors=F)
  colnames(res0) = c('OR', 'lower', 'upper')
  res0 = res0[4:8,]
  res0$var = c('')  # rename
  res0$var = factor(res0$var, levels=res0$var)
  res0$mark = c(T, F, F, F, F)
  pd = position_dodge(0.1)
  tiff(paste0('./multi OR loop/', colnames(dat)[i], '.tiff'), width=1800, height=1200, res=300)
  print(ggplot(res0, aes(x=var, y=OR, color=mark)) + geom_point(position=pd, size=3, shape=15) + 
          geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
          ggtitle(colnames(dat)[i]) + xlab('') + ylab('Odds Ratio (OR)') +
          theme_bw() + theme(panel.background=element_blank(), legend.position='none',
                             axis.text.x = element_text(size = 15),
                             axis.text.y = element_text(size = 15),
                             title = element_text(size = 18, hjust=0.5)) +
          geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) + 
          scale_color_manual(values=c('black', 'tomato')))
  dev.off()
}

library(ggpubr)
label = levels(factor(dat$dys_bdkq))
comps = combn(label, 2)
mycomp = list(comps[,1], comps[,2], comps[,3], comps[,4], comps[,5], comps[,6])
for(i in series) {
  tiff(paste0('./boxplot loop/', colnames(dat)[i], '.tiff'), width=1800, height=1200, res=300)
  print(ggplot(dat, aes(x=dys_bdkq, y=dat[,i])) + geom_boxplot() + ggtitle(colnames(dat)[i]) +
          stat_compare_means(method='anova', label.y=0.9*min(dat[,i])) + xlab('') +
          stat_compare_means(method='kruskal', label.x=2, label.y=0.9*min(dat[,i])) + 
          stat_compare_means(comparisons=mycomp, method='wilcox'))
  dev.off()
}

#### density ####
ggplot(dat, aes(x=score, color=dys_bdkq)) + geom_density(size=3) + theme_bw() +
  theme(panel.background=element_blank(), title=element_text(size=18), 
        axis.text=element_text(size=15), plot.title=element_text(hjust=0.5)) + 
  ggtitle('Density') + xlab('Score') + ylab('')
