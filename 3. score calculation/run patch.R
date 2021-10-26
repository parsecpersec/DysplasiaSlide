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
OR$var[seq(1, 2*length(filename), 2)] = OR$V2[seq(2, 2*length(filename), 2)]
OR = OR[-seq(2, 2*length(filename), 2),]
rownames(OR) = 1:length(filename)
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

#### 2. CV ####
name = 'CV'
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
dat = dat[dat$dys_bdkq != 'Cancer', ]
cv <- function(x) {return(sd(x)/mean(x))}
for(f in filename) {
  if(grepl('large', f)) {
    temp = aggregate(soft1[, f], FUN=cv, by=list(path=soft1$path))
    dat = cbind(dat, temp$x)
  } else {
    temp = aggregate(soft2[, f], FUN=cv, by=list(path=soft2$path))
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
OR$var[seq(1, 2*length(filename), 2)] = OR$V2[seq(2, 2*length(filename), 2)]
OR = OR[-seq(2, 2*length(filename), 2),]
rownames(OR) = 1:length(filename)
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
  theme_bw() + xlim(c(0 * min(dat$score), 1.25 * max(dat$score))) +
  theme(panel.background=element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle(name) + xlab('Score') + ylab('Density')

#### 3. mean + SD ####
name = 'mean + sd'
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
dat = dat[dat$dys_bdkq != 'Cancer', ]
ms <- function(x) {return(sd(x) + mean(x))}
for(f in filename) {
  if(grepl('large', f)) {
    temp = aggregate(soft1[, f], FUN=ms, by=list(path=soft1$path))
    dat = cbind(dat, temp$x)
  } else {
    temp = aggregate(soft2[, f], FUN=ms, by=list(path=soft2$path))
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
OR$var[seq(1, 2*length(filename), 2)] = OR$V2[seq(2, 2*length(filename), 2)]
OR = OR[-seq(2, 2*length(filename), 2),]
rownames(OR) = 1:length(filename)
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

#### 4. top(75) - top(25) ####
name = 'top75 - top25'
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
dat = dat[dat$dys_bdkq != 'Cancer', ]
q <- function(x) {return(quantile(x, 0.75) - quantile(x, 0.25))}
for(f in filename) {
  if(grepl('large', f)) {
    temp = aggregate(soft1[, f], FUN=q, by=list(path=soft1$path))
    dat = cbind(dat, temp$x)
  } else {
    temp = aggregate(soft2[, f], FUN=q, by=list(path=soft2$path))
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
OR$var[seq(1, 2*length(filename), 2)] = OR$V2[seq(2, 2*length(filename), 2)]
OR = OR[-seq(2, 2*length(filename), 2),]
rownames(OR) = 1:length(filename)
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
  stat_compare_means(method='anova', label.y=40) +
  stat_compare_means(method='kruskal', label.x=2, label.y=40) + 
  stat_compare_means(comparisons=mycomp, method='wilcox') + theme_bw() +
  theme(panel.background=element_blank(), legend.position='none',
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle(name) + xlab('') + ylab('Score')

ggplot(dat, aes(x=score, fill=dys_bdkq, color=dys_bdkq)) + geom_density(size=1, alpha=0.1) + 
  theme_bw() + xlim(c(0 * min(dat$score), 1.25 * max(dat$score))) +
  theme(panel.background=element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle(name) + xlab('Score') + ylab('Density')

#### 5. mean( > 0.5) ####
name = 'mean( > 0.5)'
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
dat = dat[dat$dys_bdkq != 'Cancer', ]
mean2 <- function(x) {return(ifelse(max(x) > 0.5, mean(x[x > 0.5]), 0))}
for(f in filename) {
  if(grepl('large', f)) {
    temp = aggregate(soft1[, f], FUN=mean2, by=list(path=soft1$path))
    dat = cbind(dat, temp$x)
  } else {
    temp = aggregate(soft2[, f], FUN=mean2, by=list(path=soft2$path))
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
  if(mean(dat[,(8+i)]) != 0) {
    fit = polr(dat$dys_bdkq ~ dat[,(i+8)], Hess=T)
    res = exp(cbind(OR=coef(fit), confint.lm(fit, level=0.95)[1,]))
    res = as.data.frame(res, stringsAsFactors=F)
    res$var = c(paste0('lower_', i), paste0('upper_', i))
    OR = rbind(OR, res)
  } else {
    res = data.frame(OR=1, V2=1, var=c(paste0('lower_', i), paste0('upper_', i)))
    OR = rbind(OR, res)
  }
}

OR$var = 1
OR$var[seq(1, 2*length(filename), 2)] = OR$V2[seq(2, 2*length(filename), 2)]
OR = OR[-seq(2, 2*length(filename), 2),]
rownames(OR) = 1:length(filename)
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
  theme_bw() + xlim(c(0 * min(dat$score), 1.25 * max(dat$score))) +
  theme(panel.background=element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle(name) + xlab('Score') + ylab('Density')

#### 6. mean( > median) ####
name = 'mean( > median)'
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
dat = dat[dat$dys_bdkq != 'Cancer', ]
mean3 <- function(x) {return(mean(x[x > median(x)]))}
for(f in filename) {
  if(grepl('large', f)) {
    temp = aggregate(soft1[, f], FUN=mean3, by=list(path=soft1$path))
    dat = cbind(dat, temp$x)
  } else {
    temp = aggregate(soft2[, f], FUN=mean3, by=list(path=soft2$path))
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
  if(mean(dat[,(8+i)]) != 0) {
    fit = polr(dat$dys_bdkq ~ dat[,(i+8)], Hess=T)
    res = exp(cbind(OR=coef(fit), confint.lm(fit, level=0.95)[1,]))
    res = as.data.frame(res, stringsAsFactors=F)
    res$var = c(paste0('lower_', i), paste0('upper_', i))
    OR = rbind(OR, res)
  } else {
    res = data.frame(OR=1, V2=1, var=c(paste0('lower_', i), paste0('upper_', i)))
    OR = rbind(OR, res)
  }
}

OR$var = 1
OR$var[seq(1, 2*length(filename), 2)] = OR$V2[seq(2, 2*length(filename), 2)]
OR = OR[-seq(2, 2*length(filename), 2),]
rownames(OR) = 1:length(filename)
colnames(OR) = c('OR', 'lower', 'upper')
OR$include = ifelse((OR$lower-1)*(OR$upper-1)>0, T, F)
OR$var = colnames(dat)[series]

# OR sorted by size
OR2 = plyr::arrange(OR, var, decreasing=F)
OR2$var = factor(OR2$var, levels=rev(OR2$var))
OR2[,1:3] = log(OR2[,1:3])  # if OR too large

pd = position_dodge(0.1)
ggplot(OR2, aes(x=var, y=OR, color=include)) + coord_flip() +
  ggtitle(name) + 
  geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  geom_line(position=pd, size = 0.75) + xlab('') + ylab('Odds Ratio (log scale)') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) + theme_bw() +
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
  stat_compare_means(method='anova', label.y=105) +
  stat_compare_means(method='kruskal', label.x=2, label.y=105) + 
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

#### model: case mean ####
library(MASS)
library(pROC)
library(ggplot2)

dat0 = dat  # backup

# uni
summary(fit)
prob = predict.glm(fit, type=c('response'))
dat = cbind(dat, prob)
dat$OED1 = ifelse(dat$dys_bdkq == label[1], label[1], 'Others')
dat$OED2 = ifelse(dat$dys_bdkq == label[2], label[2], 'Others')
dat$OED3 = ifelse(dat$dys_bdkq == label[3], label[3], 'Others')
dat$OED4 = ifelse(dat$dys_bdkq == label[4], label[4], 'Others')
a = roc(response=dat$OED1, predictor=dat$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat$OED2, predictor=dat$`Mild dysplasia`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat$OED3, predictor=dat$`Moderate dysplasia`, levels=c(label[3], 'Others'), auc=T, ci=T)
d = roc(response=dat$OED4, predictor=dat$`Severe dysplasia`, levels=c(label[4], 'Others'), auc=T, ci=T)
roc_data = data.frame(sen=a$sensitivities, spe=a$specificities, Category='Hyperplasia')
roc_data = rbind(roc_data, data.frame(sen=b$sensitivities, spe=b$specificities, 
                                      Category='Mild dysplasia'))
roc_data = rbind(roc_data, data.frame(sen=c$sensitivities, spe=c$specificities, 
                                      Category='Moderate dysplasia'))
roc_data = rbind(roc_data, data.frame(sen=d$sensitivities, spe=d$specificities, 
                                      Category='Severe dysplasia'))

AUC = data.frame(hyper=a$ci, mild=b$ci, moderate=c$ci, severe=d$ci)
rownames(AUC) = c('lower', 'area', 'upper')
AUC$mean = rowMeans(AUC)
AUC = round(AUC, 3)
txt = paste0('AUC of Hyperplasia: ', AUC$hyper[2], ' (', AUC$hyper[1], ' - ', AUC$hyper[3], ')\n',
             'AUC of Mild dysplasia: ', AUC$mild[2], ' (', AUC$mild[1], ' - ', AUC$mild[3], ')\n',
             'AUC of Moderate dysplasia: ', AUC$moderate[2], 
             ' (', AUC$moderate[1], ' - ', AUC$moderate[3], ')\n',
             'AUC of Severe dysplasia: ', AUC$severe[2], 
             ' (', AUC$severe[1], ' - ', AUC$severe[3], ')\n',
             'Average AUC: ', AUC$mean[2], ' (', AUC$mean[1], ' - ', AUC$mean[3], ')\n')

ggplot(data=roc_data, aes(x=1-spe, y=sen, color=Category)) + geom_path(size=1.5) + 
  geom_abline(linetype=2, size=1) + theme_bw() + xlab('1 - Specificity') + ylab('Sensitivity') +
  theme(panel.background=element_blank(), title=element_text(size=18), 
        axis.text=element_text(size=15), plot.title=element_text(hjust=0.5),
        legend.text=element_text(size=11, face='bold')) +
  ggtitle('ROC Curve of the Model') +
  annotate(geom='text', label=txt, x=0.8, y=0.2)

# multi
dat = dat0
summary(fit2)
prob = predict.glm(fit2, type=c('response'))
dat = cbind(dat, prob)
dat$OED1 = ifelse(dat$dys_bdkq == label[1], label[1], 'Others')
dat$OED2 = ifelse(dat$dys_bdkq == label[2], label[2], 'Others')
dat$OED3 = ifelse(dat$dys_bdkq == label[3], label[3], 'Others')
dat$OED4 = ifelse(dat$dys_bdkq == label[4], label[4], 'Others')
a = roc(response=dat$OED1, predictor=dat$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat$OED2, predictor=dat$`Mild dysplasia`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat$OED3, predictor=dat$`Moderate dysplasia`, levels=c(label[3], 'Others'), auc=T, ci=T)
d = roc(response=dat$OED4, predictor=dat$`Severe dysplasia`, levels=c(label[4], 'Others'), auc=T, ci=T)
roc_data = data.frame(sen=a$sensitivities, spe=a$specificities, Category='Hyperplasia')
roc_data = rbind(roc_data, data.frame(sen=b$sensitivities, spe=b$specificities, 
                                      Category='Mild dysplasia'))
roc_data = rbind(roc_data, data.frame(sen=c$sensitivities, spe=c$specificities, 
                                      Category='Moderate dysplasia'))
roc_data = rbind(roc_data, data.frame(sen=d$sensitivities, spe=d$specificities, 
                                      Category='Severe dysplasia'))

AUC = data.frame(hyper=a$ci, mild=b$ci, moderate=c$ci, severe=d$ci)
rownames(AUC) = c('lower', 'area', 'upper')
AUC$mean = rowMeans(AUC)
AUC = round(AUC, 3)
txt = paste0('AUC of Hyperplasia: ', AUC$hyper[2], ' (', AUC$hyper[1], ' - ', AUC$hyper[3], ')\n',
             'AUC of Mild dysplasia: ', AUC$mild[2], ' (', AUC$mild[1], ' - ', AUC$mild[3], ')\n',
             'AUC of Moderate dysplasia: ', AUC$moderate[2], 
             ' (', AUC$moderate[1], ' - ', AUC$moderate[3], ')\n',
             'AUC of Severe dysplasia: ', AUC$severe[2], 
             ' (', AUC$severe[1], ' - ', AUC$severe[3], ')\n',
             'Average AUC: ', AUC$mean[2], ' (', AUC$mean[1], ' - ', AUC$mean[3], ')\n')

ggplot(data=roc_data, aes(x=1-spe, y=sen, color=Category)) + geom_path(size=1.5) + 
  geom_abline(linetype=2, size=1) + theme_bw() + xlab('1 - Specificity') + ylab('Sensitivity') +
  theme(panel.background=element_blank(), title=element_text(size=18), 
        axis.text=element_text(size=15), plot.title=element_text(hjust=0.5),
        legend.text=element_text(size=11, face='bold')) +
  ggtitle('ROC Curve of the Model') +
  annotate(geom='text', label=txt, x=0.8, y=0.2)

#### try 1(23)4 ####
# uni
dat = dat0
dat$dys_bdkq = as.character(dat$dys_bdkq)
dat$dys_bdkq = ifelse(dat$dys_bdkq == 'Mild dysplasia' | dat$dys_bdkq == 'Moderate dysplasia', 
                      'Low grade', dat$dys_bdkq)
dat$dys_bdkq[dat$dys_bdkq == 'Severe dysplasia'] = 'High grade'
label = c('Hyperplasia', 'Low grade', 'High grade')
dat$dys_bdkq = factor(dat$dys_bdkq, levels=label)
fit3 = polr(dys_bdkq~score, Hess=T, data=dat)
summary(fit3)
prob = predict.glm(fit3, type=c('response'))
dat = cbind(dat, prob)
dat$OED1 = ifelse(dat$dys_bdkq == label[1], label[1], 'Others')
dat$OED2 = ifelse(dat$dys_bdkq == label[2], label[2], 'Others')
dat$OED3 = ifelse(dat$dys_bdkq == label[3], label[3], 'Others')
a = roc(response=dat$OED1, predictor=dat$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat$OED2, predictor=dat$`Low grade`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat$OED3, predictor=dat$`High grade`, levels=c(label[3], 'Others'), auc=T, ci=T)
roc_data = data.frame(sen=a$sensitivities, spe=a$specificities, Category=label[1])
roc_data = rbind(roc_data, data.frame(sen=b$sensitivities, spe=b$specificities, Category=label[2]))
roc_data = rbind(roc_data, data.frame(sen=c$sensitivities, spe=c$specificities, Category=label[3]))

AUC = data.frame(hyper=a$ci, low=b$ci, high=c$ci)
rownames(AUC) = c('lower', 'area', 'upper')
AUC$mean = rowMeans(AUC)
AUC = round(AUC, 3)
txt = paste0('AUC of Hyperplasia: ', AUC$hyper[2], ' (', AUC$hyper[1], ' - ', AUC$hyper[3], ')\n',
             'AUC of Low grade: ', AUC$low[2], ' (', AUC$low[1], ' - ', AUC$low[3], ')\n',
             'AUC of High grade: ', AUC$high[2], ' (', AUC$high[1], ' - ', AUC$high[3], ')\n',
             'Average AUC: ', AUC$mean[2], ' (', AUC$mean[1], ' - ', AUC$mean[3], ')\n')

ggplot(data=roc_data, aes(x=1-spe, y=sen, color=Category)) + geom_path(size=1.5) + 
  geom_abline(linetype=2, size=1) + theme_bw() + xlab('1 - Specificity') + ylab('Sensitivity') +
  theme(panel.background=element_blank(), title=element_text(size=18), 
        axis.text=element_text(size=15), plot.title=element_text(hjust=0.5),
        legend.text=element_text(size=11, face='bold')) +
  ggtitle('ROC Curve of the Model') +
  annotate(geom='text', label=txt, x=0.8, y=0.2)

# multi
dat = dat0
dat$dys_bdkq = as.character(dat$dys_bdkq)
dat$dys_bdkq = ifelse(dat$dys_bdkq == 'Mild dysplasia' | dat$dys_bdkq == 'Moderate dysplasia', 
                      'Low grade', dat$dys_bdkq)
dat$dys_bdkq[dat$dys_bdkq == 'Severe dysplasia'] = 'High grade'
label = c('Hyperplasia', 'Low grade', 'High grade')
dat$dys_bdkq = factor(dat$dys_bdkq, levels=label)
fit3 = polr(dys_bdkq~score+age+sex+smoke+drink, Hess=T, data=dat)
summary(fit3)
prob = predict.glm(fit3, type=c('response'))
dat = cbind(dat, prob)
dat$OED1 = ifelse(dat$dys_bdkq == label[1], label[1], 'Others')
dat$OED2 = ifelse(dat$dys_bdkq == label[2], label[2], 'Others')
dat$OED3 = ifelse(dat$dys_bdkq == label[3], label[3], 'Others')
a = roc(response=dat$OED1, predictor=dat$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat$OED2, predictor=dat$`Low grade`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat$OED3, predictor=dat$`High grade`, levels=c(label[3], 'Others'), auc=T, ci=T)
roc_data = data.frame(sen=a$sensitivities, spe=a$specificities, Category=label[1])
roc_data = rbind(roc_data, data.frame(sen=b$sensitivities, spe=b$specificities, Category=label[2]))
roc_data = rbind(roc_data, data.frame(sen=c$sensitivities, spe=c$specificities, Category=label[3]))

AUC = data.frame(hyper=a$ci, low=b$ci, high=c$ci)
rownames(AUC) = c('lower', 'area', 'upper')
AUC$mean = rowMeans(AUC)
AUC = round(AUC, 3)
txt = paste0('AUC of Hyperplasia: ', AUC$hyper[2], ' (', AUC$hyper[1], ' - ', AUC$hyper[3], ')\n',
             'AUC of Low grade: ', AUC$low[2], ' (', AUC$low[1], ' - ', AUC$low[3], ')\n',
             'AUC of High grade: ', AUC$high[2], ' (', AUC$high[1], ' - ', AUC$high[3], ')\n',
             'Average AUC: ', AUC$mean[2], ' (', AUC$mean[1], ' - ', AUC$mean[3], ')\n')

ggplot(data=roc_data, aes(x=1-spe, y=sen, color=Category)) + geom_path(size=1.5) + 
  geom_abline(linetype=2, size=1) + theme_bw() + xlab('1 - Specificity') + ylab('Sensitivity') +
  theme(panel.background=element_blank(), title=element_text(size=18), 
        axis.text=element_text(size=15), plot.title=element_text(hjust=0.5),
        legend.text=element_text(size=11, face='bold')) +
  ggtitle('ROC Curve of the Model') +
  annotate(geom='text', label=txt, x=0.8, y=0.2)

#### end ####
