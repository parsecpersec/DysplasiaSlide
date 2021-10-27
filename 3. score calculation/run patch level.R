setwd('')
library(ggplot2)
library(MASS)

files = list.files('./500/')
for(f in files) {
  temp = read.csv(paste0('./500/', f), F, stringsAsFactors=F)
  if(grepl('large', f)) {
    if(f == files[1]) {
      large = temp
      colnames(large)[1] = 'file'
    } else {
      large = cbind(large, temp$V2)
    }
    colnames(large)[ncol(large)] = gsub('[.]csv', '', f)
  } else {
    if(f == files[2]) {  # 10_small.csv
      small = temp
      colnames(small)[1] = 'file'
    } else {
      small = cbind(small, temp$V2)
    }
    colnames(small)[ncol(small)] = gsub('[.]csv', '', f)
  }
}
write.csv(large, './500/large.csv', row.names=F, quote=F)
write.csv(small, './500/small.csv', row.names=F, quote=F)

#### start here ####
dat1 = read.csv('../raw data/patch-C.csv', T, stringsAsFactors=F)
dat2 = read.csv('../raw data/patch-D.csv', T, stringsAsFactors=F)
dat = read.csv('../raw data/label.csv', T, stringsAsFactors=F)
dat = dat[dat$dysplasia != 'cancer',]

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
chi = colnames(dat1)[2:17]

rm(dat1, dat2)
dat_c = dat[substr(dat$number, 1, 1) == 'C',]
dat_d = dat[substr(dat$number, 1, 1) == 'D',]
file1 = gsub('X', '', colnames(large)[2:7])
file2 = gsub('X', '', colnames(small)[2:7])

large = large[match(dat_d$new, large$file),]
small = small[match(dat_c$new, small$file),]
large[,2:7] = 100 * large[,2:7]
small[,2:7] = 100 * small[,2:7]
dat_d = cbind(dat_d, large[,2:7])
dat_c = cbind(dat_c, small[,2:7])
label = unique(dat_c$dysplasia)
dat_c$dysplasia = factor(dat_c$dysplasia, levels=label)
dat_d$dysplasia = factor(dat_d$dysplasia, levels=label)

#### large OR ####
OR_large = data.frame()
for(i in 1:length(file1)) {
  fit = polr(dat_d$dysplasia ~ dat_d[,(21+i)], Hess=T)
  res = exp(cbind(OR_large=coef(fit), confint.lm(fit, level=0.95)[1,]))
  res = as.data.frame(res, stringsAsFactors=F)
  res$var = c(paste0('lower_', i), paste0('upper_', i))
  OR_large = rbind(OR_large, res)
}
OR_large$var = 1
OR_large$var[seq(1, 2*length(file1), 2)] = OR_large$V2[seq(2, 2*length(file1), 2)]
OR_large = OR_large[-seq(2, 2*length(file1), 2),]
rownames(OR_large) = 1:length(file1)
colnames(OR_large) = c('OR', 'lower', 'upper')
OR_large$include = ifelse((OR_large$lower-1)*(OR_large$upper-1)>0, T, F)
OR_large$var = file1

OR_large2 = plyr::arrange(OR_large, var, decreasing=F)
OR_large2$var = factor(OR_large2$var, levels=rev(OR_large2$var))
# OR_large2[,1:3] = log(OR_large2[,1:3])  # if OR_large too large

pd = position_dodge(0.1)
ggplot(OR_large2, aes(x=var, y=OR, color=include)) + coord_flip() +
  ggtitle('Large') + 
  geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  geom_line(position=pd, size = 0.75) + xlab('') + ylab('Odds Ratio') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) + theme_bw() +
  theme(panel.background=element_blank()) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        title = element_text(size = 18, hjust=0.5), legend.position='none',
        plot.title=element_text(hjust=0.5)) +
  scale_color_manual(values=c('tomato'))

#### small OR ####
OR_small = data.frame()
for(i in 1:length(file2)) {
  fit = polr(dat_c$dysplasia ~ dat_c[,(21+i)], Hess=T)
  res = exp(cbind(OR_small=coef(fit), confint.lm(fit, level=0.95)[1,]))
  res = as.data.frame(res, stringsAsFactors=F)
  res$var = c(paste0('lower_', i), paste0('upper_', i))
  OR_small = rbind(OR_small, res)
}
OR_small$var = 1
OR_small$var[seq(1, 2*length(file2), 2)] = OR_small$V2[seq(2, 2*length(file2), 2)]
OR_small = OR_small[-seq(2, 2*length(file2), 2),]
rownames(OR_small) = 1:length(file2)
colnames(OR_small) = c('OR', 'lower', 'upper')
OR_small$include = ifelse((OR_small$lower-1)*(OR_small$upper-1)>0, T, F)
OR_small$var = file2

OR_small2 = plyr::arrange(OR_small, var, decreasing=F)
OR_small2$var = factor(OR_small2$var, levels=rev(OR_small2$var))
# OR_small2[,1:3] = log(OR_small2[,1:3])  # if OR_small too large

pd = position_dodge(0.1)
ggplot(OR_small2, aes(x=var, y=OR, color=include)) + coord_flip() +
  ggtitle('Small') + 
  geom_point(position=pd, size=3, shape=15) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, position=pd, size=1.2) +
  geom_line(position=pd, size = 0.75) + xlab('') + ylab('Odds Ratio') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = "blue", size = 1) + theme_bw() +
  theme(panel.background=element_blank()) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        title = element_text(size = 18, hjust=0.5), legend.position='none',
        plot.title=element_text(hjust=0.5)) +
  scale_color_manual(values=c('grey', 'tomato'))

#### Score ####
OR = rbind(OR_large, OR_small)
rownames(OR) = OR$var

dat0 = read.csv('./saved.csv', T, stringsAsFactors=F)
dat0 = dat0[dat0$dys_bdkq != 'Cancer',]
filename = c('1_large', '2_large', '3_large', '4_large',
             '6_large', '8_small', '9_small', '10_small',
             '11_small', '12_small', '13_large', '16_small')
filename = paste0('X', filename)
series = 9:(8+length(filename))
soft1 = read.csv('./soft_patch/large.csv', T, stringsAsFactors=F)
soft2 = read.csv('./soft_patch/small.csv', T, stringsAsFactors=F)

sum = rep(0, nrow(soft1))
n = 0
for(f in file1) {
  t = log(OR[f, 1]) * soft1[, paste0('X', f)]
  sum = sum + t
  n = n + 1
}
soft1$score_large = 100 * sum / n

sum = rep(0, nrow(soft2))
n = 0
for(f in file2) {
  t = log(OR[f, 1]) * soft2[, paste0('X', f)]
  sum = sum + t
  n = n + 1
}
soft2$score_small = 100 * sum / n

# mean(OR weighted softmax)
dat0$score_large = aggregate(soft1$score_large, FUN=mean, by=list(path=soft1$path))$x
dat0$score_small = aggregate(soft2$score_small, FUN=mean, by=list(path=soft2$path))$x

dat0$dys_bdkq = factor(dat0$dys_bdkq, levels=c("Hyperplasia",
                                               "Mild dysplasia",
                                               "Moderate dysplasia",
                                               "Severe dysplasia"))
library(rms)
fit = orm(dys_bdkq~score_large, data=dat0)
summary(fit)
res = exp(cbind(OR=coef(fit), confint(fit)))
res = as.data.frame(res, stringsAsFactors=F)
colnames(res) = c('OR', 'lower', 'upper')

fit2 = orm(dys_bdkq~score_small, data=dat0)
summary(fit2)
res2 = exp(cbind(OR=coef(fit2), confint(fit2)))
res2 = as.data.frame(res2, stringsAsFactors=F)
colnames(res2) = c('OR', 'lower', 'upper')

fit3 = orm(dys_bdkq~score_large+score_small, data=dat0)
summary(fit3)
res3 = exp(cbind(OR=coef(fit3), confint(fit3)))
res3 = as.data.frame(res3, stringsAsFactors=F)
colnames(res3) = c('OR', 'lower', 'upper')

fit4 = orm(dys_bdkq~score_large+score_small+age+sex+smoke+drink, data=dat0)
summary(fit4)
res4 = exp(cbind(OR=coef(fit4), confint(fit4)))
res4 = as.data.frame(res4, stringsAsFactors=F)
colnames(res4) = c('OR', 'lower', 'upper')

#### ROC 4 - uni ####
dat00 = dat0  # backup
fit3 = polr(dys_bdkq~score_large+score_small, data=dat0, Hess=T)
prob = predict.glm(fit3, type=c('response'))
dat0 = cbind(dat0, prob)
label = levels(dat0$dys_bdkq)
dat0$OED1 = ifelse(dat0$dys_bdkq == label[1], label[1], 'Others')
dat0$OED2 = ifelse(dat0$dys_bdkq == label[2], label[2], 'Others')
dat0$OED3 = ifelse(dat0$dys_bdkq == label[3], label[3], 'Others')
dat0$OED4 = ifelse(dat0$dys_bdkq == label[4], label[4], 'Others')
a = roc(response=dat0$OED1, predictor=dat0$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat0$OED2, predictor=dat0$`Mild dysplasia`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat0$OED3, predictor=dat0$`Moderate dysplasia`, levels=c(label[3], 'Others'), auc=T, ci=T)
d = roc(response=dat0$OED4, predictor=dat0$`Severe dysplasia`, levels=c(label[4], 'Others'), auc=T, ci=T)
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

#### ROC 4 multi ####
dat0 = dat00
fit3 = polr(dys_bdkq~score_large+score_small+age+sex+smoke+drink, data=dat0, Hess=T)
prob = predict.glm(fit3, type=c('response'))
dat0 = cbind(dat0, prob)
label = levels(dat0$dys_bdkq)
dat0$OED1 = ifelse(dat0$dys_bdkq == label[1], label[1], 'Others')
dat0$OED2 = ifelse(dat0$dys_bdkq == label[2], label[2], 'Others')
dat0$OED3 = ifelse(dat0$dys_bdkq == label[3], label[3], 'Others')
dat0$OED4 = ifelse(dat0$dys_bdkq == label[4], label[4], 'Others')
a = roc(response=dat0$OED1, predictor=dat0$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat0$OED2, predictor=dat0$`Mild dysplasia`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat0$OED3, predictor=dat0$`Moderate dysplasia`, levels=c(label[3], 'Others'), auc=T, ci=T)
d = roc(response=dat0$OED4, predictor=dat0$`Severe dysplasia`, levels=c(label[4], 'Others'), auc=T, ci=T)
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

#### ROC 3 - uni ####
dat0 = dat00
dat0$dys_bdkq = as.character(dat0$dys_bdkq)
dat0$dys_bdkq = ifelse(dat0$dys_bdkq == 'Mild dysplasia' | dat0$dys_bdkq == 'Moderate dysplasia', 
                       'Low grade', dat0$dys_bdkq)
dat0$dys_bdkq[dat0$dys_bdkq == 'Severe dysplasia'] = 'High grade'
label = c('Hyperplasia', 'Low grade', 'High grade')
dat0$dys_bdkq = factor(dat0$dys_bdkq, levels=label)
fit3 = polr(dys_bdkq~score_large+score_small, data=dat0, Hess=T)
prob = predict.glm(fit3, type=c('response'))
dat0 = cbind(dat0, prob)
dat0$OED1 = ifelse(dat0$dys_bdkq == label[1], label[1], 'Others')
dat0$OED2 = ifelse(dat0$dys_bdkq == label[2], label[2], 'Others')
dat0$OED3 = ifelse(dat0$dys_bdkq == label[3], label[3], 'Others')
a = roc(response=dat0$OED1, predictor=dat0$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat0$OED2, predictor=dat0$`Low grade`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat0$OED3, predictor=dat0$`High grade`, levels=c(label[3], 'Others'), auc=T, ci=T)
roc_data = data.frame(sen=a$sensitivities, spe=a$specificities, Category='Hyperplasia')
roc_data = rbind(roc_data, data.frame(sen=b$sensitivities, spe=b$specificities, 
                                      Category='Low grade'))
roc_data = rbind(roc_data, data.frame(sen=c$sensitivities, spe=c$specificities, 
                                      Category='High grade'))

AUC = data.frame(hyper=a$ci, low=b$ci, high=c$ci)
rownames(AUC) = c('lower', 'area', 'upper')
AUC$mean = rowMeans(AUC)
AUC = round(AUC, 3)
txt = paste0('AUC of Hyperplasia: ', AUC$hyper[2], ' (', AUC$hyper[1], ' - ', AUC$hyper[3], ')\n',
             'AUC of Low Grade: ', AUC$low[2], ' (', AUC$low[1], ' - ', AUC$low[3], ')\n',
             'AUC of High Grade: ', AUC$high[2], ' (', AUC$high[1], ' - ', AUC$high[3], ')\n',
             'Average AUC: ', AUC$mean[2], ' (', AUC$mean[1], ' - ', AUC$mean[3], ')\n')

ggplot(data=roc_data, aes(x=1-spe, y=sen, color=Category)) + geom_path(size=1.5) + 
  geom_abline(linetype=2, size=1) + theme_bw() + xlab('1 - Specificity') + ylab('Sensitivity') +
  theme(panel.background=element_blank(), title=element_text(size=18), 
        axis.text=element_text(size=15), plot.title=element_text(hjust=0.5),
        legend.text=element_text(size=11, face='bold')) +
  ggtitle('ROC Curve of the Model') +
  annotate(geom='text', label=txt, x=0.8, y=0.2)

#### ROC 3 multi ####
dat0 = dat00
dat0$dys_bdkq = as.character(dat0$dys_bdkq)
dat0$dys_bdkq = ifelse(dat0$dys_bdkq == 'Mild dysplasia' | dat0$dys_bdkq == 'Moderate dysplasia', 
                       'Low grade', dat0$dys_bdkq)
dat0$dys_bdkq[dat0$dys_bdkq == 'Severe dysplasia'] = 'High grade'
label = c('Hyperplasia', 'Low grade', 'High grade')
dat0$dys_bdkq = factor(dat0$dys_bdkq, levels=label)
fit3 = polr(dys_bdkq~score_large+score_small+age+sex+smoke+drink, data=dat0, Hess=T)
prob = predict.glm(fit3, type=c('response'))
dat0 = cbind(dat0, prob)
dat0$OED1 = ifelse(dat0$dys_bdkq == label[1], label[1], 'Others')
dat0$OED2 = ifelse(dat0$dys_bdkq == label[2], label[2], 'Others')
dat0$OED3 = ifelse(dat0$dys_bdkq == label[3], label[3], 'Others')
a = roc(response=dat0$OED1, predictor=dat0$Hyperplasia, levels=c(label[1], 'Others'), auc=T, ci=T)
b = roc(response=dat0$OED2, predictor=dat0$`Low grade`, levels=c(label[2], 'Others'), auc=T, ci=T)
c = roc(response=dat0$OED3, predictor=dat0$`High grade`, levels=c(label[3], 'Others'), auc=T, ci=T)
roc_data = data.frame(sen=a$sensitivities, spe=a$specificities, Category='Hyperplasia')
roc_data = rbind(roc_data, data.frame(sen=b$sensitivities, spe=b$specificities, 
                                      Category='Low grade'))
roc_data = rbind(roc_data, data.frame(sen=c$sensitivities, spe=c$specificities, 
                                      Category='High grade'))

AUC = data.frame(hyper=a$ci, low=b$ci, high=c$ci)
rownames(AUC) = c('lower', 'area', 'upper')
AUC$mean = rowMeans(AUC)
AUC = round(AUC, 3)
txt = paste0('AUC of Hyperplasia: ', AUC$hyper[2], ' (', AUC$hyper[1], ' - ', AUC$hyper[3], ')\n',
             'AUC of Low Grade: ', AUC$low[2], ' (', AUC$low[1], ' - ', AUC$low[3], ')\n',
             'AUC of High Grade: ', AUC$high[2], ' (', AUC$high[1], ' - ', AUC$high[3], ')\n',
             'Average AUC: ', AUC$mean[2], ' (', AUC$mean[1], ' - ', AUC$mean[3], ')\n')

ggplot(data=roc_data, aes(x=1-spe, y=sen, color=Category)) + geom_path(size=1.5) + 
  geom_abline(linetype=2, size=1) + theme_bw() + xlab('1 - Specificity') + ylab('Sensitivity') +
  theme(panel.background=element_blank(), title=element_text(size=18), 
        axis.text=element_text(size=15), plot.title=element_text(hjust=0.5),
        legend.text=element_text(size=11, face='bold')) +
  ggtitle('ROC Curve of the Model') +
  annotate(geom='text', label=txt, x=0.8, y=0.2)

#### density ####
ggplot(dat00, aes(x=score_large, fill=dys_bdkq, color=dys_bdkq)) + geom_density(size=1, alpha=0.1) + 
  theme_bw() + xlim(c(0 * min(dat00$score_large), 1.5 * max(dat00$score_large))) +
  theme(panel.background=element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle('Large') + xlab('Score') + ylab('Density')

ggplot(dat00, aes(x=score_small, fill=dys_bdkq, color=dys_bdkq)) + geom_density(size=1, alpha=0.1) + 
  theme_bw() + xlim(c(0.95 * min(dat00$score_small), 1.05 * max(dat00$score_small))) +
  theme(panel.background=element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle('Small') + xlab('Score') + ylab('Density')

#### boxplot ####
label = levels(factor(dat00$dys_bdkq))
comps = combn(label, 2)
mycomp = list(comps[,1], comps[,2], comps[,3], comps[,4], comps[,5], comps[,6])

ggplot(dat00, aes(x=dys_bdkq, y=score_large, fill=dys_bdkq)) + geom_boxplot() + 
  stat_compare_means(method='anova', label.y=5) +
  stat_compare_means(method='kruskal', label.x=2, label.y=5) + 
  stat_compare_means(comparisons=mycomp, method='wilcox') + theme_bw() +
  theme(panel.background=element_blank(), legend.position='none',
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle('Large') + xlab('') + ylab('Score')

ggplot(dat00, aes(x=dys_bdkq, y=score_small, fill=dys_bdkq)) + geom_boxplot() + 
  stat_compare_means(method='anova', label.y=11) +
  stat_compare_means(method='kruskal', label.x=2, label.y=11) + 
  stat_compare_means(comparisons=mycomp, method='wilcox') + theme_bw() +
  theme(panel.background=element_blank(), legend.position='none',
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title=element_text(hjust=0.5), title = element_text(size = 18, hjust=0.5)) +
  ggtitle('Small') + xlab('') + ylab('Score')
