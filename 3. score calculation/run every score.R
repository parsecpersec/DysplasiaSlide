setwd('')
library(MASS)
library(ggplot2)
library(ggpubr)
library(pROC)

#### start ####
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
filename = paste0('X', filename)
series = 9:(8+length(filename))

soft1 = read.csv('./soft_patch/large.csv', T, stringsAsFactors=F)
soft2 = read.csv('./soft_patch/small.csv', T, stringsAsFactors=F)
dat = read.csv('./saved.csv', T, stringsAsFactors=F)
dat = dat[dat$dys_bdkq != 'Cancer', ]
dat$dys_bdkq = factor(dat$dys_bdkq, levels=c("Hyperplasia",
                                             "Mild dysplasia",
                                             "Moderate dysplasia",
                                             "Severe dysplasia"))
label = levels(dat$dys_bdkq)

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

fit = polr(data=dat, formula=dys_bdkq ~ X1_large+X2_large+X3_large+X4_large+
             X6_large+X8_small+X9_small+X10_small+X11_small+X12_small+
             X13_large+X16_small+age+sex+smoke+drink)
summary(fit)
res = exp(cbind(OR=coef(fit), confint.lm(fit, level=0.95)[1:(nrow(confint.lm(fit, level=0.95))-3),]))
res = as.data.frame(res, stringsAsFactors=F)
res$var = rownames(res)

# back up
dat0 = dat
# 4
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

# 3
dat = dat0
dat$dys_bdkq = as.character(dat$dys_bdkq)
dat$dys_bdkq = ifelse(dat$dys_bdkq == 'Mild dysplasia' | dat$dys_bdkq == 'Moderate dysplasia', 
                      'Low grade', dat$dys_bdkq)
dat$dys_bdkq[dat$dys_bdkq == 'Severe dysplasia'] = 'High grade'
label = c('Hyperplasia', 'Low grade', 'High grade')
dat$dys_bdkq = factor(dat$dys_bdkq, levels=label)
fit = polr(data=dat, formula=dys_bdkq ~ X1_large+X2_large+X3_large+X4_large+
             X6_large+X8_small+X9_small+X10_small+X11_small+X12_small+
             X13_large+X16_small+age+sex+smoke+drink)
summary(fit)
prob = predict.glm(fit, type=c('response'))
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

