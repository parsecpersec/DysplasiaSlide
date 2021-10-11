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
colnames(dat)[series] = gsub('./soft/Soft_', '', files)
for(i in series) {
  dat[,i] = 10*dat[,i]    # softmax (0-1) -> (0-10)
}
dat = dat[dat$dys_bdkq != 'Cancer',]
dat$dys_bdkq = factor(dat$dys_bdkq, levels=c("Hyperplasia",
                                             "Mild dysplasia",
                                             "Moderate dysplasia",
                                             "Severe dysplasia"))
#### view softmax ####
library(ggplot2)
ggplot()

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
OR2 = plyr::arrange(OR, include, OR, decreasing=T)
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
  
