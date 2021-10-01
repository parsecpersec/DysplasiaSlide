setwd('')

dat1 = read.csv('../raw data/patch-C.csv', T, stringsAsFactors=F)
dat2 = read.csv('../raw data/patch-D.csv', T, stringsAsFactors=F)
dat = read.csv('../raw data/label.csv', T, stringsAsFactors=F)
dat = dat[dat$dysplasia != 'cancer',]

# n_feature = 13

# dat$X5 = NULL
# dat$X7 = NULL
# dat$X14 = NULL

# C-224 please remove X3

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
# eng13 = eng[c(1:4, 6, 8:13, 15, 16)]
# chi13 = chi[c(1:4, 6, 8:13, 15, 16)]
rm(dat1, dat2)
dat_c = dat[substr(dat$number, 1, 1) == 'C',]
dat_d = dat[substr(dat$number, 1, 1) == 'D',]

table_OR = function(tab, res_df, size) {
  a = tab[1,1]
  b = tab[2,1]
  c = tab[1,2]
  d = tab[2,2]
  n = a+b+c+d
  # correction?
  if(T) {
    a = a + 0.5
    b = b + 0.5
    c = c + 0.5
    d = d + 0.5
  }
  OR = a*d/b/c
  modify = exp(1.96*sqrt(1/a+1/b+1/c+1/d))
  output = data.frame(feature=paste0(chi[i], ' - ', size),
                      number=n, positive=a+c-1, negative=b+d-1, 
                      odds_ratio=OR, lower=OR/modify, upper=OR*modify)
  res_df = rbind(res_df, output)
  return(res_df)
}

#### OR-1 dysplasia vs nondysplasia ####
dat_c1 = dat_c
dat_d1 = dat_d
dat_c1$dysplasia = ifelse(dat_c1$dysplasia == 'hyperplasia', '0_nondysplasia', '1_dysplasia')
dat_d1$dysplasia = ifelse(dat_d1$dysplasia == 'hyperplasia', '0_nondysplasia', '1_dysplasia')
OR_c1 = data.frame()
for (i in 1:16) {
  temp = dat_c1[,c((i+4), 21)]
  temp = temp[temp[,1] %in% c('YY', 'NN'),]
  temp[,1] = factor(temp[,1], levels=c('YY', 'NN'))
  temp[,2] = factor(temp[,2], levels=c('1_dysplasia', '0_nondysplasia'))
  t_tab = table(temp[,1], temp[,2])
  OR_c1 = table_OR(tab=t_tab, res_df=OR_c1, size=224)
}
OR_d1 = data.frame()
for (i in 1:16) {
  temp = dat_d1[,c((i+4), 21)]
  temp = temp[temp[,1] %in% c('YY', 'NN'),]
  temp[,1] = factor(temp[,1], levels=c('YY', 'NN'))
  temp[,2] = factor(temp[,2], levels=c('1_dysplasia', '0_nondysplasia'))
  t_tab = table(temp[,1], temp[,2])
  OR_d1 = table_OR(tab=t_tab, res_df=OR_d1, size=1024)
}
OR1 = rbind(OR_c1, OR_d1)

#### OR-2 high grade vs low grade ####
dat_c2 = dat_c
dat_d2 = dat_d
dat_c2$dysplasia = ifelse(dat_c2$dysplasia == 'hyperplasia' | dat_c2$dysplasia == 'mild dysplasia', 
                          '0_low_grade', '1_high_grade')
dat_d2$dysplasia = ifelse(dat_d2$dysplasia == 'hyperplasia' | dat_d2$dysplasia == 'mild dysplasia', 
                          '0_low_grade', '1_high_grade')
OR_c2 = data.frame()
for (i in 1:16) {
  temp = dat_c2[,c((i+4), 21)]
  temp = temp[temp[,1] %in% c('YY', 'NN'),]
  temp[,1] = factor(temp[,1], levels=c('YY', 'NN'))
  temp[,2] = factor(temp[,2], levels=c('1_high_grade', '0_low_grade'))
  t_tab = table(temp[,1], temp[,2])
  OR_c2 = table_OR(tab=t_tab, res_df=OR_c2, size=224)
}
OR_d2 = data.frame()
for (i in 1:16) {
  temp = dat_d2[,c((i+4), 21)]
  temp = temp[temp[,1] %in% c('YY', 'NN'),]
  temp[,1] = factor(temp[,1], levels=c('YY', 'NN'))
  temp[,2] = factor(temp[,2], levels=c('1_high_grade', '0_low_grade'))
  t_tab = table(temp[,1], temp[,2])
  OR_d2 = table_OR(tab=t_tab, res_df=OR_d2, size=1024)
}
OR2 = rbind(OR_c2, OR_d2)

#### OR-3 severe (high risk) vs others ####
dat_c3 = dat_c
dat_d3 = dat_d
dat_c3$dysplasia = ifelse(dat_c3$dysplasia == 'severe dysplasia', '1_high_risk', '0_low_risk')
dat_d3$dysplasia = ifelse(dat_d3$dysplasia == 'severe dysplasia', '1_high_risk', '0_low_risk')
OR_c3 = data.frame()
for (i in 1:16) {
  temp = dat_c3[,c((i+4), 21)]
  temp = temp[temp[,1] %in% c('YY', 'NN'),]
  temp[,1] = factor(temp[,1], levels=c('YY', 'NN'))
  temp[,2] = factor(temp[,2], levels=c('1_high_risk', '0_low_risk'))
  t_tab = table(temp[,1], temp[,2])
  OR_c3 = table_OR(tab=t_tab, res_df=OR_c3, size=224)
}
OR_d3 = data.frame()
for (i in 1:16) {
  temp = dat_d3[,c((i+4), 21)]
  temp = temp[temp[,1] %in% c('YY', 'NN'),]
  temp[,1] = factor(temp[,1], levels=c('YY', 'NN'))
  temp[,2] = factor(temp[,2], levels=c('1_high_risk', '0_low_risk'))
  t_tab = table(temp[,1], temp[,2])
  OR_d3 = table_OR(tab=t_tab, res_df=OR_d3, size=1024)
}
OR3 = rbind(OR_c3, OR_d3)

#### save OR ####
title = c(1,1,1,2,2,2,3,3,3)
ORs = cbind(OR1, cbind(OR2[,5:7], OR3[,5:7]))
colnames(ORs)[5:13] = paste0(colnames(ORs)[5:13], '_', title)
write.csv(ORs, './OR.csv', quote=F, row.names=F)

#### OR analysis ####
setwd('')
library(ggplot2)
ORs = read.csv('./OR.csv', T, stringsAsFactors=F)
ORs$size = c(rep('small', 16), rep('large', 16))
ORs$feature = rep(paste0(sprintf('%02d', 1:16), ' - ', eng), 2)
ORs$feature = paste0(ORs$feature, c(rep(' - 224', 16), rep(' - 1024', 16)))
ORs$include_1 = ifelse((ORs$lower_1-1)*(ORs$upper_1-1)>0, T, F)  # sig
ORs$include_2 = ifelse((ORs$lower_2-1)*(ORs$upper_2-1)>0, T, F)
ORs$include_3 = ifelse((ORs$lower_3-1)*(ORs$upper_3-1)>0, T, F)
for(i in 5:13) {
  ORs[,i] = log(ORs[,i])  # upper too large, use log scale
}
# ORs$feature = factor(ORs$feature, levels=paste0(1:16, ' - ', eng))
ORs$size = factor(ORs$size, levels=c('small', 'large'))
ORs = ORs[order(ORs$feature),,drop=F]
ORs$feature = factor(ORs$feature, levels=rev(ORs$feature))
ORs$num = substr(ORs$feature, 1, 2)

plot1 = ORs[,c(1,5:7,14,15,18)]
plot2 = ORs[,c(1,8:10,14,16,18)]
plot3 = ORs[,c(1,11:13,14,17,18)]

pd = position_dodge(0.1)

ggplot(plot1, aes(x=feature, y=odds_ratio_1, color=size, size=include_1)) + 
  geom_point(shape=15, size=4) + coord_flip() +
  ggtitle('Dysplasia over Nondysplasia') +
  geom_errorbar(aes(ymin=lower_1, ymax=upper_1, linetype=include_1), position=pd, width=0) +
  theme_bw() + theme(panel.background=element_blank(), text=element_text(face='bold',size=14),
                     plot.title=element_text(hjust=0.5, size=20)) +
  ylab('Odds Ratio (log scale)') + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) +
  scale_size_manual(values=c(0.8, 2.5)) + scale_linetype_manual(values=c(6,1)) +
  geom_text(aes(label=num, y=lower_1-0.5), position=position_dodge(0.1), color='black', size=5)

ggplot(plot2, aes(x=feature, y=odds_ratio_2, color=size, size=include_2)) + 
  geom_point(shape=15, size=4) + coord_flip() +
  ggtitle('High Grade over Low Grade') +
  geom_errorbar(aes(ymin=lower_2, ymax=upper_2, linetype=include_2), position=pd, width=0) +
  theme_bw() + theme(panel.background=element_blank(), text=element_text(face='bold',size=14),
                     plot.title=element_text(hjust=0.5, size=20)) +
  ylab('Odds Ratio (log scale)') + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) +
  scale_size_manual(values=c(0.8, 2.5)) + scale_linetype_manual(values=c(6,1)) +
  geom_text(aes(label=num, y=lower_2-0.5), position=position_dodge(0.1), color='black', size=5)

ggplot(plot3, aes(x=feature, y=odds_ratio_3, color=size, size=include_3)) + 
  geom_point(shape=15, size=4) + coord_flip() +
  ggtitle('High Risk over Others') +
  geom_errorbar(aes(ymin=lower_3, ymax=upper_3, linetype=include_3), position=pd, width=0) +
  theme_bw() + theme(panel.background=element_blank(), text=element_text(face='bold',size=14),
                     plot.title=element_text(hjust=0.5, size=20)) +
  ylab('Odds Ratio (log scale)') + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) +
  scale_size_manual(values=c(0.8, 2.5)) + scale_linetype_manual(values=c(6,1)) +
  geom_text(aes(label=num, y=lower_3-0.5), position=position_dodge(0.1), color='black', size=5)

#### final selection ####
ORs$include = ORs$include_1 & ORs$include_2 & ORs$include_3
write.csv(ORs, './log OR - sig.csv', quote=F, row.names=F)

ORs = read.csv('./log OR - sig.csv', T, stringsAsFactors=F)
ORs[,grepl('include_', colnames(ORs))] = NULL
ORs$changes = c(rep('Architectural', 16), rep('Cytological', 16))
ORs$freq = T
ORs$freq[ORs$num == 5 | ORs$num == 7 | ORs$num == 14] = F
ORs$freq[ORs$size == 'small' & ORs$num == 3] = F
ORs$freq[ORs$size == 'small' & ORs$num == 15] = F
ORs$sig = ORs$include
ORs$include = NULL
# large: 5 #
# small: 7, 14 #
ORs$size_compare = F
ORs$size_compare[ORs$num == 6 | ORs$num == 11 | ORs$num == 12] = T
ORs$size_compare[ORs$size == 'large' & ORs$num %in% c(1,2,3,4,5,9,10,13,16)] = T
ORs$size_compare[ORs$size == 'small' & ORs$num %in% c(7,8,14,15)] = T

ORs$final_selection = ORs$freq & ORs$size_compare
# key sentence #
ORs$final_selection[ORs$freq & ORs$sig] = T

ORs$to_train = 'not trained'
ORs$to_train[ORs$size == 'large' & ORs$num %in% c(1,2,3,4,6,8,9,10,13,15,16)] = 'trained'
ORs$to_train[ORs$size == 'small' & ORs$num %in% c(1,10,11,12,13)] = 'trained'
ORs$need = ORs$final_selection & ORs$to_train == 'not trained'
ORs$test_accuracy = NA
ORs$test_AUC =NA

write.csv(ORs, './log OR - models.csv', quote=F, row.names=F)

library(ggplot2)

ORs$feature = factor(ORs$feature, levels=rev(ORs$feature))
ORs$size = factor(ORs$size, levels=c('small', 'large'))
ORs$num = substr(ORs$feature, 1, 2)

plot1 = ORs[,c(1,5:7,14:16,20)]
plot2 = ORs[,c(1,8:10,14:16,20)]
plot3 = ORs[,c(1,11:13,14:16,20)]

pd = position_dodge(0.1)

ggplot(plot1, aes(x=feature, y=odds_ratio_1, color=size, size=final_selection)) + 
  geom_point(shape=15, size=4) + coord_flip() +
  ggtitle('Dysplasia over Nondysplasia') +
  geom_errorbar(aes(ymin=lower_1, ymax=upper_1, linetype=final_selection), position=pd, width=0) +
  theme_bw() + theme(panel.background=element_blank(), text=element_text(face='bold',size=14),
                     plot.title=element_text(hjust=0.5, size=20)) +
  ylab('Odds Ratio (log scale)') + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) +
  scale_size_manual(values=c(0.8, 2.5)) + scale_linetype_manual(values=c(6,1)) +
  geom_text(aes(label=num, y=lower_1-0.5), position=position_dodge(0.1), color='black', size=5)

ggplot(plot2, aes(x=feature, y=odds_ratio_2, color=size, size=final_selection)) + 
  geom_point(shape=15, size=4) + coord_flip() +
  ggtitle('High Grade over Low Grade') +
  geom_errorbar(aes(ymin=lower_2, ymax=upper_2, linetype=final_selection), position=pd, width=0) +
  theme_bw() + theme(panel.background=element_blank(), text=element_text(face='bold',size=14),
                     plot.title=element_text(hjust=0.5, size=20)) +
  ylab('Odds Ratio (log scale)') + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) +
  scale_size_manual(values=c(0.8, 2.5)) + scale_linetype_manual(values=c(6,1)) +
  geom_text(aes(label=num, y=lower_2-0.5), position=position_dodge(0.1), color='black', size=5)

ggplot(plot3, aes(x=feature, y=odds_ratio_3, color=size, size=final_selection)) + 
  geom_point(shape=15, size=4) + coord_flip() +
  ggtitle('High Risk over Others') +
  geom_errorbar(aes(ymin=lower_3, ymax=upper_3, linetype=final_selection), position=pd, width=0) +
  theme_bw() + theme(panel.background=element_blank(), text=element_text(face='bold',size=14),
                     plot.title=element_text(hjust=0.5, size=20)) +
  ylab('Odds Ratio (log scale)') + xlab('') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "blue", size = 1) +
  scale_size_manual(values=c(0.8, 2.5)) + scale_linetype_manual(values=c(6,1)) +
  geom_text(aes(label=num, y=lower_3-0.5), position=position_dodge(0.1), color='black', size=5)
