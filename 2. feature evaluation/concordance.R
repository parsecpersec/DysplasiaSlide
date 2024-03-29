setwd('D:/work/OLK_B/chi/')
library(ggplot2)
library(irr)
library(extrafont)

#### record ####
for (i in 2:17) {
  condition = (dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes') |
    (dat1[,i] == 'No' & dat2[,i] == 'No' & dat3[,i] == 'No')
  dat[,i][condition] = 'Yes'
  dat[,i][!condition] = 'No'
}

kc = data.frame()
for (i in 2:17) {
  kc = rbind(kc, 
             data.frame(feature=i-1,
                        k1=kappa2(data.frame(P1=dat1[,i], P2=dat2[,i]))$value,
                        k2=kappa2(data.frame(P1=dat1[,i], P2=dat3[,i]))$value,
                        k3=kappa2(data.frame(P1=dat2[,i], P2=dat3[,i]))$value))
}
kc$k = rowMeans(kc[,2:4])
write.csv(kc, 'kappa_C.csv', quote=F, row.names=F)

for (i in 2:17) {
  condition = (dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes') |
    (dat1[,i] == 'No' & dat2[,i] == 'No' & dat3[,i] == 'No')
  dat[,i][condition] = 'Yes'
  dat[,i][!condition] = 'No'
}

kd = data.frame()
for (i in 2:17) {
  kd = rbind(kd, 
             data.frame(feature=i-1,
                        k1=kappa2(data.frame(P1=dat1[,i], P2=dat2[,i]))$value,
                        k2=kappa2(data.frame(P1=dat1[,i], P2=dat3[,i]))$value,
                        k3=kappa2(data.frame(P1=dat2[,i], P2=dat3[,i]))$value))
}
kd$k = rowMeans(kd[,2:4])
write.csv(kd, 'kappa_D.csv', quote=F, row.names=F)

#### C+D 3 old ####
for(i in 1:16) {
  tiff(filename=paste0('./plot/C_positive_', i, '.tiff'), res=300, height=1000, width=1000)
  print(ggvenn(data=list(A=dat1[,1][dat1[,(i+1)] == 'Yes'],
                         B=dat2[,1][dat2[,(i+1)] == 'Yes'],
                         C=dat3[,1][dat3[,(i+1)] == 'Yes']),
               show_percentage=F, fill_color=c('dodgerblue', 'indianred1', 'purple')))
  dev.off()
}
for(i in 1:16) {
  tiff(filename=paste0('./plot/C_negative_', i, '.tiff'), res=300, height=1000, width=1000)
  print(ggvenn(data=list(A=dat1[,1][dat1[,(i+1)] == 'No'],
                         B=dat2[,1][dat2[,(i+1)] == 'No'],
                         C=dat3[,1][dat3[,(i+1)] == 'No']),
               show_percentage=F, fill_color=c('dodgerblue', 'indianred1', 'purple')))
  dev.off()
}

for(i in 1:16) {
  tiff(filename=paste0('./plot/D_positive_', i, '.tiff'), res=300, height=1000, width=1000)
  print(ggvenn(data=list(A=dat1[,1][dat1[,(i+1)] == 'Yes'],
                         B=dat2[,1][dat2[,(i+1)] == 'Yes'],
                         C=dat3[,1][dat3[,(i+1)] == 'Yes']),
               show_percentage=F, fill_color=c('dodgerblue', 'indianred1', 'purple')))
  dev.off()
}
for(i in 1:16) {
  tiff(filename=paste0('./plot/D_negative_', i, '.tiff'), res=300, height=1000, width=1000)
  print(ggvenn(data=list(A=dat1[,1][dat1[,(i+1)] == 'No'],
                         B=dat2[,1][dat2[,(i+1)] == 'No'],
                         C=dat3[,1][dat3[,(i+1)] == 'No']),
               show_percentage=F, fill_color=c('dodgerblue', 'indianred1', 'purple')))
  dev.off()
}

#### C+D 3 ####
dat1 = readxl::read_xlsx('sheet.xlsx', sheet='C1')
dat2 = readxl::read_xlsx('sheet.xlsx', sheet='C2')
dat3 = readxl::read_xlsx('sheet.xlsx', sheet='C3')
dat = dat1
for (i in 2:17) {
  condition = (dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes') |
    (dat1[,i] == 'No' & dat2[,i] == 'No' & dat3[,i] == 'No')
  dat[,i][condition] = 'Yes'
  dat[,i][!condition] = 'No'
}
datp1 = data.frame()
for(i in 2:17) {
  datp1 = rbind(datp1, data.frame(feature=gsub('[0-9]*', '', colnames(dat)[i]),
                                  num=length(dat[,i][dat[,i] == 'Yes']),
                                  common='Yes',
                                  size='C - 224'))
}
for(i in 2:17) {
  datp1 = rbind(datp1, data.frame(feature=gsub('[0-9]*', '', colnames(dat)[i]),
                                  num=length(dat[,i][dat[,i] == 'No']),
                                  common='No',
                                  size='C - 224'))
}

dat1 = readxl::read_xlsx('sheet.xlsx', sheet='D1')
dat2 = readxl::read_xlsx('sheet.xlsx', sheet='D2')
dat3 = readxl::read_xlsx('sheet.xlsx', sheet='D3')
dat = dat1
for (i in 2:17) {
  condition = (dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes') |
    (dat1[,i] == 'No' & dat2[,i] == 'No' & dat3[,i] == 'No')
  dat[,i][condition] = 'Yes'
  dat[,i][!condition] = 'No'
}
datp2 = data.frame()
for(i in 2:17) {
  datp2 = rbind(datp2, data.frame(feature=gsub('[0-9]*', '', colnames(dat)[i]),
                                  num=length(dat[,i][dat[,i] == 'Yes']),
                                  common='Yes',
                                  size='D - 1024'))
}
for(i in 2:17) {
  datp2 = rbind(datp2, data.frame(feature=gsub('[0-9]*', '', colnames(dat)[i]),
                                  num=length(dat[,i][dat[,i] == 'No']),
                                  common='No',
                                  size='D - 1024'))
}
datp = rbind(datp1, datp2)
datp$common = factor(datp$common, levels=c('No', 'Yes'))
datp$feature = factor(datp$feature, levels=rev(gsub('[0-9]*', '', colnames(dat))))
ggplot(datp, aes(x=feature, y=num, fill=common)) + geom_bar(stat='identity', position='stack') +
  scale_fill_manual(values=c('indianred1', 'dodgerblue'), name=c(''), labels=c('No', 'Yes')) + 
  facet_grid(~size) + theme_bw() + coord_flip() +
  theme(panel.background=element_blank(), axis.text.x=element_text(family='Times New Roman'),
        axis.text.y=element_text(size=12, face='bold', family='SimHei'),
        strip.text=element_text(size=15, face='bold', family='Times New Roman'), 
        legend.text=element_text(size=15, face='bold', family='SimHei')) +
  xlab('') + ylab('') + ggtitle('')

#### A+B 3 ####
dat1 = dat[1:60,]
dat2 = dat[61:299,]
datp1 = as.data.frame(table(dat1$common))
datp2 = as.data.frame(table(dat2$common))
datp1$prop = datp1$Freq / sum(datp1$Freq)
datp2$prop = datp2$Freq / sum(datp2$Freq)
datp1 = datp1[c(2,1),]
datp2 = datp2[c(2,1),]
datp1$ypos = cumsum(datp1$prop) - 0.5 * datp1$prop
datp2$ypos = cumsum(datp2$prop) - 0.5 * datp2$prop
datp = rbind(datp1, datp2)
datp$group = c('A: kappa = ka', 'A: kappa = ka', 'B: kappa = kb', 'B: kappa = kb')
datp$text = paste0(as.character(round(datp$prop * 100, 1)), '%')

# plot
ggplot(datp, aes(x='', y=prop, fill=Var1)) + geom_bar(stat='identity') +
  scale_fill_manual(values=c('indianred1', 'dodgerblue'), name='', labels=c('no', 'yes')) +
  xlab('') + ylab('') + coord_polar('y', start=-pi/2, direction=-1) + facet_grid(~group) + theme_void() +
  theme(strip.text=element_text(face='bold', size=24, family='Times New Roman', vjust=1), 
        legend.text=element_text(size=15, family='SimHei')) +
  geom_text(aes(y=ypos, label=text), color="black", size=6, family='Times New Roman')
