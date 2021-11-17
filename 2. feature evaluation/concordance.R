

setwd('D:/work/OLK_B/chi/')
library(ggplot2)
library(irr)
library(extrafont)

#### record ####
for (i in 2:17) {
  dat[,i][dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes'] = 'Yes'
  dat[,i][!(dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes')] = 'No'
}

kc = data.frame()
for (i in 2:17) {
  kc = rbind(kc, 
             data.frame(feature=i-1,
                        kappa=kappam.fleiss(data.frame(P1=dat1[,i], P2=dat2[,i], P3=dat3[,i]))$value))
}
write.csv(kc, 'kappa_C.csv', quote=F, row.names=F)

for (i in 2:17) {
  dat[,i][dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes'] = 'Yes'
  dat[,i][!(dat1[,i] == 'Yes' & dat2[,i] == 'Yes' & dat3[,i] == 'Yes')] = 'No'
}

kd = data.frame()
for (i in 2:17) {
  kd = rbind(kd, 
             data.frame(feature=i-1,
                        kappa=kappam.fleiss(data.frame(P1=dat1[,i], P2=dat2[,i], P3=dat3[,i]))$value))
}
write.csv(kd, 'kappa_D.csv', quote=F, row.names=F)

#### C+D 3 ####


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
datp$group = c('A: kappa = 0.29', 'A: kappa = 0.29', 'B: kappa = 0.28', 'B: kappa = 0.28')
datp$text = paste0(as.character(round(datp$prop * 100, 1)), '%')

# plot
ggplot(datp, aes(x='', y=prop, fill=Var1)) + geom_bar(stat='identity') +
  scale_fill_manual(values=c('indianred1', 'dodgerblue'), name='', labels=c('no', 'yes')) +
  xlab('') + ylab('') + coord_polar('y', start=-pi/2, direction=-1) + facet_grid(~group) + theme_void() +
  theme(strip.text=element_text(face='bold', size=24, family='Times New Roman', vjust=1), 
        legend.text=element_text(size=15, family='SimHei')) +
  geom_text(aes(y=ypos, label=text), color="black", size=6, family='Times New Roman')