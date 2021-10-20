setwd('')

#### generate patch file ####
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
