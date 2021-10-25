setwd('')
library(ggplot2)

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
