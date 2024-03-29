setwd('')

dat1 = read.csv('./patch-C.csv', T, stringsAsFactors=F)
dat2 = read.csv('./patch-D.csv', T, stringsAsFactors=F)
dir1 = '../../Annotation/patch224_rename/'
dir2 = '../../Annotation/patch1024_rename/'

# files = c(paste0(1:16, '_large'), paste0(1:16, '_small'))
# for(f in files) {
#    dir.create(paste0('./questionable/', f))
# }
# small = files[17:32]
# large = files[1:16]

for(i in 2:ncol(dat1)) {
  files = dat1[,1][dat1[,i] == 'YN' | dat1[,i] == 'NY']
  for(f in files) {
    file.copy(from=paste0(dir1, f, '.jpeg'), 
              to=paste0('./questionable/', small[i-1], '/', f, '_', dat1[,i][dat1[,1] == f], '.jpeg'))
  }
}

for(i in 2:ncol(dat2)) {
  files = dat2[,1][dat2[,i] == 'YN' | dat2[,i] == 'NY']
  for(f in files) {
    file.copy(from=paste0(dir2, f, '.jpeg'), 
              to=paste0('./questionable/', large[i-1], '/', f, '_', dat2[,i][dat2[,1] == f], '.jpeg'))
  }
}

img = list.files(path='./questionable/', pattern='.jpeg', recursive=T)
img_file = gsub('^.*/', '', img)
img_file = gsub('_.*$', '', img_file)
img_file = unique(img_file)
img_C = img_file[grepl('C',img_file)]
img_D = img_file[grepl('D',img_file)]

for(f in img_C) {
  file.copy(from=paste0(dir1, f, '.jpeg'), 
            to=paste0('../../Annotation2/image/C/', f, '.jpeg'))
}

for(f in img_D) {
  file.copy(from=paste0(dir2, f, '.jpeg'), 
            to=paste0('../../Annotation2/image/D/', f, '.jpeg'))
}

img_C = img_C[order(img_C)]
img_D = img_D[order(img_D)]
write.csv(c(img_C, img_D), 'CD.csv', quote=F, row.names=F)

#### add previous anno ####
dat_C = dat1[dat1$patch %in% img_C,]
dat_C[dat_C == 'NN' | dat_C == 'YY'] = ''

dat_D = dat2[dat2$patch %in% img_D,]
dat_D[dat_D == 'NN' | dat_D == 'YY'] = ''

write.csv(dat_C, 'C.csv', quote=F, row.names=F)
write.csv(dat_D, 'D.csv', quote=F, row.names=F)

