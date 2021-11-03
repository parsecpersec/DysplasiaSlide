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
