# view #
for(p in path) {
  num = c()
  for(s in soft) {
    filename = paste0('./soft_patch/temp_small/', p, '_', s, '.csv')
    temp = read.csv(filename, header=F, stringsAsFactors=F)
    num = c(num, nrow(temp))
  }
  if(max(num) != min(num)) {print(paste0(num)); print(p)}
}

# remove #
for(p in path) {
  for(s in soft[6:7]) {
    filename = paste0('./soft_patch/temp_small/', p, '_', s, '.csv')
    temp = read.csv(filename, header=F, stringsAsFactors=F)
    temp2 = temp[1:(nrow(temp)/2),]
    write.table(temp2, filename, sep=',', quote=F, row.names=F, col.names=F)
  }
}
