#!/bin/Rscript

#read library

library(data.table)

#read files and filenames

files <- (Sys.glob("~/Área de Trabalho/BDS/rfmix/*.bed"))
filenames <- list.files(path = "~/Área de Trabalho/BDS/rfmix/", pattern = ".bed")
names(files) <- c(filenames)

#store all files

df1 <- lapply(files, function(x) read.table(x, header = F))

#the journey by genes and snps

####################################APOE#####################################################

#rs440446
chr19 <- lapply(df1, function(x) x[x$V1 == 19,])
chr19a <- lapply(chr19, function(x) x[x$V2 <= 45409167 & x$V3 >= 45409167,])
for(i in 1:length(chr19a)){
  write.table(data.frame(chr19a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
}
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs440446.txt', row.names = F)

#rs769450
chr19 <- lapply(df1, function(x) x[x$V1 == 19,])
chr19a <- lapply(chr19, function(x) x[x$V2 <= 45410444 & x$V3 >= 	45410444,])
for(i in 1:length(chr19a)){
  write.table(data.frame(chr19a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs769450.txt', row.names = F)

#rs429358
chr19 <- lapply(df1, function(x) x[x$V1 == 19,])
chr19a <- lapply(chr19, function(x) x[x$V2 <= 45411941 & x$V3 >= 45411941,])
for(i in 1:length(chr19a)){
  write.table(data.frame(chr19a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs429358.txt', row.names = F)

#rs7412
chr19 <- lapply(df1, function(x) x[x$V1 == 19,])
chr19a <- lapply(chr19, function(x) x[x$V2 <= 45412079 & x$V3 >= 45412079,])
for(i in 1:length(chr19a)){
  write.table(data.frame(chr19a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs7412.txt', row.names = F)

##########################################ACE####################################
#rs1799752
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 61565904 & x$V3 >=	61565904,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs1799752.txt', row.names = F)

#################################FOXO3###############################################
#rs2764264
  chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
  chr6a <- lapply(chr6, function(x) x[x$V2 <= 108934461 & x$V3 >=	108934461,])
  for(i in 1:length(chr6a)){
    write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
    cat('\n',file='test.txt',append = TRUE)  
  } 
  test <- fread("~/Área de Trabalho/test.txt")
  result <- data.frame(filenames, test)
  write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2764264.txt', row.names = F)

#rs13217795
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108974098 & x$V3 >=	108974098,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs13217795.txt', row.names = F)

#rs13220810
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108913201 & x$V3 >= 108913201,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs13220810.txt', row.names = F)

#rs2802292
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108908518 & x$V3 >= 108908518,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2802292.txt', row.names = F)

#rs74998990
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108887366 & x$V3 >= 108887366,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs74998990.txt', row.names = F)

#rs9384681
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108890069 & x$V3 >= 108890069,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs9384681.txt', row.names = F)

#rs12055603
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108892235 & x$V3 >= 108892235,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs12055603.txt', row.names = F)

#rs2883881
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108897386 & x$V3 >= 108897386,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2883881.txt', row.names = F)

#rs12206094
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108906200 & x$V3 >= 108906200,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs12206094.txt', row.names = F)

#rs117831673
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108940824 & x$V3 >= 108940824,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs117831673.txt', row.names = F)

#rs117831673
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108940824 & x$V3 >= 108940824,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs117831673.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs7762395
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108945107 & x$V3 >= 108945107,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs7762395.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs62428909
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108949277 & x$V3 >= 108949277,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs62428909.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs13207511
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108953975 & x$V3 >= 108953975,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs13207511.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs12207868
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108964885 & x$V3 >= 108964885,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs12207868.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs117420572
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108971213 & x$V3 >= 108971213,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs117420572.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs77061140
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108979306 & x$V3 >= 108979306,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs77061140.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs4946933
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108981117 & x$V3 >= 108981117,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs4946933.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs9398171
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108983527 & x$V3 >= 108983527,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs9398171.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs9374040
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108997435 & x$V3 >= 108997435,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs9374040.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3800232
chr6 <- lapply(df1, function(x) x[x$V1 == 6,])
chr6a <- lapply(chr6, function(x) x[x$V2 <= 108998953 & x$V3 >= 108998953,])
for(i in 1:length(chr6a)){
  write.table(data.frame(chr6a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3800232.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#################################TP53###############################################

#rs2909430
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 7578645 & x$V3 >= 7578645,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2909430.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs9895829
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 7578679 & x$V3 >= 7578679,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs9895829.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs1042522
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 7579472 & x$V3 >= 7579472,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs1042522.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs8078476
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 7581228 & x$V3 >= 7581228,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs8078476.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2078486
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 7583083 & x$V3 >= 7583083,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2078486.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs11652704
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 7584400 & x$V3 >= 7584400,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs11652704.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs17882854
chr17 <- lapply(df1, function(x) x[x$V1 == 17,])
chr17a <- lapply(chr17, function(x) x[x$V2 <= 7588680 & x$V3 >= 7588680,])
for(i in 1:length(chr17a)){
  write.table(data.frame(chr17a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs17882854.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#################################PON1###############################################

#rs854551
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94927677 & x$V3 >= 94927677,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs854551.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917577
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94927708 & x$V3 >= 94927708,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917577.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917576
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94927865 & x$V3 >= 94927865,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917576.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs854552
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94927924 & x$V3 >= 94927924,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs854552.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917572
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94929190 & x$V3 >= 94929190,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917572.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917558
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94932904 & x$V3 >= 94932904,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917558.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917551
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94934455 & x$V3 >= 94934455,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917551.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917549
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94935200 & x$V3 >= 94935200,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917549.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs662
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94937446 & x$V3 >= 94937446,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs662.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs397841819 (rs3917525)
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94940480 & x$V3 >= 94940480,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs397841819.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2299257
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94942765 & x$V3 >= 94942765,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2299257.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917507
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94943884 & x$V3 >= 94943884,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917507.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs62467349
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94946795 & x$V3 >= 94946795,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs62467349.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2074351
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94947799 & x$V3 >= 94947799,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2074351.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs854562
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94947969 & x$V3 >= 94947969,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs854562.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3917490
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94948841 & x$V3 >= 94948841,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3917490.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2049649
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94949329 & x$V3 >= 94949329,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2049649.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2299260
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94949537 & x$V3 >= 94949537,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2299260.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2299261
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94949663 & x$V3 >= 94949663,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2299261.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs854569
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94950055 & x$V3 >= 94950055,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs854569.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2237584
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94950837 & x$V3 >= 94950837,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2237584.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs854570
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 94952692 & x$V3 >= 94952692,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs854570.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#################################IL-6###############################################

#rs1800795
chr7 <- lapply(df1, function(x) x[x$V1 == 7,])
chr7a <- lapply(chr7, function(x) x[x$V2 <= 22766645 & x$V3 >= 22766645,])
for(i in 1:length(chr7a)){
  write.table(data.frame(chr7a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs1800795.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#################################IL-10###############################################

#rs3024496
chr1 <- lapply(df1, function(x) x[x$V1 == 1,])
chr1a <- lapply(chr1, function(x) x[x$V2 <= 206941864 & x$V3 >= 206941864,])
for(i in 1:length(chr1a)){
  write.table(data.frame(chr1a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3024496.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs3024495
chr1 <- lapply(df1, function(x) x[x$V1 == 1,])
chr1a <- lapply(chr1, function(x) x[x$V2 <= 206942413 & x$V3 >= 206942413,])
for(i in 1:length(chr1a)){
  write.table(data.frame(chr1a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs3024495.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs1554286
chr1 <- lapply(df1, function(x) x[x$V1 == 1,])
chr1a <- lapply(chr1, function(x) x[x$V2 <= 206944233 & x$V3 >= 206944233,])
for(i in 1:length(chr1a)){
  write.table(data.frame(chr1a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs1554286.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs1518111
chr1 <- lapply(df1, function(x) x[x$V1 == 1,])
chr1a <- lapply(chr1, function(x) x[x$V2 <= 206944645 & x$V3 >= 206944645,])
for(i in 1:length(chr1a)){
  write.table(data.frame(chr1a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs1518111.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#################################FOXO1A###############################################

#rs17592236
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41131932 & x$V3 >= 41131932,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs17592236.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2755211
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41142223 & x$V3 >= 41142223,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2755211.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs75700692
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41145211 & x$V3 >= 41145211,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs75700692.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2755213
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41146301 & x$V3 >= 41146301,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2755213.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs2995991
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41151134 & x$V3 >= 41151134,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs2995991.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs79123773
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41186382 & x$V3 >= 41186382,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs79123773.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs12876443
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41196877 & x$V3 >= 41196877,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs12876443.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs17531346
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41208607 & x$V3 >= 41208607,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs17531346.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs9532571
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41218670 & x$V3 >= 41218670,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs9532571.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs1334241
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41223110 & x$V3 >= 41223110,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs1334241.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#rs9603777
chr13 <- lapply(df1, function(x) x[x$V1 == 1,])
chr13a <- lapply(chr13, function(x) x[x$V2 <= 41230062 & x$V3 >= 41230062,])
for(i in 1:length(chr13a)){
  write.table(data.frame(chr13a[[i]]),file='~/Área de Trabalho/test.txt',append=TRUE, row.names=F, col.names = F)
  cat('\n',file='test.txt',append = TRUE)  
} 
test <- fread("~/Área de Trabalho/test.txt")
result <- data.frame(filenames, test)
write.table(result, file='~/Área de Trabalho/BDS/Genes/rs9603777.txt', row.names = F)
file.remove("~/Área de Trabalho/test.txt")

#bind all snps from multiple files in a specific gene directory

##################APOE_results####################################

files <- (Sys.glob("~/Área de Trabalho/BDS/Genes/APOE/*.txt"))
filenames <- list.files(path="~/Área de Trabalho/BDS/Genes/APOE/", pattern = ".txt")
names(files) <- c(filenames)
df <- lapply(files, function(x) read.table(x, header = F))
df1 <- do.call(cbind, df)
write.table(df1, file='~/Área de Trabalho/BDS/Genes/APOE/apoe_results.txt', row.names = F)

##################FOXO1_results####################################

files <- (Sys.glob("~/Área de Trabalho/BDS/Genes/FOXO1/*.txt"))
filenames <- list.files(path="~/Área de Trabalho/BDS/Genes/FOXO1/", pattern = ".txt")
names(files) <- c(filenames)
df <- lapply(files, function(x) read.table(x, header = F))
df1 <- do.call(cbind, df)
write.table(df1, file='~/Área de Trabalho/BDS/Genes/FOXO1/foxo1_results.txt', row.names = F)

##################FOXO3_results####################################

files <- (Sys.glob("~/Área de Trabalho/BDS/Genes/FOXO3/*.txt"))
filenames <- list.files(path="~/Área de Trabalho/BDS/Genes/FOXO3/", pattern = ".txt")
names(files) <- c(filenames)
df <- lapply(files, function(x) read.table(x, header = F))
df1 <- do.call(cbind, df)
write.table(df1, file='~/Área de Trabalho/BDS/Genes/FOXO3/foxo3_results.txt', row.names = F)

##################IL-10_results####################################

files <- (Sys.glob("~/Área de Trabalho/BDS/Genes/IL-10/*.txt"))
filenames <- list.files(path="~/Área de Trabalho/BDS/Genes/IL-10/", pattern = ".txt")
names(files) <- c(filenames)
df <- lapply(files, function(x) read.table(x, header = F))
df1 <- do.call(cbind, df)
write.table(df1, file='~/Área de Trabalho/BDS/Genes/IL-10/il10_results.txt', row.names = F)

##################PON1_results####################################

files <- (Sys.glob("~/Área de Trabalho/BDS/Genes/PON1/*.txt"))
filenames <- list.files(path="~/Área de Trabalho/BDS/Genes/PON1/", pattern = ".txt")
names(files) <- c(filenames)
df <- lapply(files, function(x) read.table(x, header = F))
df1 <- do.call(cbind, df)
write.table(df1, file='~/Área de Trabalho/BDS/Genes/PON1/pon1_results.txt', row.names = F)

##################TP53_results####################################

files <- (Sys.glob("~/Área de Trabalho/BDS/Genes/TP53/*.txt"))
filenames <- list.files(path="~/Área de Trabalho/BDS/Genes/TP53/", pattern = ".txt")
names(files) <- c(filenames)
df <- lapply(files, function(x) read.table(x, header = F))
df1 <- do.call(cbind, df)
write.table(df1, file='~/Área de Trabalho/BDS/Genes/TP53/tp53_results.txt', row.names = F)
