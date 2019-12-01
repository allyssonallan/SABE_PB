#read the library
library(data.table)

#input (map + ped)
apoea <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_APOE_29_10_2019.map")
apoeb <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_APOE_29_10_2019.ped")
foxo1a <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_FOXO1_30_10_2019.map")
foxo1b <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_FOXO1_30_10_2019.ped")
foxo3a <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_FOXO3_30_10_2019.map")
foxo3b <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_FOXO3_30_10_2019.ped")
gpx1a <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_GPX1_30_10_2019.map")
gpx1b <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_GPX1_30_10_2019.ped")
il10a <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_IL10_30_10_2019.map")
il10b <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_IL10_30_10_2019.ped")
pon1a <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_PON1_30_10_2019.map")
pon1b <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_PON1_30_10_2019.ped")
tp53a <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_TP53_30_10_2019.map")
tp53b <- fread("~/Área de Trabalho/BDS/Genes/73_IDOSOS_TP53_30_10_2019.ped")

#reservando os snps dos arquivos map de todos os genes
apoea1 <- t(apoea$V2)
foxo1a1 <- t(foxo1a$V2)
foxo3a1 <- t(foxo3a$V2)
gpx1a1 <- t(gpx1a$V2)
il10a1 <- t(il10a$V2)
pon1a1 <- t(pon1a$V2)
tp53a1 <- t(tp53a$V2)

#fusionando os arquivos dos snps
all <- cbind(apoea1, foxo1a1, foxo3a1, gpx1a1, il10a1, pon1a1, tp53a1)

#retirando os genotipos dos arquivos ped
apoeb1 <- apoeb[,-c(1:6)]
foxo1b1 <- foxo1b[,-c(1:6)]
foxo3b1 <- foxo3b[,-c(1:6)] 
gpx1b1 <- gpx1b[,-c(1:6)]
il10b1 <- il10b[,-c(1:6)]
pon1b1 <- pon1b[,-c(1:6)]
tp53b1 <- tp53b[,-c(1:6)]

#juntando os genotipos extraidos
all1 <- cbind(apoeb1, foxo1b1, foxo3b1, gpx1b1, il10b1, pon1b1, tp53b1)

#juntando rs com os genotipos
all2 <- rbind(all, all1, use.names=F)
all3 <- setnames(all1, all)
all4 <- cbind(apoeb$V2, all3)

library(Hmisc)
results <- describe(all3)

install.packages("readr")
library(readr)
fwrite(results1, "~/Área de Trabalho/all.txt")

descriptives <- as.data.table(table(unlist(all3)))

write.table(results1, "~/Área de Trabalho/results.txt", col.names=TRUE)
              