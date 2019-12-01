#!/bin/Rscript

library(data.table)
setwd("~/Área de Trabalho/BDS/rfmix/")
filelist <- list.files("~/Área de Trabalho/BDS/rfmix/", pattern = "*.bed")
datalist <- lapply(filelist, read.table)
a <- fread("~/Downloads/TABELA ALYSSON.csv")

