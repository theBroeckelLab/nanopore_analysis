###TODO
###Adress sex so we can use X and Y if needed


library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

flist <- list.files(path = "/data/svi/prom/analysis/coverage_analysis/downsample_30x",full.names = TRUE, pattern = ".*\\_100k\\_coverage\\.txt")
flist

coverageResults <- read_tsv("/data/svi/prom/analysis/coverage_analysis/downsample_30x/GRCh38.100k.bed", col_names = c("chr","start","end"))

for (f in flist){
  samplename <- sub("/data/svi/prom/analysis/coverage_analysis/downsample_30x/","",f)
  samplename <- sub("_100k_coverage.txt","",samplename)
  sampleCov <- read_tsv(f)
  sampleCov <- sampleCov %>% select(rname,endpos,meandepth)
  
  sampleCov <- sampleCov %>% rename(!!samplename := "meandepth")
  
  coverageResults <- left_join(coverageResults,sampleCov, by = c("chr"="rname","end"="endpos"))
  
}

coverageResults <- coverageResults %>% 
  rowwise() %>% 
  mutate(
    mean = mean(c_across(4:ncol(coverageResults))),
    sd = sd(c_across(4:ncol(coverageResults)))
  )
chr.toplot = coverageResults[which(coverageResults$chr=="chr2"),]
chr.toplot = chr.toplot[which(chr.toplot$mean<100),]
ggplot(data = chr.toplot, aes(x=start, y=sd))+
  geom_point()
pairs(chr.toplot[,4:7])
