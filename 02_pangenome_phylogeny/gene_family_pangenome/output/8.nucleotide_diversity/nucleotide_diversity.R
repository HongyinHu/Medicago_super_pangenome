# 不同类型的泛基因组家族类型的核酸多样性计算
library(tidyverse)
library(ape)
library(pegas)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
row_num <- args[2]
outfile <- args[3]

nuc_div <- read_lines(file_name,skip = 1) %>%
  matrix(ncol=2,byrow = TRUE) %>%
  as.data.frame() %>%
  pull(V2) %>%
  str_split(pattern = "") %>% 
  unlist() %>%
  matrix(nrow = as.numeric(row_num), byrow = TRUE) %>%
  as.DNAbin() %>%
  nuc.div()

write(nuc_div,file=outfile)
# 
# 
# data <- read_lines("GF_000150.fasta.msa.nuc",skip = 1)
# data <- matrix(data,ncol=2,byrow = TRUE)
# data <- as.data.frame(data)
# data <- pull(data,V2)
# data <- str_split(data,pattern = "") %>% unlist()
# data <- subset(data,data != "-")
# data <- matrix(data,nrow = 5,byrow = TRUE)
# as.DNAbin(data) %>% nuc.div()
# 
# c("A","T","C","G","A","A","C","G") %>% 
#   matrix(nrow = 2,byrow = TRUE) %>% 
#   ape::as.DNAbin() %>% 
#   pegas::nuc.div()
