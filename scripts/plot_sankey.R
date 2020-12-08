#!/usr/bin/env Rscript 

args <- commandArgs(trailingOnly=TRUE)

basedir <- args[1]
outfile <- args[2]

library(tidyverse)

index <- length(str_split(basedir, "/")[[1]])

df <- list.files(path=basedir, pattern="clusters.tsv", recursive=T) %>%
      enframe %>%
      mutate(filename = paste0(basedir, "/", value)) %>%
      select(-name, -value) %>%
      mutate(df = map(filename, read_tsv)) %>%
      unnest() %>%
      mutate(ftokens = str_split(filename, "/")) %>%
      mutate(str_clusters = unlist(lapply(ftokens, `[[`, index+1))) %>%
      select(-ftokens) %>%
      mutate(group = paste(str_clusters,label,sep="_")) %>%
      mutate(num_clusters = str_sub(str_clusters, 2, length(str_clusters))) %>%
      select(-str_clusters, -label) %>% 
      mutate(num_clusters = as.numeric(num_clusters)) %>%
      arrange(num_clusters) %>%
      select(-filename) %>%
      spread(num_clusters, group)

numsamples <- nrow(df %>% distinct(sample))

library(networkD3)

links <- tibble()
for (k in 2:4) {
    lnks <- tibble(source=pull(df[,k]), target=pull(df[,k+1]), value=rep(1,numsamples))
    links <- bind_rows(links, lnks)
}

nodes <- data.frame(
    name = c(as.character(links$source),
    as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

p1 <- sankeyNetwork(Links = links, Nodes=nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value="value", NodeID="name", fontSize=15,
                   sinksRight=F)

saveNetwork(p1, file = outfile)
