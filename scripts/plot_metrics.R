#!/usr/bin/env Rscript 

args <- commandArgs(trailingOnly=TRUE)

basedir <- args[1]
outfile <- args[2]

library(tidyverse)
library(purrr)
library(reticulate)

np <- import("numpy")
get_metrics <- function(x) {
    f <- np$load(x$filename2)
    pac <- f$f[["pac"]]
    ccc <- f$f[["cophenet"]]
    return(tibble(PAC = pac, CCC = ccc))
}

df <- list.files(basedir, "sumo_results.npz", recursive=T) %>%
    enframe() %>%
    mutate(filename = paste0(basedir, "/", value)) %>%
    select(-value) %>%
    mutate(filename2 = filename) %>%
    group_by(name, filename) %>%
    nest() %>%
    mutate(tbl = purrr::map(data, get_metrics)) %>%
    select(-data) %>%
    unnest() %>%
    ungroup() %>%
    select(-name) %>%
    gather(metric, tbl, -filename) %>%
    group_by(filename, metric) %>%
    summarise(ymin=min(tbl), ymax=max(tbl), ymed=median(tbl))

dirlist <- str_split(df$filename, "/")
nummems <- length(dirlist[[1]])
df$name <- unlist(lapply(dirlist, `[[`, nummems - 2))
df$numclusters <- as.numeric(str_sub(unlist(lapply(dirlist, `[[`, nummems - 1)),2))


# I like drawing 16 PAC curves in a 7x7, but it can get crowded beyond that
side <- 7 * length(unique(df$name)) / 16
if (side < 7) {
    side <- 7
}

pdf(outfile, width=side, height=side)
df %>% 
    ggplot(aes(x=numclusters, y=ymed, ymin=ymin, ymax=ymax, 
               color=metric, fill=metric)) + 
        geom_ribbon(alpha=0.3) + 
        geom_line() + 
        geom_point() + 
        geom_hline(yintercept=0.1, color="grey") + 
        geom_hline(yintercept=0.95, color="grey") + 
        facet_wrap(name~.) + 
        theme_bw() +
        xlab("Number of clusters") + 
        ylab("metric")
dev.off()
