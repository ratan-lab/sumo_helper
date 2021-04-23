#!/usr/bin/env Rscript

load_packages <- function(){
  library(tidyverse)
  library(ComplexHeatmap)
  # reticulate::use_python(Sys.which('python3'), required = TRUE) 
  # NOTE run above command before attaching the reticulate library to spcify correct python version/path if unable to open .npz files
  library(reticulate)
  # NOTE when using reticulate you may encounter an error in initialize_python which prevents python bindings being loaded,
  # to solve this issue make sure "--enable-shared" option was used when building python from source  
  library(ggpubr)
  library(viridis)
  library(latex2exp)
}

run_diagnostics <- function(dir_name){
  stopifnot(dir.exists(dir_name))
  outname <- paste0("sumo_diagnostics.", dir_name, ".pdf")
  k_dirs <- list.files(dir_name)[grepl('^k',list.files(dir_name))]
  np <- import("numpy")
  quality_metrics <- NULL
  factorizations <- NULL
  debug_flag <- FALSE
  max_iter <- c()
  costs <- NULL

  for (k_dir in k_dirs){
    k <- as.numeric(strsplit(k_dir, 'k')[[1]][-1])
    fname <- file.path(dir_name, k_dir, "sumo_results.npz")
    npz <- np$load(fname, allow_pickle = TRUE)
    config <- npz$f[['config']]
    max_iter <- c(max_iter, unlist(config[,2][config[,1] == "max_iter"]) %>% as.numeric())
    selected_eta <- unlist(config[,2][config[,1] == "sparsity"]) %>% as.numeric()
    debug_flag <- any(grepl("cost", npz$files))
    
    k_metrics <- tibble(pac=npz$f[['pac']], ccc=npz$f[['cophenet']], k=k)
    quality_metrics <- if (is.null(quality_metrics)) k_metrics else quality_metrics %>% full_join(k_metrics, by = c("pac", "ccc", "k"))
    steps <- c(npz$f[['steps']])
    
    if (debug_flag){
      cost_ncol <- length(npz$f[['cost0']])
      second_term <- sapply(0:(length(steps)-1), function(x){npz$f[[paste0('cost',x)]][cost_ncol]})
      first_term <- sapply(0:(length(steps)-1), function(x){npz$f[[paste0('cost',x)]][cost_ncol-1]})
      cost_tib <- tibble(step=steps, first_term=first_term, second_term=second_term, k=k, rep=1:length(steps))
      factorizations <- if (is.null(factorizations)) cost_tib else factorizations %>% 
        full_join(cost_tib, by = c("step", "first_term", "second_term", "k", "rep"))
    } else {
      step_tib <- tibble(step=steps, k=k, rep=1:length(steps))
      factorizations <- if (is.null(factorizations)) step_tib else factorizations %>% full_join(step_tib, by = c("step", "k", "rep"))
    }
    
    log_fname <- file.path(dir_name, k_dir, paste0("eta_",selected_eta,".log"))
    log_lines <- system(paste0("cat ",log_fname), intern = T)
    n <- NULL
    for (line in log_lines){
      if (is.null(n) || n <= 5){
        if (grepl("N=", line)){
          n <- strsplit(strsplit(line, "N=")[[1]][-1], "\\)") %>% unlist() %>% as.numeric()
        } 
        if (grepl("Step", line)){
          s = strsplit(line, '\t')[[1]]
          step_tib <- tibble(step=strsplit(strsplit(s[1], 'Step\\(')[[1]][-1], "\\)")[[1]][1] %>% as.numeric(),
                             cost=as.numeric(strsplit(s[2], ":")[[1]][-1]), n=n, k=k)
          costs <- if (is.null(costs)) step_tib else costs %>% full_join(step_tib, by = c("step", "cost", "n", "k"))
        }
      } else {
        break
      }
    }
  }
  
  print(paste0("Saving results to: ", outname))
  pdf(file=outname, width = 13, height = 7)
  
  metric_p <- quality_metrics %>%
    gather(metric, val, -k) %>%
    group_by(k, metric) %>%
    summarise(minval=min(val), maxval=max(val), medval=median(val), .groups = 'drop') %>%
    mutate(metric=as.factor(metric), k=as.factor(k)) %>%
    ggplot() +
    geom_line(aes(x=k, y=medval, group=metric, color=metric), size=1) +
    geom_point(aes(x=k, y=medval, color=metric)) +
    geom_ribbon(aes(x=k, ymin=minval, ymax=maxval, group=metric, fill=metric), alpha=0.4) +
    theme_bw() +
    labs(y="metric value") +
    theme(legend.position = "bottom")
  
  p <- ggarrange(metric_p + ylim(0,1) + ggtitle("K selection"), metric_p + facet_wrap(metric~., scales="free", ncol=1), 
            nrow=1, ncol=2, common.legend = TRUE, legend = "bottom")
  print(p)
  
  if (debug_flag){
    c1_p <- factorizations %>%
      mutate(k=as.factor(k)) %>%
      ggplot() +
      geom_line(aes(x=rep, y=first_term, group=k, color=k), size=1) + 
      labs(x="factorization repetition", y="final cost function value [first term]", 
           title=TeX("$\\sum_{i=1}^{t}\\lambda_i|| W_i * (A_i - H S_i H^T)||_F^2$")) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    c2_p <- 
      factorizations %>%
      mutate(k=as.factor(k)) %>%
      ggplot() +
      geom_line(aes(x=rep, y=second_term, group=k, color=k), size=1) + 
      labs(x="factorization repetition", y="final cost function value [second term]", 
           title=TeX("$eta||H||_F^2$")) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    p2 <- ggarrange(c1_p, c2_p, nrow=1, ncol=2, common.legend = TRUE, legend = "bottom", align = "h")
    print(p2)
  }
  
  costs_p <- costs %>%
    mutate(k=as.factor(k)) %>%
    ggplot() + 
    geom_line(aes(x=step, y=cost, group=n, color=k), alpha=0.4, size=1) +
    # geom_point(aes(x=step, y=cost, group=n, color=k), alpha=0.4) +
    facet_wrap(k~., scales = "free") +
    theme_bw() +
    theme(legend.position = "null") +
    labs(y="cost function value [all terms]", x="steps", title = "Cost function value throughout the factorization (selected runs)")
  print(costs_p)
  
  steps_p <- factorizations %>%
    mutate(k=as.factor(k)) %>%
    ggplot() +
    geom_hline(yintercept = max_iter[1], color="red", size=1, alpha=0.4) +
    geom_point(aes(x=rep, y=step, group=k), size=1.5) + 
    facet_wrap(k~.) +
    labs(x="factorization repetition (red line marks set 'max_iter')", y="number of iterations/steps", 
         title="Number of iterations reached by solver") +
    theme_bw()
  print(steps_p)
  
  for (k_dir in k_dirs){
    k <- as.numeric(strsplit(k_dir, 'k')[[1]][-1])
    fname <- file.path(dir_name, k_dir, "sumo_results.npz")
    npz <- np$load(fname, allow_pickle = TRUE)
    
    ucon <- Heatmap(npz$f[['unfiltered_consensus']], cluster_rows = F, cluster_columns = F, 
                    name='Consensus unfiltered', heatmap_legend_param = list(direction = "horizontal"))
    con <- Heatmap(npz$f[['consensus']], cluster_rows = F, cluster_columns = F, name='Consensus filtered',
                   heatmap_legend_param = list(direction = "horizontal"))
    draw(ucon + con, heatmap_legend_side = "bottom", column_title=paste0("CONSENSUS ", k_dir))
    
    labels_fname <- file.path(dir_name, k_dir, "clusters.tsv")
    labels <- read_tsv(labels_fname, col_types = cols())
    
    steps <- c(npz$f[['steps']])
    rand_idx <- sample(0:(length(steps)-1), 3, replace = FALSE)
    debug_flag <- any(grepl("cost", npz$files))
    
    if(debug_flag){
      hlist <- lapply(1:3, function(x){npz$f[[paste0('h', rand_idx[x])]]})
      grid.newpage()
      for (i in 1:length(hlist)){
        pushViewport(viewport(x = (i-1)/3, width = 1/3, just = "left"))
        indices <- sapply(npz$f[[paste0('samples', rand_idx[i])]], function(x){grep(paste0(x,'$'),labels$sample)})
        ra <- rowAnnotation(labels = labels$label[indices], annotation_legend_param = list(labels = list(direction = "horizontal")),
                            col = list(labels=setNames(viridis_pal(option = "D")(k), 0:(k-1)))) 
        h <- Heatmap(hlist[[i]], name=paste0("H [rep:",rand_idx[i],"]"), heatmap_legend_param = list(direction = "horizontal"), 
                     split=labels$label[indices], cluster_row_slices = TRUE)
        draw(h+ra, newpage = FALSE, heatmap_legend_side = "bottom", column_title=paste0("Randomly selected final H matrices for ", k_dir))
        popViewport()
      }
    }
  }
  
  dev.off()
  
}

args = commandArgs(trailingOnly = T)

if (length(args) <= 0){
  print("./sumo_diagnostics.R sumo_results_dir1,sumo_results_dir2,...")
} else {
  load_packages()
  for (result_dir in args){
    print(paste0("Analyzing ", result_dir, " directory"))
    run_diagnostics(dir_name=result_dir)
  } 
}

# args <- c("miter1000_tol1e-09", "miter10000_tol1e-09", "miter15000_tol1e-08", "miter15000_tol1e-09",  
#   "miter15000_tol1e-10", "miter500_tol1e-08", "miter500_tol1e-09", "miter500_tol1e-10", "miter5000_tol1e-09" )
