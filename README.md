# sumo_helper

### Helper scripts

The following scripts are included in the scripts folder

1. plot_metrics.R : Frequently we want to run SUMO using various parameters, and then subsequently inspect the PAC and CCC curves. This scripts takes a directory name as input and plots the PAC and CCC curves for those directories in a single plot. Run it as
```
plot_metrics.R /Users/ar7jq/benchmark_analysis output.pdf
```

2. plot_sankey.R : This plot is helpful in inspecting how the cluster memberships change as the number of clusters are changed. Run it as 

```
./plot_sankey.R /Users/ar7jq/benchmark_analysis output.pdf
```

3. cluster_layer: Cluster the samples in a particular data type using Symmetric NMF.

```
./cluster_layer prepared.npz 0 2
```

Here "0" is the label for the layer, and you would like to see how the clustering performs with k=2.

4. examine_layers: Cluster each layer separately using Symmetric NMF from k=2:15. Also calculates adjusted NMI for these various clusterings.

```
./examine_layers -c clusters.txt prepare.npz
```

5. Create plot summary of SUMO run results, based on various metrics extracted from the output directory (for SUMO v0.2.7+).

```
./sumo_diagnostics.R sumo_results_dir1,sumo_results_dir2,...
```
Created plots illustrate metrics used for K-selection, cost function distribution throughout the factorization, number of iterations reached by solver in each factorization repetition and heatmaps of cosnsensus matrices.
If -log DEBUG flag was used for SUMO run, additional plots of final cost function values (separated by terms) and heatmaps of selected final H matrices are created.
