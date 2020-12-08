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

