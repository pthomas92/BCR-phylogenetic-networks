# BCR-phylogenetic-networks

Script written by Peter Thomas, Laura McCoy's Lab, UCL IIT, 2024
Manuscript currently in preparation. A selection data will be available upon manuscript submission/acceptance.

Required packages:
- `ape` (script uses v5.8)
- `igraph` (script uses v2.0.3)
- `alakazam` (script uses v1.3.0)
- `adephylo` (script uses v1.1.16)
- `tidyverse` (script uses v2.0.0)
- `foreach` (script uses v1.5.2)
- `doParallel` (script uses v1.0.17)
- `datatable` (script uses v1.15.4)
- `ggpubr` (script uses v0.6.0)

By default, the script will check if a version of the above packages has been installed, and install the latest version if none exists. Functionality should be retained across later package versions but this is not guaranteed.


## Outline:
Script creates B cell phylogenetic networks from tip to tip distances in phylogenetic trees. Briefly, it measures the tip to tip distances between nodes in igphyml generated trees, aggregates these distances to produce an overall threshold per set of trees, and creates a directed B cell network based on these values. From these networks, various node and network centrality scores are computed. The aim is to understand B cell evolutionary trajectory based on its importance in the network (i.e., how many close neighbours it has).

## Requirements:
1. A group of source files (metadata, can be either csv or tsv/tab) with:
  * a unique sequence identifier column, named 'sequence_id'.
  * a lineage definitions column, named 'clone_id'                                                             
2. A group of lineage tree files, each one corresponding to the data in each source file                          
  * For example, I have a data frame comprising lineages with all tested MBC subsets (e.g., 'all-subsets-igphyml.tsv'), and a matching igphyml output file ('all-subsets-igphyml_igphyml-pass.tab')     
  * 'all-subsets-igphyml.tsv' was the source data frame input to the BuildTrees.py script to generate igphyml  
     trees.
  * sequence_id column should be identical between data tables to allow data joining and analysis              
## Outputs:                                                                                                         
1. A set of networks per tree type, as well as their summary statistics                                         
2. Statistical analysis of network descriptors (hub score and degree ratio [in vs out is favoured])
3. A selection of network plots (if plot_trees == T)
4. Bootstrap resampling statistics to check if comparisons are robust                                           
