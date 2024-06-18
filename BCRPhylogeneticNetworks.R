
####################################################################################################################
#                                                                                                                  #
# Script written by Peter Thomas, Laura McCoy's Lab, UCL IIT, 2024                                                 #
#                                                                                                                  #
# Manuscript currently in preparation                                                                              #
#                                                                                                                  #
# Outline:                                                                                                         #
#   Script creates B cell phylogenetic networks from tip to tip distances in phylogenetic trees.                   #
#   Briefly, it measures the tip to tip distances between nodes in igphyml generated trees, aggregates these       #
#   distances to produce an overall threshold per set of trees, and creates a directed B cell network based on     #
#   these values. From these networks, various node and network centrality scores are computed. The aim is to      #
#   understand B cell evolutionary trajectory based on its importance in the network (i.e., how many close         #
#   neighbours it has.                                                                                             #
#                                                                                                                  #
# Requirements:                                                                                                    #
#   A group of source files (metadata, can be either csv or tsv/tab) with:                                         #
#     * a unique sequence identifier column, named 'sequence_id'.                                                  #
#     * a lineage definitions column, named 'clone_id'                                                             #
#   A group of lineage tree files, each one corresponding to the data in each source file                          #
#     * For example, I have a data frame comprising lineages with all tested MBC subsets                           #
#     ('all-subsets-igphyml.tsv'), and a matching igphyml output file ('all-subsets-igphyml_igphyml-pass.tab')     #
#     * 'all-subsets-igphyml.tsv' was the source data frame input to the BuildTrees.py script to generate igphyml  #
#     trees.                                                                                                       #
#     * sequence_id column should be identical between data tables to allow data joining and analysis              #
#                                                                                                                  #
# Outputs:                                                                                                         #
#   * A set of networks per tree type, as well as their summary statistics                                         #
#   * Statistical analysis of network descriptors (hub score and degree ratio [in vs out is favoured])             #
#   * A selection of network plots (if plot_trees == T)                                                            #
#   * Bootstrap resampling statistics to check if comparisons are robust                                           #
#                                                                                                                  #
####################################################################################################################

load_installPackages <- function(){
  
  ### vector for packages needed
  requiredPackages <- c('ape', 'igraph', 'alakazam', 'adephylo', 'tidyverse', 'foreach', 'doParallel', 'data.table', 'ggpubr')
  versions <- c("5.8", "2.0.3", "1.3.0", "1.1.16", "2.0.0", "1.5.2", "1.0.17", "1.15.4", "0.6.0")
  
  ### get all packages installed
  packages <- installed.packages()
  
  ### find packages not installed
  packages_to_install <- requiredPackages %in% packages == F
  
  ### loop through vector and install ones marked as T
  for(i in 1:length(packages_to_install)){
    
    if(packages_to_install[i] == T){
      
      cat(paste('>Package', requiredPackages[i], 'not installed. Installing...'))
      
      install.packages(requiredPackages[i])
      
    }
    
  }
  
  cat('>Loading packages...\n')
  
  ### loop through packages and load all sequentially
  for(i in 1:length(packages_to_install)){
    
    library(requiredPackages[i], character.only = TRUE)
    
  }
  
  message(paste('>Original script uses ',
                paste(paste(requiredPackages, paste('v', versions, sep = '')), collapse = ', '),
                '. Compatibility is not assessed for other package versions.', sep = ''))
  
  cat('>Packages loaded\n')
  
}

getThreshold = function(distances_){
  
  t_ = sort(distances_) # arrange the distances in ascending order
  y_ = density(t_) # create a density curve
  
  differentiate_ =  diff(y_$y) / diff(y_$x) # work out the lagged differences between points
  
  cutoff_ = pracma::findpeaks(differentiate_) # find the peaks in this dataset (where the data changes most rapidly)
  cutoff_pos = cutoff_[which.max(cutoff_[,1]), 2] # take the midpoint of the tallest peak (where most data is concentrated)
  cutoff_ = y_$x[cutoff_pos] # select this point from the density curve
  
  plot(differentiate_) # visualise the plot
  abline(v = cutoff_pos) # draw the cut off
  
  # visualise the plot, but with the density curve rather than differenced curve
  print(ggplot()+
          geom_density(aes(x = t_), linewidth = 1)+
          geom_vline(xintercept = cutoff_, colour = 'firebrick', linewidth = 0.8)+
          theme_bw()+
          labs(y = 'Density',
               x = 'Phylogenetic distance')+
          theme(text = element_text(size = 15)))
  
  return(cutoff_)
  
}

networkStatistics <- function(g){
  
  node_in_degree = degree(g, mode = 'in')
  node_out_degree = degree(g, mode = 'out')
  nodes_ = length(V(g)$name)
  
  eigen_ = eigen_centrality(g, directed = T, scale = T)
  
  close_in = closeness(g, mode = 'in', normalized = T)
  close_out = closeness(g, mode = 'out', normalized = T)
  
  pr_ = page_rank(g, directed = T, damping = page_rank_damping)
  
  
  scaled_node_in_degree = node_in_degree / nodes_
  scaled_node_out_degree = node_out_degree / nodes_
  average_in_degree_network = mean(scaled_node_in_degree)
  average_out_degree_network = mean(scaled_node_out_degree)
  
  return(list(in_degree = node_in_degree,
              out_degree = node_out_degree,
              in_degree.scaled = scaled_node_in_degree,
              out_degree.scaled = scaled_node_out_degree,
              in_degree.scaled.network = average_in_degree_network,
              out_degree.scaled.network = average_out_degree_network,
              eigenvector = eigen_$vector,
              eigenvalue = eigen_$value,
              closeness_in = close_in,
              closeness_out = close_out,
              page_rank = pr_$vector,
              authority_ = igraph::authority.score(g)$vector,
              hubs_ = igraph::hub.score(g)$vector))
  
}

plotNetworks <- function(network_list, tree_id,
                         density_plot = 'hubs_', node_size = 7, node_colour = 'steelblue'){
  
  g_ = network_list[[tree_id]]$graph
  
  layout(mat = matrix(c(2, 1), 
                      nrow = 2, 
                      ncol = 1),
         heights = c(0.4, 1),    # Heights of the two rows
         widths = c(1, 1))
  
  # Plot 1: Network
  par(mar=c(0,0,0,0))
  plot.igraph(g_, vertex.label = NA, vertex.size = node_size, edge.arrow.size = 0.5,
              vertex.color = node_colour, edge.color = 'black')
  
  # Plot 2: Top (height) boxplot
  par(mar=c(0,0,0,0)+2)
  boxplot(network_list[[tree_id]][[density_plot]],
          horizontal = T,
          ylim = c(-0.05, 1.05),
          frame = F,
          col = 'steelblue')
  
  
}

repeatedTestParallel <- function(network_stats,
                                 source_data_frame_list,
                                 mode = c('Hub', 'Degree.ratio'),
                                 stats_group = 'isotype',
                                 nproc = 6,
                                 repetitions = 1000,
                                 sample_size = 100,
                                 replacement = F,
                                 seed = 1234){
  
  cl = makeCluster(nproc)
  registerDoParallel(cl)
  
  stats_ = foreach(i=1:repetitions) %dopar% {
    
    library(tidyverse)
    
    tmp = do.call('rbind',
                  lapply(network_stats, function(x) x %>% 
                           group_by(Type) %>% 
                           slice_sample(n = sample_size, replace = replacement)))
    tmp = left_join(tmp,
                    do.call('rbind', source_data_frame_list) %>%
                      select(sequence_id, clone_id, isotype))
    
    formula = paste(mode, Type, sep = ' ~ ')
    
    stats = compare_means(as.formula(formula), tmp, group.by = stats_group, p.adjust.method = 'BH')
    stats$comparison = paste(stats$group1, stats$group2, sep = '--')
    return(stats)
    
  }
  
  stopImplicitCluster()
  
  return(do.call('rbind', stats_))
  
}

##### load libraries and define file paths and constants (e.g., colour palettes) #####

load_installPackages()

# rename the below as you need

tree_path_ = '/path/to/igphyml/files/'
file_pattern_ = 'igphyml_igphyml' # a string variable to allow specific processing of only
metadata_path_ = '/path/to/metadata/files/'
metadata_delim_pattern= '\\tsv$'

plot_trees = F

palette_ = c('#fa800e', '#1ca5a0', '#8000FF')

#####

##### read data #####

trees_ = lapply(list.files(tree_path_,
                           pattern = file_pattern_,
                           full.names = T),
                readIgphyml, format = 'phylo', branches = 'distance')

names(trees_) = gsub('-igphyml.*?$', '',
                     list.files(tree_path_, pattern = file_pattern_))

data_file = lapply(list.files(metadata_path_,
                              pattern = metadata_delim_pattern,
                              full.names = T),
                   fread)

data_file = lapply(data_file, as.data.frame)

names(data_file) = gsub('-igphyml.*?$', '',
                        list.files(metadata_path_,
                                   pattern = metadata_delim_pattern))

#####

##### calculate the network statistics #####

trees_dist_ = lapply(1:length(data_file), function(i){
  
  # data file is currently a list of data frames, each of which contains lots of unsuitable lineages
  # select only those which are in the trees object
  tmp_df = data_file[[i]][data_file[[i]]$clone_id %in% names(trees_[[i]]$trees),]
  
  # split the tmp_df by into data frames per clone
  tmp_df = split(tmp_df, tmp_df$clone_id)
  # use this line to reorder the list of data frames to be in the same order as the trees object
  tmp_df = tmp_df[names(trees_[[i]]$trees)]
  
  # calculate the tip to tip distance for each tree
  trees_dist_ =  lapply(trees_[[i]]$trees, distTips)
  
  # look through the tip distances and select the closest connection
  trees_dist_flatten = lapply(trees_dist_[[i]], function(x){
    
    x = as.matrix(x)
    out_ = rep(NA, ncol(x))
    for(ii in 1:ncol(x)){
      
      mod_ = x[,ii]
      mod_[ii] = 1000 # set the self-connection high to remove it from comparisons
      
      pos_ = which.min(mod_) # select the closest connection
      out_[ii] = x[pos_, ii] # store in the output vector
      
    }
    
    return(out_)
    
  })
  
  # combine everything into a single vector per tree group
  trees_dist_flatten = do.call('c', trees_dist_flatten)
  
  # apply the threshold finding function to cut the tip distances
  threshold = getThreshold(trees_dist_flatten)
  
  # choose a value for page rank damping (0.4 is selected arbitrarily here)
  page_rank_damping = 0.4
  
  # create an adjacency matrix based on the calculated threshold
  trees_dist_ = lapply(trees_dist_, function(x, t_ = threshold){
    
    x = as.matrix(x) # convert to distance object to square matrix
    
    adj_matrix = apply(x, 2, function(y, threshold = t_){
      
      return(ifelse(y <= threshold, 1, 0)) # return 1 if below threshold, else 0
      
    })
    
    for(i in 1:ncol(adj_matrix)){
      adj_matrix[i,i] = 0 # zero the diagonal
    }
    
    germ_ = grep('GERM', colnames(x)) # find the position of the germline sequence
    
    for(i in 1:ncol(x)){
      
      links_ = which(adj_matrix[,i] == 1) # find the edges for each node
      
      if(is_empty(links_)){ 
        
        next
        
      } else { # if the clone has edges...
        
        links_ = links_[links_ != i] # remove self edges
        germ_dist_source = x[germ_, i] # find distance from germline to the query sequence
        germ_dist_rest = x[links_, germ_] # get the germline distance of each linked clone
        closer_ = ifelse(germ_dist_source < germ_dist_rest, 1, 0) # if the query is closer than the rest, assign 1, else 0
        
        adj_matrix[links_,i] = closer_ # input into the table
        
      }
      
    }
    
    # remove the germline sequence from the matrix
    adj_matrix = adj_matrix[1:(germ_-1), 1:(germ_-1)] 
    
    if(length(unique(as.numeric(adj_matrix))) > 1){
      
      # https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.link_analysis.hits_alg.hits.html
      hits_score = networkR::hits(adj_matrix)
      
    } else {
      
      # create an empty table if the adjacency matrix is size 1
      hits_score = list(authorities = rep(0, ncol(adj_matrix)),
                        hubs = rep(0, ncol(adj_matrix)),
                        iterations = ncol(adj_matrix))
      
    }
    
    # create a graph from this adjacency matrix
    g = graph_from_adjacency_matrix(adj_matrix, mode = 'directed') 
    
    # take values for the final output
    edgelist = as_edgelist(g) # convert the graph to a list of edges
    node_in_degree = degree(g, mode = 'in') # calculate the in degree, as a directed network
    node_out_degree = degree(g, mode = 'out') # calculate the out degree, as a directed network
    nodes_ = length(V(g)$name) # get the names of the nodes
    
    close_in = closeness(g, mode = 'in', normalized = T) # calculate closeness centrality of paths to a node
    close_out = closeness(g, mode = 'out', normalized = T) # calculate closeness centrality of paths from a node
    
    pr_ = page_rank(g, directed = T, damping = page_rank_damping) # calculate the page rank score per node
                                                                  # can adjust the damping factor as required
                                                                  # 0.4 is a randomly selected value
    
    scaled_node_in_degree = node_in_degree / nodes_ # normalise the in degree per node
    scaled_node_out_degree = node_out_degree / nodes_ # normalise the out degree per node
    average_in_degree_network = mean(scaled_node_in_degree) # average the in degree over the network
    average_out_degree_network = mean(scaled_node_out_degree) # average the out degree over the network
    
    # return network descriptors
    return(list(graph = g,
                in_degree = node_in_degree,
                out_degree = node_out_degree,
                in_degree.scaled = scaled_node_in_degree,
                out_degree.scaled = scaled_node_out_degree,
                in_out.ratio.scaled = scaled_node_in_degree / scaled_node_out_degree,
                in_degree.scaled.network = average_in_degree_network,
                out_degree.scaled.network = average_out_degree_network,
                closeness_in = close_in,
                closeness_out = close_out,
                page_rank_raw = pr_$vector,
                page_rank_minmax = (pr_$vector - min(pr_$vector)) / max(pr_$vector - min(pr_$vector)),
                authority_ = hits_score$authorities,
                hubs_ = hits_score$hubs))
    
    
  })
  
  # ensure consistent ordering of the data frame to the network statistics
  tmp_df = tmp_df[names(trees_dist_)]
  
  # loop through the network statistics list
  for(ii in 1:length(trees_dist_)){
    
    # reorder the data frame based on the names (sequence_id) of the node names, and take the subset
    cell.type = tmp_df[[ii]][match(V(trees_dist_[[ii]]$graph)$name,
                                   table = tmp_df[[ii]]$sequence_id), 'subset']
    
    # add the cell type as a vertex attribute
    V(trees_dist_[[ii]]$graph)$cell.type = cell.type
    
  }
  
  # return the labelled network statistics
  return(trees_dist_)
  
})

### reformat the network statistics output
network_statistics_ = list()
for(i in 1:length(trees_dist_)){
  
  tmp = lapply(trees_dist_[[i]], function(y){
    
    return(data.frame(sequence_id = V(y$graph)$name,
                      Type = V(y$graph)$cell.type,
                      In.Degree = y$in_degree.scaled,
                      Out.Degree = y$out_degree.scaled,
                      Degree.ratio = y$out_degree / y$in_degree,
                      In.Closeness = y$closeness_in,
                      Out.Closeness = y$closeness_out,
                      Page.Rank = y$page_rank_raw,
                      Page.Rank.Min_Max = y$page_rank_minmax,
                      Authority = y$authority_,
                      Hub = y$hubs_))
    
  })
  
  network_statistics_[[i]] = do.call('rbind', tmp)
  
}
rm(tmp)

#####

##### plot first 10 networks per set of trees (default is false) #####

if(plot_trees == T){

  for(i in 1:length(trees_dist_)){
    
    for(j in 1:length(trees_dist_[[i]])){
      
      plotNetworks(trees_dist_[[i]],
                   tree_id = j)
      
      if(j > 10){
        break
      }
      
    }
    
  }
}

layout(mat = matrix(c(1, 1), 
                    nrow = 1, 
                    ncol = 1),
       heights = c(1),    # Heights of the two rows
       widths = c(1))

#####

##### plot all comparisons #####

network_statistics_ = do.call('rbind', network_statistics_)
network_statistics_ = left_join(network_statistics_,
                     do.call('rbind', data_file) %>% 
                       mutate(sequence_id = gsub(',', '-', sequence_id)) %>%
                       select(sequence_id, clone_id, isotype))

# fix degree ratios, as infinite values can be produced
network_statistics_$Degree.ratio = ifelse(is.infinite(network_statistics_$Degree.ratio), 0, network_statistics_$Degree.ratio)

# plot the hub score
ggplot(network_statistics_,
       aes(x = Type, y = Hub, colour = Type))+
  geom_boxplot()+
  theme_bw()+
  scale_y_continuous(limits = c(NA, 1.5), breaks = seq(0, 1, 0.2))+
  facet_grid(rows = vars(isotype))+
  labs(y = 'Hub score',
       x = 'Subset',
       colour = 'Subset')+
  scale_colour_manual(values = palette_)+
  theme(text = element_text(size = 15),
        legend.position = 'none',
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'))

# Statistically compare the result, using benjamini-hochberg p-value adjustment
set.seed(42)
compare_means(Hub ~ Type,
              network_statistics_ %>% 
                mutate(donor = gsub('^(.*?)_\\d.*?$', '\\1', sequence_id)) %>% 
                group_by(Type, donor) %>% 
                slice_sample(n = 100),
              group.by = 'isotype',
              p.adjust.method = 'BH')

# plot the degree ratio (add 0.0001 to make log scale display easier)
ggplot(network_statistics_,
       aes(x = Type, y = Degree.ratio+0.0001, colour = Type))+
  geom_hline(yintercept = 1)+
  geom_boxplot()+
  scale_y_log10(limits = c(NA, 5000))+
  theme_bw()+
  labs(y = 'Degree Ratio',
       x = 'Subset')+
  facet_grid(rows = vars(isotype))+
  scale_colour_manual(values = palette_)+
  theme(text = element_text(size = 15),
        legend.position = 'none',
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'))

# Statistically compare the result, using benjamini-hochberg p-value adjustment
set.seed(42)
compare_means(Degree.ratio ~ Type,
              network_statistics_ %>% 
                ungroup() %>% 
                filter(!is.nan(Degree.ratio)) %>% 
                mutate(donor = gsub('^(.*?)_\\d.*?$', '\\1', sequence_id)) %>% 
                group_by(Type, donor) %>% 
                slice_sample(n = 100),
              group.by = 'isotype',
              p.adjust.method = 'BH')

#####


##### repeated stats measures #####

# bootstrap resampling statistics to ensure comparisons are valid
# (i.e., check that there is a significant p-value shift during multiple comparisons)
# hub score and degree.ratio are calculated since these are the values displayed in the manuscript

hub_stats = repeatedTestParallel(network_stats = network_statistics_,
                                 source_data_frame_list = data_file,
                                 stats_group = 'isotype',
                                 mode = 'Hub')

degree_stats = repeatedTestParallel(network_stats = network_statistics_,
                                    source_data_frame_list = data_file,
                                    stats_group = 'isotype',
                                    mode = 'Degree.ratio')

#####



