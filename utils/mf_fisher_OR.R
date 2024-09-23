test_fisher <- function(count.dist) {
  library(purrr)
  library(dplyr)
  library(tidyr)
  
  count.dist.df <- as.data.frame.matrix(count.dist)
  count.dist.df$cluster <- rownames(count.dist.df)
  count.dist.long.df <- count.dist.df %>%
    pivot_longer(cols = -cluster, names_to = "type", values_to = "count")
  
  fisher_results <- purrr::map_dfr(seq_len(nrow(count.dist.long.df)), function(i) {
    row <- count.dist.long.df[i, ]
    this_cluster <- row$cluster
    this_type <- row$type
    this_count <- as.numeric(row$count)
    
    other_clusters <- sum(as.numeric(count.dist.long.df$count[count.dist.long.df$cluster != this_cluster & count.dist.long.df$type == this_type]))
    other_cluster_counts <- sum(as.numeric(count.dist.long.df$count[count.dist.long.df$cluster == this_cluster & count.dist.long.df$type != this_type]))
    
    contingency_table <- matrix(
      c(
        this_count,
        sum(as.numeric(count.dist.long.df$count[count.dist.long.df$cluster == this_cluster])) - this_count,
        other_clusters,
        sum(as.numeric(count.dist.long.df$count[count.dist.long.df$type == this_type])) - other_clusters
      ),
      ncol = 2
    )
    
    fisher_test <- fisher.test(contingency_table)
    
    data.frame(
      cluster = this_cluster,
      type = this_type,
      p.value = fisher_test$p.value,
      OR = fisher_test$estimate
    )
  })
  
  fisher_results <- fisher_results %>% 
    mutate(adj.p.value = p.adjust(p.value, "BH"))
  
  OR_matrix <- fisher_results %>%
    select(cluster, type, OR) %>%
    pivot_wider(names_from = type, values_from = OR)
  
  p_value_matrix <- fisher_results %>%
    select(cluster, type, adj.p.value) %>%
    pivot_wider(names_from = type, values_from = adj.p.value)
  
  return(list(p.value = p_value_matrix, OR = OR_matrix, raw.results = fisher_results))
}
