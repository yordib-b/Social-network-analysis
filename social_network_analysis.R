# Load libraries
library(igraph)
library(ggplot2)
library(dplyr)
library(scales)

# Load data (skip the header line that starts with %)
df <- read.table("C:/Users/dell/Downloads/CollegeMsg.txt/arenas-jazz/out.arenas-jazz", header = FALSE, comment.char = "%", col.names = c("source", "target"))
cat("Shape:", nrow(df), "rows x", ncol(df), "columns\n")
head(df, 10)

# Data Cleaning
cat("Missing values:", sum(is.na(df)), "\n")
cat("Duplicates:", sum(duplicated(df)), "\n")
df_clean <- df[!duplicated(df), ]
df_clean <- df_clean[df_clean$source != df_clean$target, ]  # Remove self-loops
cat("Total edges:", nrow(df_clean), "\n")
cat("Unique musicians:", length(unique(c(df_clean$source, df_clean$target))), "\n")
head(df_clean)

# Create Network (undirected since collaborations are mutual)
G <- graph_from_data_frame(df_clean, directed = FALSE)
G <- simplify(G)  # Remove any duplicate edges

cat("Nodes:", vcount(G), "\n")
cat("Edges:", ecount(G), "\n")
#Only 14% of all possible connections exist
cat("Density:", round(edge_density(G), 4), "\n")
#Each musician worked with ~28 others on average
cat("Degree - Mean:", round(mean(degree(G)), 2), "Max:", max(degree(G)), "Min:", min(degree(G)), "\n")

# Network Visualization
#Central dense region: Core musicians who collaborated frequently with many others
#Peripheral nodes: Musicians with fewer connections
degrees <- degree(G)

set.seed(42)
plot(G, 
     layout = layout_with_fr(G),
     vertex.size = rescale(degrees, to = c(2, 12)),
     vertex.color = heat.colors(100)[100 - rescale(degrees, to = c(1, 99))],
     vertex.label = NA,
     edge.width = 0.5,
     edge.color = rgb(0.5, 0.5, 0.5, 0.3),
     main = "Jazz Musicians Collaboration Network")

# Degree Centrality
degree_cent <- degree(G)
top_degree <- sort(degree_cent, decreasing = TRUE)[1:10]
cat("Top 10 Musicians by Degree Centrality:\n")
print(top_degree)


# Degree Distribution
hist(degree_cent, breaks = 30, col = "steelblue", main = "Degree Distribution", xlab = "Degree")

# Betweenness Centrality
betweenness_cent <- betweenness(G, normalized = TRUE)
top_betweenness <- sort(betweenness_cent, decreasing = TRUE)[1:10]

cat("Top 10 Musicians by Betweenness Centrality:\n")
print(round(top_betweenness, 4))

# Closeness Centrality
closeness_cent <- closeness(G, normalized = TRUE)
top_closeness <- sort(closeness_cent, decreasing = TRUE)[1:10]

cat("Top 10 Musicians by Closeness Centrality:\n")
print(round(top_closeness, 4))

# Eigenvector Centrality
eigen_cent <- eigen_centrality(G)$vector
top_eigen <- sort(eigen_cent, decreasing = TRUE)[1:10]

cat("Top 10 Musicians by Eigenvector Centrality:\n")
print(round(top_eigen, 4))

# Create centrality summary dataframe
centrality_df <- data.frame(
  musician = V(G)$name,
  degree = degree_cent,
  betweenness = betweenness_cent,
  closeness = closeness_cent,
  eigenvector = eigen_cent
)

# Show top 10 by degree
top10 <- centrality_df[order(-centrality_df$degree), ][1:10, ]
print(top10)

# Connected Components
components <- components(G)
cat("Number of components:", components$no, "\n")
cat("Sizes of components:", sort(components$csize, decreasing = TRUE), "\n")
cat("Largest component size:", max(components$csize), "(", round(100*max(components$csize)/vcount(G), 1), "% of network)\n")


# Community Detection using Louvain algorithm
communities <- cluster_louvain(G)
cat("Number of communities:", length(communities), "\n")
cat("Modularity:", round(modularity(communities), 4), "\n")
cat("Community sizes:", sort(sizes(communities), decreasing = TRUE), "\n")

# Visualize Louvain Communities
# OUTPUT: Left = network colored by community; Right = bar chart of community sizes
set.seed(42)
par(mfrow = c(1, 2))

# Network plot with communities
plot(communities, G,
     layout = layout_with_fr(G),
     vertex.size = rescale(degrees, to = c(3, 10)),
     vertex.label = NA,
     edge.width = 0.3,
     main = "Louvain Community Structure")

# Community size distribution
comm_sizes <- sizes(communities)
barplot(sort(comm_sizes, decreasing = TRUE),
        col = rainbow(length(comm_sizes)),
        main = "Community Size Distribution",
        xlab = "Community", ylab = "Number of Musicians",
        names.arg = 1:length(comm_sizes))

par(mfrow = c(1, 1))

# Walktrap Community Detection (using random walks of length 4 and 5)
walktrap_4 <- cluster_walktrap(G, steps = 4)
walktrap_5 <- cluster_walktrap(G, steps = 5)

cat("=== Walktrap with 4 steps ===\n")
cat("Number of communities:", length(walktrap_4), "\n")
cat("Modularity:", round(modularity(walktrap_4), 4), "\n")
cat("Community sizes:", sort(sizes(walktrap_4), decreasing = TRUE), "\n")
cat("\n=== Walktrap with 5 steps ===\n")
cat("Number of communities:", length(walktrap_5), "\n")
cat("Modularity:", round(modularity(walktrap_5), 4), "\n")
cat("Community sizes:", sort(sizes(walktrap_5), decreasing = TRUE), "\n")


# Visualize Walktrap communities (using 4 steps)
set.seed(42)
par(mfrow = c(1, 2))

# Walktrap 4 steps
plot(walktrap_4, G,
     layout = layout_with_fr(G),
     vertex.size = rescale(degrees, to = c(3, 10)),
     vertex.label = NA,
     edge.width = 0.3,
     main = "Walktrap (4 steps)")

# Walktrap 5 steps
plot(walktrap_5, G,
     layout = layout_with_fr(G),
     vertex.size = rescale(degrees, to = c(3, 10)),
     vertex.label = NA,
     edge.width = 0.3,
     main = "Walktrap (5 steps)")

par(mfrow = c(1, 1))

# Edge Betweenness Community Detection (Girvan-Newman)
edge_betweenness_comm <- cluster_edge_betweenness(G)

cat("=== Edge Betweenness (Girvan-Newman) ===\n")
cat("Number of communities:", length(edge_betweenness_comm), "\n")
cat("Modularity:", round(modularity(edge_betweenness_comm), 4), "\n")
cat("Community sizes:", sort(sizes(edge_betweenness_comm), decreasing = TRUE), "\n")

# Plot Dendrogram for Edge Betweenness Community Detection
dendro <- as.dendrogram(edge_betweenness_comm)
plot(dendro, 
     main = "Dendrogram - Edge Betweenness (Girvan-Newman)",
     xlab = "Nodes", 
     ylab = "Height",
     leaflab = "none")  # Hide leaf labels for clarity since network is large

# Compare all community detection methods

cat("=== Community Detection Comparison ===\n\n")
cat(sprintf("%-20s %12s %12s\n", "Algorithm", "Communities", "Modularity"))
cat(rep("-", 46), "\n", sep="")
cat(sprintf("%-20s %12d %12.4f\n", "Louvain", length(communities), modularity(communities)))
cat(sprintf("%-20s %12d %12.4f\n", "Walktrap (4 steps)", length(walktrap_4), modularity(walktrap_4)))
cat(sprintf("%-20s %12d %12.4f\n", "Walktrap (5 steps)", length(walktrap_5), modularity(walktrap_5)))
cat(sprintf("%-20s %12d %12.4f\n", "Edge Betweenness", length(edge_betweenness_comm), modularity(edge_betweenness_comm)))

# Clustering Coefficients
global_cc <- transitivity(G, type = "global")
local_cc <- transitivity(G, type = "local")
avg_local_cc <- mean(local_cc, na.rm = TRUE)

cat("Global Clustering Coefficient:", round(global_cc, 4), "\n")
cat("Average Local Clustering Coefficient:", round(avg_local_cc, 4), "\n")

# Assortativity
# OUTPUT: Range -1 to +1; positive = popular connect with popular; negative = hub-spoke
degree_assortativity <- assortativity_degree(G)
cat("=== Assortativity Analysis ===\n")
cat("Degree Assortativity Coefficient:", round(degree_assortativity, 4), "\n\n")

if (degree_assortativity > 0.1) {
  cat("Interpretation: ASSORTATIVE - high-degree nodes connect with high-degree nodes\n")
} else if (degree_assortativity < -0.1) {
  cat("Interpretation: DISASSORTATIVE - high-degree nodes connect with low-degree nodes\n")
} else {
  cat("Interpretation: NEUTRAL - no degree preference\n")
}
# Path Analysis
diameter_val <- diameter(G)
avg_path <- mean_distance(G)

cat("Diameter (longest shortest path):", diameter_val, "\n")
cat("Average Path Length:", round(avg_path, 4), "\n")
# Summary Statistics
cat("=" %>% rep(50) %>% paste(collapse=""), "\n")
cat("JAZZ MUSICIANS NETWORK SUMMARY\n")
cat("=" %>% rep(50) %>% paste(collapse=""), "\n")
cat("Nodes (Musicians):", vcount(G), "\n")
cat("Edges (Collaborations):", ecount(G), "\n")
cat("Density:", round(edge_density(G), 4), "\n")
cat("Average Degree:", round(mean(degree(G)), 2), "\n")
cat("Diameter:", diameter_val, "\n")
cat("Average Path Length:", round(avg_path, 4), "\n")
cat("Global Clustering Coefficient:", round(global_cc, 4), "\n")
cat("Number of Communities:", length(communities), "\n")
cat("Modularity:", round(modularity(communities), 4), "\n")
cat("=" %>% rep(50) %>% paste(collapse=""), "\n")

