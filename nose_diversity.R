library(dplyr)
library(tidyr)
library(ggplot2)

source("C:/Users/Marcelo/Documents/Github/microbiome-help/diversity_data_helper_functions.R")

source("C:/Users/Marcelo/Documents/Github/microbiome-help/table_importers.R")

source("C:/Users/Marcelo/Documents/Github/microbiome-help/graphs.R")

setwd("C:/Users/Marcelo/Desktop/HMP_nose_data_analysis/qiime_analyses/qiime_analyses_asv")

# ASVs
nose_biom_path <- "./3_resultados/6_table_w_tax_strain.biom"

asv_table_nose <- load_biom(biom_path = nose_biom_path, tax_rank = "Strain", order_table = TRUE)

# Select only the 30 more abundant species.
asv_table_nose30 <- asv_table_nose[1:30,]

# Barplot
barplot_from_feature_table(feature_table = asv_table_nose30)

write.table(asv_table_nose, "./3_resultados/nose_asv_table.csv", sep = ",", col.names = FALSE, quote = FALSE)

write.table(asv_table_nose30, "./3_resultados/nose_asv_table30.csv", sep = ",", col.names = FALSE, quote = FALSE)

####################################################################################
# Graph of means for each ASV/OTU

asv_table_nose2 <- asv_table_nose30

# Calculate means for each ASV/OTU
asv_table_nose2$Mean<- rowMeans(asv_table_nose2)

# Reduce table only to Means, discard all sample values
asv_table_nose2 <- select(asv_table_nose2, Mean)

# Barplot
barplot_from_feature_table(asv_table_nose2)

# Species co-occurrence analyses
####################################################################################

library(cooccur)

# Transforming abundance data to presence/abscence
otu_table_pa <- vegan::decostand(asv_table_nose30, method = "pa")

# Infering co-ocurrences
cooccur.otus <- cooccur(otu_table_pa,
                        type = "spp_site",
                        spp_names = TRUE)

summary(cooccur.otus)
plot(cooccur.otus)


# Getting only the significant interactions
co <- print(cooccur.otus)

# Create a data frame of the nodes in the network. 
nodes <- data.frame(id = 1:nrow(otu_table_pa),
                    label = rownames(otu_table_pa),
                    color = "#606482",
                    shadow = TRUE) 

# Create an edges dataframe from the significant pairwise co-occurrences.
edges <- data.frame(from = co$sp1, to = co$sp2,
                    color = ifelse(co$p_lt <= 0.05,
                                   "#B0B2C1", "#3C3F51"),
                    dashes = ifelse(co$p_lt <= 0.05, TRUE, FALSE))

# Plotting network
library(visNetwork)
visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(layout = "layout_with_kk")

# Export networka as edges list.
write.csv(edges, "C:/Users/marce/Desktop/coocur_network.csv")

# Export graph in "graphml" format
library(igraph)

g <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
print(g, e=TRUE, v=TRUE)
plot(g)
write.graph(g, "C:/Users/marce/Desktop/coocur_network.graphml", format = "graphml")

####################################################################################
# Heatmaps

col_fun = circlize::colorRamp2(c(-0.7, 2, 5.5), c("#5F8D4E", "white", "#FF6464"))

# Scaled by colum (i.e. by sample).
asv_table30_scaled_by_sample <- scale(asv_table_nose30, center = TRUE, scale = TRUE)

# Graph simple heatmap
heatmap(asv_table30_scaled_by_sample, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")#heatmap(data.matrix(asv_proportions), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")

# complex heatmap
htmp <- ComplexHeatmap::Heatmap(asv_table30_scaled_by_sample,
                                #name = "Scaled ASV abundance",
                                show_column_names = FALSE,
                                col = col_fun,
                                row_names_gp = grid::gpar(fontsize = 8),
                                heatmap_legend_param = list(
                                  title = "Scaled abundance",
                                  labels_rot = 0,
                                  direction = "horizontal"
                                ))

ComplexHeatmap::draw(htmp, heatmap_legend_side="bottom")

##### Scaling by ASV
asv_table30_scaled_by_asv <- t(scale(t(asv_table_nose30), center = TRUE, scale = TRUE))

col_fun = circlize::colorRamp2(c(-2, 2, 4), c("#5F8D4E", "white", "#FF6464"))

htmp_by_asv <- ComplexHeatmap::Heatmap(asv_table30_scaled_by_asv,
                                #name = "Scaled ASV abundance",
                                show_column_names = FALSE,
                                col = col_fun,
                                heatmap_legend_param = list(
                                  title = "Scaled abundance",
                                  labels_rot = 0,
                                  direction = "horizontal"
                                ))
ComplexHeatmap::draw(htmp_by_asv, heatmap_legend_side="bottom")

##### Normalization by sample

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply Min-Max normalization
asv_normalized <- as.data.frame(lapply(asv_table_nose30, min_max_norm), row.names =  row.names(asv_table_nose30))

col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("#5F8D4E", "white", "#FF6464"))
htmp_norm <- ComplexHeatmap::Heatmap(as.matrix(asv_normalized),
                                 #name = "Scaled ASV abundance",
                                 show_column_names = FALSE,
                                 col = col_fun,
                                 heatmap_legend_param = list(
                                   title = "Normalized abundance",
                                   labels_rot = 0,
                                   direction = "horizontal"
                                 ))
ComplexHeatmap::draw(htmp_norm, heatmap_legend_side="bottom")

### Proportions heatmap

asv_proportions <- t(t(asv_table_nose30)/rowSums(t(asv_table_nose30)))

col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("#5F8D4E", "white", "#FF6464"))
htmp_prop <- ComplexHeatmap::Heatmap(as.matrix(asv_proportions),
                                     #name = "Scaled ASV abundance",
                                     show_column_names = FALSE,
                                     col = col_fun,
                                     heatmap_legend_param = list(
                                       title = "Normalized abundance",
                                       labels_rot = 0,
                                       direction = "horizontal"
                                     ))
ComplexHeatmap::draw(htmp_prop, heatmap_legend_side="bottom")
