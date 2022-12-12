library(dplyr)
library(tidyr)
library(ggplot2)

source("C:/Users/Marcelo/Documents/Github/microbiome-help/diversity_data_helper_functions.R")

setwd("C:/Users/Marcelo/Desktop/HMP_nose_data_analysis/qiime_analyses/qiime_analyses_asv")

# ASVs
otu_table_path <- "./3_resultados/7_table.from_biom_w_taxonomy_strain.txt"

biom_path <- "./3_resultados/6_table_w_tax_strain.biom"

load_test <- load_biom(biom_path, tax_rank = "Strain")
load_test2 <- load_biom(biom_path, tax_rank = "Family")

# Read table
otu_table <- read_qiime_otu_table2(otu_table_path)

# Order table by larger to lower mean abundance of bacteria (rows)
otu_table_ordered_means <- otu_table[order(rowMeans(otu_table), decreasing = TRUE),]

# Remove unassigned counts
row_names_df_to_remove<-c("k__Bacteria","Unassigned")
otu_table_ordered_means <- otu_table_ordered_means[!(row.names(otu_table_ordered_means) %in% row_names_df_to_remove),]

# Remove ASVs without NCBI refseq confident assignation.
row_names_df_to_remove2<-c("Neisseriaceae sp","Streptophyta sp")

row_names_df_to_remove2<-c("Neisseriaceae sp","Streptophyta sp", "Corynebacterium sp", "Streptococcus sp", "Bacilli sp", "Rothia mucilaginosa", "Staphylococcus epidermidis", "Anaerococcus sp", "Lachnospiraceae sp", "Actinomyces sp")

otu_table_ordered_means <- otu_table_ordered_means[!(row.names(otu_table_ordered_means) %in% row_names_df_to_remove2),]

# Fix names of samples that do NOT begin with a letter.
otu_table_ordered_means <- otu_table_ordered_means %>% dplyr::rename_all(make.names)

# Select only the 30 more abundant species.
otu_table2 <- otu_table_ordered_means[1:30,]

# Generate a column with the names of ASVs/OTUs using rownames.
otu_table2["bacteria"] <- row.names(otu_table2)

otu_g <- gather(otu_table2, X5a950f27980b5d93e4c16da1244c9c15:d57eb430d669de8329be1769d4e8f962 , key = "sample", value = "counts")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2","brown1", "#CC79A7", "olivedrab3", "rosybrown",
                "darkorange3", "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
                "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
                "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
                "#e3ae78", "#a23f3f", "#290f76", "#ce7e00", "#386857",
                "#738564", "#e89d56", "#cd541d", "#1a3a46", "#ffe599",
                "#583E26", "#A78B71", "#F7C815", "#EC9704", "#9C4A1A",
                "firebrick2", "#C8D2D1", "#14471E", "#EE9B01", "#DA6A00",
                "#4B1E19", "#C0587E", "#FC8B5E", "#EA592A", "#FEF4C0")

ggplot(otu_g, aes(x=sample, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

write.table(otu_table_ordered, "C:/Users/marce/Desktop/otu_table.csv", sep = ",", col.names = FALSE, quote = FALSE)

#write.table(otu_table_ordered, "C:/Users/marce/Desktop/otu_table_100.csv", sep = ",", col.names = FALSE, quote = FALSE)

####################################################################################
# Graph of means for each ASV/OTU

otu_table3 <- otu_table_ordered_means

# Calculate means for each ASV/OTU
otu_table3$Mean<- rowMeans(otu_table3)

# Reduce table only to Means, discard all sample values
otu_table3 <- select(otu_table3, Mean)

# Generate a column with the names of ASVs/OTUs using rownames.
otu_table3["bacteria"] <- row.names(otu_table3)

otu_table3["sample"] <- "mean"

#otu_table3<- otu_table3[!(row.names(otu_table3) %in% row_names_df_to_remove),]

ggplot(otu_table3[1:30,], aes(x=sample, y=Mean, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)



# Species co-occurrence analyses
####################################################################################

library(cooccur)

# Transforming abundance data to presence/abscence
otu_table_pa <- vegan::decostand(otu_table_ordered_means[1:30,], method = "pa")

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

write.csv(edges, "C:/Users/marce/Desktop/coocur_network.csv")

library(igraph)

g <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
print(g, e=TRUE, v=TRUE)
plot(g)
write.graph(g, "C:/Users/marce/Desktop/coocur_network.graphml", format = "graphml")

####################################################################################

row_names_df_to_remove<-c("k__Bacteria","Unassigned")
otu_table_ordered_means<- otu_table_ordered_means[!(row.names(otu_table_ordered_means) %in% row_names_df_to_remove),]

otu_table_100 <- otu_table_ordered_means[1:100,]

otu_table_100 <- otu_table_100 %>% dplyr::rename_all(make.names)

write.table(otu_table_100, "C:/Users/marce/Desktop/otu_table.csv", sep = ",", col.names = FALSE, quote = FALSE)

####################################################################################

#scaled by column
otu_table_scaled <- scale(otu_table_ordered_means[1:30,])

# Graph heatmap
heatmap(otu_table_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")

species_scaled_df <- scale(t(otu_table_ordered_means[1:30,]))
heatmap(t(species_scaled_df), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale ="none")

heatmap(data.matrix(otu_table2), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale = "column")
