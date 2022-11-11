library(dplyr)
library(tidyr)
library(ggplot2)

source("C:/Users/marce/Documents/Repos/microbiome-help/diversity_data_helper_functions.R")

otu_table_path <- "D:/hmp/qiime_analysis/7_table.from_biom_w_taxonomy_strain_level.txt"

greengenes_parser_strain_level <- function(string){
  result_string <- ""
  pieces <- strsplit(string, split  = ";")[[1]]
  len_pieces <- length(pieces)
  
  pieces_counter <- len_pieces
  # while result string is empty.
  while(nchar(result_string) == 0){
    # if we arrived to the end of the categories and not reached taxonomy return the original string, e.g. "Undetermined". To make sure loop stops.
    if(pieces_counter < 2){
      result_string <- string
    }else{
      # lets analyze the last two pieces of the string
      last_piece <- strsplit(pieces[pieces_counter], split  = "__")[[1]]
      ap_piece <- strsplit(pieces[pieces_counter - 1], split  = "__")[[1]]
      if (length(last_piece) == 2 && (last_piece[1] == " n" )) {
        genus_piece <- strsplit(pieces[pieces_counter - 2], split  = "__")[[1]]
        result_string <- paste(genus_piece[2], ap_piece[2], last_piece[2])
      }
      # if the last piece has length two (means is has a name on it), "s" means species. Ideal taxonomy resolution case.
      else if(length(last_piece) == 2 && (last_piece[1] == " s" || last_piece[1] == "s")){
        result_string <- paste(ap_piece[2], last_piece[2])
      }else if(length(last_piece) == 2 && last_piece[1] != " s"){ #if last piece has a name but it is not species, add sp
        result_string <- paste(last_piece[2], "sp")
      }else if(length(ap_piece) == 2){ #if antepenultimate piece has a name, add sp
        result_string = paste(ap_piece[2], "sp")
      }
    }
    # we go to the next taxonomy level if we did not find a resolved taxonomy in these levels.
    pieces_counter <- pieces_counter - 1
  }
  return(result_string)
}

read_qiime_otu_table2 <- function(table_path){
  if (!"readr" %in% installed.packages()) install.packages("readr")
  if (!"collections" %in% installed.packages()) install.packages("collections")
  if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
  if (!"tibble" %in% installed.packages()) install.packages("tibble")
  
  # read otu_table "as is"
  otu_table <- readr::read_delim(table_path, skip = 1 ,delim = "\t")
  
  # getting a vector of parsed taxonomy
  
  tax_col <- apply(otu_table["taxonomy"], 1, greengenes_parser_strain_level)
  
  # renaming species in taxonomy.
  # 
  
  # assigning taxonomy to column parsed_taxonomy
  otu_table["taxonomy"] <- tax_col
  
  # moving tax column to the first column
  otu_table <- cbind(otu_table[, ncol(otu_table)], otu_table[1:nrow(otu_table), 2:(ncol(otu_table)-1)])
  # renaming tax to taxonomy. rename() is a dplyr function.
  #otu_table <- dplyr::rename(otu_table, taxonomy = parsed_taxonomy)
  
  otu_table <- otu_table %>%
    group_by(taxonomy) %>%
    summarise_all(sum)
  # setting row names and dropping rownames column
  otu_table <- tibble::column_to_rownames(otu_table, var = "taxonomy")
  return(otu_table)
}

otu_table <- read_qiime_otu_table2(otu_table_path)

otu_table_ordered_sums <- otu_table[order(rowSums(otu_table), decreasing = TRUE),]

otu_table_ordered_means <- otu_table[order(rowMeans(otu_table), decreasing = TRUE),]

row_names_df_to_remove<-c("k__Bacteria","Unassigned")

otu_table_ordered_means <- otu_table_ordered_means[!(row.names(otu_table_ordered_means) %in% row_names_df_to_remove),]

otu_table_ordered_means <- otu_table_ordered_means %>% dplyr::rename_all(make.names)



otu_table2 <- otu_table_ordered_means[1:50,]

otu_table2["bacteria"] <- row.names(otu_table2)

otu_g <- gather(otu_table2, X5a950f27980b5d93e4c16da1244ee6c4:d57eb430d669de8329be1769d4e8d74a , key = "sample", value = "counts")

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
  scale_fill_manual(values=cbbPalette)

write.table(otu_table_ordered, "C:/Users/marce/Desktop/otu_table.csv", sep = ",", col.names = FALSE, quote = FALSE)

write.table(otu_table_ordered, "C:/Users/marce/Desktop/otu_table_100.csv", sep = ",", col.names = FALSE, quote = FALSE)

####################################################################################

otu_table3 <- otu_table_ordered_means

otu_table3$Mean<- rowMeans(otu_table3)

otu_table3 <- select(otu_table3, Mean)

otu_table3["bacteria"] <- row.names(otu_table3)

otu_table3["sample"] <- "mean"

otu_table3<- otu_table3[!(row.names(otu_table3) %in% row_names_df_to_remove),]

ggplot(otu_table3[1:50,], aes(x=sample, y=Mean, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)



# Species co-occurrence analyses
####################################################################################

library(cooccur)

# Transforming abundance data to presence/abscence
otu_table_pa <- vegan::decostand(otu_table_ordered[1:100,], method = "pa")

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
heatmap(otu_table_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"))


row_scaled_df <- scale(t(otu_table_ordered_means[1:30,]))

heatmap(row_scaled_df, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"))
