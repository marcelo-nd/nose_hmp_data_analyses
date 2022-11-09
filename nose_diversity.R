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

otu_table_ordered <- otu_table[order(rowSums(otu_table), decreasing = TRUE),]


otu_table2 <- otu_table_ordered[,1:50]

otu_table2 <- otu_table2 %>% dplyr::rename_all(make.names)

otu_table2["bacteria"] <- row.names(otu_table_ordered)

otu_g <- gather(otu_table2, X5a950f27980b5d93e4c16da1244ee6c4:d57eb430d669de8329be1769d4dd9678 , key = "sample", value = "counts")

ggplot(otu_g, aes(x=sample, y=counts, fill=bacteria)) + 
  geom_bar(position="stack", stat="identity")

# Species co-occurrence analyses
####################################################################################

library(cooccur)

# Transforming abundance data to presence/abscence
otu_table_pa <- vegan::decostand(otu_table_ordered, method = "pa")

# Infering co-ocurrences
cooccur.otus <- cooccur(otu_table_pa,
                        type = "spp_site",
                        spp_names = TRUE)

summary(cooccur.otus)
plot(cooccur.otus)