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

otu_table <- read_qiime_otu_table(otu_table_path, level = "Strain")


