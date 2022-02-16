#import(ape)

parse_input_file <- function(file_path, genetrees = TRUE){

  ### Reads input file, containing
  ### gene trees in Newick format and
  ### their frequencues

  if (genetrees){

    iftree = FALSE
    tree_list <- vector(mode = "list")
    freqs_vector <- vector()
    input_file = file(file_path, "r")

    i <- 1
    j <- 1

    while (TRUE) {

      line = readLines(input_file, n = 1)
      if (length(line) == 0){ #if no more lines, break out of loop
        break
      }

      if (startsWith(line, "(")){ #line containing a newick tree
        tree <- ape::read.tree(text = line)
        tree_list[[i]] <- tree
        i <- i + 1
      }

      if (startsWith(line, "0") | startsWith(line, ".")){ #line containing tree frequency
        freq <- as.double(line)
        freqs_vector[j] <- freq
        j <- j + 1
      }

      tree_list[[i]] <- freqs_vector #add vector of tree freqs to end of list

    }

    close(input_file)

    return(tree_list)

  } else{

    input_file = file(file_path, "r")

    text_line <- readLines(input_file, n = 1)
    tree <- ape::read.tree(text = text_line)

    close(input_file)
    return(tree)
  }

}
