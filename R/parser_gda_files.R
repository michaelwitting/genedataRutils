#' @title Importing Genedata Expressionist for MS .gda files
#'
#' @description
#'
#' `readGda` Genedata Expressionist for MS can export obtained data in multiple file formats.
#' One is the .gda format, which can be read by the Analyst module. This file
#' contains the actual data (e.g. peak intensities, areas etc...) and additional
#' annotations.
#' This functions reads directly a .gda file and returns a list with three 
#' distinct data frames. The first element contains the actual MS data, the 
#' second contains annotation of the rows and the thrid the annotation of the
#' columns.
#'
#' @param file `character` Path to .gda file that shall be read.
#'
#' @return A list of three data frames, the first contains the actual MS data,
#'    the second the row annotation and the third the column annotations.
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
#' file <- system.file("extdata/20191108_Pesticides_test_Cluster.gda", package = "genedataRutils")
#' pest_cluster <- readGda(file)
#' 
readGda <- function(file) {
  
  # create data frames for data
  row_anno <- data.frame(stringsAsFactors = FALSE)
  ms_data <- data.frame()
  clusterNames <- NULL
  
  # read file line by line
  con <- file(file, "r")
  
  while (TRUE) {
    line <- readLines(con, n = 1)
    
    # exit loop
    if (length(line) == 0) {
      break
    }
    
    # check for different line content
    if (grepl("^Name", line)) {
      
      # get content splitted
      content <- strsplit(line, "\t")[[1]]
      
      # isolate sample names and add to ColumnDf
      samples_names <- as.character(content[2:(length(content) - noOfRowAnno)])
      col_anno <- data.frame(File = character(length(content) - noOfRowAnno - 1))
      col_anno["File"] <- samples_names
      
      # isolate row annotation names
      row_anno_names <- as.character(content[(length(content) - noOfRowAnno + 1):length(content)])
    } else if (grepl("^# Row Annotations:", line)) {
      
      # split and get number of row annotations
      noOfRowAnno <- as.numeric(regmatches(line, regexpr("\\d+", line)))
      
    } else if (grepl("\\[[A-Z]\\]", line)) {
      
      # get content splitted
      content <- strsplit(line, "\t")[[1]]
      
      # get new annotation
      colAnnotationName <- content[[1]][1]
      colAnnotation <- content[2:(length(samples_names) + 1)]
      
      # add to annotation data frame
      col_anno[paste(colAnnotationName)] <- colAnnotation
      
    } else if (grepl("^(Peak|Cluster)", line)) {
      
      # get content splitted
      content <- strsplit(line, "\t")[[1]]
      
      # get cluster names
      clusterNames <- c(clusterNames, content[1])
      
      # get peak values
      peakValues <- as.numeric(content[2:(length(content) - noOfRowAnno)])
      rowAnnotations <- content[(length(content) - noOfRowAnno + 1):length(content)]
      
      # add to row_anno and ms_data
      row_anno <- rbind.data.frame(row_anno,
                                    rowAnnotations,
                                    stringsAsFactors = FALSE
      )
      
      ms_data <- rbind.data.frame(
        ms_data,
        peakValues
      )
    }
  }
  
  # close connection to file
  close(con)
  
  # adjust column headers
  colnames(ms_data) <- samples_names
  colnames(row_anno) <- row_anno_names
  
  # adjust row names
  row.names(ms_data) <- clusterNames
  row.names(row_anno) <- clusterNames
  
  # adjust data type
  for (name in colnames(col_anno)) {
    if (grepl("\\[N|n\\]", name)) {
      
      col_anno[name] <- as.numeric(unlist(col_anno[name]))
      
    } else if (grepl("\\[C|c\\]", name)) {
      
      col_anno[name] <- as.factor(unlist(col_anno[name]))
      
    }
  }
  
  row_anno <- type.convert(row_anno)
  col_anno <- type.convert(col_anno)
  ms_data <- type.convert(ms_data)
  
  # return
  list(ms_data, row_anno, col_anno)
}

#' @title Get feature data from read .gda file
#'
#' @description
#'
#' `getMsData` returns the actual table with features from the read .gda file
#'
#' @param x `list` List with data read from .gda file.
#'
#' @return A data frame with the feature data
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
#' file <- system.file("extdata/20191108_Pesticides_test_Cluster.gda", package = "genedataRutils")
#' pest_cluster <- readGda(file)
#' getMsData(pest_cluster)
#'
getMsData <- function(x) {
  x[[1]]
}

#' @title Get row annotation data from read .gda file
#'
#' @description
#'
#' `getRowData` returns the actual table with features from the read .gda file
#'
#' @param x `list` List with data read from .gda file.
#'
#' @return A data frame wit row annotations.
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
#' file <- system.file("extdata/20191108_Pesticides_test_Cluster.gda", package = "genedataRutils")
#' pest_cluster <- readGda(file)
#' getRowAnno(pest_cluster)
#'
getRowAnno <- function(x) {
  x[[2]]
}

#' @title Get column annotation data from read .gda file
#'
#' @description
#'
#' `getColData` returns the actual table with features from the read .gda file
#'
#' @param x `list` List with data read from .gda file.
#'
#' @return A data frame with the column annotations.
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
#' file <- system.file("extdata/20191108_Pesticides_test_Cluster.gda", package = "genedataRutils")
#' pest_cluster <- readGda(file)
#' getColAnno(pest_cluster)
#'
getColAnno <- function(x) {
  x[[3]]
}
