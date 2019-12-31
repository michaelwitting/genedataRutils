# this function reads a .gda file and creates three data frames
#' @title Importing Genedata Expressionist for MS .gda files
#' 
#' Genedata Expressionist for MS can export obtained data in multiple file formats.
#' One is the .gda format, which can be read by the Analyst module. This file
#' contains the actual data (e.g. peak intensities, areas etc...) and additional
#' annotations.
#' This functions read directly a .gda file and creates three distinct data frames.
#' One contains the actual MS data (ms_data), the other two the annotations of the columns
#' and rows (row_anno and col_anno).
#' 
#' @param path_to_file File path to .gda file
#' @param prefix Name prefix for the generated data frames. Should be used if multiple .gda files
#' shall be read. The prefix is added before each data frame name, e.g. prefix = "cluster" generates
#' cluster_ms_data as name of data frame variable instead of ms_data. Data frame are automatically assigned
#' to the current working environment.
#' 
#' @examples 
#' file <- system.file('')
#' read_gda_file(file, prefix = 'pesticide')
#'
#' @import stringr
#'
#' @export
read_gda_file <- function(path_to_file, prefix = "") {
  
  require(stringr)
  
  # create data frames for data
  rowAnnoDf <- data.frame(stringsAsFactors = FALSE)
  dataDf <- data.frame()
  clusterNames <- NULL
  
  # read file line by line
  con <- file(path_to_file, "r")
  
  while (TRUE) {
    line <- readLines(con, n = 1)

    # exit loop
    if (length(line) == 0) {
      break
    }
    
    # check for different line content
    if (grepl("^Name", line)) {
      
      # get content splitted
      content <- stringr::str_split(line, "\t")[[1]]
      
      # isolate sample names and add to ColumnDf
      SampleNames <- as.character(content[2:(length(content) - noOfRowAnno)])
      colAnnoDf <- data.frame(SampleNames = character(length(content) - noOfRowAnno - 1))
      colAnnoDf["SampleNames"] <- SampleNames
      
      # isolate row annotation names
      rowAnnotationNames <- as.character(content[(length(content) - noOfRowAnno + 1):length(content)])
    } else if (grepl("^# Row Annotations:", line)) {
      
      # split and get number of row annotations
      noOfRowAnno <- as.numeric(stringr::str_extract(stringr::str_split(line, ":")[[1]][2], "\\d+"))
      print(as.numeric(stringr::str_extract(stringr::str_split(line, ":")[[1]][2], "\\d+")))
      
      # print(line)
    } else if (grepl("\\[[A-Z]\\]", line)) {
      
      # get content splitted
      content <- stringr::str_split(line, "\t")[[1]]
      
      # get new annotation
      colAnnotationName <- content[[1]][1]
      colAnnotation <- content[2:(length(SampleNames) + 1)]
      
      # add to annotation data frame
      colAnnoDf[paste(colAnnotationName)] <- colAnnotation
    } else if (grepl("^(Peak|Cluster)", line)) {
      
      # get content splitted
      content <- stringr::str_split(line, "\t")[[1]]
      
      # get cluster names
      clusterNames <- c(clusterNames, content[1])
      
      # get peak values
      peakValues <- as.numeric(content[2:(length(content) - noOfRowAnno)])
      rowAnnotations <- content[(length(content) - noOfRowAnno + 1):length(content)]
      
      # add to rowAnnoDf and dataDf
      rowAnnoDf <- rbind.data.frame(rowAnnoDf,
                                    rowAnnotations,
                                    stringsAsFactors = FALSE
      )
      
      dataDf <- rbind.data.frame(
        dataDf,
        peakValues
      )
    }
  }
  
  # close connection to file
  close(con)
  
  # adjust column headers
  colnames(dataDf) <- SampleNames
  colnames(rowAnnoDf) <- rowAnnotationNames
  
  row.names(dataDf) <- clusterNames
  row.names(rowAnnoDf) <- clusterNames
  
  # adjust data type
  for (name in colnames(colAnnoDf)) {
    if (grepl("\\[N|n\\]", name)) {
      colAnnoDf[name] <- as.numeric(unlist(colAnnoDf[name]))
    } else if (grepl("\\[C|c\\]", name)) {
      colAnnoDf[name] <- as.factor(unlist(colAnnoDf[name]))
    }
  }

  rowAnnoDf <- type.convert(rowAnnoDf)
  colAnnoDf <- type.convert(colAnnoDf)
  dataDf <- type.convert(dataDf)
  
  # assign data frames
  assign(paste0(prefix, "ms_data"), dataDf, envir = parent.frame())
  assign(paste0(prefix, "col_anno"), colAnnoDf, envir = parent.frame())
  assign(paste0(prefix, "row_anno"), rowAnnoDf, envir = parent.frame())
}
