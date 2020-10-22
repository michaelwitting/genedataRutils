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
#' @importFrom utils type.convert
#'
#' @export
#'
#' @examples
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
    if (grepl("^(Name|Column)", line)) {
      
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
  
  row_anno <- type.convert(row_anno, as.is = TRUE)
  col_anno <- type.convert(col_anno, as.is = TRUE)
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
getMsData <- function(x) {
  x[[1]]
}

#' @title Sets feature data from read .gda file
#'
#' @description
#'
#' `setMsData` Sets the MS data, data is checked if row names match
#'
#' @param x `list` List with data read from .gda file.
#' @param msData `data.frame` Replacement for MS data in x, rownames must match
#'
#' @return A list with MS data, row annotations, column annotations
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
setMsData <- function(x, msData) {
  
  # isolate data
  row_anno <- getRowAnno(x)
  col_anno <- getColAnno(x)
  
  # check that sample names are the same
  if(!identical(sort(rownames(row_anno)),
                sort(rownames(msData)))) {
    
    stop("Row names are not matching!")
    
  }
  
  # check that sample names are the same
  if(!identical(sort(col_anno$File),
                sort(colnames(msData)))) {
    
    stop("Column names are not matching!")
    
  }
  
  x[[1]] <- msData
  
  x
  
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
getRowAnno <- function(x) {
  x[[2]]
}

#' @title Sets row annotation data from read .gda file
#'
#' @description
#'
#' `setRowAnno` Sets the row annotation data, data is checked if row names match
#'
#' @param x `list` List with data read from .gda file.
#' @param rowAnno `data.frame` Replacement for row annotation data in x,
#'      rownames must match
#'
#' @return A list with MS data, row annotations, column annotations
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
setRowAnno <- function(x, rowAnno) {
  
  # isolate data
  row_anno <- getRowAnno(x)
  
  # check that sample names are the same
  if(!identical(sort(rownames(row_anno)),
                sort(rownames(rowAnno)))) {
    
    stop("Row names are not matching!")
    
  }
  
  x[[2]] <- rowAnno
  
  x
  
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
getColAnno <- function(x) {
  x[[3]]
}

#' @title Sets feature data from read .gda file
#'
#' @description
#'
#' `setColAnno` Sets the column annotation data, data is checked if row names match
#'
#' @param x `list` List with data read from .gda file.
#' @param colAnno `data.frame` Replacement for column annotation data in x,
#'      rownames must match
#'
#' @return A list with MS data, row annotations, column annotations
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
setColAnno <- function(x, colAnno) {
  
  # isolate data
  col_anno <- getColAnno(x)

  # check that sample names are the same
  if(!identical(sort(col_anno$File),
                sort(colnames(colAnno)))) {
    
    stop("Column names are not matching!")
    
  }
  
  x[[3]] <- colAnno
  
  x
  
}

#' @title Write .gda file
#'
#' @description
#'
#' `writeGda` Writes a list of row and column annotation and ms data to a file
#' 
#' @param x `list` List with data read from .gda file. 
#' @param file `character` Name of file to write to
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
writeGda <- function(x, file = "export.gda") {
  
  # sanity checks
  if(!length(x) == 3) {
    
    stop("Input is not of length 3. Sure it contains data from a .gda file?")
    
  }
  
  # isolate data and sort
  ms_data <- getMsData(x)
  row_anno <- getRowAnno(x)
  col_anno <- getColAnno(x)
  
  
  ms_data <- ms_data[,order(names(ms_data))]
  col_anno <- col_anno[order(col_anno$File),]
  
  # check names are in same order
  if(!all(names(ms_data) == col_anno$File)) {
    
    stop("Name not in same order")
    
  }
  
  ms_data_full <- merge(ms_data, row_anno, by = "row.names")
  
  if(file.exists(file)) {
    
    file.remove(file)
    
  }
  
  # helper function for writing file
  con <- file(description = file, open = "at")
  
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  #write meta data block
  .cat("# Namespace: RExport\n")
  .cat("# Rowtype: Metabolites\n")
  .cat("# Observable: Max. Intensity\n")
  .cat("# Row Annotations: ", ncol(row_anno), "\n")
  .cat("# Column Annotations: ", ncol(col_anno) - 1, "\n")
  .cat("# Transformation: LOG\n")
  .cat("# Version: 1\n")
  .cat("# Author: RExport\n")
  
  # write header
  header <- colnames(ms_data_full)
  header[1] <- "Name"
  
  .cat(paste0(header, collapse = "\t"))
  .cat("\n")
  
  # write column annotation
  col_anno_names <- colnames(col_anno)
  col_anno_names <- col_anno_names[col_anno_names != "File"]
  
  for(col_anno_name in col_anno_names) {
    
    .cat(paste0(c(col_anno_name, as.vector(col_anno[,col_anno_name])), collapse = "\t"))
    .cat("\n")
    
  }
  
  close(con)
  
  write.table(ms_data_full, file = file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}


#' @title Remove specific samples from data
#'
#' @description
#'
#' `removeSamples` Removes specific samples from the data
#'
#' @param x `list` List with data read from .gda file.#' 
#' @param samples `character` Vector with the sample names to be removed
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
removeSamples <- function(x, samples) {
  
  # sanity checks
  if(!length(x) == 3) {
    
    stop("Input is not of length 3. Sure it contains data from a .gda file?")
    
  }
  
  # isolate data and sort
  ms_data <- getMsData(x)
  row_anno <- getRowAnno(x)
  col_anno <- getColAnno(x)
  
  # test if samples are contained in x
  if(!samples %in% colnames(ms_data)) {
    
    stop("sample not contained in data set")
    
  }
  
  ms_data <- ms_data[,!colnames(ms_data) %in% samples]
  col_anno <- col_anno[which(!col_anno$File %in% samples),]
  
  list(ms_data, row_anno, col_anno)
  
}

#' @title Remove specific cluster from data
#'
#' @description
#'
#' `removeCluster` Removes specific clusterfrom the data
#'
#' @param x `list` List with data read from .gda file.
#' 
#' @param cluster `character` Vector with the cluster to be removed
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
removeCluster <- function(x, cluster) {
  
  # sanity checks
  if(!length(x) == 3) {
    
    stop("Input is not of length 3. Sure it contains data from a .gda file?")
    
  }
  
  # isolate data and sort
  ms_data <- getMsData(x)
  row_anno <- getRowAnno(x)
  col_anno <- getColAnno(x)
  
  # test if samples are contained in x
  if(!cluster %in% row.names(ms_data)) {
    
    stop("sample not contained in data set")
    
  }
  
  ms_data <- ms_data[!(row.names(ms_data) %in% cluster),]
  row_anno <- row_anno[!(row.names(row_anno) %in% cluster),]
  
  list(ms_data, row_anno, col_anno)
  
}
