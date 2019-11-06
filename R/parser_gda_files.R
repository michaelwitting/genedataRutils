# this function reads a .gda file and creates three data frames
read_gda_file <- function(path_to_file, prefix = "") {
  
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
      print(name)
      colAnnoDf[name] <- as.numeric(unlist(colAnnoDf[name]))
    } else if (grepl("\\[C|c\\]", name)) {
      colAnnoDf[name] <- as.factor(unlist(colAnnoDf[name]))
    }
  }
  
  # correct the format of the columns
  rowAnnoDf$m.z <- as.numeric(rowAnnoDf$`m/z`)
  rowAnnoDf$RT <- as.numeric(rowAnnoDf$RT)
  rowAnnoDf$`RT Min` <- as.numeric(rowAnnoDf$`RT Min`)
  rowAnnoDf$`RT Max` <- as.numeric(rowAnnoDf$`RT Max`)
  rowAnnoDf$`m/z Min` <- as.numeric(rowAnnoDf$`m/z Min`)
  rowAnnoDf$`m/z Max` <- as.numeric(rowAnnoDf$`m/z Max`)
  
  
  # assign data frames
  assign(paste0(prefix, "ms_data"), dataDf, envir = parent.frame())
  assign(paste0(prefix, "col_anno"), colAnnoDf, envir = parent.frame())
  assign(paste0(prefix, "row_anno"), rowAnnoDf, envir = parent.frame())
}