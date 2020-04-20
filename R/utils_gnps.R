#' @title Importing Genedata Expressionist for MS .gda files
#'
#' @description
#'
#' `createGnpsFiles` Genedata Expressionist for MS can export obtained data in multiple file formats.
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
createGnpsFiles <- function(x, spectra, path = "") {
  
  ms_data <- getMsData(x)
  row_anno <- getRowAnno(x)

  # sanity checks
  # does the MS2 have correct columns
  # if not perform add ID first
  if(!"CLUSTER_ID" %in% colnames(spectra_id@elementMetadata)) {
    stop("No column CLUSTER_ID found in spectra. Perform ms2_add_id() first")
  }
  
  # create feature table
  # create new DF with feature table in XCMS3 format
  feature_table <- data.frame(Row.names = row.names(row_anno),
                              mzmed = row_anno$`m/z`,
                              mzmin = row_anno$`m/z Min`,
                              mzmax = row_anno$`m/z Max`,
                              rtmed = row_anno$RT * 60,
                              rtmin = row_anno$`RT Min` * 60,
                              rtmax = row_anno$`RT Max` * 60,
                              npeaks = ncol(ms_data),
                              sample = ncol(ms_data))
  
  # rename features according to XCMS and add intensities
  feature_table$Row.names <- gsub("Cluster_", "FT", feature_table$Row.names)
  feature_table <- cbind.data.frame(feature_table, ms_data)
  
  # convert names in MS2 metadata to XCMS3 format
  # create fields required
  mcols(spectra)$SCANS <- as.numeric(str_extract(mcols(spectra)$FEATURE_ID, "\\d+"))
  mcols(spectra)$FEATURE_ID <- gsub("Cluster_", "FT", mcols(spectra)$FEATURE_ID)
  mcols(spectra)$PEAK_ID <- mcols(spectra)$FEATURE_ID
  mcols(spectra)$COMPOUND <- mcols(spectra)$FEATURE_ID
  
  # filter spectra
  spectra_filtered <- spectra[which(peaksCount(spectra) > 0)]
  spectra_filtered <- spectra_filtered[which(!is.na(mcols(spectra_filtered)$FEATURE_ID))]
  
  # adjust scan index and acquisition number
  for(i in 1:length(spectra_filtered)) {
    
    print(paste0(i, " of ", length(spectra_filtered), " spectra processed"))
    
    spectra_filtered[[i]]@scanIndex <- as.integer(mcols(spectra_filtered[i])$SCANS)
    spectra_filtered[[i]]@acquisitionNum <- as.integer(mcols(spectra_filtered[i])$SCANS)
    
  }
  
  # drop scans column
  mcols(spectra_filtered)$SCANS <- NULL
  
  # write spectra and MS1 data
  writeMgfData(spectra_filtered,
               paste0(path, "/spectra_gnps.mgf"))
  
  write.table(feature_table,
              paste0(path, "/ms1_feature_table_gnps.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  
}

# ==============================================================================
# DEPRECATED FUNCTIONS
# ==============================================================================
#' @title Create input for GNPS feature based molecular networking
#'
#' This functions takes data from MS1 cluster as well as MS2 spectra that 
#' contain CLUSTER_ID in the metadata. Two output files are generated. The first
#'  one a .mgf file with all the metadata required for linking the MS1 and MS2 
#'  data. The second file contains the MS1 table in a format similar to XCMS3.
#'  Both the .mgf file and the feature table can be used with the XCMS3 branch 
#'  of the GNPS feature based molecular networking
#'
#' @param ms_data Data frame with the intensity or peak area values
#' @param row_anno Data frame with the row annotations
#' @param spectra Spectra object with the MS2 spectra, CLUSTER_ID column in metadata is required
#' @param path_to_save Path where the results files shall be saved
#' 
#' @export
create_gnps_files <- function(ms_data, row_anno, spectra, path_to_save) {
  
  .Deprecated("createGnpsFiles")
  
  require(MSnbase)
  require(tidyverse)
  
  # sanity checks
  # does the MS2 have correct columns
  # if not perform add ID first
  if(!"CLUSTER_ID" %in% colnames(spectra_id@elementMetadata)) {
    stop("No column CLUSTER_ID found in spectra. Perform ms2_add_id() first")
  }
  
  # create feature table
  # create new DF with feature table in XCMS3 format
  feature_table <- data.frame(Row.names = row.names(row_anno),
                              mzmed = row_anno$`m/z`,
                              mzmin = row_anno$`m/z Min`,
                              mzmax = row_anno$`m/z Max`,
                              rtmed = row_anno$RT * 60,
                              rtmin = row_anno$`RT Min` * 60,
                              rtmax = row_anno$`RT Max` * 60,
                              npeaks = ncol(ms_data),
                              sample = ncol(ms_data))
  
  # rename features according to XCMS and add intensities
  feature_table$Row.names <- str_replace(feature_table$Row.names, "Cluster_", "FT")
  feature_table <- bind_cols(feature_table, ms_data)
  
  # convert names in MS2 metadata to XCMS3 format
  # create fields required
  mcols(spectra)$SCANS <- as.numeric(str_extract(mcols(spectra)$FEATURE_ID, "\\d+"))
  mcols(spectra)$FEATURE_ID <- str_replace(mcols(spectra)$FEATURE_ID, "Cluster_", "FT")
  mcols(spectra)$PEAK_ID <- mcols(spectra)$FEATURE_ID
  mcols(spectra)$COMPOUND <- mcols(spectra)$FEATURE_ID
  
  # filter spectra
  spectra_filtered <- spectra[which(peaksCount(spectra) > 0)]
  spectra_filtered <- spectra_filtered[which(!is.na(mcols(spectra_filtered)$FEATURE_ID))]
  
  # adjust scan index and acquisition number
  for(i in 1:length(spectra_filtered)) {

    print(paste0(i, " of ", length(spectra_filtered), " spectra processed"))
    
    spectra_filtered[[i]]@scanIndex <- as.integer(mcols(spectra_filtered[i])$SCANS)
    spectra_filtered[[i]]@acquisitionNum <- as.integer(mcols(spectra_filtered[i])$SCANS)
    
  }
  
  # drop scans column
  mcols(spectra_filtered)$SCANS <- NULL
  
  # write spectra and MS1 data
  writeMgfData(spectra_filtered, paste0(path_to_save, "/spectra_gnps.mgf"))
  write_tsv(feature_table, paste0(path_to_save, "/ms1_feature_table_gnps.tsv"))
  
}