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
#' @param ms2_spectra_id Spectra object with the MS2 spectra, CLUSTER_ID column in metadata is required
#' @param path_to_save Path where the results files shall be saved
#' 
#' @export
create_gnps_files <- function(ms_data, row_anno, ms2_spectra_id, path_to_save) {
  
  require(MSnbase)
  require(tidyverse)
  
  # sanity checks
  # does the MS2 have correct columns
  # if not perform add ID first
  if(!"CLUSTER_ID" %in% colnames(ms2_spectra_id_id@elementMetadata)) {
    stop("No column CLUSTER_ID found in ms2_spectra_id. Perform ms2_add_id() first")
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
  mcols(ms2_spectra_id)$SCANS <- as.numeric(str_extract(mcols(ms2_spectra_id)$FEATURE_ID, "\\d+"))
  mcols(ms2_spectra_id)$FEATURE_ID <- str_replace(mcols(ms2_spectra_id)$FEATURE_ID, "Cluster_", "FT")
  mcols(ms2_spectra_id)$PEAK_ID <- mcols(ms2_spectra_id)$FEATURE_ID
  mcols(ms2_spectra_id)$COMPOUND <- mcols(ms2_spectra_id)$FEATURE_ID
  
  # filter spectra
  ms2_spectra_id_filtered <- ms2_spectra_id[which(peaksCount(ms2_spectra_id) > 0)]
  ms2_spectra_id_filtered <- ms2_spectra_id_filtered[which(!is.na(mcols(ms2_spectra_id_filtered)$FEATURE_ID))]
  
  # adjust scan index and acquisition number
  for(i in 1:length(ms2_spectra_id_filtered)) {

    print(paste0(i, " of ", length(ms2_spectra_id_filtered), " spectra processed"))
    
    ms2_spectra_id_filtered[[i]]@scanIndex <- as.integer(mcols(ms2_spectra_id_filtered[i])$SCANS)
    ms2_spectra_id_filtered[[i]]@acquisitionNum <- as.integer(mcols(ms2_spectra_id_filtered[i])$SCANS)
    
  }
  
  # drop scans column
  mcols(ms2_spectra_id_filtered)$SCANS <- NULL
  
  # write spectra and MS1 data
  writeMgfData(ms2_spectra_id_filtered, paste0(path_to_save, "/ms2_spectra_id_gnps.mgf"))
  write_tsv(feature_table, paste0(path_to_save, "/ms1_feature_table_gnps.tsv"))
  
}