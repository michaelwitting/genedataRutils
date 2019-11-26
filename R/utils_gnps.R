create_gnps_files <- function(ms_data, row_anno, ms2_spectra, path_to_save) {
  
  # sanity checks
  # does the MS2 have correct columns
  # if not perform add ID first
  if(!"CLUSTER_ID" %in% colnames(ms2_spectra_id@elementMetadata)) {
    stop("No column CLUSTER_ID found in ms2_spectra. Perform ms2_add_id() first")
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
  mcols(ms2_spectra)$SCANS <- as.numeric(str_extract(mcols(ms2_spectra)$FEATURE_ID, "\\d+"))
  mcols(ms2_spectra)$FEATURE_ID <- str_replace(mcols(ms2_spectra)$FEATURE_ID, "Cluster_", "FT")
  mcols(ms2_spectra)$PEAK_ID <- mcols(ms2_spectra)$FEATURE_ID
  mcols(ms2_spectra)$COMPOUND <- mcols(ms2_spectra)$FEATURE_ID
  
  # filter spectra
  ms2_spectra_filtered <- ms2_spectra[which(peaksCount(ms2_spectra) > 0)]
  ms2_spectra_filtered <- ms2_spectra_filtered[which(!is.na(mcols(ms2_spectra_filtered)$FEATURE_ID))]
  
  # adjust scan index and acquisition number
  for(i in 1:length(ms2_spectra_filtered)) {

    print(paste0(i, " of ", length(ms2_spectra_filtered), " spectra processed"))
    
    ms2_spectra_filtered[[i]]@scanIndex <- as.integer(mcols(ms2_spectra_filtered[i])$SCANS)
    ms2_spectra_filtered[[i]]@acquisitionNum <- as.integer(mcols(ms2_spectra_filtered[i])$SCANS)
    
  }
  
  # drop scans column
  mcols(ms2_spectra_filtered)$SCANS <- NULL
  
  # write spectra and MS1 data
  writeMgfData(ms2_spectra_filtered, paste0(path_to_save, "/ms2_spectra_gnps.mgf"))
  write_tsv(feature_table, paste0(path_to_save, "/ms1_feature_table_gnps.mgf"))
  
}