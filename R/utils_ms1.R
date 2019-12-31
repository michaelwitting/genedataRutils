# Genedata does not allow to export MS1 isotope pattern directly. They have to be reconstructed from the Peak and Cluster data
#'
#'
#' @import tidyverse
#' @export
recon_iso_pattern <- function(peak_row_anno, cluster_row_anno, peak_ms_data) {
  
  # add here function to estimate the cluster format
  cluster_ids <- row.names(cluster_row_anno)
  cluster_ids <- stringr::str_remove(cluster_ids, "Cluster_")
  
  cluster_id_pad <- max(nchar(cluster_ids))
  
  # reshape peak ms data
  peak_ms_data_melt <- reshape2::melt(rownames_to_column(peak_ms_data), id = c("rowname"))
  
  # reshape peak row anno
  peak_row_anno_melt <- rownames_to_column(peak_row_anno) %>% select(rowname, `Cluster [C]`, `m/z`)
  
  # get all cluster idss
  unique_cluster_id <- unique(peak_row_anno_melt$`Cluster [C]`)
  unique_sample_names <- unique(peak_ms_data_melt$variable)
  
  # create empty Spectra object
  ms1_spectra <- Spectra()
  
  # iterate over all clusters
  for(cluster in unique_cluster_id) {
    
    # extract number of cluster to work on
    cluster_num <- str_extract(cluster, "\\d+")
    
    # get peak ids
    unique_peak_ids <- peak_row_anno_melt %>% filter(`Cluster [C]` == cluster_num)
    #print(unique_peak_ids)
    
    # innerloop with foreach?
    ms1_spectra_clipboard_list <- foreach(i = 1:length(unique_sample_names),
                                          .packages=c("tidyverse", "MSnbase"),
                                          .combine = 'c') %dopar% {
      
      # get sample to work on
      sample <- unique_sample_names[[i]]
      
      # get peak data
      peak_ms_data_melt_filter <- peak_ms_data_melt %>% filter(rowname %in% unique_peak_ids$rowname,
                                                               variable == sample)
      
      # get values for spectrum
      mz <- unique_peak_ids$`m/z`
      int <- peak_ms_data_melt_filter$value
      
      # replace NA by 0
      int[is.na(int)] <- 0
      
      # create Spectrum1 object
      ms1_spectrum <- new("Spectrum1",
                          mz = mz,
                          intensity = int,
                          centroided = TRUE)

      # return reconstructed MS1 spectrum
      ms1_spectrum
    }
    
    # convert
    ms1_spectra_clipboard_spectra <- Spectra(ms1_spectra_clipboard_list)
    mcols(ms1_spectra_clipboard_spectra)$CLUSTER_ID <- paste0("Cluster_", stringr::str_pad(cluster, cluster_id_pad, pad = "0"))
    mcols(ms1_spectra_clipboard_spectra)$RAWFILE <- unique_sample_names
    
    # add to previous list
    ms1_spectra <- append(ms1_spectra, ms1_spectra_clipboard_spectra)
    
  }
  
  return(ms1_spectra)
  
}

filter_mz_ms1_data <- function(row_anno, ms_data, mz, mzTol = 0.005, mzTolType = "abs") {
  
  # filter row_anno DF based on values that are in the mz list
  
  # filter ms_data based on the row names of row_anno and ms_data
  
  
}