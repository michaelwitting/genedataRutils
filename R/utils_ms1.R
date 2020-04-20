#' @title Importing Genedata Expressionist for MS .gda files
#'
#' @description
#'
#' `reconstructIsoPattern` Genedata does not allow to export MS1 isotope pattern
#'      directly. They have to be reconstructed from the Peak and Cluster data
#'
#' @param peaks `list` List with data read from .gda file containing individual 
#'     MS1 peaks
#' @param cluster `list`List with data read from .gda file containing grouped 
#'     MS1 cluster
#'
#' @return A list of three data frames, the first contains the actual MS data,
#'    the second the row annotation and the third the column annotations.
#'
#' @author Michael Witting
#' 
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#'
#' 
reconstructIsoPattern <- function(peaks, cluster) {
  
  # isolate data from lists
  peak_row_anno <- getRowAnno(peaks)
  peak_ms_data <- getMsData(peaks)

  cluster_row_anno <- getRowAnno(cluster)
  
  # add here function to estimate the cluster format
  cluster_ids <- row.names(cluster_row_anno)
  cluster_ids <- gsub("Cluster_", "", cluster_ids)
  
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
    mcols(ms1_spectra_clipboard_spectra)$CLUSTER_ID <- paste0("Cluster_", stringr::str_pad(cluster, cluster_id_pad, pad = "0")) #sprintf(eval(paste0("%0",105, "d")), 5)
    mcols(ms1_spectra_clipboard_spectra)$RAWFILE <- unique_sample_names
    
    # add to previous list
    ms1_spectra <- append(ms1_spectra, ms1_spectra_clipboard_spectra)
    
  }
  
  return(ms1_spectra)
  
}

# ==============================================================================
# DEPRECATED FUNCTIONS
# ==============================================================================
# Genedata does not allow to export MS1 isotope pattern directly. They have to be reconstructed from the Peak and Cluster data
#'
#'
#' @import tidyverse
#' @export
recon_iso_pattern <- function(peak_row_anno, cluster_row_anno, peak_ms_data) {
  
  .Deprecated("reconstructIsoPattern")
  
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