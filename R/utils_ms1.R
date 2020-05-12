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
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @import Spectra
#'
#' @export
#'
#' @examples
#'
#' 
reconstructIsoPattern <- function(peaks, cluster) {
  
  # isolate data
  ms_data_cluster <- getMsData(cluster)
  row_anno_cluster <- getRowAnno(cluster)
  ms_data_peaks <- getMsData(peaks)
  row_anno_peaks <- getRowAnno(peaks)
  
  # check that sample names are the same
  if(!identical(sort(colnames(ms_data_cluster)),
                sort(colnames(ms_data_peaks)))) {
    
    stop("peaks and cluster are not from the same experiment")
    
  }
  
  clusterLength <- length(gsub("Cluster_", "", row.names(ms_data_cluster)))
  
  clusterIds <- row.names(ms_data_cluster)
  clusterIds <- as.integer(gsub("Cluster_", "", clusterIds))
  
  all(clusterIds %in% row_anno_peaks$`Cluster [C]`)
  
  peaks_full <- merge(row_anno_peaks, ms_data_peaks, by = "row.names")
  
  samples <- colnames(ms_data_peaks)
  
  ms1_spectra <- Spectra()
  
  for(clusterId in clusterIds) {
    
    peaks_full_filter <- peaks_full[which(peaks_full$`Cluster [C]` == clusterId),]
    mz <- peaks_full_filter$`m/z`
    
    for(sample in samples) {
      
      int <- peaks_full_filter[,sample]
      
      # create Spectra object
      ms1_spectrum <- DataFrame(
        msLevel = 1L,
        mz = mz,
        intensity = int,
        sample = sample,
        CLUSTER_ID = paste0("Cluster_", sprintf(eval(paste0("%0", clusterLength, "d")), clusterId))
      )
      
      ms1_spectrum$mz <- IRanges::NumericList(ms1_spectrum$mz)
      ms1_spectrum$intensity <- IRanges::NumericList(ms1_spectrum$intensity)
      
      ms1_spectrum <- Spectra(ms1_spectrum)
      
      ms1_spectra <- c(ms1_spectra, ms1_spectrum)
      
      
    }
    
  }
  
  return(ms1_spectra)
  
}

