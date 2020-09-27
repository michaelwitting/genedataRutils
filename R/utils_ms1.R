#' @title Reconstructing Isotope Pattern from MS1 data
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
#' @return `Spectra` Returns a Spectra object with the reconstructed MS1 isotope
#'     pattern
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
  
  # determine length of number to correctly create string at the end
  clusterLength <- nchar(gsub("Cluster_", "", row.names(row_anno_cluster)[1]))
  
  clusterIds <- row.names(row_anno_cluster)
  clusterIds <- as.integer(gsub("Cluster_", "", clusterIds))
  
  if(!all(clusterIds %in% row_anno_peaks$`Cluster [C]`)) {
    
    warning("Not all Cluster IDs in peak table found")
    
  }
  
  # prepare data
  peaks_full <- merge(ms_data_peaks, row_anno_peaks,by = "row.names")
  samples <- colnames(ms_data_peaks)
  
  # creacte empty Spectra object to add data
  ms1_spectra <- Spectra()
  
  for(clusterId in clusterIds) {
    
    # get data relevant for cluster
    peaks_full_filter <- peaks_full[which(peaks_full$`Cluster [C]` == clusterId),]
    mz <- peaks_full_filter$`m/z`
    adduct <- .adduct_conversion(row_anno_cluster[paste0("Cluster_", sprintf(eval(paste0("%0", clusterLength, "d")), clusterId)), "Adduct [C]"])
    rtime <- as.numeric(row_anno_cluster[paste0("Cluster_", sprintf(eval(paste0("%0", clusterLength, "d")), clusterId)), "RT"]) * 60
    
    spd <- DataFrame(
      msLevel = rep(1L, length(samples)),
      sample = samples,
      rtime = rep(rtime, length(samples)),
      CLUSTER_ID = rep(paste0("Cluster_", sprintf(eval(paste0("%0", clusterLength, "d")), clusterId)), length(samples)),
      ADDUCT = rep(adduct, length(samples))
    )
    
    spd$mz <- lapply(vector(mode = "list", length = length(samples)), function(x) mz)
    spd$intensity <- unname(as.list(peaks_full_filter[,2:(1+length(samples))]))
    
    ms1_spectra <- c(ms1_spectra, Spectra(spd))
    
    # for(sample in samples) {
    #   
    #   int <- as.vector(peaks_full_filter[,sample])
    # 
    #   int[is.na(int)] <- 0
    #   
    #   # create Spectra object
    #   ms1_spectrum <- DataFrame(
    #     msLevel = 1L,
    #     mz = mz,
    #     intensity = int,
    #     sample = sample,
    #     rtime = rtime,
    #     CLUSTER_ID = paste0("Cluster_", sprintf(eval(paste0("%0", clusterLength, "d")), clusterId)),
    #     ADDUCT = adduct
    #   )
    #   
    #   ms1_spectrum$mz <- IRanges::NumericList(ms1_spectrum$mz)
    #   ms1_spectrum$intensity <- IRanges::NumericList(ms1_spectrum$intensity)
    #   
    #   ms1_spectrum <- Spectra(ms1_spectrum)
    #   
    #   ms1_spectra <- c(ms1_spectra, ms1_spectrum)
    #   
    #   
    # }
    
  }
  
  return(ms1_spectra)
  
}


.adduct_conversion <- function(x) {
  
  if(x == "M+H") {
    
    return("[M+H]+")
    
  } else if(x == "M+Na") {
    
    return("[M+Na]+")
    
  } else if(x == "M+NH4") {
    
    return("[M+NH4]+")
    
  } else if(x == "M-H") {
    
    return("[M-H]-")
    
  } else if(x == "M+FA") {
    
    return("[M+CH2O2-H]-")
    
  } else if(x == "M+HAc") {
    
    return("[M+C2H4O2-H]-")
    
  }
  
  NA_character_
  
}

#' @title Fill with empty MS2 spectra
#' 
#' @description 
#' 
#' `fillMs1Spectra` adds empty spectra for cluster that have no MS2 data
#'     
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra` Spectra object containing MS1 spectra for 
#'     which the RT shall be corrected
#'
#'
#' @export
#' 
fillMs1Spectra <- function(x, spectra) {
  
  # check if spectra have CLUSTER_ID
  if(!"CLUSTER_ID" %in% spectraVariables(spectra)) {
    
    stop("No column CLUSTER_ID found in MS1 spectra.")
    
  }
  
  row_anno <- getRowAnno(x)
  
  row_anno_cluster <- sort(rownames(row_anno))
  spectra_cluster <- sort(unique(spectra$CLUSTER_ID))
  
  diff_cluster <- setdiff(row_anno_cluster, spectra_cluster)
  
  for(diff in diff_cluster) {
    
    # get information
    precusorMz <- row_anno[diff, "m/z"]
    rtime <- row_anno[diff, "RT"] * 600
    
    spd <- DataFrame(
      msLevel = c(2L),
      polarity = c(NA_integer_),
      precursorMz = precursorMz,
      rtime = rtime,
      CLUSTER_ID = diff)
    
    ## Assign m/z and intensity values.
    spd$mz <- list(c())
    spd$intensity <- list(c())
    
    sps <- Spectra(spd)
    
    spectra <- c(spectra, sps)
    
  }
  
  spectra
  
}
