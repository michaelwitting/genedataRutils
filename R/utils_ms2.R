#' @title Add MS1 ID to MS2 spectra
#'
#' @description
#'
#' `ms2AddId` allows to add MS1 IDs, e.g. "Cluster_0001" to the respective MS2 
#'      data present in a <code>Spectra</code> object. This allows to link MS1 
#'      and MS2 data.
#'
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra` Spectra object containing MS2 spectra to be
#'     matched with the MS1 data
#'
#' @return `Spectra` with the additional column CLUSTER_ID.
#'
#' @author Michael Witting
#' 
#' @importFrom BiocGenerics lapply
#'
#' @export
#'
#' @examples
#'
#' 
ms2AddId <- function(x, spectra) {
  
  row_anno <- getRowAnno(x)

  # add cluster ID
  spectra$CLUSTER_ID <- unlist(spectrapply(spectra, function(x, df) {
    
    id <- row.names(df[which(df$`RT Min` * 60 < x$rtime &
                               df$`RT Max`* 60 > x$rtime &
                               df$`m/z Min` < x$precursorMz &
                               df$`m/z Max` > x$precursorMz),])

    id[1]
    
  }, df = row_anno))
  
  # return spectra
  spectra
  
}

#' @title Correction of retention time in Spectra objects
#' 
#' @description 
#' 
#' `correctRtime` allows to correct the retention time in <code>Spectra</code> 
#'     objects. It uses the RT from the MS1 data supplied from a .gda read data. 
#'     <code>Spectra</code> can be of either MS1 or MS2 level
#'     
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra` Spectra object containing MS1 or MS2 spectra for 
#'     which the RT shall be corrected
#'
#'
#' @export
#' 
correctRtime <- function(x, spectra) {
  
  # check if spectra have CLUSTER_ID
  if(!"CLUSTER_ID" %in% spectraVariables(spectra)) {
    
    stop("No column CLUSTER_ID found in MS1 spectra.")
    
  }
  
  row_anno <- getRowAnno(x)
  
  for(i in 1:length(spectra)) {
    
    clusterid <- spectra[i]$CLUSTER_ID
    rtime <- row_anno[clusterid, "RT"] * 60
    spectra$rtime[i] <- rtime
    
  }
  
  spectra
  
}

#' @title Fill with empty MS2 spectra
#' 
#' @description 
#' 
#' `fillMs2Spectra` adds empty spectra for cluster that have no MS2 data
#'     
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra` Spectra object containing MS2 spectra for 
#'     which the RT shall be corrected
#'
#'
#' @export
#' 
fillMs2Spectra <- function(x, spectra) {
  
  # check if spectra have CLUSTER_ID
  if(!"CLUSTER_ID" %in% spectraVariables(spectra)) {
    
    stop("No column CLUSTER_ID found in MS2 spectra.")
    
  }
  
  row_anno <- getRowAnno(x)
  
  row_anno_cluster <- sort(rownames(row_anno))
  spectra_cluster <- sort(unique(spectra$CLUSTER_ID))
  
  diff_cluster <- setdiff(row_anno_cluster, spectra_cluster)
  
  for(diff in diff_cluster) {
    
    # get information
    precursorMz <- row_anno[diff, "m/z"]
    rtime <- row_anno[diff, "RT"] * 600
    
    spd <- DataFrame(
      msLevel = c(2L),
      polarity = c(NA_integer_),
      precursorMz = as.numeric(precursorMz),
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
