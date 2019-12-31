# this functions read spectra in mgf files and adds peak or cluster id
#'
#'
#' @import MSnbase
#'
#' @export
ms2_add_id <- function(ms2_spectra, row_anno) {

  # add cluster ID
  mcols(ms2_spectra)$CLUSTER_ID <- unlist(lapply(ms2_spectra, function(x, df) {
    
    id <- row.names(df[which(df$`RT Min` * 60 < rtime(x) &
                               df$`RT Max`* 60 > rtime(x) &
                               df$`m/z Min` < precursorMz(x) &
                               df$`m/z Max` > precursorMz(x)),])

    id[1]
    
  }, df = row_anno))
  
  # return values
  return(ms2_spectra)
  
}

# this function converts the absolute scale to relative scale (scaled to 999 for base peak)
#'
#' @import MSnbase
#' 
#' @export
ms2_abs_to_rel <- function(ms2_spectra, abs_cutoff = 0, rel_cutoff = 10) {
  
  # iterate over Spectra object
  for(i in 1:length(ms2_spectra)) {
    
    # get spectrum
    mz_old <- mz(ms2_spectra[[i]])
    int_old <- intensity(ms2_spectra[[i]])
    
    # filter based on absolute threshold
    mz_old <- mz_old[int_old > abs_cutoff]
    int_old <- int_old[int_old > abs_cutoff]
    
    # convert to rel scale
    mz_new <- mz_old
    int_new <- .abs_to_rel(int_old)
    
    # filter based on relative treshold
    mz_new <- mz_new[int_new > rel_cutoff]
    int_new <- int_new[int_new > rel_cutoff]
    
    # update Spectrum2 object within Spectra
    ms2_spectra[[i]]@mz <- mz_new
    ms2_spectra[[i]]@intensity <- int_new
    ms2_spectra[[i]]@peaksCount <- length(mz_new)
    
  }
  
  return(ms2_spectra)
  
}

# this function converts the absolute scale to relative scale (scaled to 999 for base peak)
#'
#' @import MSnbase
#' @export
ms2_round_mz <- function(ms2_spectra, digits = 4) {
  
  # iterate over Spectra object
  for(i in 1:length(ms2_spectra)) {
    
    # get spectrum
    mz_old <- mz(ms2_spectra[[i]])
    
    # convert to rel scale
    mz_new <- round(mz_old, digits)
    
    # update Spectrum2 object within Spectra
    ms2_spectra[[i]]@mz <- mz_new
    
  }
  
  return(ms2_spectra)
  
}

.abs_to_rel <- function(abs_int) {
  
  rel <- abs_int / max(abs_int) * 999
  
  return(rel)
  
}