# this functions read spectra in mgf files and adds peak or cluster id
#'
#'
#' @import MSnbase
#' @import masstrixR
#'
ms2_add_id <- function(path_to_mgf, row_anno) {
  
  # get list with all .mgf files
  mgf_files <- list.files(path_to_mgf,
                          pattern = ".mgf$",
                          full.names = TRUE)
  
  # create empty Spectra ojbect
  ms2_spectra <- Spectra()
  
  # read all .mgf files and add peak or cluster id
  for(i in 1:length(mgf_files)) {
    
    # read spectra (requires masstrixR)
    ms2_clipboard <- readAnnotatedMgfData(mgf_files[[i]])
    
    # add cluster ID
    mcols(ms2_clipboard)$id <- unlist(lapply(ms2_clipboard, function(x, df) {
      
      id <- row.names(df[which(df$`RT Min` < rtime(x) &
                       df$`RT Max` > rtime(x) &
                       df$`m/z Min` < precursorMz(x) &
                       df$`m/z Max` > precursorMz(x)),])
      
      id
      
    }, df = row_ann))
    
    # append to previous MS2 spectra
    ms2_spectra <- append(ms2_spectra, ms2_clipboard)
    
  }
  
  # return values
  return(ms2_spectra)
  
}

# this function converts the absolute scale to relative scale (scaled to 999 for base peak)
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
    
    
  }
  
}

.abs_to_rel <- function(abs_int) {
  
  rel <- abs_int / max(abs_int) * 999
  
  return(rel)
  
}