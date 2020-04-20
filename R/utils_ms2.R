#' @title Importing Genedata Expressionist for MS .gda files
#'
#' @description
#'
#' `reconstructIsoPattern` Genedata does not allow to export MS1 isotope pattern
#'      directly. They have to be reconstructed from the Peak and Cluster data
#'
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra`MSnbase Spectra object containing MS2 spectra to be
#'     matched with the MS1 data
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
#' 
ms2AddId <- function(x, spectra) {
  
  row_anno <- getRowAnno(x)

  # add cluster ID
  mcols(spectra)$CLUSTER_ID <- unlist(lapply(spectra, function(x, df) {
    
    id <- row.names(df[which(df$`RT Min` * 60 < rtime(x) &
                               df$`RT Max`* 60 > rtime(x) &
                               df$`m/z Min` < precursorMz(x) &
                               df$`m/z Max` > precursorMz(x)),])

    id[1]
    
  }, df = row_anno))
  
  # return spectra
  spectra
  
}

# ==============================================================================
# DEPRECATED FUNCTIONS
# ==============================================================================
# this functions read spectra in mgf files and adds peak or cluster id
#'
#'
#' @import MSnbase
#'
#' @export
ms2_add_id <- function(ms2_spectra, row_anno) {
  
  .Deprecated("ms2AddId")
  
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
