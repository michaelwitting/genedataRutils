#' @title Importing Genedata Expressionist for MS .gda files
#'
#' @description
#'
#' `ms2AddId` Genedata does not allow to export MS1 isotope pattern
#'      directly. They have to be reconstructed from the Peak and Cluster data
#'
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra`MSnbase Spectra object containing MS2 spectra to be
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
  spectra$CLUSTER_ID <- unlist(lapply(spectra, function(x, df) {
    
    id <- row.names(df[which(df$`RT Min` * 60 < x$rtime &
                               df$`RT Max`* 60 > x$rtime &
                               df$`m/z Min` < x$precursorMz &
                               df$`m/z Max` > x$precursorMz),])

    id[1]
    
  }, df = row_anno))
  
  # return spectra
  spectra
  
}
