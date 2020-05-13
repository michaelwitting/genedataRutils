#' @title Importing Genedata Expressionist for MS .gda files
#'
#' @description
#'
#' `writeSiriusFile` Genedata does not allow to export MS1 isotope pattern
#'      directly. They have to be reconstructed from the Peak and Cluster data
#'
#' @param ms1_spectrum `Spectra`
#' @param ms2_spectra `Spectra`
#' @param adduct `character`
#' @param con `character`
#'
#' @author Michael Witting
#' 
#' @import Spectra
#'
#' @export
#'
#' @examples
#'
#' 
writeSiriusFile <- function(ms1_spectrum, ms2_spectra = Spectra(), adduct, con = "sirius.ms") {
  
  # some sanity checks here
  if(!"CLUSTER_ID" %in% spectraVariables(ms1_spectrum)) {
    stop("No column CLUSTER_ID found in MS1 spectra.")
  }
  
  if(length(ms2_spectra) > 0 & !"CLUSTER_ID" %in% spectraVariables(ms2_spectra)) {
    stop("No column CLUSTER_ID found in MS2 spectra. Perform ms2AddId() first")
  }
  
  # helper function for writing file
  con <- file(description = con, open = "at")
  on.exit(close(con))
  
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  # sanity check for adduct definitions
  # [M+H]+
  # [M+Na]+
  # [M+ADDUCT]+
  # [M+ADDUCT]-

  .cat(paste0(">compound ", ms1_spectrum$CLUSTER_ID, "\n"))
  .cat(paste0(">ionization ", adduct, "\n"))
  .cat(paste0(">parentmass ", min(ms1_spectrum$mz), "\n"))
  .cat("\n>ms1\n")
  .cat(paste(unlist(ms1_spectrum$mz),
             unlist(ms1_spectrum$intensity),
             collapse = "\n"))
  

  if(length(ms2_spectra) > 0) {
    # iterate over ms2 spectra
    for(i in 1:length(ms2_spectra)) {
      
      if("collisionEnergy" %in% spectraVariables(ms2_spectra[i]) && !is.na(ms2_spectra[i]$collisionEnergy)) {
        
        .cat(paste0("\n\n>collision ", ms2_spectra[i]$collisionEnergy, "\n"))
        .cat(paste(unlist(ms2_spectra[i]$mz),
                   unlist(ms2_spectra[i]$intensity),
                   collapse = "\n"))
        
      } else {
        
        .cat("\n\n>ms2\n")
        .cat(paste(unlist(ms2_spectra[i]$mz),
                   unlist(ms2_spectra[i]$intensity),
                   collapse = "\n"))
        
      }
    }
  }
}
