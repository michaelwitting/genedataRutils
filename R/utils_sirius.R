#' @title Importing Genedata Expressionist for MS .gda files
#'
#' @description
#'
#' `writeSiriusFile` Genedata does not allow to export MS1 isotope pattern
#'      directly. They have to be reconstructed from the Peak and Cluster data
#'
#' @param ms1_spectrum `Spectra` <code>Spectra</code> containing isotope pattern
#' @param ms2_spectra `Spectra` <code>Spectra</code> containing MS2 data
#' @param folder `character` path to folder where .ms files shall be stored
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
writeSiriusFile <- function(ms1_spectra = Spectra(), ms2_spectra = Spectra(), folder = ".") {
  
  # some sanity checks here
  if(!"CLUSTER_ID" %in% Spectra::spectraVariables(ms1_spectra)) {
    
    stop("No column CLUSTER_ID found in MS1 spectra.")
    
  }
  
  if(length(ms2_spectra) > 0 && !"CLUSTER_ID" %in% Spectra::spectraVariables(ms2_spectra)) {
    
    stop("No column CLUSTER_ID found in MS2 spectra. Perform ms2AddId() first")
    
  }
  
  if(length(ms1_spectra) > length(unique(ms1_spectra$CLUSTER_ID))) {
    
    stop("Only single MS1 spectra per CLUSTER_ID are allowed, perform consolidation first.")
    
  }
  
  # get unique cluster ids
  clusters <- unique(ms1_spectra$CLUSTER_ID)
  
  for(cluster in clusters) {

    file_path <- paste0(folder, "/", cluster, ".ms")
    
    if(file.exists(file_path)) {
      
      file.remove(file_path)
      
    }
    
    # helper function for writing file
    con <- file(description = file_path, open = "at")
    #on.exit(close(con))
    
    .cat <- function(..., file = con, sep = "", append = TRUE) {
      cat(..., file = file, sep = sep, append = append)
    }
    
    # select correct spectra
    ms1_spectrum <- ms1_spectra[which(ms1_spectra$CLUSTER_ID == cluster)]
    ms2_spectra_filter <- ms2_spectra[which(ms2_spectra$CLUSTER_ID == cluster)]
    
    # sanity check for adduct definitions
    # [M+H]+
    # [M+Na]+
    # [M+ADDUCT]+
    # [M+ADDUCT]-
    # if(is.na(ms1_spectrum$ADDUCT) && ms1_spectrum$polarity == 0) {
    #   
    #   adduct <- "[M+ADDUCT]-"
    #   
    # } else if(is.na(ms1_spectrum$ADDUCT) && ms1_spectrum$polarity == 1) {
    #   
    #   adduct <- "[M+ADDUCT]+"
    #   
    # } else {
    #   
    #   adduct <- ms1_spectrum$ADDUCT
    #   
    # }
    
    if(!is.null(ms1_spectrum$ADDUCT) && !is.na(ms1_spectrum$ADDUCT)) {
      
      adduct <- ms1_spectrum$ADDUCT
      
    } else {
      
      adduct <- ""
      
    }
    
    .cat(paste0(">compound ", cluster, "\n"))
    .cat(paste0(">ionization ", adduct, "\n"))
    .cat(paste0(">parentmass ", min(ms1_spectrum$mz), "\n"))
    .cat(paste0(">rt ", ms1_spectrum$rtime, "\n"))
    .cat("\n>ms1\n")
    .cat(paste(unlist(ms1_spectrum$mz),
               unlist(ms1_spectrum$intensity),
               collapse = "\n"))
    
    
    if(length(ms2_spectra_filter) > 0) {
      # iterate over ms2 spectra
      for(i in 1:length(ms2_spectra_filter)) {
        
        if("collisionEnergy" %in% spectraVariables(ms2_spectra_filter[i]) && !is.na(ms2_spectra_filter[i]$collisionEnergy)) {
          
          .cat(paste0("\n\n>collision ", ms2_spectra_filter[i]$collisionEnergy, "\n"))
          .cat(paste(unlist(ms2_spectra_filter[i]$mz),
                     unlist(ms2_spectra_filter[i]$intensity),
                     collapse = "\n"))
          
        } else {
          
          .cat("\n\n>ms2\n")
          .cat(paste(unlist(ms2_spectra_filter[i]$mz),
                     unlist(ms2_spectra_filter[i]$intensity),
                     collapse = "\n"))
          
        }
      }
    }
    
    close(con)
    
  }
  
  #on.exit(close(con))
}
