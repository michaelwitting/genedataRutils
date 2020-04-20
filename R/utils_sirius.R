#'
#'
#'
#' @export
write_ms_file <- function(ms1_spectrum, ms2_spectra, precursor, adduct, con) {
  
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

  .cat(paste0(">compound ", ms1_spectrum[1]@elementMetadata$CLUSTER_ID, "\n"))
  .cat(paste0(">ionization ", adduct, "\n"))
  .cat(paste0(">parentmass ", precursor, "\n"))
  .cat("\n>ms1\n")
  .cat(paste(mz(ms1_spectrum[[1]]), intensity(ms1_spectrum[[1]]), collapse = "\n"))
  
  # iterate over ms2 spectra
  for(i in 1:length(ms2_spectra)) {
    
    if(!length(collisionEnergy(ms2_spectra[[i]])) == 0) {
      
      .cat(paste0("\n\n>collision ", collisionEnergy(ms2_spectra[[i]]), "\n"))
      .cat(paste(mz(ms2_spectra[[i]]), intensity(ms2_spectra[[i]]), collapse = "\n"))
      
    } else {
      
      .cat("\n\n>ms2\n")
      .cat(paste(mz(ms2_spectra[[i]]), intensity(ms2_spectra[[i]]), collapse = "\n"))
      
    }
  }
}


#'
#'
#'
#' @export
write_ms_file_ms1only <- function(ms1_spectrum, precursor, adduct, con) {
  
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
  
  .cat(paste0(">compound ", ms1_spectrum[1]@elementMetadata$CLUSTER_ID, "\n"))
  .cat(paste0(">ionization ", adduct, "\n"))
  .cat(paste0(">parentmass ", precursor, "\n"))
  .cat("\n>ms1\n")
  .cat(paste(mz(ms1_spectrum[[1]]), intensity(ms1_spectrum[[1]]), collapse = "\n"))
  
}
