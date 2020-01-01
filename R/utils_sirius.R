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
#' This file reads the results from SIRIUS back into a data frame
#'
#' @export
read_sirius_result <- function(x) {
  
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

#'
#' This file reads the results from SIRIUS back into a data frame
#'
#' @export
read_sirius_result_ms1only <- function(x) {
  
  con <- file(x)
  sirius_result <- readLines(con)
  close(con)
  
  sirius_result_df <- data.frame()
  
  for(i in 1:length(sirius_result)) {
    
    #print(sirius_result[i])
    
    line_split <- unlist(stringr::str_split(sirius_result[i], "\t"))
    
    info <- line_split[1:3]
    results <- line_split[-c(1,2,3)]
    
    dim(results) <- c(2, length(results) / 2)
    
    results <- as.data.frame(t(results))
    colnames(results) <- c("formula", "isotope.score")
    
    
    clipboard <- data.frame(CLUSTER_ID = info[1],
                            exact_mass = info[2],
                            adduct = info[3],
                            formula = results$formula,
                            isotope.score = results$isotope.score)
    
    sirius_result_df <- rbind.data.frame(sirius_result_df, clipboard)
    
  }
  
  return(sirius_result_df)
  
}
