#' @title Exporting Spectra as .mgf file
#'
#' @description
#'
#' `exportMgf` Genedata does not allow to export MS1 isotope pattern
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
exportMgf <- function(splist, con, COM = NULL, TITLE = NULL) {
  
  if (class(con) == "character" && file.exists(con)) {
    
    message("Overwriting ", con, "!")
    unlink(con)
    
  }

  if (class(con)[1] == "character") {
    con <- file(description = con, open = "at")
    on.exit(close(con))
  }
  
  if (is.null(COM)) {
    COM <- paste0(ifelse(length(splist) <= 1, "Spectrum", "Experiment"),
                  "exported by gnedataRutils on ", date())
  }
  
  cat(paste0("COM=",COM), file = con, sep = "")
  

  for (i in seq(along=splist)) {

    .writeMgfContent(splist[i], TITLE = TITLE, con = con)
    
  }

}

.writeMgfContent <- function(sp, TITLE = NULL, con) {
  
  # custom cat function for writing of content
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  .cat("\nBEGIN IONS\n",
       "SCANS=",sp$acquisitionNum)
  
  if (is.null(TITLE)) {
    .cat("\nTITLE=msLevel ", sp$msLevel,
         "; retentionTime ", sp$rtime,
         "; scanNum ", sp$acquisitionNum)
    
    if (length(sp$scanIndex)) {
      .cat("; scanIndex ", sp$scanIndex)
    }
    
    if (sp$msLevel > 1) {
      .cat("; precMz ", precursorMz(sp),
           "; precCharge ", precursorCharge(sp))
    }
  } else {
    .cat("\nTITLE=", TITLE)
  }
  
  if (sp$msLevel > 1) {
    .cat("\nRTINSECONDS=", sp$rtime,
         "\nPEPMASS=", sp$precursorMz)
    
    if (length(sp$precursorCharge) && !is.na(sp$precursorCharge)) {
      
      .cat("\nCHARGE=", sp$precursorCharge, "+")
      
    }
    
  } else {
    
    .cat("\nRTINSECONDS=", sp$rtime)
    
  }
  
  if(length(sp$CLUSTER_ID) && !is.na(sp$CLUSTER_ID)) {
    
    .cat("\nCLUSTERID=", sp$CLUSTER_ID)
    
  }

  .cat("\n", paste(sp$mz, sp$intensity, collapse = "\n"))
  .cat("\nEND IONS\n")
}
