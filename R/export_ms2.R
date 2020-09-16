#' @title Exporting Spectra as .mgf file
#'
#' @description
#'
#' `exportMgf` This function allows to export <code>Spectra</code> object to .mgf
#'  files. If the <code>Spectra</code> object contains the CLUSTER_ID and RAWFILE variables 
#'  this is additionally exported as additional metadata fields.
#'
#' @param splist `Spectra` Spectra object containing spectra to be exported
#' 
#' @param con `character` File name for export to .mgf
#' 
#' @param COM `character` Optional comment
#' 
#' @param TITLE `character` Optional title
#'
#' @author Michael Witting
#' 
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

#' private function for writing content
#' @noRd
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
  
  if("RAWFILE" %in% spectraVariables(sp) && length(sp$RAWFILE) && !is.na(sp$RAWFILE)) {
    
    .cat("\nRAWFILE=", sp$RAWFILE)
    
  }
  
  if("CLUSTER_ID" %in% spectraVariables(sp) && length(sp$CLUSTER_ID) && !is.na(sp$CLUSTER_ID)) {
    
    .cat("\nCLUSTER_ID=", sp$CLUSTER_ID)
    
  }

  .cat("\n", paste(sp$mz[[1]], sp$intensity[[1]], collapse = "\n"))
  .cat("\nEND IONS\n")
}
