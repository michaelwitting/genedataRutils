#' @title Reformat .gda and MS2 spectra for GNPS FBMN
#'
#' @description
#'
#' `createGnpsFiles` reformats the Gendata data matrix and spectra to be consistent
#'     with the requirement of GNPS FBMN. The resulting data can be fed into the
#'     XCMS3 workflow.
#'
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra`MSnbase Spectra object containing MS2 spectra to be
#'     matched with the MS1 data
#'
#' @return A `list` with the reformated feature table in XCMS3 format as first 
#'     element and a reformated `Spectra` as second element.
#'
#' @author Michael Witting
#'
#' @export
#'
#' @examples
#'
#' 
createGnpsFiles <- function(x, spectra) {
  
  ms_data <- getMsData(x)
  row_anno <- getRowAnno(x)
  
  # sanity checks
  # does the MS2 have correct columns
  # if not perform add ID first
  if(!"CLUSTER_ID" %in% spectraVariables(spectra)) {
    stop("No column CLUSTER_ID found in spectra. Perform ms2_add_id() first")
  }
  
  # check if data contains NAs and remove
  if(any(is.na(spectra$CLUSTER_ID))) {
    # filter spectra
    spectra <- spectra[which(!is.na(spectra$CLUSTER_ID))]
  }
  
  # create feature table
  # create new DF with feature table in XCMS3 format
  feature_table <- data.frame(Row.names = row.names(row_anno),
                              mzmed = row_anno$`m/z`,
                              mzmin = row_anno$`m/z Min`,
                              mzmax = row_anno$`m/z Max`,
                              rtmed = row_anno$RT * 60,
                              rtmin = row_anno$`RT Min` * 60,
                              rtmax = row_anno$`RT Max` * 60,
                              npeaks = ncol(ms_data),
                              sample = ncol(ms_data))
  
  # rename features according to XCMS and add intensities
  feature_table$Row.names <- gsub("Cluster_", "FT", feature_table$Row.names)
  feature_table <- cbind.data.frame(feature_table, ms_data)
  
  # convert names in MS2 metadata to XCMS3 format
  # create fields required
  spectra$scanIndex <- as.integer(regmatches(spectra$CLUSTER_ID, regexpr("\\d+", spectra$CLUSTER_ID)))
  spectra$FEATURE_ID <- gsub("Cluster_", "FT", spectra$CLUSTER_ID)
  spectra$PEAK_ID <- spectra$FEATURE_ID
  spectra$COMPOUND <- spectra$FEATURE_ID

  list(feature_table, spectra)
  
}
