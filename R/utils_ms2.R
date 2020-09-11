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

#'
#'
#'
#' @export
#' 
correctRtime <- function(x, spectra) {
  
  # check if spectra have CLUSTER_ID
  if(!"CLUSTER_ID" %in% spectraVariables(spectra)) {
    
    stop("No column CLUSTER_ID found in MS1 spectra.")
    
  }
  
  row_anno <- getRowAnno(x)
  
  for(i in 1:length(spectra)) {
    
    clusterid <- spectra[i]$CLUSTER_ID
    rtime <- row_anno[clusterid, "RT"] * 60
    spectra$rtime[i] <- rtime
    
  }
  
  spectra
  
}

#'
#' @param x `Spectra` Spectra for which library matching shall be performed
#' @param librarySpectra `Spectra` Spectra object containing library
#' @param treshold `numeric` minimum spectral similarity, default is 0.7 for ndotproduct
#' @param tolerance `numeric` absolute tolerance for matching of peaks
#' @param ppm `numeric` relative tolerance for matching of peaks
#' @param rtOffset `numeric` known offset between measurement and library RT
#' @param rtimeTolerance `numeric` tolerance for retention time search
#'
#' @export
compareSpectraInHouse <- function(x, librarySpectra, treshold = 0.7,
                                  tolerance = 0, ppm = 0,
                                  rtOffset = 0, 
                                  rtimeTolerance = 10, ...) {
  
  # create empty data frame for results
  results <- data.frame()
  
  for(i in 1:length(x)) {
    
    # get information on precursor in data
    mz <- x[i]$precursorMz
    rtime <- x[i]$rtime
    clusterid <- x[i]$CLUSTER_ID
    
    # filter spectra according to precursor information
    # first based on precursorMz
    librarySpectra_filter <- librarySpectra[which(
      abs(librarySpectra$precursorMz - mz) < tolerance + ppm(mz, ppm))]
    
    # second on rtime
    librarySpectra_filter <- librarySpectra_filter[which(
      abs(librarySpectra_filter$rtime + rtOffset - rtime) < rtimeTolerance)]
    
    # get results and format
    if(length(librarySpectra_filter)) {
      
      res <- compareSpectra(x[i],
                            librarySpectra_filter,
                            tolerance = tolerance,
                            ppm = ppm,
                            ...)
      
      res_df <- .isolate_metadata(librarySpectra_filter[res > treshold])
      
      if(nrow(res_df) > 0) {
        
        results <- rbind.data.frame(results ,cbind.data.frame(data.frame(CLUSTER_ID = rep(clusterid, nrow(res_df)),
                                                                         precursorMz = rep(mz, nrow(res_df)),
                                                                         rtime = rep(rtime, nrow(res_df)),
                                                                         similarity = res[res > treshold]),
                                                              res_df))

      }
    }
  }

  results
  
}

#'
#' @param x `Spectra` Spectra for which library matching shall be performed
#' @param librarySpectra `Spectra` Spectra object containing library
#' @param treshold `numeric` minimum spectral similarity, default is 0.7 for ndotproduct
#' @param tolerance `numeric` absolute tolerance for matching of peaks
#' @param ppm `numeric` relative tolerance for matching of peaks
#'
#' @export
compareSpectraExternal <- function(x, librarySpectra, treshold = 0.7,
                                   tolerance = 0, ppm = 0, ...) {
  
  # create empty data frame for results
  results <- data.frame()
  
  for(i in 1:length(x)) {
    
    # get information on precursor in data
    mz <- x[i]$precursorMz
    rtime <- x[i]$rtime
    clusterid <- x[i]$CLUSTER_ID
    
    # filter spectra according to precursor information
    librarySpectra_filter <- librarySpectra[which(
      abs(librarySpectra$precursorMz - mz) < tolerance + ppm(mz, ppm))]
    
    # get results and format
    if(length(librarySpectra_filter)) {
      
      res <- compareSpectra(x[i],
                            librarySpectra_filter,
                            tolerance = tolerance,
                            ppm = ppm,
                            ...)
      
      res_df <- .isolate_metadata(librarySpectra_filter[res > treshold])
      
      if(nrow(res_df) > 0) {
        
        results <- rbind.data.frame(results ,cbind.data.frame(data.frame(CLUSTER_ID = rep(clusterid, nrow(res_df)),
                                                                         precursorMz = rep(mz, nrow(res_df)),
                                                                         rtime = rep(rtime, nrow(res_df)),
                                                                         similarity = res[res > treshold]),
                                                              res_df))
        
      }
    }
  }
  
  results
  
}

#'
#'
#' helper function to isolate metadata from different backends
.isolate_metadata <- function(x) {
  
  if(class(x@backend) == "MsBackendMassbank") {
    
    lib_id <- x$accession
    lib_names <- unlist(lapply(x$name, function(x) {paste0(x, collapse = ";")}))
    lib_formula <- x$formula
    lib_precursorMz <- x$precursorMz
    lib_rtime <- x$rtime
    lib_adduct <- x$focus_precursor_type
    lib_inchikey <- x$link_inchikey
    lib_inchi <- x$iupac
    lib_smiles <- x$smiles
    
    
  } else if(class(x@backend) == "MsBackendMsp") {
    
    lib_id <- x$accession
    lib_names <- unlist(lapply(x$name, function(x) {paste0(x, collapse = ";")}))
    lib_formula <- x$formula
    lib_precursorMz <- x$precursorMz
    lib_rtime <- x$rtime
    lib_adduct <- x$focus_precursor_type
    lib_inchikey <- x$link_inchikey
    lib_inchi <- x$iupac
    lib_smiles <- x$smiles
    
  } else {
    
  }
  
  data.frame(lib_id = lib_id,
             lib_names = lib_names,
             lib_formula = lib_formula,
             lib_libMz = lib_precursorMz,
             lib_libRtime = lib_rtime,
             lib_adduct = lib_adduct,
             lib_inchikey = lib_inchikey,
             lib_inchi = lib_inchi,
             lib_smiles = lib_smiles)
  
}
