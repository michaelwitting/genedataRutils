#' @title Add MS1 ID to MS2 spectra
#'
#' @description
#'
#' `ms2AddId` allows to add MS1 IDs, e.g. "Cluster_0001" to the respective MS2 
#'      data present in a <code>Spectra</code> object. This allows to link MS1 
#'      and MS2 data.
#'
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra` Spectra object containing MS2 spectra to be
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
  spectra$CLUSTER_ID <- unlist(spectrapply(spectra, function(x, df) {
    
    id <- row.names(df[which(df$`RT Min` * 60 < x$rtime &
                               df$`RT Max`* 60 > x$rtime &
                               df$`m/z Min` < x$precursorMz &
                               df$`m/z Max` > x$precursorMz),])

    id[1]
    
  }, df = row_anno))
  
  # return spectra
  spectra
  
}

#' @title Correction of retention time in Spectra objects
#' 
#' @description 
#' 
#' `correctRtime` allows to correct the retention time in <code>Spectra</code> 
#'     objects. It uses the RT from the MS1 data supplied from a .gda read data. 
#'     <code>Spectra</code> can be of either MS1 or MS2 level
#'     
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param spectra `Spectra` Spectra object containing MS1 or MS2 spectra for 
#'     which the RT shall be corrected
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

#' @title  Compare against a in-house MS2 library
#' 
#' @description 
#' 
#' `compareSpectraInHouse` compares measured spectra stored in a <code>Spectra
#'     </code> object against a library also stored in a <code>Spectra</code> 
#'     object. Difference to `compareSpectraExternal` is that RT is used for 
#'     additional filtering.
#' 
#' @param x `Spectra` Spectra for which library matching shall be performed
#' @param librarySpectra `Spectra` Spectra object containing library
#' @param treshold `numeric` minimum spectral similarity, default is 0.7 for ndotproduct
#' @param tolerance `numeric` absolute tolerance for matching of peaks
#' @param ppm `numeric` relative tolerance for matching of peaks
#' @param rtOffset `numeric` known offset between measurement and library RT
#' @param rtimeTolerance `numeric` tolerance for retention time search
#' 
#' @return `data.frame` with the results
#'
#' @export
#' 
#' @examples 
#' 
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

#' @title  Compare against a MS2 library
#' 
#' @description 
#' 
#' `compareSpectraExternal` compares measured spectra stored in a <code>Spectra
#'     </code> object against a library also stored in a <code>Spectra</code> 
#'     object. Difference to `compareSpectraInHouse` is that RT is ignored.
#' 
#' @param x `Spectra` Spectra for which library matching shall be performed
#' @param librarySpectra `Spectra` Spectra object containing library
#' @param treshold `numeric` minimum spectral similarity, default is 0.7 for ndotproduct
#' @param tolerance `numeric` absolute tolerance for matching of peaks
#' @param ppm `numeric` relative tolerance for matching of peaks
#' 
#' @return `data.frame` with the results
#'
#' @export
#' 
#' @examples 
#' 
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
