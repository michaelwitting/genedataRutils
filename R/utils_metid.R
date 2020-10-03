#' @title  Annotate MS1 data
#' 
#' @description 
#' 
#' `annotateMz` compares measured m/z values with a MS1 library
#' 
#' @param x `list` List with data read from .gda file containing grouped 
#'     MS1 cluster
#' @param ms1library `data.frame` Data frame containing MS1 library, minimal 
#'     columns are name, adduct and mz
#' @param tolerance `numeric` absolute tolerance for matching of peaks
#' @param ppm `numeric` relative tolerance for matching of peaks
#' @param matchAdducts `boolean` Indicates if MS1 library shall be first 
#'     filtered based on the values in `adducts`
#' @param adducts `character` Vector with adducts to be used. Names of adducts 
#'     between `ms1library` and `adducts` have to match!
#' 
#' @return `data.frame` with the results
#'
#' @export
#' 
#' @examples 
#' 
annotateMz <- function(x, ms1Library,
                       tolerance = 0, ppm = 0,
                       rtimeTolerance = Inf,
                       matchAdduct = FALSE,
                       adducts = c("[M+H]+")) {

  # sanity checks on ms1Library
  if(!all(c("name", "adduct", "mz") %in% colnames(ms1Library))) {
    
    stop("The columns name, adduct and mz are required")
    
  }

  # filter based on the adducts  
  if(matchAdduct) {
    
    ms1library <- ms1library[which(ms1library$adduct %in% adducts),]
    
  }
  
  row_anno <- getRowAnno(x)
  
  # order library according to mz for closest function
  ms1library <- ms1library[order(ms1library$mz),]
  
  results <- data.frame()
  
  for(i in 1:nrow(row_anno)) {
    
    # use MsCoreUtils closest() to get closests hits
    matches <- closest(ms1library$mz, row_anno$`m/z`[i],
                       tolerance = tolerance, ppm = ppm,
                       duplicates = "keep") == 1
    matches[is.na(matches)] <- FALSE
    
    if(any(matches)) {
      
      results <- rbind.data.frame(results,
                                  cbind.data.frame(row_anno[i,],
                                                   ms1library[matches,]))
      
    }
  }
  
  results

}


#' @title  Compare against a MS2 library
#' 
#' @description 
#' 
#' `compareSpectraLibrary` compares measured spectra stored in a <code>Spectra
#'     </code> object against a library also stored in a <code>Spectra</code> 
#'     object. By default rtimeTolerance is set to Inf, if RT matching shall be 
#'     performed a different value in seconds needs to be set.
#' 
#' @param x `Spectra` Spectra for which library matching shall be performed
#' @param ms2library `Spectra` Spectra object containing library
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
compareSpectraLibrary <- function(x,
                                  ms2library,
                                  tolerance = 0,
                                  ppm = 0,
                                  rtOffset = 0, 
                                  rtimeTolerance = Inf,
                                  plot = FALSE) {
  
  # sanity checks
  if(!all(c("accession", "name", "exactmass", "adduct", "precursorMz") %in%
     spectraVariables(ms2library))) {
    
    stop("Missing basic library information: accession, name, exactmass, adduct and precursorMz")
    
  }

  results_list <- spectrapply(x,
                              FUN = .compareLibrary,
                              ms2library = ms2library,
                              tolerance = tolerance,
                              ppm = ppm,
                              rtOffset = rtOffset,
                              rtimeTolerance = rtimeTolerance,
                              plot = plot)

  do.call(rbind, results_list)

}

#'
#'
#' @noRd
.compareLibrary <- function(x,
                            ms2library,
                            tolerance = 0,
                            ppm = 0,
                            rtOffset = 0,
                            rtimeTolerance = Inf,
                            plot = FALSE) {

  # get precursor information
  mz <- x$precursorMz
  rtime <- x$rtime
  cluster_id <- x$CLUSTER_ID
  
  # filter ms2library according to precursor information
  ms2library <- ms2library[which(abs(ms2library$precursorMz - mz) < tolerance + ppm(mz, ppm))]
  ms2library <- ms2library[which(abs(ms2library$rtime + rtOffset - rtime) < rtimeTolerance)]
  
  # costum matching function
  specNrow <-  function(x, y, ...) nrow(x)
  mydotproduct <- function(x, y, type = NA, ...) MsCoreUtils::ndotproduct(x, y, ...)
  
  if(length(ms2library)) {
    
    # type = "right" keeps only matching peaks = reverse dot product
    # type = "outer" keeps all peaks = forward dot product
    # use specNrow only with type = "inner"
    
    res_forward <- compareSpectra(x,
                                  ms2library,
                                  FUN = mydotproduct,
                                  tolerance = tolerance,
                                  ppm = ppm,
                                  type = "outer")
    
    res_backward <- compareSpectra(x,
                                   ms2library,
                                   FUN = mydotproduct,
                                   tolerance = tolerance,
                                   ppm = ppm,
                                   type = "right")
    
    res_count <- compareSpectra(x,
                                ms2library,
                                FUN = specNrow,
                                tolerance = tolerance,
                                ppm = ppm,
                                type = "inner")
    
    res_df <- data.frame(CLUSTER_ID = rep(cluster_id, length(ms2library)),
                         precursorMz = rep(mz, length(ms2library)),
                         rtime = rep(rtime, length(ms2library)),
                         collisionEnergy = rep(x$collisionEnergy, length(ms2library)),
                         forward = res_forward,
                         backward = res_backward,
                         count = res_count,
                         lib_accession = ms2library$accession,
                         lib_name = paste0(unlist(ms2library$name), collapse = ";"),
                         lib_exactmass = ms2library$exactmass,
                         lib_adduct = ms2library$adduct,
                         lib_precursorMz = ms2library$precursorMz,
                         lib_collisionEnergy = ms2library$collisionEnergy)
    
    
    return(res_df)
    
  } else {
    
  }
}