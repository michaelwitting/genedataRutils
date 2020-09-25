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
    
    # res_df <- data.frame(CLUSTER_ID = NA,
    #                            precursorMz = NA,
    #                            rtime = NA,
    #                            forward = NA,
    #                            backward = NA,
    #                            count = NA,
    #                            lib_accession = NA,
    #                            lib_name = NA,
    #                            lib_exactmass = NA,
    #                            lib_adduct = NA,
    #                            lib_precursorMz = NA)
    # 
    # return(res_df)
    
  }
}