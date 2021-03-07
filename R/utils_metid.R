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
#' @param rtOffset `numeric` Known offset between data and database
#' @param rtimeTolerance `numeric` tolerance for matching with RT
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
annotateMz <- function(x, ms1library,
                       tolerance = 0,
                       ppm = 0,
                       rtOffset = 0,
                       rtimeTolerance = Inf,
                       matchAdduct = FALSE,
                       adducts = c("[M+H]+")) {

  # sanity checks on ms1Library
  if(!all(c("name", "adduct", "mz") %in% colnames(ms1library))) {
    
    stop("The columns name, adduct and mz are required")
    
  }

  if(!c("rtime") %in% colnames(ms1library)) {
    
    ms1library$rtime <- 0
    
  }

  # filter based on the adducts  
  if(matchAdduct) {
    
    ms1library <- ms1library[which(ms1library$adduct %in% adducts),]
    
  }
  
  row_anno <- getRowAnno(x)
  row_anno$id <- row.names(row_anno)
  
  # order library according to mz for closest function
  ms1library <- ms1library[order(ms1library$mz),]
  
  results <- data.frame()
  
  for(i in 1:nrow(row_anno)) {
    
    # filter based on RT
    ms1libraryfilter <- ms1library[which(abs(ms1library$rtime - row_anno$RT[i]) < rtOffset + rtimeTolerance),]
    
    # use MsCoreUtils closest() to get closests hits
    matches <- closest(ms1libraryfilter$mz, row_anno$`m/z`[i],
                       tolerance = tolerance, ppm = ppm,
                       duplicates = "keep") == 1
    
    matches[is.na(matches)] <- FALSE
    
    if(any(matches)) {
      
      results <- rbind.data.frame(results,
                                  cbind.data.frame(row_anno[i,],
                                                   ms1libraryfilter[matches,]))
      
    }
  }
  
  # TODO add filtering on RT and CCS here
  
  # return result values
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
  
  if("CLUSTER_ID" %in% spectraVariables(x)) {
    
    cluster_id <- x$CLUSTER_ID
    
  } else {
    
    cluster_id <- paste0("M",
                         mz,
                         "T",
                         as.integer(rtime))
    
  }
  
  
  # filter ms2library according to precursor information
  ms2library <- ms2library[which(abs(ms2library$precursorMz - mz) < tolerance + ppm(mz, ppm))]
  ms2library <- ms2library[which(abs(ms2library$rtime - rtime) < rtOffset + rtimeTolerance)]
  
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
                         lib_name = unlist(lapply(ms2library$name, function(x) {paste0(x, collapse = ";")})),
                         #lib_name = paste0(unlist(ms2library$name), collapse = ";"),
                         lib_exactmass = ms2library$exactmass,
                         lib_adduct = ms2library$adduct,
                         lib_precursorMz = ms2library$precursorMz,
                         lib_collisionEnergy = ms2library$collisionEnergy)
    
    
    return(res_df)
    
  } else {
    
  }
}

#' @title  Compare against a MS2 library
#' 
#' @description 
#' 
#' `matchIonMode` allows to matches Cluster from different ion modes against each
#'     other. Cluster are matched by retention time and then checked if they share
#'     a common neutral mass.
#' 
#' @param pos `list` List with data read from .gda file containing grouped 
#'     MS1 cluster from positive mode
#' @param pos `list` List with data read from .gda file containing grouped 
#'     MS1 cluster from negative mode
#' @param pos_adducts `character` vector with adducts in positive mode to check
#' @param neg_adducts `character` vector with adducts in negative mode to check
#' @param tolerance `numeric` absolute tolerance for matching of peaks
#' @param ppm `numeric` relative tolerance for matching of peaks
#' @param rtOffset `numeric` known offset between measurement and library RT
#' @param rtimeTolerance `numeric` tolerance for retention time search
#' 
#' @return `data.frame` with the results
#'
#' @export
matchIonMode <- function(pos,
                         neg,
                         pos_adducts = c("[M+H]+", "[M+Na]+"),
                         neg_adducts = c("[M-H]-"),
                         tolerance = 0,
                         ppm = 0,
                         rtOffset = 0,
                         rtimeTolerance = Inf) {

  # get data
  row_anno_pos <- getRowAnno(pos)
  row_anno_neg <- getRowAnno(neg)
  row_anno_pos$id <- row.names(row_anno_pos)
  row_anno_neg$id <- row.names(row_anno_neg)
  
  match_df <- expand.grid(id.pos = row_anno_pos$id,
                          id.neg = row_anno_neg$id,
                          stringsAsFactors = FALSE)
  
  match_df <- merge(match_df, row_anno_pos, by.x = "id.pos", by.y = "id")
  match_df <- merge(match_df, row_anno_neg, by.x = "id.neg", by.y = "id", suffixes = c(".pos", ".neg"))
  
  match_df <- match_df[which(abs(match_df$RT.pos - match_df$RT.neg) < rtOffset + rtimeTolerance),]
  
  for(i in 1:nrow(match_df)) {
    
    adduct_df <- expand.grid(adduct.pos = pos_adducts,
                             adduct.neg = neg_adducts,
                             stringsAsFactors = FALSE)
    
    adduct_df$mz.pos <- match_df$`m/z.pos`[i]
    adduct_df$mz.neg <- match_df$`m/z.neg`[i]

    adduct_df <- cbind(adduct_df, match = mapply(.adduct_match,
                                                            adduct_df$adduct.pos,
                                                            adduct_df$adduct.neg,
                                                            adduct_df$mz.pos,
                                                            adduct_df$mz.neg,
                                                            tolerance,
                                                            ppm))
    
    
    
    adduct_df <- adduct_df[which(adduct_df$match),]
    
    if(nrow(adduct_df) > 0) {

      match_df$match[i] <- paste0(adduct_df$adduct.pos[1], "<->",adduct_df$adduct.neg[1])
      
    } else {
      
      match_df$match[i] <- NA
      
    }
    
  }
  
  match_df[which(!is.na(match_df$match)),]
  
}

#' @importFrom MetaboCoreUtils mz2mass
#' @importFrom MsCoreUtils ppm
.adduct_match <- function(pos_adduct,
                          neg_adduct,
                          pos_mz,
                          neg_mz,
                          tolerance,
                          ppm) {
  
  # calculate neutral masses
  pos_m <- mz2mass(pos_mz, pos_adduct)
  neg_m <- mz2mass(neg_mz, neg_adduct)
  
  # calculate matching error
  if(neg_m > pos_m) {
    tolerance_new <- tolerance + ppm(neg_m, ppm)
  } else {
    tolerance_new <- tolerance + ppm(pos_m, ppm)
  }
  
  if(abs(pos_m - neg_m) < tolerance_new) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' @title  Compare against a MS2 library
#' 
#' @description 
#' 
#' `addMetid` allows to matches Cluster from different ion modes against each
#'     other. Cluster are matched by retention time and then checked if they share
#'     a common neutral mass.
#' 
#' @param x `list` List with data read from .gda file.
#' @param metids `data.frame`
#' 
#' @return `list` Same as x, but additionally contains columns from metids
#'
#' @import MetaboCoreUtils
#'
#' @export
addMetid <- function(x, metids) {
  
  # sanity checks
  if(!length(x) == 3) {
    
    stop("Input is not of length 3. Sure it contains data from a .gda file?")
    
  }
  
  if(!all(c("DBID", "Formula", "SMILES", "InChI", "Name", "MSI") %in% colnames(metids))) {
    
    stop("metids must contain minimal: DBID, Formula, SMILES, InChI, Name, MSI")
    
  }
  
  if(!all(row.names(metids) %in% rownames(getRowAnno(x)))) {
    
    stop("IDs in metids not found in data")
    
  }
  
  row_anno <- getRowAnno(x)
  
  row_anno <- merge(row_anno, metids, by = "row.names", all = TRUE)
  row.names(row_anno) <- row_anno$Row.names
  row_anno <- row_anno[, !(colnames(row_anno) %in% c("Row.names"))]

  x[[2]] <- row_anno
  
  x
  
}
