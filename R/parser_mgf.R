#' This function allows to read mgf files which contain additional fields to standard
#' ased on the code contributed by Guangchuang Yu \email{guangchuangyu@@gmail.com}
#' Modified by Sebastian Gibb \email{mail@@sebastiangibb.de}
#' Modified by Michael Witting \email{michael.witting@@helmholtz-muenchen.de}
#'
#'
#' @import MSnbase
#' @export
readAnnotatedMgfData <- function(filename,
                                 filterAddFields = NULL,
                                 centroided = TRUE,
                                 smoothed = FALSE,
                                 cache = 1) {
  
  mgf <- scan(file = filename, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)
  
  ## From http://www.matrixscience.com/help/data_file_help.html#GEN
  ## Comment lines beginning with one of the symbols #;!/ can be included,
  ## but only outside of the BEGIN IONS and END IONS statements that delimit an MS/MS dataset.
  cmts <- grep("^[#;!/]", mgf)
  if (length(cmts))
    mgf <- mgf[-cmts]
  
  # find begin and end of spectra
  begin <- grep("BEGIN IONS", mgf) + 1L
  end <- grep("END IONS", mgf) - 1L
  
  n <- length(begin)
  
  # prepare listes
  spectra <- vector("list", length = n)
  fdata <- vector("list", length = n)
  addFieldsVec <- vector("list", length = n)
  
  # iterate over all spectra in .mgfi file
  for (i in seq(along = spectra)) {
    
    specInfo <- extractMgfSpectrum2Info(mgf[begin[i]:end[i]],
                                        centroided = centroided,
                                        filterAddFields =  filterAddFields)
    spectra[[i]] <- specInfo$spectrum
    fdata[[i]] <- specInfo$fdata
    addFieldsVec[[i]] <- specInfo$addData
    
  }
  
  # convert to data frame
  fdata <- do.call(rbind, fdata)
  addFieldsVec <- do.call(rbind, addFieldsVec)
  
  addFieldsDf <- as.data.frame(addFieldsVec)
  
  # check if RAWFILE exists, otherwise add
  if(!"RAWFILE" %in% colnames(addFieldsDf))
  {
    addFieldsDf$RAWFILE <- basename(filename)
  }
  
  # convert data types (based on tidyverse)
  addFieldsDf <- type.convert(addFieldsDf)
  
  # create Spectra object and add metadata
  ms2Spectra <- Spectra(spectra)
  mcols(ms2Spectra) <- as.data.frame(addFieldsVec)
  
  # return SpectraList
  return(ms2Spectra)
  
}

#'
#'
#'
#' @import MSnbase
extractMgfSpectrum2Info <- function(mgf, centroided, filterAddFields = NULL) {
  
  # grep description
  desc.idx <- grep("=", mgf)
  desc <- mgf[desc.idx]
  spec <- mgf[-desc.idx]
  
  ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
  mode(ms) <- "double"
  
  if(!length(ms)) {
    ms <- matrix(numeric(), ncol = 2L)
  }
  
  r <- regexpr("=", desc, fixed = TRUE)
  desc <- setNames(substring(desc, r + 1L, nchar(desc)), substring(desc, 1L, r - 1L))
  fdata <- desc
  
  # some data prep
  desc[c("PEPMASSMZ", "PEPMASSINT")] <- strsplit(desc["PEPMASS"], "[[:space:]]+")[[1L]][1:2]
  desc["CHARGE"] <- sub("[+-]", "", desc["CHARGE"])
  
  # select only values of interest and convert to numeric (base fields)
  voi <- c("RTINSECONDS", "CHARGE", "SCANS", "PEPMASSMZ", "PEPMASSINT")
  desc.base <- setNames(as.numeric(desc[voi]), voi)
  desc.base[is.na(desc.base[voi])] <- 0L
  cat(".")
  
  # select additional values
  fieldNames <- names(desc)
  addFieldNames <- fieldNames[!fieldNames %in% c(voi, "PEPMASS", "TITLE")]
  
  # filter only fields wanted by user
  if(!is.null(filterAddFields)) {
    addFieldNames <- filterAddFields
  }
  
  # get additional descriptors
  desc.add <- setNames(desc[addFieldNames], addFieldNames)
  
  
  # create spectrum
  sp <- new("Spectrum2",
            rt = as.numeric(unname(desc["RTINSECONDS"])),
            scanIndex = as.integer(unname(as.integer(desc["SCANS"]))),
            precursorMz = as.numeric(unname(desc["PEPMASSMZ"])),
            precursorIntensity = as.numeric(unname(desc["PEPMASSINT"])),
            precursorCharge = as.integer(unname(as.integer(desc["CHARGE"]))),
            mz = ms[, 1L],
            intensity = ms[, 2L],
            fromFile = 1L,
            centroided = centroided)
  
  # return values
  return(list(spectrum = sp, fdata = fdata, addData = desc.add))
}
