test_that("Util functions for working with MS2", {
  
  # load MS1 data
  gda_file <- system.file("extdata/20191108_Pesticides_test_Cluster.gda",
                          package = "genedataRutils")
  
  gda_result <- readGda(gda_file)
  
  # load MS2 data
  mgf_files <- dir(system.file("extdata",
                               package = "genedataRutils"),
                   pattern = "mgf$",
                   full.names = TRUE)
  
  spectra <- Spectra::Spectra(mgf_files, source = MsBackendMgf())
  
  # add MS1 id
  spectra_new <- ms2AddId(gda_result, spectra)
  
  gnps <- createGnpsFiles(gda_result, spectra_new)
  
  gnps_feature <- gnps[[1]]
  gnps_spectra <- gnps[[2]]
  
  # perform tests
  expect_equal(length(spectra), 512)
  expect_equal(length(spectra_new), 512)
  expect_equal(length(gnps_spectra), 320)
  expect_equal(unique(spectra_new$CLUSTER_ID), c("Cluster_0615", "Cluster_0331", NA, "Cluster_0444", "Cluster_0311"))
  expect_equal(unique(gnps_spectra$CLUSTER_ID), c("Cluster_0615", "Cluster_0331", "Cluster_0444", "Cluster_0311"))
  expect_equal(unique(gnps_spectra$scanIndex), c(615, 331, 444, 311))
  
})


