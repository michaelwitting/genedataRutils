test_that("Util functions for working with MS2", {
  
  # read MS1 data
  gda_cluster <- readGda(system.file("extdata/20191108_Pesticides_test_Cluster.gda",
                                     package = "genedataRutils"))
  
  gda_peaks <- readGda(system.file("extdata/20191108_Pesticides_test_Peaks.gda",
                                   package = "genedataRutils"))
  
  spectra <- reconstructIsoPattern(gda_peaks, gda_cluster)
  
  # perform tests
  expect_equal(length(spectra), 16)
  expect_equal(unique(spectra$CLUSTER_ID), c("Cluster_0311", "Cluster_0331", "Cluster_0444", "Cluster_0615"))
  
})