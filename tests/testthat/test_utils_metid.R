test_that("Util functions for working with MS1", {
  
  # read MS1 data
  gda_cluster <- readGda(system.file("extdata/20191108_Pesticides_test_Cluster.gda",
                                     package = "genedataRutils"))
  
  ms1library_file <- system.file("extdata/pesticideLibMs1.tsv", package = "genedataRutils")
  
  ms1library <- read.table(ms1library_file, sep = "\t", header = TRUE, comment.char = "")
  
  adducts <- c("[M+H]+")
  
  result <- annotateMz(gda_cluster,
             ms1library,
             tolerance = 0.005,
             adducts = adducts)
  
  expect_equal(nrow(result), 4)
  
})