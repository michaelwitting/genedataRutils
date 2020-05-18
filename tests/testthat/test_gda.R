test_that("parsing of .gda files", {
  
  # load example file (Cluster)
  gda_file <- system.file("extdata/20191108_Pesticides_test_Cluster.gda",
                     package = "genedataRutils")
  
  # read example file
  gda_result <- readGda(gda_file)
  ms_data <- getMsData(gda_result)
  row_anno <- getRowAnno(gda_result)
  col_anno <- getColAnno(gda_result)
  
  # perform tests
  expect_equal(length(gda_result), 3)
  expect_equal(ncol(ms_data), 4)
  expect_equal(nrow(ms_data), 4)
  expect_equal(ncol(row_anno), 15)
  expect_equal(nrow(row_anno), 4)
  expect_equal(ncol(col_anno), 7)
  expect_equal(nrow(col_anno), 4)
  
  # load example file (Peaks)
  gda_file <- system.file("extdata/20191108_Pesticides_test_Peaks.gda",
                          package = "genedataRutils")
  
  # read example file
  gda_result <- readGda(gda_file)
  ms_data <- getMsData(gda_result)
  row_anno <- getRowAnno(gda_result)
  col_anno <- getColAnno(gda_result)
  
  # perform tests
  expect_equal(length(gda_result), 3)
  expect_equal(ncol(ms_data), 4)
  expect_equal(nrow(ms_data), 15)
  expect_equal(ncol(row_anno), 11)
  expect_equal(nrow(row_anno), 15)
  expect_equal(ncol(col_anno), 7)
  expect_equal(nrow(col_anno), 4)
  
})
