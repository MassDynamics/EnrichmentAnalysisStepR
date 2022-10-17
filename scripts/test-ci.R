.libPaths("/usr/local/lib64/R/library/")

devtools::load_all()
results <- devtools::test()

for (i in 1:length(results)) {
  if (!is(results[[i]]$result[[1]], "expectation_success")) {
    stop("There were test failures")
  }
}
