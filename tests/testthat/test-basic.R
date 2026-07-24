test_that("the installed pipeline exists", {
  expect_true(dir.exists(dona_pipeline_path()))
  expect_true(file.exists(file.path(dona_pipeline_path(), "config.R")))
  expect_true(file.exists(file.path(dona_pipeline_path(), "run_pipeline.R")))
})

test_that("the pipeline can be copied", {
  destination <- tempfile("dona-")
  expect_message(path <- dona_copy_pipeline(destination),"pipeline")
  expect_true(dir.exists(path))
  expect_true(file.exists(file.path(path, "config.R")))
})
