test_that("Ol_Reliable returns a ggplot2 theme object", {
  th <- Ol_Reliable()
  expect_s3_class(th, "theme")
  expect_s3_class(th, "gg")
})
