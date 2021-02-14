test_that("all_equal() tests equaly", {
  a <- rep(0, 10)
  b <- rep("0", 10)
  c <- rep(TRUE, 10)
  d<- rep(factor(b))
  
  expect_identical(all_equal(a), TRUE)
  expect_identical(all_equal(b), TRUE)
  expect_identical(all_equal(c), TRUE)
  expect_identical(all_equal(d), TRUE)
  expect_identical(all_equal(c(a,b)), TRUE)
  expect_identical(all_equal(c(c,d)), TRUE)
  expect_identical(all_equal(c(a,b,c,d)), FALSE)
})
