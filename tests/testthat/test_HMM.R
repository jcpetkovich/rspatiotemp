test_that("Viterbi Function",{
  expect_equal(viterbi(matrix(c(0.1,0.9,0.9,0.1),2,2), matrix(c(0.3,0.7,0.7,0.3),2,2), c(0.2,0.8), c(1,0,0,0,1)), c(0,1,0,1,0))
})
