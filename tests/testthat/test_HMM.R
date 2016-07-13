transProb = matrix(c(0.1,0.2,0.5,0.3,0.3,0.1,0.1,0.2,0.2,0.6,0.2,0.4,0.4,0.1,0.2,0.1),4,4)
emisProb = matrix(c(0.2,0.1,0.5,0.3,0.1,0.3,0.1,0.2,0.6,0.2,0.2,0.4,0.1,0.4,0.2,0.1),4,4)
initProb = c(0.25,0.25,0.25,0.25)

test_that("Viterbi Function",{
  expect_equal(viterbi(matrix(c(0.1,0.9,0.9,0.1),2,2), matrix(c(0.3,0.7,0.7,0.3),2,2), c(0.2,0.8), c(1,0,0,0,1)), c(0,1,0,1,0))
})

test_that("Forward Algorithm",{
  expect_equal(forward(transProb,emisProb,initProb,c(1,3),2,2),0.02)
  expect_equal(forward(transProb,emisProb,initProb,c(1,3,2,1,0),3,5),8.301125e-05)
})

test_that("Viterbi Probability Values",{
  expect_equal(viterbiProbVal(transProb,emisProb,initProb,c(1,3),c(2,1)),c(0.025,0.001))
})
