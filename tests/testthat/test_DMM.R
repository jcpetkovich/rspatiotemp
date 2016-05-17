s1 = c(0,1,2,3,3,2,1,0,0,3,2,0,1,1,2,3)
s2 = c(1,1,0,2,4,0,3,1,2,3,4,0,0,1,3,4)

test_that("Generate Matrix",{
  resultTransProb = matrix(c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,0.5,0,0,0.5,0,0,0,0,0,0,1,0,0,0,0,0,0.5,0,0,0.5,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0),16,4)
  resultTransCount = matrix(c(1,2,0,1,1,1,2,0,1,1,0,1,0,0,2,1),16,1)
  resultEmisProb = matrix(c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0.5,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0.5,0.5,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0.5,0,0,0,0,0),25,4)
  resultEmisCount = matrix(c(1,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,0,0,2,2,0,0,0,0),25,1)

  resultTrans = list(Probability = resultTransProb,Counter = resultTransCount,AlphabetSize = 4, GroupSize = 2)
  resultEmis = list(Probability = resultEmisProb, Counter = resultEmisCount, AlphabetSize = 5, GroupSize = 2)
  result = list(Transition = resultTrans, Emission = resultEmis)
  expect_equal(createProbMatX(s2,s1,2,5,4),result)
})

test_that("Update Matrix",{
  resultTransProb = matrix(c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,0.5,0,0,0.5,0,0,0,0,0,0,1,0,0,0,0,0,0.5,0,0,0.5,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0),16,4)
  resultTransCount = matrix(c(1,2,0,1,1,1,2,0,1,1,0,1,0,0,2,1),16,1)*2
  resultEmisProb = matrix(c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0.5,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0.5,0.5,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0.5,0,0,0,0,0),25,4)
  resultEmisCount = matrix(c(1,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,0,0,2,2,0,0,0,0),25,1)*2

  resultTrans = list(Probability = resultTransProb,Counter = resultTransCount,AlphabetSize = 4, GroupSize = 2)
  resultEmis = list(Probability = resultEmisProb, Counter = resultEmisCount, AlphabetSize = 5, GroupSize = 2)
  result = list(Transition = resultTrans, Emission = resultEmis)
  expect_equal(updateProb(createProbMatX(s2,s1,2,5,4),s2,s1,2,5,4),result)
})
