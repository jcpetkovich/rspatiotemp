s1 = c(0,1,2,3,3,2,1,0,0,3,2,0,1,1,2,3)
s2 = c(1,1,0,2,4,0,3,1,2,3,4,0,0,1,3,4)

test_that("Generate Matrix",{
  resultTransProb = matrix(c(0.25,0.5,0,0.25,0.25,0.25,0.5,0,0.25,0.25,0,0.5,0,0,2/3,1/3),4,4)
  resultTransCount = c(4,4,4,3)
  resultEmisProb = matrix(c(rep(0,7),1/3,rep(0,8),1/3,0,0,0,1/3,0,0,0,0,0.25,0.25,0,0.25,0,0,0.25,rep(0,23),0.25,0,0,0.25,rep(0,10),0.25,0.25,rep(0,6),0.25,rep(0,10),0.25,0.25,0,0,0,0,0.25,0,0,0,0,0),25,4)
  resultEmisCount = c(3,4,4,4)
  resultRevEmisProb = matrix(c(0,1,0,0 ,0,1,0,0, 0,0,0,1, 0,1,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0, 1,0,0,0, 0,0,1,0,rep(0, 19),1, 0,0,0,1, 0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0.5,0.5, 0.5,0,0.5,0,rep(0, 16)), 4, 25)
  resultRevEmisCount = c(1,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,0,0,2,2,0,0,0,0)

  result = list(TransProb = resultTransProb, TransCount = resultTransCount,EmisProb = resultEmisProb,EmisCount = resultEmisCount, RevEmisProb = resultRevEmisProb, RevEmisCount = resultRevEmisCount)
  expect_equal(createProbMatX(s2,s1,2,5,4),result)
})

test_that("Update Matrix",{
  resultTransProb = matrix(c(0.25,0.5,0,0.25,0.25,0.25,0.5,0,0.25,0.25,0,0.5,0,0,2/3,1/3),4,4)
  resultTransCount = c(4,4,4,3)*2
  resultEmisProb = matrix(c(rep(0,7),1/3,rep(0,8),1/3,0,0,0,1/3,0,0,0,0,0.25,0.25,0,0.25,0,0,0.25,rep(0,23),0.25,0,0,0.25,rep(0,10),0.25,0.25,rep(0,6),0.25,rep(0,10),0.25,0.25,0,0,0,0,0.25,0,0,0,0,0),25,4)
  resultEmisCount = c(3,4,4,4)*2
  resultRevEmisProb = matrix(c(0,1,0,0 ,0,1,0,0, 0,0,0,1, 0,1,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0, 1,0,0,0, 0,0,1,0,rep(0, 19),1, 0,0,0,1, 0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0.5,0.5, 0.5,0,0.5,0,rep(0, 16)), 4, 25)
  resultRevEmisCount = c(1,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,0,0,2,2,0,0,0,0)*2

  result = list(TransProb = resultTransProb, TransCount = resultTransCount,EmisProb = resultEmisProb,EmisCount = resultEmisCount, RevEmisProb = resultRevEmisProb, RevEmisCount = resultRevEmisCount)
  expect_equal(updateProbMatX(createProbMatX(s2,s1,2,5,4),s2,s1,2,5,4),result)
})
