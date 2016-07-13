test_that("Get Energy Function",{
  expect_equal(getEnergy(matrix(c(1,2,3,4,5,6),2,3)), c(1.58113883,3.5355339059,5.522680509))
})

test_that("Format RUL HMM Function",{
  expect_equal(formatRULHMM(c(1,1,1,2,2,2,3,3,3,1,2,1,1)),list(redSeq = c(1,2,3,1,2,1), repSeq = c(3,3,3,1,1,2)))
})

test_that("Mean and Standard Deviation Function",{
  expect_equal(meanStdDev(c(1,2,3,1,2,1),c(3,3,3,1,1,2),3,T),list(mean = c(2,2,3),stdDev = c(0.8164965809,0.99999999995,0),tab = c(3,2,1)))
})

test_that("Critical Path - Single Path greater than One Step",{
  expect_equal(criticalPath(matrix(c(0.05,0.05,0.5,0.9,0.05,0.25,0.05,0.9,0.25),3,3),0,2,0:2),c(2,1,0))
})

test_that("Critical Path - One Step path",{
  expect_equal(criticalPath(matrix(c(0.1,0.3,0.3,0.1,0.3,0.3,0.8,0.4,0.4),3,3),0,2,0:2),c(0,2))
})

test_that("Critical Path - Two Paths greater than One Step",{
  expect_equal(criticalPath(matrix(c(0.1,0.1,0.25,0,0.85,0.1,0.25,0.05,0,0.1,0.25,0.9,0.05,0.7,0.25,0.05),4,4),0,2,0:3),c(2,3,1,0))
})

test_that("Compute RUL Bounds Function",{
  expect_equal(computeRULBounds(c(0,1,2),c(1,1.5,1.25),c(0.5,0.75,1),0.25),list(Lower = 3.1875,Mean = 3.75, Upper = 4.3125))
})

test_that("Viterbi Probability Depmix",{
  transProb = matrix(c(0.1,0.3,0.25,0.6,0.4,0.1,0.25,0.1,0.3,0.4,0.3,0.2,0.2,0.2,0.2,0.1),4,4)
  emisProb = matrix(c(0,1,0.25,1,0.5,1,1,1),2,4)
  hidSeq = c(1,1,2,3,2)
  obsSeq = c(0.1,0.2,0.6,0.7,0.3)
  expect_equal(viterbiProbDepmix(transProb,emisProb,obsSeq,hidSeq),-6.68682442)
})

test_that("Viterbi Continuous",{
  transProb = matrix(c(0.1,0.3,0.25,0.6,0.4,0.1,0.25,0.1,0.3,0.4,0.3,0.2,0.2,0.2,0.2,0.1),4,4)
  emisProb = matrix(c(0,1,0.25,1,0.5,1,1,1),2,4)
  seq = c(0.1,0.2,0.6,0.7,0.3)
  expect_equal(viterbiCont(transProb,emisProb,seq), c(0,1,2,3,0))
})
