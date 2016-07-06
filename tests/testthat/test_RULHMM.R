test_that("Critical Path - Single Path greater than One Step",{
  expect_equal(criticalPath(matrix(c(0.05,0.05,0.5,0.9,0.05,0.25,0.05,0.9,0.25),3,3),0,2,0:2),c(2,1,0))
})

test_that("Critical Path - One Step path",{
    expect_equal(criticalPath(matrix(c(0.1,0.3,0.3,0.1,0.3,0.3,0.8,0.4,0.4),3,3),0,2,0:2),c(0,2))
})

test_that("Critical Path - Two Paths greater than One Step",{
    expect_equal(criticalPath(matrix(c(0.1,0.1,0.25,0,0.85,0.1,0.25,0.05,0,0.1,0.25,0.9,0.05,0.7,0.25,0.05),4,4),0,2,0:3),c(2,3,1,0))
})