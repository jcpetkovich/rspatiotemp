testSTPN <- function(tabTrained, testObs,testHid,testOSize,testHSize){
  for(i in 1:10){
    ob = testObs[0:(50*i)]
    hi = testHid[0:(50*i)]
    .Call('rspatiotemp_cosMeasure', PACKAGE = 'rspatiotemp', tabTrained, ob, hi, testOSize, testHSize)
    }
}
