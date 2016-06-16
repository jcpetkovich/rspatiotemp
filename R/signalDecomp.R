#' Decompose a time series signal into several subseries using band-pass filtering based decomposition
#' @title Signal Decomposition (decomp)
#' @param data The time series signal as a vector
#' @param peakThreshold The chosen threshold for which to identify peaks in the spectrum
#' @param gridColNum The number of columns when the displaying the subseries plots. The number of rows is also calculated from this value.
#' @return None. Computes and plots the decompsed subseries of a time series signal
#' @export
decomp <- function(data, peakThreshold, gridColNum){
  f = data.frame(coef = fft(data), freqindex = c(1:length(data)))
  peaks<- Mod(f$coef) > peakThreshold

  midindex <- ceiling((length(f$coef)-1)/ 2) + 1
  peakind <- f[abs(f$coef) > peakThreshold & f$freqindex > 1 & f$freqindex < midindex,]
  lindex <- length(f$coef)

  lowerind <- 1

  subsignals <- lapply(c(peakind$freqindex, midindex+1), function(x){
    upperind <- x
    fsub <- f
    notnullind <- ((fsub$freqindex >= lowerind
                    & fsub$freqindex < upperind)
                   |
                     (fsub$freqindex >  (lindex - upperind + 2)
                      & fsub$freqindex <= (lindex - lowerind + 2)))
    fsub[!notnullind,"coef"] <- 0
    lowerind <<- upperind
    Re(fft(fsub$coef, inverse=TRUE)/length(fsub$coef))
  })

  grid::grid.newpage()
  rowNum = length(subsignals)/gridColNum
  if(!is.integer(rowNum))
    rowNum = as.integer(rowNum)+1
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(rowNum,gridColNum)))

  rowCount = 1
  colCount = 1

  for(x in 1:length(subsignals)){
    #browser()
    h <- data.frame(index = c(1:length(subsignals[[x]])),  orders = subsignals[[x]])
    lab <- paste("Subseries ", as.character(x), sep="")
    #print in viewport. Note: crashes when too many plots. ~6 plots
    #print(ggplot2::qplot(index, orders, data = h, geom = "line", main=lab), vp = grid::viewport(layout.pos.row = rowCount, layout.pos.col = colCount))
    filename = paste(x,".png",sep="")

    #save to disk. Alternative to print in viewport
    ggplot2::qplot(index, orders, data = h, geom = "line", main=lab)
    ggplot2::ggsave(filename)

    colCount = colCount+1
    if(colCount > gridColNum){
      colCount = 1
      rowCount = rowCount+1
    }
  }
}
