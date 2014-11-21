##' @rdname rawgenoBin
##' @title Locate bins using the RawGeno algorithm
##'
##' @description \code{fsaRGbin} identifies bins in a peak table
##'
##' @details TBA
##'
##' @references
##'
##' This function uses the `RawGeno' binning algorithm of
##' \href{http://www.biomedcentral.com/1471-2105/10/33}{Arrigo et al.
##' 2009}. I reimplemented this following their flow-chart, rather than
##' directly copying the code in the RawGeno package. As a consequence,
##' there may be some small numerical differences between the two
##' implementations. I elected not to copy the code as the style of the two
##' packages are not easily reconciled.
##'
##' Arrigo, N., Tuszynski, J. W., Ehrich, D., Gerdes, T., & Alvarez, N.
##' (2009). Evaluating the impact of scoring parameters on the structure of
##' intra-specific genetic variation using RawGeno, an R package for
##' automating AFLP scoring. BMC bioinformatics, 10(1), 33.
##'
##' @param pt a peak table, as produced by \code{fsa2PeakTab}
##'
##' @author Tyler Smith
##'
##' @export
##'
##' @return
##'
##' TBA
##'
##' @examples
##'
##' TBA
##'
##' @keywords aflp binning
fsaRGbin <- function(pt, start = 49.999, end = 495, mxbin = 1.5, mnbin = 1,
                     verbose = FALSE){
  breaks <- c(start, end)
  pt <- pt[which(pt$bp > start & pt$bp < end), ]
  breaks <- splitBin(breaks, pt, start, end, mxbin, mnbin, verbose)
  return(matrix(breaks, byrow = TRUE, ncol = 2))
}

splitBin <- function(breaks, pt, clow, chigh, mxbin, mnbin, verbose = FALSE){
  ptsub <- pt[pt$bp > clow & pt$bp < chigh, ]
  width <- diff(range(ptsub$bp))        # calculate width on actual
                                        # fragment distribution, not on bin
                                        # boundaries!!
  ## debugging a deeply recursive function is tricky. Hence all my messages.
  if(verbose) {
    message("")
    message("*******************")
    message("width = ", width)  
    message("chigh = ", chigh)
    message("clow = ", clow)
    message("mxbin = ", mxbin)
  }
  if(width > mxbin) {
    ## message("-- big bin")
    big <- which.max(diff(ptsub$bp))
    mp <- mean(ptsub$bp[c(big, big + 1)]) #mp = midpoint
    ## print(ptsub$bp)
    mp.low <- mp - 0.05
    ## message("- splitting low half")
    ## message("- clow ", clow)
    ## message("- mp.low ", mp.low)
    breaks <- splitBin(sort(c(breaks, mp.low)), pt, clow = clow, 
                       chigh = mp.low, mxbin, mnbin)
    ## message("- splitting high half")    
    breaks <- splitBin(sort(c(breaks, mp)), pt, clow = mp, 
                       chigh = chigh, mxbin, mnbin)
  } else if(width > mnbin) {
    ## message("-- greater than mnbin")
    if (nrow(ptsub) > length(unique(ptsub$sample.name))) {
      ## message("- multiple peaks")
      ## multiple fragments from the same individual present!
      big <- which.max(diff(ptsub$bp))
      mp <- mean(ptsub$bp[c(big, big + 1)]) #mp = midpoint
      ## WARNING - this magic number, 0.01, can cause problems when the
      ## algorithm creates a bin boundary between two very close
      ## fragments!! Something more sophisticated is necessary to avoid
      ## problematic edge cases.
      mp.low <- mp - 0.01                   
      ## message("head(ptsub$bp): ", paste(round(head(ptsub$bp), 4), collapse = " "))
      ## message("big ", big)
      ## message("mp ", mp)
      ## message("mp.low ", mp.low)
      breaks <- splitBin(sort(c(breaks, mp.low)), pt, clow = clow, 
                         chigh = mp.low, mxbin, mnbin)
      breaks <- splitBin(sort(c(breaks, mp)), pt, clow = mp, 
                         chigh = chigh, mxbin, mnbin)
    }
  }
  return(breaks)
}

