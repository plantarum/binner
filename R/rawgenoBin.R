##' @rdname rawgenoBin
##' @title Locate bins using the RawGeno algorithm
##'
##' @description \code{fsaRGbin} identifies bins in a peak table
##'
##' @details
##'
##' fsaRGbin calculates bins using all peaks in a peak table, including
##' very small ones. As a consequence, some bins may include only peaks <
##' 50 rfus, and will therefore appear to be 'empty' in \code{scanGel()}.
##' These empty bins should be filtered out before analysis, and are not
##' otherwise of any concern.
##'
##' Note that the way the algorithm is coded, the largest bin will always
##' be bounded on the upper side by 495. Consequently, this bin may appear
##' quite large. It will only ever include fragments that are collectively
##' within mxbin base pairs of each other, but may include a sizeable empty
##' section above the actual fragments. Again, this is not a concern for
##' subsequent analysis.
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
  bins <- matrix(breaks, byrow = TRUE, ncol = 2)
  return(bins)
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
  if((width > mxbin) |                  # bin too wide
     ((width > mnbin) & (nrow(ptsub) >  # bin has multiple peaks from same
                                        # sample 
                           length(unique(ptsub$sample.name))))) {
    ## message("-- big bin")
    big <- which.max(diff(ptsub$bp))
    mp <- max(mean(ptsub$bp[c(big, big + 1)]), # mp = midpoint
              ptsub$bp[big + 1] - 0.2) # just below lowest frag in upper
                                        # bin 
    ## print(ptsub$bp)
    mp.low <- min(mp - 0.01, # just below midpoint
                  ptsub$bp[big] + 0.2) # just above highest frag in lower bin
    ## message("- splitting low half")
    ## message("- clow ", clow)
    ## message("- mp.low ", mp.low)
    breaks <- splitBin(sort(c(breaks, mp.low)), pt, clow = clow, 
                       chigh = mp.low, mxbin, mnbin)
    ## message("- splitting high half")    
    breaks <- splitBin(sort(c(breaks, mp)), pt, clow = mp, 
                       chigh = chigh, mxbin, mnbin)
  }  return(breaks)
}

