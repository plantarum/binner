##' @rdname editFSA
##' @title View and manipulate fsa files
##'
##' @description Functions to view electropherograms, rename or remove
##' samples, combine fsa objects and normalize data.
##'
##' @details
##'
##' \code{fsaDrop} removes a sample from an \code{fsa} object.
##'
##' \code{fsaRename} renames a sample in an \code{fsa} object.
##'
##' \code{fsaCombine} combines two \code{fsa} objects
##'
##' \code{fsaNormalize} normalizes the rfu values across each sample in
##' the fsa. 
##'
##' \code{fsa2PeakTab} converts the fsa to a peakscanner-style peaks
##' table. 
##'
##' @param fsa An fsa object, as produced by \code{\link{readFSA}}.
##'
##' @param epn The index of one or more samples. Can be an integer or a
##' character string. Integers refer to samples in the order they are
##' stored in an \code{fsa} object; character strings refer to them by
##' name.
##'
##' @param newname The new name that will be applied to sample \code{epn}.
##'
##' @param fsa1 An fsa object, as produced by \code{readFSA}
##'
##' @param fsa2 An fsa object, as produced by \code{readFSA}
##'
##' @param dye Which dye to read when converting to a Peakscanner table.
##' Valid values include "d6.FAM", "VIC, "NED" and "PET".
##'
##' @author Tyler Smith
##'
##' @return
##'
##' \code{fsaDrop} returns an \code{fsa} object with the indicated samples
##' removed.
##'
##' \code{fsaRename} returns an \code{fsa} object with the indicated sample
##' renamed.
##' 
##' \code{fsaNormalize} returns an \code{fsa} with the rfu values
##' normalized. The values are scaled such that the sum of the rfus for
##' each sample will equal the mean of the sum for all samples. Different
##' dyes are normalized separately.
##'
##' \code{plot.fsa} is the most convenient way to display an
##' electropherogram. The \code{epn} argument allows you to pick which
##' sample you wish to plot, and the \code{raw} argument determines whether
##' you scale the lines by the raw reads (\code{raw = TRUE}), or by `base
##' pairs' (\code{raw = FALSE}, the default).
##' 
##' \code{fsa2PeakTab} converts the \code{fsa} to a peak table
##' \code{data.frame} in a layout similar to that produced by peakscanner.
##' The data.frame has three columns: sample.name, bp, and rfu. Each row
##' includes the data for a single peak detected for a single sample. The
##' rest of the chromatogram data is dropped.
##' 
##' @examples
##'
##' \dontrun{
##' ## The following examples assume fsa.data is an fsa object created by a
##' ## call to readFSA()
##' 
##' ## Remove the tenth sample
##' fsa.data.B <- fsaDrop(fsa.data, 10)
##'
##' ## Remove sample `unit3'
##' fsa.data.C <- fsaDrop(fsa.data, "unit3")
##'
##' ## Rename sample `d.x11' to `unit4'
##' fsa.data.D <- fsaRename(fsa.data, epn = "d.x11", newname = "unit4")
##'
##' ## Combine fsa objects from two different file sets
##' fsa.combined <- fsaCombine(fsa1 = block1, fsa2 = block2)
##'
##' ## Normalize the RFU data across all samples in an fsa object
##' fsa.norm <- fsaNormalize(fsa.data)
##'
##' }
##' @keywords aflp fsa genemapper peakscanner microsatellite ssr
##' 

##' @rdname editFSA
##' @export
fsaDrop <- function(fsa, epn){
  if(is.numeric(epn)){
    fsa$area <- fsa$area[-epn, ]
    fsa$ep <- fsa$ep[-epn]
  } else {
    fsa$area <- fsa$area[fsa$area$sample != epn, ]
    fsa$ep <- fsa$ep[names(fsa$ep) != epn]
  }
  return(fsa)
}

##' @rdname editFSA
##' @export
fsaRename <- function(fsa, epn, newname){
  fsa$area[epn, "sample"] <- newname
  if(is.numeric(epn)){
    names(fsa$ep)[epn] <- newname
  } else {
    names(fsa$ep)[which(names(fsa$ep) == epn)] <- newname 
  }
  return(fsa)
}

##' @rdname editFSA
##' @export
fsaCombine <- function(fsa1, fsa2){
  fsa1$area <- rbind(fsa1$area, fsa2$area)
  fsa1$ep <- c(fsa1$ep, fsa2$ep)
  return(fsa1)
}

##' @rdname editFSA
##' @export 
fsaNormalize <- function(fsa) {
  for(i in fsa$dyes){
    rfcorrection <- mean(fsa$area[,i]) / fsa$area[,i]
    for(j in seq_along(fsa$ep)){
      fsa$ep[[j]]$scans[,i] <- fsa$ep[[j]]$scans[,i] * rfcorrection[j]
    }
  }
  return(fsa)
}

## Not sure this is necessary/useful?
## ' @rdname editFSA
## ' @export
## fsaTrim <- function(fsa) {
##   fsa$dat <- fsa$dat[which(fsa$dat$peak & ! is.na(fsa$dat$bp)), ]
##   fsa
## }

##' @rdname editFSA 
##' @export
fsa2PeakTab <- function(fsa, dye){
  pt <- data.frame(sample.name = character(), bp = numeric(),
                   height = numeric())
  for(ep in fsa$ep){
    dat <- ep$scans[, c("bp", dye)]
    dat <- cbind(rep(clean.label(ep$sample), nrow(dat)), dat)
    dat <- dat[ep$peaks[[dye]], ]
    dat <- dat[complete.cases(dat), ]
    pt <- rbind(pt, dat)
  }
  pt <- pt[order(pt$bp),]
  ## I'm not sure why these names don't stick??
  names(pt) <- c("sample.name", "bp", "height")
  return(pt)
}

#######################
## Generic Functions ##
#######################

##' @rdname editFSA
##' @export
print.fsa <- function(fsa){
  cat(paste("fsa list with", length(fsa$ep), "samples\n"))
}

##' @rdname editFSA
##' @export
summary.fsa <- function(fsa){
  ## To be expanded when I decide what info I want summarized. 
  print(fsa)
  print(names(fsa$ep))
}

##' @rdname editFSA
##' @export
plot.fsa <- function(fsa, epn, ...){
  ep1 <- fsa$ep[[epn]]

  if(is.numeric(epn))
    ptitle <- names(fsa$ep)[epn]
  else
    ptitle <- epn
  
  plot(ep1, main = ptitle, ...)
}

##' @rdname editFSA
##' @export
plot.electropherogram <- function(ep, raw = FALSE, ...){
  if(raw){
    xvals <- 1:nrow(ep$scans)
    ylims <- range(ep$scans[, colnames(ep$scans) != "bp"], na.rm = TRUE)
    xlabs <- "time"
  } else {
    xvals <- ep$scans$bp
    xlabs <- "base pairs"
    ylims <-
      range(ep$scans[! is.na(ep$scans$bp), colnames(ep$scans) != "bp"],
            na.rm = TRUE) 
  }

  plot(x = xvals, y = ep$scans$standard, col = 'red', type = 'l', xlab =
         xlabs, ylim = ylims, ylab = "RFU", ...)
  for(dye in 2 + seq_along(colnames(ep$scans)[3:(ncol(ep$scans))])){
    lines(x = xvals, y = ep$scans[[dye]], col = dye)
    points(x = xvals[ ep$peaks[[colnames(ep$scans)[dye]]] ],
           y = ep$scans[[dye]][ ep$peaks[[colnames(ep$scans)[dye]]] ],
           col = dye)
  }
}
