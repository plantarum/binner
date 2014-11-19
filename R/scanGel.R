##' @rdname scanGel
##' @title Gel-like visualization of electropherograms, and tools to
##' interactively manipulate bin boundaries
##' 
##' @description \code{scanGel} plots the raw data from an fsa file.
##'
##' @details This is a placeholder
##' @param pt A peak table, as produced by \code{fsa2PeakTab}
##'
##' @param dye Which dye to read when converting to a Peakscanner table.
##' Valid values include "FAM", "VIC, "NED" and "PET".
##' 
##' @param lower Lower limit, in base pairs, of the plot
##' 
##' @param upper Upper limit, in base pairs, of the plot
##'
##' @param step The size of the step, in base pairs, to make when scrolling
##' the plot
##'
##' @param bin.lines The bin boundaries to use. Can be a two column table,
##' with the first column being the lower boundary and the second the upper
##' boundary; or a vector, with odd elements lower boundaries and even
##' elements upper boundaries. Omit if you want to define bins completely
##' manually, otherwise use the output of \code{fsaRGbin}.
##'
##' @param use.bins 
##'
##' @param fsa 
##'
##' @param maxsamples
##'
##' @param ylim
##'
##' @return Another placeholder
##' @export
##' @author Tyler Smith
scanGel <- function(pt, dye, lower = 50, upper = lower + 12, step = 10, bin.lines = numeric(),
                     use.bins = FALSE, fsa = NULL, maxsamples = 50, ylim = c(0, 3500)){
  if(class(bin.lines) %in% c("matrix", "data.frame"))
    bin.lines = sort(c(bin.lines[,1], bin.lines[,2]))
  
  tmp <- plot.gel(pt, dye, lower, upper, step, bin.lines, use.bins)
  geldev <- dev.cur()
  ## if(double.plot) {
  ##   dev.new()
  ##   fsa.plot(fsa, dye, sample = sample(1:length(levels(fsa$dat$tag)), min(maxsamples,
  ##                      length(levels(fsa$dat$tag)))), bp =
  ##            c(lower, upper), ylim = ylim)
  ##   abline(v = bin.lines, col = c("orange", "purple"))
  ##   sampledev <- dev.cur()
  ##   dev.set(geldev)
  ## }
  oldlow <- lower
  res <- place.bins(tmp)
  while(!is.numeric(res)){
    ## if(double.plot & oldlow != res$lower) {
    ##   dev.set(sampledev)
    ##   fsa.plot(fsa, dye, sample = sample(1:length(levels(fsa$dat$tag)), min(maxsamples,
    ##                      length(levels(fsa$dat$tag)))), bp =
    ##            c(res$lower, res$upper), ylim = ylim)
    ##   abline(v = res$bin.lines, col = c("orange", "purple"))
    ##   dev.set(geldev)
    ##   oldlow <- res$lower
    ## }
    res <- place.bins(res)
  }
  
  return(matrix(res, ncol = 2, byrow = TRUE))
}

plot.gel <- function(pt, dye, lower = 50, upper = lower + 6, step = 10, bin.lines = numeric(),
                     use.bins = FALSE){

  if("bin" %in% colnames(pt)) {
    bins <- subset(pt[,3:4], !duplicated(bin) | c(!duplicated(bin)[-1], FALSE))
  } else {
    bins <- NA
  }
     
  if (use.bins & (length(bin.lines) == 0)) {
    bin.lines <- bins$Size
  }

  samples <- length(levels(pt$sample.name))
  id <- sub("REP", "", levels(pt$sample.name))
  reps <- id[duplicated(id) | duplicated(id, fromLast = TRUE)]

  par(mar = c(2, 2, 2, 2))
  plot(x = 0, ylim = c(lower, upper), xlim = c(0, samples + 1), type = 'n', yaxp =
       c(round(lower, 0), round(upper, 0), round(upper, 0) - round(lower, 0)))
  abline(v = 0:500, col = "grey", lty = 2)
  abline(h = 0:500, col = "black", lty = 1)
  abline(h = lower, lwd = 3)
  abline(h = upper, lwd = 3)

  if(! is.na(bins))
    abline(h = bins$Size, lty = c(2, 3), col = c("blue", "purple"))

  abline(v = which(id %in% reps), col = 'red', lty = 4)
  
  heights <- pt[, dye]
  hcol <- heights
  hcol[heights < 50] <- 0
  hcol[heights < 100 & heights >= 50] <- 2
  hcol[heights < 150 & heights >= 100] <- 3
  hcol[heights >= 150 & heights < 4000] <- 1
  hcol[heights >= 4000] <- 4

  hlen <- 0.2 * (log(heights) - 3)

  ## hcol <- 1
  ## heights <- log(pt[, dye]) * 0.5
  ## hlen <- 1/12
  
  segments(x0 = unclass(pt$sample.name), x1 = unclass(pt$sample.name) + 1,
           y0 = pt$bp, y1 = pt$bp, col = hcol, lwd = 12 * hlen)
  
  abline(h = bin.lines, col = c("orange", "purple"))
  
  invisible(list(pt = pt, dye = dye, lower = lower, upper = upper,
                 step = step, bin.lines = bin.lines)) 
}

place.bins <- function(bv){
  command <- locator(1)
  if(grconvertX(command$x, "user", "npc") < 0){
    if(length(bv$bin.lines) %% 2 != 0)
      message("odd number of bin boundaries!!")
    return(sort(bv$bin.lines))
  }
  if(grconvertX(command$x, "user", "npc") > 1){
    if(grconvertY(command$y, "user", "npc") > 0.5) {
      if(bv$lower < bv$upper - 3){
        if(bv$step < 3)
          bv$step <- 3
        tmp <- plot.gel(bv$pt, bv$dye, bv$lower + 1, bv$upper - 1, bv$step - 2, bv$bin.lines)
        return(tmp)
      } else {
        tmp <- plot.gel(bv$pt, bv$dye, bv$lower, bv$upper, bv$step, bv$bin.lines)
        return(tmp)
      }
    }
    else if (grconvertY(command$y, "user", "npc") <= 0.5) {
      tmp <- plot.gel(bv$pt, bv$dye, bv$lower - 1, bv$upper + 1, bv$step + 2, bv$bin.lines)
      return(tmp)
    }
  }
  else if(command$x < 0 & length(bv$bin.lines) > 0) {
    diffs <- abs(bv$bin.lines - command$y)
    if (min(diffs) < 0.25) {
      abline(h = bv$bin.lines[which(diffs == min(diffs))], col = "red", lwd = 2)
      bv$bin.lines <- bv$bin.lines[-which(diffs == min(diffs))]
    }
      return(bv)
  }
  else if(command$y < bv$lower){
    tmp <- plot.gel(bv$pt, bv$dye, bv$lower - bv$step, bv$upper - bv$step, bv$step, bv$bin.lines)
    return(tmp)
  }
  else if(command$y > bv$upper){
    tmp <- plot.gel(bv$pt, bv$dye, bv$lower + bv$step, bv$upper + bv$step, bv$step, bv$bin.lines)
    return(tmp)
  }
  else {
    abline(h = command$y, col = "green")
    text(x = command$x, y = command$y, labels = round(command$y, 3))
    bv$bin.lines <- sort(c(bv$bin.lines, command$y))
    return(bv)
  }
}  

