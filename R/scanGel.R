##' @rdname scanGel
##' @title Gel-like visualization of electropherograms, and tools to
##' interactively manipulate bin boundaries
##' 
##' @description \code{scanGel} plots the raw data from an fsa file.
##'
##' @details This is a placeholder
##' @param pt A peak table, as produced by \code{fsa2PeakTab}
##'
##' @param bin.lines The bin boundaries to use. Can be a two column table,
##' with the first column being the lower boundary and the second the upper
##' boundary; or a vector, with odd elements lower boundaries and even
##' elements upper boundaries. Omit if you want to define bins completely
##' manually, otherwise use the output of \code{fsaRGbin}.
##'
##' @return Another placeholder
##' @export
##' @author Tyler Smith
##' 
scanGel <- function(pt, bin.lines = numeric()){
  if(class(bin.lines) %in% c("matrix", "data.frame"))
    bin.lines = sort(c(bin.lines[,1], bin.lines[,2]))

  win <- gwindow("gel window", visible = FALSE, expand = TRUE, fill = TRUE)
  topgrp <- ggroup(horizontal = FALSE, container = win)
  buttongrp <- ggroup(horizontal = TRUE, container = topgrp)
  scrollUp <- gbutton("scroll up", container = buttongrp)
  scrollDown <- gbutton("scroll down", container = buttongrp)
  plotwin <- ggraphics(container = topgrp)

  .plower <- 50
  .pupper <- 62
  .pstep <- 10
  
  addHandlerChanged(scrollUp, handler = function(h, ...) {
    .plower <<- .plower - .pstep
    .pupper <<- .pupper - .pstep
    plot.gel(pt, bin.lines) 
  })

  addHandlerChanged(scrollDown, handler = function(h, ...) {
    .plower <<- .plower + .pstep
    .pupper <<- .pupper + .pstep
    plot.gel(pt, bin.lines) 
  })

  ID <- addHandlerExpose(plotwin, handler = function(h, ...){
    plot.gel(pt, bin.lines)
  })

  plot.gel <- function(pt, bin.lines = numeric()){

    if("bin" %in% colnames(pt)) {
      bins <- subset(pt[,3:4], !duplicated(bin) |
                       c(!duplicated(bin)[-1], FALSE))
    } else {
      bins <- NA
    }
    
    samples <- length(levels(pt$sample.name))
    id <- sub("REP", "", levels(pt$sample.name))
    reps <- id[duplicated(id) | duplicated(id, fromLast = TRUE)]

    par(mar = c(2, 2, 2, 2))
    plot(x = 1, ylim = c(.plower, .pupper), xlim = c(1, samples + 1), type = 'n',
         yaxp = c(round(.plower, 0), round(.pupper, 0),
           round(.pupper, 0) - round(.plower, 0)))
    abline(v = 0:500, col = "grey", lty = 2)
    abline(h = 0:500, col = "black", lty = 1)
    abline(h = .plower, lwd = 3)
    abline(h = .pupper, lwd = 3)

    if(! is.na(bins))
      abline(h = bins$Size, lty = c(2, 3), col = c("blue", "purple"))

    abline(v = which(id %in% reps), col = 'red', lty = 4)
    
    heights <- pt[, "height"]
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
    
    invisible(list(pt = pt, bin.lines = bin.lines)) 
  }
  
  visible(win) <- TRUE
  ## geldev <- dev.cur()
  ## if(double.plot) {
  ##   dev.new()
  ##   fsa.plot(fsa, dye, sample = sample(1:length(levels(fsa$dat$tag)), min(maxsamples,
  ##                      length(levels(fsa$dat$tag)))), bp =
  ##            c(lower, upper), ylim = ylim)
  ##   abline(v = bin.lines, col = c("orange", "purple"))
  ##   sampledev <- dev.cur()
  ##   dev.set(geldev)
  ## }
  ## oldlow <- lower
  ## res <- place.bins(tmp)
  ## while(!is.numeric(res)){
  ##   ## if(double.plot & oldlow != res$lower) {
  ##   ##   dev.set(sampledev)
  ##   ##   fsa.plot(fsa, dye, sample = sample(1:length(levels(fsa$dat$tag)), min(maxsamples,
  ##   ##                      length(levels(fsa$dat$tag)))), bp =
  ##   ##            c(res$lower, res$upper), ylim = ylim)
  ##   ##   abline(v = res$bin.lines, col = c("orange", "purple"))
  ##   ##   dev.set(geldev)
  ##   ##   oldlow <- res$lower
  ##   ## }
  ##   res <- place.bins(res)
  ## }
  
  ## return(matrix(res, ncol = 2, byrow = TRUE))
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
        tmp <- plot.gel(bv$pt, bv$lower + 1, bv$upper - 1, bv$step - 2, bv$bin.lines)
        return(tmp)
      } else {
        tmp <- plot.gel(bv$pt, bv$lower, bv$upper, bv$step, bv$bin.lines)
        return(tmp)
      }
    }
    else if (grconvertY(command$y, "user", "npc") <= 0.5) {
      tmp <- plot.gel(bv$pt, bv$lower - 1, bv$upper + 1, bv$step + 2, bv$bin.lines)
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
    tmp <- plot.gel(bv$pt, bv$lower - bv$step, bv$upper - bv$step, bv$step, bv$bin.lines)
    return(tmp)
  }
  else if(command$y > bv$upper){
    tmp <- plot.gel(bv$pt, bv$lower + bv$step, bv$upper + bv$step, bv$step, bv$bin.lines)
    return(tmp)
  }
  else {
    abline(h = command$y, col = "green")
    text(x = command$x, y = command$y, labels = round(command$y, 3))
    bv$bin.lines <- sort(c(bv$bin.lines, command$y))
    return(bv)
  }
}  

