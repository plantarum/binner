##' @rdname scanGel
##' @title Gel-like visualization of electropherograms, and tools to
##' interactively manipulate bin boundaries
##' 
##' @description \code{scanGel} plots the raw data from an fsa file.
##'
##' @details This is a placeholder
##' @param pt A peak table, as produced by \code{fsa2PeakTab}
##'
##' @param bin.lines The bin boundaries to use. Can be a two column matrix
##' or data.frame, with the first column being the lower boundary and the
##' second the upper boundary. Omit if you want to define bins completely
##' manually, otherwise use the output of \code{fsaRGbin}.
##'
##' @return Another placeholder
##' @export
##' @author Tyler Smith
##'
##' @example
##'
##' \dontrun{
##' mybins <- fsaRGbin(oxyPT)
##' scanGel(oxyPT, mybins)
##'
##' ## Any changes made to the bins in the GUI will be stored in the
##' ## variable editedBins:
##'
##' dim(mybins)
##' dim(editedBins)
##'
##' ## If you want to store your changes, make sure to save editedBins to a
##' ## csv file or .Rdata object!!
##'
##' }
scanGel <- function(pt, bin.lines = numeric()){

  editedBins <- bin.lines
  
  win <- gwindow("gel window", visible = FALSE, expand = TRUE, fill = TRUE)
  topgrp <- ggroup(horizontal = FALSE, container = win)
  plotwin <- ggraphics(container = topgrp)

  pageUpAction <- function(h, ...) {
    .plower <<- .plower + .pstep
    .pupper <<- .pupper + .pstep
    plot.gel(pt, bin.lines) 
  }

  pageDownAction <- function(h, ...) {
    .plower <<- .plower - .pstep
    .pupper <<- .pupper - .pstep
    plot.gel(pt, bin.lines)
  }

  zoomInAction <- function(h, ...) {
    inc <- .pstep * 0.2
    .pstep <<- .pstep - inc
    .plower <<- .plower + 0.5 * inc
    .pupper <<- .pupper - 0.5 * inc
    plot.gel(pt, bin.lines)
  }

  zoomOutAction <- function(h, ...) {
    inc <- .pstep * 1.25
    .pstep <<- .pstep + inc
    .plower <<- .plower - 0.5 * inc
    .pupper <<- .pupper + 0.5 * inc
    plot.gel(pt, bin.lines)
  }

  resetZoomAction <- function(h, ...) {
    .pstep <<- 10
    .pupper <<- .plower + 12
    plot.gel(pt, bin.lines)
  }

  refreshAction <- function(h, ...) {
    plot.gel(pt, bin.lines)
  }
  
  deleteBinAction <- function(h, ...){
    sel <- locator(1)
    selbin <- which(bin.lines[,1] < sel$y & bin.lines[,2] > sel$y)
    rect(ybottom = bin.lines[selbin, 1], ytop = bin.lines[selbin, 2],
         xleft = 1, xright = length(levels(pt$sample.name)) + 1,
         col = "#00000055")
    bin.lines <<- bin.lines[ -selbin, ]
    assign("editedBins", bin.lines, envir = .GlobalEnv)
  }

  addBinAction <- function(h, ...){
    unblockHandler(plotwin, addBin)
  }

  addBin <- addHandlerChanged(plotwin, handler = function(h, ...){
    bintop <- h$y[2]
    binbottom <- h$y[1]
    if(sum(bin.lines <= bintop & bin.lines >= binbottom) > 0){
      galert("invalid bin, pick again!")
    } else {
      rect(ybottom = binbottom, ytop = bintop,
           xleft = 1, xright = length(levels(pt$sample.name)) + 1,
           col = "#99000055")
      bin.lines <<- rbind(bin.lines, h$y)
      bin.lines <<- bin.lines[order(bin.lines[,1]), ]
      assign("editedBins", bin.lines, envir = .GlobalEnv)
      blockHandler(plotwin, addBin)       # return to blocked when done
    }
  })

  blockHandler(plotwin, addBin)         # block until button clicked!
  
  actionList <- list(
    page_up = gaction(label = "Page Up", icon = "Page Up", handler =
                        pageUpAction, parent = win),
    page_down = gaction(label = "Page Down", icon = "Page Down", handler =
                        pageDownAction, parent = win),
    zoom_in = gaction(label = "Zoom In", icon = "Zoom In", handler =
                        zoomInAction, parent = win),
    zoom_out = gaction(label = "Zoom Out", icon = "Zoom Out", handler =
                        zoomOutAction, parent = win),
    reset_zoom = gaction(label = "Reset Zoom", icon = "Reset Zoom",
      handler = resetZoomAction, parent = win), 
    sep1 = gseparator(),
    delete_bin = gaction(label = "Delete Bin", icon = "Delete Bin",
      handler = deleteBinAction),
    add_bin = gaction(label = "Add Bin", icon = "Add Bin",
      handler = addBinAction),
    sep2 = gseparator(),
    refresh_action = gaction(label = "Refresh", icon = "Refresh",
      handler = refreshAction)
  )
  
  tool_bar <- gtoolbar(actionList, cont = win)
  
  .plower <- 50
  .pupper <- 62
  .pstep <- 10
  
  ## ID <- addHandlerExpose(plotwin, handler = function(h, ...){
  ##   message("exposed")
  ##   plot.gel(pt, bin.lines)
  ## })

  plot.gel <- function(pt, bin.lines = NA){
    
    samples <- length(levels(pt$sample.name))
    ## These two lines necessary only for coloring the replicate lanes,
    ## which needs more thought.
    ## id <- sub("REP", "", levels(pt$sample.name))
    ## reps <- id[duplicated(id) | duplicated(id, fromLast = TRUE)]

    par(mar = c(2, 2, 2, 2))
    plot(x = 1, ylim = c(.plower, .pupper), xlim = c(1, samples + 1), type = 'n',
         yaxp = c(round(.plower, 0), round(.pupper, 0),
           round(.pupper, 0) - round(.plower, 0)))
    abline(v = 0:500, col = "grey", lty = 2)
    abline(h = 0:500, col = "black", lty = 1)
    abline(h = .plower, lwd = 3)
    abline(h = .pupper, lwd = 3)

    ## Not sure this makes sense, highlighting lanes that contain replicate
    ## samples?
    
    ## abline(v = which(id %in% reps), col = 'red', lty = 4)
    
    heights <- pt[, "height"]
    hcol <- heights
    hcol[heights < 50] <- 5
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
    
    if((is.matrix(bin.lines) | is.data.frame(bin.lines))
       & nrow(bin.lines) > 0){
      abline(h = bin.lines[, 1], col = "orange", lty = 2, lwd = 2)
      abline(h = bin.lines[, 2], col = "purple", lty = 3, lwd = 2)
      segments(x0 = 0.9, x1 = 0.9, col = "purple", lwd = 4,
               y0 = bin.lines[, 1], y1 = bin.lines[, 2])
      segments(x0 = samples + 1.1, x1 = samples + 1.1, lwd = 4, col = "purple",
               y0 = bin.lines[, 1], y1 = bin.lines[, 2])
    }
    
    enabled(actionList$zoom_in) <- .pstep > 5
    enabled(actionList$zoom_out) <- .pstep < 200
    enabled(actionList$page_down) <- .plower > 0
    enabled(actionList$page_up) <- .pupper < 600
    enabled(actionList$delete_bin) <- ! is.na(bin.lines)
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

