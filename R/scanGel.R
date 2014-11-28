##' @rdname scanGel
##' @title Gel-like visualization of electropherograms, and GUI tools to
##' interactively manage bins
##' 
##' @description \code{scanGel} plots the raw data from a set of fsa files
##' in a slab-gel format, optionally including a set of bins. The GUI
##' interface provides options for users to interactively add and delete
##' bins.
##'
##' @details Electropherograms from multiple samples are collated into a
##' single gel-like image, such as would have been produced on a slab-gel
##' style instrument. If bin lines are provided, they will be overlaid on
##' the gel image. The user may then add or delete bins with the mouse.
##'
##' @param pt A peak table, as produced by \code{fsa2PeakTab}
##'
##' @param bin.lines The bin boundaries to use. Can be a two column matrix
##' or data.frame, with the first column being the lower boundary and the
##' second the upper boundary. Omit if you want to define bins completely
##' manually, otherwise use the output of \code{fsaRGbin}.
##'
##' @return This function does not return anything, but is used to launch
##' an interactive graphical viewer displaying the data in a peak table.
##' However, see below for a discussion of the object \code{editedBins},
##' which is placed directly in the global environment to 'return' the
##' modified bins to the user for further use.
##'
##' Once launched, the viewer displays each sample in a column, and the
##' y-axis indicates the size in base pairs of individual fragments. (If
##' you see an empty box, click the \code{refresh} button at the top to
##' replot the image.) Samples aren't labelled, by design, as when
##' evaluating scoring you should not be considering the identify of
##' individuals, but rather the quality of the data (you knew that,
##' right?).
##'
##' The height of the electropherogram peak for each fragment is captured
##' by both the color and width of the plotted band. Line width increases
##' with the log of the peak height, so fatter bands are higher peaks.
##' Further, color is used to indicate absolute ranges:
##'
##' \describe{
##' \item{< 50 RFU: }{cyan. Note that some very low peaks are
##' effectively invisible, although they may have been used to define bins.
##' This causes some bins to appear empty. Such peaks are filtered out
##' before data analysis, so not a reason for concern.}
##' \item{50-100 RFU: }{red}
##' \item{100-150 RFU: }{green}
##' \item{150-4000 RFU: }{black}
##' \item{4000+ RFU: }{dark blue}
##' }
##' 
##' Later in the analysis you will be able to set a threshold to eliminate
##' short peaks. Most studies reject peaks less than 100 or 150 RFU,
##' corresponding to the cyan, red and possibly green bands on the image.
##'
##' If you have provided the \code{bin.lines} argument, your bins will be
##' displayed as well. The lower boundary of a bin is indicated with a
##' dashed orange line, and the upper boundary with a dotted purple line.
##' In addition, a purple bar is placed in the margin indicating the size
##' of the bin.
##'
##' The initial view starts at 50-62 bp. You can page up and down through
##' the data with the toolbar buttons, as well as zoom the view in and out.
##'
##' @section Editing Bins:
##' 
##' Why would you want to edit the bins by hand? In short, because no
##' algorithm will provide a perfect binning of your data. Rather than
##' trust to the objectivity of an arbitrary algorithm, I choose to remove
##' bins that look to have inaccurately grouped fragments.
##'
##' More rarely, I may split a single large bin into two, or otherwise
##' manipulate the boundaries to better reflect the fragment pattern on the
##' gel. Comparing the results of my hand-edited binning vs the automated
##' binning, there is rarely a noticeable difference. But if you're
##' obsessive about your data, you may appreciate having the ability to
##' exclude the worst, least defined groups of fragments.
##'
##' In the example, scroll up to base pairs 82-85. If you feel comfortable
##' with the algorithmic binning of this section, no problem. However, you
##' may prefer to exclude this section, as the boundaries between
##' successive bins are not clear.
##' 
##' If you click the \code{Delete Bin} button on the tool bar, you will be
##' able to remove one bin from the view. Simply click anywhere in the area
##' of the bin. The entire bin will be shaded out. As you do this, the
##' modified bin list will be written to your global environment with the
##' name \code{editedBins}. If you click the \code{Refresh} button, you
##' will see the revised data set with the deleted bin removed.
##'
##' \code{Add Bin} works similarly, except that you mark your new bin by
##' dragging the mouse across the gel to define the upper and lower
##' boundaries with a temporary rectangle.
##'
##' Note: there is no undo! However, the variable you passed in to
##' \code{scanGel} as the \code{bin.lines} argument remains unchanged. So
##' you can easily return to the initial bins if you don't like the changes
##' you've made. For this reason, you should *not* use \code{editedBins} as
##' the argument for \code{scanGel}!
##'
##' On the other hand, there is no automatic saving of your new bins
##' either. So if you want to keep you modified bins for future use, be
##' sure to save \code{editedBins} to a file for safe-keeping.
##'
##' @export
##' @author Tyler Smith
##'
##' @examples
##'
##' \dontrun{
##' mybins <- fsaRGbin(oxyPT)
##' scanGel(oxyPT, mybins)
##'
##' ## You may have to click 'Refresh' to see the gel!
##' 
##' ## Any changes made to the bins in the GUI will be stored in the
##' ## variable editedBins:
##'
##' dim(mybins)
##' dim(editedBins)
##'
##' ## Delete some bins
##'
##' dim(mybins) ## unchanged
##'
##' dim(editedBins) ## deleted bins are removed
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
    if(length(selbin) > 0){
      rect(ybottom = bin.lines[selbin, 1], ytop = bin.lines[selbin, 2],
           xleft = 1, xright = length(levels(pt$sample.name)) + 1,
           col = "#00000055")
      bin.lines <<- bin.lines[ -selbin, ]
      assign("editedBins", bin.lines, envir = .GlobalEnv)
    }
    galert("There's no bin there!")
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

