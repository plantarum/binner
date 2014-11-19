##' @rdname readFSA
##' @title Read and size .fsa files
##'
##' @description \code{readFSA} reads and processes raw .fsa files into R.
##'
##' @details \code{pretrim} and \code{posttrim} are regexps, passed to
##' grep. The substring at the front of each rowname matching
##' \code{pretrim} (or the end for \code{posttrim}) is removed. To cancel
##' trimming, set these to \code{NA}.
##'
##' \code{ladder.check} In the standard ladder GS500, the 250bp fragment
##' commonly migrates at an odd rate, making it inappropriate for use in
##' sizing. Setting \code{ladder.check = 250}, which is the default, will
##' exclude this fragment from the sizing process. Set \code{ladder.check =
##' NA} if you want to use all the peaks in \code{ladder} in sizing the
##' data.
##'
##' @param files A list of fsa files to read. If NULL (the default), all
##' .fsa files in the directory specified by \code{path} will be read.
##' 
##' @param path The directory to search for \code{files}. The default is
##' the current directory.
##'
##' @param dye A vector of dyes to include when reading data. Valid values
##' include: "FAM", "VIC", "NED", "PET".
##' 
##' @param lad.channel Which .fsa data channel has the size standard ladder
##' data. The default is 105, which is the value for our system.
##' 
##' @param pretrim A regexp - text to trim off the front of the sample names.
##' 
##' @param posttrim A regexp - text to trim off the end of the sample names.
##' 
##' @param thresh For fsa files, lower limit of peak data to read. If this
##' value is less than -10, all data for the selected channels will be
##' read. It's probably not necessary to set this unless you have very
##' large sample sizes (> 100 individuals), and if you do it will interfere
##' with normalization later on.
##'
##' @param ladder A vector with the fragments present in the ladder, in
##' order. The default is the standard GS500(-250)LIZ ladder.
##'
##' @param SNR This is a cut-off value, used to exclude the primer-dimer
##' spike at the beginning of the run from being erroneously interpreted as
##' a ladder fragment. This spike is usually > 6000 rfus, and the true
##' ladder peaks are usually (always?) well below this cut-off. Not setting
##' this value may lead to slower, and poorer ladder-fitting.
##'
##' @param ladder.check If not null, the size of a ladder fragment that is
##' present but \emph{not} used for sizing. This size of this fragment will
##' be estimated, and the estimate reported during scanning. Otherwise, it
##' will be ignored. See below.
##'
##' @param sizing Currently two options are supported, "local" and "cubic".
##' "local" provides the local Southern method, identical to the one used
##' in PeakScanner et al., and recommended. "cubic" uses a cubic spline
##' function.
##'
##' @param bin.width The width in basepairs of each bin. Used to tune the
##' peak-finding algorithm of the internal function \code{get.peaks}.
##'
##' @param min.peak.height The minimum rfu value to consider a true peak,
##' passed to \code{get.peaks}. Note that you can exclude low peaks later
##' on in the process. This is preferable, because normalization isn't done
##' inside \code{readFSA}; the height of some peaks may be increased by
##' normalization, if they aren't excluded here.
##'
##' @param baseline.width The width of the window to use when 'correcting'
##' the rfu intensity. Each rfu value will be corrected by having the
##' running minimum from a window \code{baseline.width} units wide
##' subtracted from it. Without this correction, the rfu values on some
##' runs will gradually decline over the course of the run. As I understand
##' it, this is identical to the implementation in PeakScanner.
##'
##' @param verbose Do you want to see all the details scroll by or not?
##' \code{readFSA} can take a while, so this gives you something to watch
##' while you wait.
##'
##' @param smoothing This is a tuning value. If smoothing is > 1, the rfu
##' values will be converted to the running mean of the actual values, with
##' a window width of of 'smoothing'. 3 seems to work nicely and is the
##' default. 1 may be fine too. Even numbers or non-integer values may
##' break the time-space continuum (untested).
##'
##' @author Tyler Smith
##'
##' @export
##' 
##' @return \code{readFSA} returns an object of class \code{fsa}. The
##' elements include:
##'
##' \describe{
##'
##' \item{ep: }{A list of electropherogram objects, each
##' corresponding to one fsa file (see below.)} 
##'
##' \item{dyes: }{A list of the data channels (fluorescent dyes) read.}
##' 
##' \item{area: }{A data frame recording the total area under the curve
##' used for each dye/sample combination, used for normalizing results.}
##'
##' \item{error: }{If present, a vector of sample names for all samples
##' that produced unsatisfactory sizing results. Most likely bad reactions
##' that should be removed.}
##' }
##'
##' \code{electropherogram} objects have three components:
##'
##' \describe{
##'
##' \item{scans: }{A data frame, the columns of which are the heights (in
##' RFUs) of each dye, including the size standard, for each time step in
##' the capillary run. The data is ordered, with the first reads at the
##' beginning of the table. There is an additional column, `bp', which
##' stores the size, in base pairs, of each row in the table.}
##'
##' \item{peaks: }{A list of vectors, each of which contains the position
##' of the peaks for each dye in the electropherogram, in base pairs.}
##'
##' \item{sample: }{The original sample name for the fsa file.}
##' 
##' }
##'
##' @seealso \code{\link{fsaNormalize}}, \code{\link{fsa2PeakTab}},
##' \code{\link{plot.fsa}}, \code{\link{fsaRGbin}}, \code{\link{binSet}}
##' 
##' @examples
##'
##' \dontrun{
##' ## A set of fsa files are included in this package, which you can read
##' ## with the following example. For your own data replace
##' ## \code{system.file(...)} with the path to your fsa files.
##'
##' ## Read the raw files:
##' ## Pretrim and postrim are optional, and serve only to remove
##' ## extraneous components of the sample name added by the sequencing
##' ## lab. 
##' fsa.data <- readFSA(path = system.file("pp5", package = "fragread"),
##'                      pretrim = "AFLP_.*AFLP_", posttrim = "-5_Frag.*",
##'                      dye = "FAM")
##'
##' ## The print function for fsa objects doesn't do much yet:
##' fsa.data
##' summary(fsa.data)
##' 
##' ## Plot the second sample
##' plot(fsa.data, 2)
##' 
##' ## Normalize the electropherograms
##' fsa.norm <- fsaNormalize(fsa.data)
##'
##' ## Plot the second sample again, note the peak heights (y-axis) have
##' ## changed
##' plot(fsa.norm, 2)
##' 
##' ## Convert the electropherograms into a peak table
##' peaktab <- fsa2PeakTab(fsa.norm, dye = "FAM")
##' head(peaktab)
##' 
##' ## Binning:
##' bins <- fsaRGbin(peaktab)
##' aflp <- binSet(peaktab, bins, pref = "A")
##'
##' ## Extract the scoring data and proceeed with analysis:
##' mydata <- aflp[, , "alleles"]
##' 
##' }
##' @keywords aflp fsa genemapper peakscanner microsatellite ssr
readFSA <- function(files = NULL, path = "./", dye, lad.channel = 105,
                     pretrim = NA, posttrim = ".fsa", thresh = -100,
                     ladder = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350,
                         400, 450, 490, 500), SNR = 6000, ladder.check = 250,
                     sizing = "local", bin.width = 1, min.peak.height = 50, 
                     baseline.width = 51, verbose = TRUE, smoothing = 3){ 

######################################################
## Wrapper to fsa.proc, where the real work is done ##
######################################################
  
    if(is.null(files))
        files <- list.files(path, pattern = "\\.fsa$", full.names = TRUE)
    else
        files <- paste(path, files, sep = "")

    .sizing.errors <- character()
    .total.areas <- data.frame(sample = character(),
                               stringsAsFactors = FALSE) 
    
    for (i in dye){
      .total.areas <- cbind(.total.areas, numeric())
      colnames(.total.areas)[ncol(.total.areas)] <- clean.label(i)
    }

    res <- list()

    ## Originally used lapply here, but it gets complicated when you have
    ## multiple results to split out.
    for(fl in seq_along(files)){
      tmp <- fsa.proc(file = files[fl], files, dye,
                      lad.channel, pretrim, posttrim,
                      thresh, ladder, SNR, ladder.check,
                      sizing, bin.width, min.peak.height,
                      baseline.width, verbose, smoothing,
                      .total.areas, .sizing.errors)
      ## split out results here
      res[[fl]] <- tmp$result
      .total.areas <- tmp$ta
      .sizing.errors <- tmp$err
    }

    names(res) <- clean.label(sapply(res, function(x) x$sample))

    fsa <- list(ep = res)
    ## Threshold won't work anymore, and probably isn't needed.
    ## if (thresh > -10) res <- subset(res, height > thresh)

    if(length(.sizing.errors) > 0){
        message(length(.sizing.errors), " sizing errors: ")
        message(paste("\"", paste(.sizing.errors, collapse = "\", \""), "\"", sep = ""))
    }
    
    fsa$area <- .total.areas
    fsa$errors <- .sizing.errors
    fsa$dyes <- dye
    class(fsa) <- "fsa"
    return(fsa)
}

fsa.proc <- function(file, files, dye, lad.channel, pretrim, posttrim, thresh,
                     ladder, SNR, ladder.check, sizing, bin.width, min.peak.height, 
                     baseline.width, verbose, smoothing, .total.areas,
                     .sizing.errors) { 
  ## This function started out as an anonymous function passed as an argument to do.call(rbind...)
  ## from inside readFSA. This made it really difficult to debug. So I've pulled it out into it's
  ## own function now. Note that it uses the variables .sizing.errors and .total.areas to pass
  ## additional information back to readFSA.
  
  abif <- read.abif(file)
  tag <- tag.trimmer(basename(file), pretrim, posttrim)
  if (verbose) message("\n", which(files == file), "/", length(files), ": ", tag)
  lad.dat <- abif$Data[[paste('DATA.', lad.channel, sep='')]]
  lad.dat <- lad.dat - runmin(lad.dat, baseline.width)
  if(smoothing > 1)
    lad.dat <- runmean(lad.dat, k = smoothing)

  scans <- data.frame(standard = lad.dat)
  
  tmp <- set.ladder(lad.dat, ladder, SNR, ladder.check, verbose)
  if(tmp$val < 0.9999 & tmp$val > 0.99){
    message("re-sizing with lower peak threshold!")
    tmp <- set.ladder(lad.dat, ladder, SNR, ladder.check, verbose, bad.size = TRUE)
  }
  
  scans$bp <- tmp$bp                   # ladder added to res1
  val <- tmp$val                       # how successful was the ladder fit?

  lad.mat <- cbind(bp = scans$bp[!is.na(scans$bp)], time = which(!is.na(scans$bp)))

  dyeNames <- c(abif$Data$DyeN.1, abif$Data$DyeN.2, abif$Data$DyeN.3,
                abif$Data$DyeN.4)

  dyeNames <- clean.dye(dyeNames)       # 6-FAM -> FAM
  
  sig.channel <-
    c("DATA.1", "DATA.2", "DATA.3", "DATA.4")[dyeNames %in% dye] 
  
  dyeSafe <- clean.label(dyeNames)

  tot.area <- data.frame(sample = clean.label(tag))
  
  for (i in seq_along(sig.channel)) {
    sc <- sig.channel[i]
    dy <- dyeSafe[i]
    chan.dat <- abif$Data[[sc]]
    chan.dat <- chan.dat - runmin(chan.dat, baseline.width)
    if(smoothing > 1)
      chan.dat <- runmean(chan.dat, k = smoothing)
    scans[[dy]] <- chan.dat
    tot.area[[dy]] <- sum(chan.dat)
  }

  .total.areas <- rbind(.total.areas, tot.area)
  
  if(! is.null(ladder.check)) {
    ctime <- lad.mat[, "time"][lad.mat[, "bp"] != ladder.check]
    cbp <- lad.mat[, "bp"][lad.mat[, "bp"] != ladder.check]
    check <- lad.mat[, "time"][lad.mat[, "bp"] == ladder.check]
    scans$standard[which(scans$bp == ladder.check)] <- NA
  } else {
    ctime <- lad.mat[, "time"]
    cbp <- lad.mat[, "bp"]
  }

  if(sizing == "cubic") {
    calibrate <- splinefun(ctime, cbp)  # cubic spline interpolation

    if(! is.null(ladder.check))
      message(paste("  ladder.check: ", ladder.check, ", cubic sized at: ",
                    round(calibrate(check), 2), sep = ""))
    
    scans$bp[is.na(scans$bp)] <- calibrate(which(is.na(scans$bp)))
  } else if (sizing == "local") {
    scans$bp[is.na(scans$bp)] <-
      size.frags((1:length(lad.dat))[is.na(scans$bp)], data.frame(time = ctime, bp = cbp))
    if(! is.null(ladder.check)) {
      message(paste("ladder.check: ", ladder.check, ", local Southern sized at: ",
                    round(size.frags(check,
                                     data.frame(time = ctime, bp = cbp)), 2),
                    sep = ""))
    }
  } else
    stop("invalid sizing selected!")

  if(val < 0.9999) {
    message("  !! sizing error: ", tag, "!!")
    .sizing.errors <- c(.sizing.errors, tag)
  }

  result <- list()
  result[["scans"]] <- scans
  result[["peaks"]] <-
    get.peaks(scans, bin.width, min.peak.height = min.peak.height)
  result[["sample"]] <- tag
  class(result) <- "electropherogram"
  return(list(result = result, ta = .total.areas, err = .sizing.errors))
}

get.peaks <- function(scans, bin.width = 1, min.peak.height){
  peakList <- list()
  for(i in names(scans)[-1:-2]){         # exclude standard and bp
    samp <- scans[, i]
    ## First determine how many units per bin width
    ## User provides the desired bin width in base pairs. Convert that to raw time steps
    bp <- scans$bp[! is.na(scans$bp)]
    bin.width <- ceiling(bin.width / ((max(bp) - min(bp)) / length(bp)))
    
    ## now identify local maxima for each bin.width size window
    lmax <- runmax(samp, k = bin.width)
            
    ## select matching peaks, dropping those below minimum threshold
    peaks <- which(samp == lmax & samp > min.peak.height)
    peakList[[i]] <- peaks
  }
  
  return(peakList)
}

clean.label <- function(x){
  ## Tidy labels for use as column names in data frames:
  ## - prepend letters to the front if they start with a number
  ## - convert all punctuation to .
  x <- gsub("[[:punct:]]", ".", x)
  numInd <- which(grepl(pattern = "[0-9]", substring(x, 0, 1)))
  x[numInd] <- paste("d", x[numInd], sep = "")
  x
}

clean.dye <- function(x){
  ## Tidy labels for use as column names in data frames:
  ## - prepend letters to the front if they start with a number
  ## - convert all punctuation to .
  x <- gsub("[[:punct:]]", "", x)
  x <- gsub("[0-9]", "", x)
  x
}

tag.trimmer <- function(x, pretrim = NA, posttrim = NA) {
  ## utility to remove pretrim and postrim from sample names
  if(! is.na(pretrim)) {
    x <- sub(paste("^", pretrim, sep = ""), "", x)
  }
  if(! is.na(posttrim)){
    x <- sub(paste(posttrim, "$", sep = ""), "", x)
  }
  x
}

set.ladder <- function(lad.dat, ladder, SNR, ladder.check = NULL, verbose = TRUE,
                       bad.size = FALSE) {
  if(verbose) message("set.ladder ->")

  ## Adapted from the AFLP R package of
  ## http://r-forge.r-project.org/R/?group_id=1027 
  n <- 10 * length(ladder)

  bp <- rep(NA, length(lad.dat))
  if(bad.size)
    n <- n + 50  # We had a bad sizing the first time through, so try again
                 # with a lowered peak threshold.
  
  ## Peak takes the vector of height readings, and divides it into sections
  ## above and below a threshold value.  Each section gets it's own unique
  ## number: 
  Peak <-
    cumsum(abs(c(0, diff(lad.dat >
                           quantile(lad.dat, 1 - n / length(lad.dat))))))
  
  ## Peak <- peakclean(Peak)               # not ready for prime time
  
  ## Next, pull out the even-numbered sections of Peak. These will be sections that have
  ## heights above the threshold. That means there will be a peak in that section.

  ## To identify the actual peak, find the value in each section that is equal to the
  ## maximum value for that section.
  
  Index <-
    seq_len(length(lad.dat))[Peak %% 2 == 1][lad.dat[Peak %% 2 == 1] == 
                                          ave(lad.dat[Peak %% 2 == 1],
                                              Peak[Peak %% 2 == 1], FUN = max)]    

  ## Added by Tyler: ignore all peaks at or below the primer peak. This is
  ## tuned by the SNR argument. SNR in this case is a height threshold -
  ## anything higher than this is not a real peak, but a primer peak. We
  ## almost always get one of these somewhere around 30 bp, and they're
  ## usually well over 6000 rfu. Real peaks are well under this value. So I
  ## ignore this high peak. I also ignore everything smaller (fewer bp).
  ## The primer peak often has a few spurious shoulder peaks, which confuse
  ## the ladder fitting algorithm.
  if(max(lad.dat) > SNR){
    if (verbose) message("removing primer peak")
    primer.peak <-
      tail(which(lad.dat[Index] == max(lad.dat)), 1) 
    Index <- Index[-1 * (1:primer.peak)] # exclude the primer peak and
                                        # everything smaller 
  }

  
  while(length(Index) < length(ladder)){
    ## From the AFLP package: if we don't have enough peaks, lower the
    ## threshold and try again. Originally they lowered the threshold by
    ## setting
    ## n <- n + 1,
    ## but I  find this works better with bigger steps, so:
    n <- n + 10
    Peak <- cumsum(abs(c(0, diff(lad.dat > quantile(lad.dat, 1 - n /
                                                             length(lad.dat))))))
    ## Peak <- peakclean(Peak) # not ready for prime time!
    Index <-
      seq_len(length(lad.dat))[Peak %% 2 == 1][lad.dat[Peak %% 2 == 1]
                                         == ave(lad.dat[Peak %% 2 == 1],
                                                        Peak[Peak %%
                                                             2 == 1],
                                                        FUN = max)]  
    if(max(lad.dat) > SNR){
      if (verbose) message("removing primer peak")
      primer.peak <- which(lad.dat[Index] == max(lad.dat))
      Index <- Index[-1 * (1:primer.peak)] # exclude the primer peak and everything smaller
    }

  }
  if(length(Index) > length(ladder)){
    ## If the there are more peaks than steps on the ladder, generate a matrix of all
    ## possible combinations of the right number of peaks, and pick the one with the best
    ## fit: 
    if(length(Index) - length(ladder) == 1){
      Index <- Index[-which.max(sapply(seq_along(Index), function(i){
        summary(lm(ladder ~ stats::poly(Index[-i], 2)))$r.squared
      }))]
    } else {
      toTry <- combn(length(Index), length(Index) - length(ladder))

      ## From the AFLP package - they used 1000 combinations as the maximum
      ## number they would try. This is a bit low - with a 15-step ladder,
      ## excluding the 250 bp fragment, that means if you get just three
      ## additional fragments you will exclude nearly 2/3 of the possible
      ## combos. It's rare, but not impossible, to get three additional
      ## peaks associated with the primer peak. I noticed this when I found
      ## a reaction that would pass the peak sizing most of the time, but
      ## occasionally would fail. Very unsettling to get different results
      ## from the same data!
      if(ncol(toTry) > 4000){
        toTry <- toTry[, sample(ncol(toTry), 4000)]
      }
      Index <- Index[-toTry[, which.max(sapply(seq_len(ncol(toTry)), function(i){
        summary(lm(ladder ~ stats::poly(Index[-toTry[, i]], 2)))$r.squared
      }))]]
    }
  }

  bp[Index] <- ladder

  if(! is.null(ladder.check)) {
    ## remove ladder check peaks before evaluating success. For the ladder
    ## we use, the 250 bp peak is known to be migrate unreliably. So I use
    ## it to check the consistency between runs, but not for actually
    ## sizing the data. In practice, a good ladder fit results in the '250'
    ## bp ladder fragment being sized around 246 bp. This is just an
    ## observation, not a rule of thumb - the size of this fragment may
    ## vary on different equipment!
    lad.ind <- which(ladder == ladder.check)
    ladder <- ladder[-lad.ind]
    Index <- Index[-lad.ind]
  }

  fit.value <- summary(lm(ladder ~ stats::poly(Index, 2)))$r.squared
  
  if(verbose) {
    message("  ladder fit r2: ", round(fit.value, 5))
    ##message(paste(Index, collapse = " "))
  }

  return(list(bp = bp, val = fit.value))
}

peakclean <- function(Peak){

  ## This is an experimental function. I'm not sure it's worth the effort to properly sort out
  ## this mess. It is only necessary for bad reactions, reactions that will likely fail even if
  ## you sort out the split peaks issue.
  
  ## Sometimes the ladder peaks will be split - two separate peaks very close together. This
  ## can cause problems if the split peaks are higher than a true peak that is below the
  ## current threshold. So we need to collapse the peaks together if they are not separated by
  ## sufficient distance. The ladders fragments are at least 9-10 bp apart, so setting the
  ## minimum distance to 50 (which is about 4-5bp) should eliminate these double peaks.

  ## Loop through all valley regions that are too short, and collapse the short valley and the
  ## following peak into the preceding peak.

  ## The vector of short valleys:
  shortvalleys <- unique(Peak[Peak%%2 == 0])[table(Peak[Peak%%2 == 0]) < 50]

  for (shortvalley in shortvalleys){
    
    ## Collapse the short valley into the previous peak:
    Peak[Peak == shortvalley] <- shortvalley - 1
    
    ## Collapse the next peak into the previous peak:
    Peak[Peak == shortvalley + 1] <- shortvalley - 1
    
    ## Move all subsequent peaks and valleys back 2:
    Peak[Peak > shortvalley + 1] <-  Peak[Peak > shortvalley + 1] - 2
  }
  return(Peak)
}


## I'll just leave this here for my use. You shouldn't need to use it. It
## converts from data from the format used in an older, unpublished version
## of binner to the S3 classes described here.
fsaConvert <- function(x){
  res <- structure(list(), class = "fsa")
  dfRows <- length(subset(x$dat, tag == levels(x$dat$tag)[1] &
                            chan == "standard")$height)
  for(tg in levels(x$dat$tag)){
    print(tg)
    res[[tg]] <- list()
    res[[tg]][["peaks"]] <- list()
    res[[tg]][["area"]] <- list()
    dat <- data.frame(standard = numeric(dfRows))
    for(ch in levels(x$dat$chan)) {
      sb <- subset(x$dat, tag == tg & chan == ch)
      dat[[ch]] <- sb$height
      if(ch != "standard") {
        res[[tg]][["peaks"]][[ch]] <-
          sb$bp[which(sb$peak == TRUE)]
        ## Jeebus, hold onto your hat
        res[[tg]][["peaks"]][[ch]] <-
          res[[tg]][["peaks"]][[ch]][!is.na(res[[tg]][["peaks"]][[ch]])]
      } else {
        res[[tg]][["peaks"]][[ch]] <-
          sb$bp[! is.na(sb$bp)]
      }
    }
    ## I apologise for the next line.
    dat[["bp"]] <- subset(x$dat, tag == tg &
                            chan == levels(x$dat$chan)[levels(x$dat$chan)
                              != "standard"][1])$bp
    ## WTF right? I needed to extract the bp values for one channel, but it
    ## can't be the standard channel. However, I don't know ahead of time
    ## which channels are available. So I select the first channel from the
    ## available channels, after filtering out the standard channel.
    res[[tg]][["scans"]] <- dat
  }
  return(res)
}

