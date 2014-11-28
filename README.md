# binner

`binner` is an R package for processing DNA fingerprint data. At present, it
provides a complete workflow for AFLP analysis:

1. reading ABI `.fsa` files
2. normalizing electropherograms
2. identifying and sizing peaks
4. viewing individual electropherograms
5. dropping, adding and renaming samples
6. Automated peak-binning, using the `RawGeno` algorithm
7. Generating presence-absence matrices for further analysis in R (or
   export for use in other programs, if you like)

## Installation

You can install `binner` using Hadley Wickham's
[devtools](https://github.com/hadley/devtools) package.

    install.packages("devtools")
    install_github("plantarum/binner")

## Help

There is no vignette for the package. I have provided a reasonably complete
example in the help file for the function `readFSA`. Open it from with R
with `?readFSA` (after you've loaded binner of course).

## Alternatives
There are two R-based alternatives you may like to consider as well:
- [RawGeno](http://sourceforge.net/projects/rawgeno/)
- [AFLP](https://r-forge.r-project.org/projects/aflp/)

### [RawGeno](http://sourceforge.net/projects/rawgeno/)

`RawGeno` currently requires the use of a second program,
[PeakScanner](https://www.lifetechnologies.com/order/catalog/product/4381867),
to read the raw fsa data. `binner` is a bit slower than `PeakScanner` for
reading fsa files, but implements the same sizing algorithm (local
southern). Parts of `binner` may make their way into the next version of
`RawGeno`, and the current version of `binner` has pilfered some key bits of
`RawGeno` (most notably the binning algorithm).

If working entirely in R is useful to you, `binner` is currently preferable
to `RawGeno`. `binner` also provides a GUI for checking and editing your bins.
Beyond that it's a matter of preference.

### [AFLP](https://r-forge.r-project.org/projects/aflp/)

`AFLP` provides a more sophisticated workflow than `binner` or `RawGeno`, and
provides a very powerful system for controlling for variation between
capillaries, sequencers, sequencing runs etc. However, it seems to use a
different sizing system, making it hard to make direct comparisons between
electropherograms processed with `AFLP` and those using the
`PeakScanner`/`binner` approach. 
