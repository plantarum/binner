## Code follows the explanation from Genographer:
## http://home.cc.umanitoba.ca/~psgendb/birchhomedir/BIRCHDEV/doc/genographer/lsouthern.html
## Equations after code below

size.frags <- function(frags, ladder) {
  ## frags is a vector of fragments to size
  ## ladder is a two-column matrix, the first column is raw fragment size, the second
  ## column is the fragment size in bp

  makeSizer <- function(n) {
    P <- lsParams(ladder[n + (0:3), ])
    sizefun <- function(f) lsSize(f, P)
    
    return(sizefun)
  }

  funindex <- numeric(length(frags))

  for (i in 2:(nrow(ladder)))
    funindex[frags > ladder[i, 1]] <- funindex[frags > ladder[i, 1]] + 1

  funindex[funindex > nrow(ladder) - 3] <- 0

  sizerlist <- list()
  
  for(i in 1:max(funindex))
    sizerlist[[i]] <- makeSizer(i)
  
  res <- rep(NA, length(frags))

  for(i in seq_along(funindex))
    if (funindex[i] > 0)
      res[i] <- sizerlist[[funindex[i]]](frags[i])

  return(res)
}

lsSize <- function(m, P) {
  if(is.null(P$slope)) {
    (P$C1/(m - P$M1) + P$L1 + P$C2/(m - P$M2) + P$L2) / 2
  } else if (is.infinite(P$C1) & ! is.infinite(P$C2)) {
    ((P$C2/(m - P$M2) + P$L2) + (P$slope * m + P$intercept)) / 2
  } else if (is.infinite(P$C2) & ! is.infinite(P$C1)) {
    ((P$C1/(m - P$M1) + P$L1) + (P$slope * m + P$intercept)) / 2
  } else if (is.infinite(P$C2) & is.infinite(P$C1)) {
    P$slope * m + P$intercept
  }
}
  
lsParams <- function(mat){
  x1 <- mat[1,1]
  x2 <- mat[2,1]
  x3 <- mat[3,1]
  x4 <- mat[4,1]

  y1 <- mat[1,2]
  y2 <- mat[2,2]
  y3 <- mat[3,2]
  y4 <- mat[4,2]  

  C1 <- -1 * ((-x3 + x2) * (-x3 + x1) * (x1 - x2) * (-y3 + y2) * (y1 - y3) * (-y2 + y1)) /
    (-y1 * x3 + y2 * x3 - y2 * x1 + x2 * y1 - x2 * y3 + x1 * y3)^2

  L1 <- -1 * (-y1 * x1 * y3 + y1 * y3 * x3 - y2 * x2 * y1 + y2 * x2 * y3 - y2 * y3 * x3 +
             y2 * y1 * x1) / (-y1 * x3 + y2 * x3 - y2 * x1 + x2 * y1 - x2 * y3 + x1 * y3)

  M1 <- (-y2 * x2 * x1 + y2 * x2 * x3 - y1 * x1 * x3 + x1 * y3 * x3 - x2 * y3 * x3 + x2 *
        y1 * x1)/ (-y1 * x3 + y2 * x3 - y2 * x1 + x2 * y1 - x2 * y3 + x1 * y3) 

  C2 <- -1 * ((-x4 + x3) * (-x4 + x2) * (x2 - x3) * (-y4 + y3) * (y2 - y4) * (-y3 + y2)) /
    (-y2 * x4 + y3 * x4 - y3 * x2 + x3 * y2 - x3 * y4 + x2 * y4)^2

  L2 <- -1 * (-y2 * x2 * y4 + y2 * y4 * x4 - y3 * x3 * y2 + y3 * x3 * y4 - y3 * y4 * x4 +
             y3 * y2 * x2) / (-y2 * x4 + y3 * x4 - y3 * x2 + x3 * y2 - x3 * y4 + x2 * y4)

  M2 <- (-y3 * x3 * x2 + y3 * x3 * x4 - y2 * x2 * x4 + x2 * y4 * x4 - x3 * y4 * x4 + x3 *
        y2 * x2)/ (-y2 * x4 + y3 * x4 - y3 * x2 + x3 * y2 - x3 * y4 + x2 * y4) 

  slope <- intercept <- NULL

  if(is.infinite(C1) & is.infinite(C2)) {
    slope <- (y4 - y1)/(x4 - x1)
    intercept <- y4 - slope * x4
  } else if(is.infinite(C1)) {
    slope <- (y3 - y1) / (x3 - x1)
    intercept <- y3 - slope * x3
  } else if(is.infinite(C2)) {
    slope <- (y4 - y2) / (x4 - x2)
    intercept <- y4 - slope * x4
  }
  
  return(list(C1 = C1, L1 = L1, M1 = M1, C2 = C2, L2 = L2, M2 = M2, slope = slope,
              intercept = intercept))
}



##         c
##   L = --------- + Lo
##        m - Mo

## L is the size, m is the mobility, c, M0 and L0 are constants

## Use three known points to calculate the constants:

##           c 
##   y1 = ------- + Lo         (1)
##        x1 - Mo
 

##           c
##   y2 = ------- + Lo         (2)
##        x2 - Mo
 

##           c
##   y3 = ------- + Lo         (3)
##        x3 - Mo


##        -(-x3 + x2) (-x3 + x1) (x1 - x2) (-y3 + y2) (y1 - y3) (-y2 + y1)
##   c =  ----------------------------------------------------------------
##                                                             2
##             (-y1 x3 + y2 x3 - y2 x1 + x2 y1 - x2 y3 + x1 y3)
 

##        - (-y1 x1 y3 + y1 y3 x3 - y2 x2 y1 + y2 x2 y3 - y2 y3 x3 + y2 y1 x1)
##  Lo =  --------------------------------------------------------------------
##             (-y1 x3 + y2 x3 - y2 x1 + x2 y1 - x2 y3 + x1 y3)
 

##          (-y2 x2 x1 + y2 x2 x3 - y1 x1 x3 + x1 y3 x3 - x2 y3 x3 + x2 y1 x1)
##  Mo =  --------------------------------------------------------------------
##             (-y1 x3 + y2 x3 - y2 x1 + x2 y1 - x2 y3 + x1 y3)

##        c1/(m - m1) + L1 + c2/(m - m2) + L2
##   s = -------------------------------------
##                        2

