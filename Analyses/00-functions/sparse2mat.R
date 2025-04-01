sparse2mat <- function (x) 
{ #From pedicure version 3.0.1
  rn <- attr(x, "rowNames")
  nrow <- max(x[, 1])
  ncol <- max(x[, 2])
  y <- rep(0, nrow * ncol)
  y[(x[, 2] - 1) * nrow + x[, 1]] <- x[, 3]
  y[(x[, 1] - 1) * nrow + x[, 2]] <- x[, 3]
  a <- matrix(y, nrow = nrow, ncol = ncol, byrow = FALSE)
  dimnames(a) <- list(rn, rn)
  a
}
