mat2sparse <- function (x, rowNames = dimnames(x)[[1]], tol = 1e-12) 
{
  which <- (abs(x) > tol & lower.tri(x, diag = TRUE))
  df <- as.matrix(data.frame(row = t(row(x))[t(which)], col = t(col(x))[t(which)], 
                             val = t(x)[t(which)]))
  if (is.null(rowNames)) 
    rowNames <- as.character(1:nrow(x))
  attr(df, "rowNames") <- rowNames
  df
}
