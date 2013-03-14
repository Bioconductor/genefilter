dist2 = function (x,
                  fun = function(a, b) mean(abs(a - b), na.rm = TRUE),
                  diagonal = 0) {

  if (!(is.numeric(diagonal) && (length(diagonal) == 1)))
    stop("'diagonal' must be a numeric scalar.")

  if (missing(fun)) {
    res = apply(x, 2, function(w) colMeans(abs(x-w), na.rm=TRUE))
  } else {
    res = matrix(diagonal, ncol = ncol(x), nrow = ncol(x))
    if (ncol(x) >= 2) {
      for (j in 2:ncol(x))
        for (i in 1:(j - 1))
          res[i, j] = res[j, i] = fun(x[, i], x[, j])
    } # if
  } # else

  colnames(res) = rownames(res) = colnames(x)
  return(res)
}
