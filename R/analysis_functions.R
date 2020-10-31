#' @title Analysis functions for EROSBIT dataset
#'
#' @description NA
#'
#' @details NA
#'
#' @author Brandon McNellis
#' @name analysis_functions
#' @rdname analysis_functions
NULL
#' @rdname analysis_functions
#' @export
AggregateTemplate <- function(BITraster, FUN) {
  require(raster) # implicit

  stopifnot(
    validObject(BITraster),
    nlayers(BITraster) > 1
  )

  zlist <- vector('list', nlayers(BITraster))

  for (i in seq_along(1:nlayers(BITraster))) {
    if (i == 1) next
    ii <- subset(BITraster, i)
    xx <- zonal(ii, subset(BITraster, 1), fun = FUN)
    zlist[[i]] <- xx
  }

  zlist <- zlist[-1]

}
