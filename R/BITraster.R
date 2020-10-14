#' @title S4 class for extending Raster* to EROS BIT dataset
#'
#' @description NA
#'
#' @details NA
#'
#' @include generics.R
#' @import sp
#' @import raster
#'
#' @author Brandon McNellis
#' @name BITraster
#' @rdname BITraster
NULL
#' An S4 class for extending Raster* to EROS BIT
#'
#' @slot variables Character vector describing the contained BIT variables
#' @slot data_directory Full path to un-zipped BIT data, defaults to current working directory
#' @slot temp_directory Optional place to put temporary raster files (they get big)
#' @slot years Years of the analysis
#'
#' @rdname BITraster
#' @export
BITraster <- setClass(
  'BITraster',
  contains = "RasterStack",
  slots = list(
    variables = "character",
    data_directory = "character",
    temp_directory = "character",
    years = "numeric"
  )
)
#' @export
setValidity('BITraster', function(object) {
  errors <- character()

  if (!dir.exists(object@data_directory)) {
    msg <- paste0('data directory does not exist')
    errors <- c(errors, msg)
  }
  #if (!dir.exists(object@temp_directory)) {
  #  msg <- paste0('temp directory does not exist or has been moved')
  #  errors <- c(errors, msg)
  #}
  if (!all(object@variables %in% VarNames())) {
    msg <- paste0('malformed variable names')
    errors <- c(errors, msg)
  }
  if (!all(object@years %in% AllYears())) {
    msg <- paste0('invalid year range')
    errors <- c(errors, msg)
  }

  # returns
  if (length(errors) == 0) {
    TRUE
  } else {
    errors
  }
})
#' @rdname BITraster
#' @export
setMethod('initialize',
          signature(.Object = 'BITraster'),
          function (.Object, ...) {
            params <- list(...)

            if ('temp_directory' %in% names(params)) {
              .Object@temp_directory <- params$temp_directory
            } else {
              .Object@temp_directory <- tempdir()
            }

            if ('data_directory' %in% names(params)) {
              .Object@data_directory <- params$data_directory
            } else {
              .Object@data_directory <- getwd()
            }

            # other defaults

            # returns
            .Object <- callNextMethod()
            mt <- validObject(.Object)
            if (isTRUE(mt)) {
              return(.Object)
            } else {
              return(mt)
            }
          }
)
setAs('RasterLayer', 'BITraster', def = function(from) {
  require(raster)
  to <- as(from, 'RasterStack')
  to <- as(to, 'BITraster')
  to
})
#' @rdname BITraster
#' @export
ImportBIT <- function(dir, variables, years) {
  # Imports BIT meta-data to the BITraster format

  # Convenience wrapper for 'new' that can have detailed error messages
  # Also includes a masking option

  require(raster)

  stopifnot(
    dir.exists(dir),
    is.numeric(years),
    is.character(variables),
    all(variables %in% EROSBIT::VarNames()),
    all(years %in% EROSBIT::AllYears())
  )

  out_BIT <- new('BITraster', data_directory = dir, years = years, variables = variables)

  return(out_BIT)

}
#' @rdname BITraster
#' @export
.subset_BIT_points <- function(BITraster, points, variables, years, parallel = FALSE) {

  require(sp)
  require(raster)

  # This function subsets 30-m EROSBIT pixels from a SpatialPoints object
  # Adapted from FDF project code

  # Returns a SpatialPointsDataFrame
  # Future options to return a raster?

  # This could be changed to a method for BITraster if need be

  # sanity checks
  stopifnot(
    inherits(BITraster, 'BITraster'),
    #inherits(points, 'SpatialPoints'),
    is.character(variables),
    is.numeric(years),
    validObject(BITraster),
    all(years %in% BITraster@years),
    all(variables %in% BITraster@variables)
  )

  # main loop
  out_df <- points[, 2, drop = FALSE]

  if (parallel) {
    # implement a foreach here with foreach package
    require(foreach)
    stop('not implemented')

  } else if (!parallel) {

    for (i in seq_along(years)) {
      ii <- years[i]

      ilayer <- subset(BITraster, which(BITraster@years == ii))
      ivals <- extract(ilayer, points[, 1])
      out_df <- cbind(out_df, ivals)
      colnames(out_df)[ncol(out_df)] <- as.character(ii)

    }

    return(out_df)

  }

}
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'matrix'),
          function(BITraster, points, variables, years, parallel = FALSE) {
            .subset_BIT_points(
              BITraster = BITraster,
              points = points,
              variables = variables,
              years = years,
              parallel = parallel)

          })
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'data.frame'),
          function(BITraster, points, variables, years, parallel = FALSE) {
            .subset_BIT_points(
              BITraster = BITraster,
              points = as.data.frame(points),
              variables = variables,
              years = years,
              parallel = parallel)

          })
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'SpatialPoints'),
          function(BITraster, points, variables, years, parallel = FALSE) {

            require(sp)
            # xyFromCell from the parent BITraster gets coordinates

            stopifnot(
              identical(crs(BITraster) == crs(points)),
              identical(bbox(BITraster), bbox(points))
            )

            points0 <- BIT_spatial@coords

            .subset_BIT_points(
              BITraster = BITraster,
              points = points0,
              variables = variables,
              years = years,
              parallel = parallel)

          })
