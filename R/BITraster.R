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
#' @slot character Character vector with filepath to a mask. See details
#'
#' @details
#'
#' `mask` must be a RasterLayer with the same extent and resolution as the BITraster
#'
#' @rdname BITraster
#' @export
BITraster <- setClass(
  'BITraster',
  contains = "RasterBrick",
  slots = list(
    variables = "character",
    data_directory = "character",
    temp_directory = "character",
    years = "numeric",
    mask = "character"
  )
)
#' @export
setValidity('BITraster', function(object) {
  require(raster)

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
  if (!identical(raster::origin(object), raster::origin(EROSBIT::sample_BIT))) {
    msg <- paste0('bad origin - should be 15, 15')
    errors <- c(errors, msg)
  }
  if (!identical(raster::res(object), raster::res(EROSBIT::sample_BIT))) {
    msg <- paste0('bad resolution - should be 30, 30')
    errors <- c(errors, msg)
  }

  mask <- object@mask
  if (length(mask) > 1) {
    msg <- paste0('mask must be length-1 character')
    errors <- c(errors, msg)
  } else if (length(mask) == 1) {
    if (!file.exists(mask)) {
      msg <- paste0('mask file does not exist')
      errors <- c(errors, msg)
    }
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
            .Object@crs <- EROSBIT::EROSBIT_CRS
            .Object@mask <- character()

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
  to <- as(from, 'RasterBrick')
  to <- as(to, 'BITraster')
  to
})
setAs('RasterStack', 'BITraster', def = function(from) {
  require(raster)
  to <- as(from, 'RasterBrick')
  to <- as(to, 'BITraster')
  to
})
#' @rdname BITraster
#' @export
.make_mask <- function(BITraster, mask, filename) {
  stopifnot(validObject(BITraster))
  require(raster)

  raster::values(mask) <- ifelse(is.na(raster::values(mask)), NA, 1)
  raster::writeRaster(x = mask, filename = filename, format = 'GTiff', overwrite = T)

  BITraster@mask <- filename
  stopifnot(CheckMaskValidity(BITraster), validObject(BITraster))
  return(BITraster)

}
#' @rdname BITraster
#' @export
setMethod('MakeMask',
          signature(mask = 'SpatialPolygonsDataFrame'),
          function(BITraster, mask, filename) {
            require(raster)
            require(sp)

            # change mask to RasterLayer format
            mask0 <- sp::spTransform(mask, EROSBIT::EROSBIT_CRS)

            # define target raster
            mask_blank <- raster::raster(mask0)
            # set resolution to 30m x 30m
            raster::res(mask_blank) <- raster::res(EROSBIT::sample_BIT)
            # adjust origin
            s0 <- raster::origin(EROSBIT::sample_BIT) - raster::origin(mask_blank)
            mask_blank <- raster::shift(mask_blank, dx = s0[1], dy = s0[2])
            # convert template to raster
            mask0 <- raster::rasterize(mask0, mask_blank)

            .make_mask(
              BITraster = BITraster,
              mask = mask0,
              filename = filename
            )
          })
#' @rdname BITraster
#' @export
.template_BIT <- function(template, dir, variables, years) {
  require(raster)

  # Input checks:
  stopifnot(
    dir.exists(dir),
    is.numeric(years),
    is.character(variables),
    all(variables %in% EROSBIT::VarNames()),
    all(years %in% EROSBIT::AllYears()),
    inherits(template, 'RasterLayer')
  )

  # Shift origin to the 15, 15 in EROSBIT
  raster::values(template) <- NA
  s0 <- origin(EROSBIT::sample_BIT) - origin(template)
  template <- raster::shift(template, dx = s0[1], dy = s0[2])

  # Convert to BITraster and add meta-data
  out_BIT <- as(template, 'BITraster')
  out_BIT@data_directory <- dir
  out_BIT@variables <- variables
  out_BIT@years <- years

  # Output checks:

  stopifnot(
    identical(crs(out_BIT), crs(template)),
    identical(bbox(out_BIT), bbox(template))
  )

  gc()
  return(out_BIT)

}
#' @rdname BITraster
#' @export
setMethod('TemplateBIT',
          signature(template = 'SpatialPolygonsDataFrame'),
          function(template, dir, variables, years) {
            require(raster)
            require(sp)

            # TODO: parallelize the rasterize? seems slow

            # Conver SPDF to RasterLayer

            # convert CRS
            template_new <- sp::spTransform(template, EROSBIT::EROSBIT_CRS)

            # define target raster
            template_rast <- raster::raster(template_new)
            # set resolution to 30m x 30m
            raster::res(template_rast) <- raster::res(EROSBIT::sample_BIT)
            # convert template to raster
            template_rast <- raster::rasterize(template_new, template_rast)

            .template_BIT(
              template = template_rast,
              dir = dir,
              variables = variables,
              years = years
            )
          })
#' @rdname BITraster
#' @export
setMethod('TemplateBIT',
          signature(template = 'RasterLayer'),
          function(template, dir, variables, years) {
            .template_BIT(
              template = template,
              dir = dir,
              variables = variables,
              years = years
            )
          })
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
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'SpatialPolygonsDataFrame'),
          function(BITraster, points, variables, years, parallel = FALSE) {

            require(sp)
            # xyFromCell from the parent BITraster gets coordinates

            stopifnot(
              identical(crs(BITraster) == crs(points)),
              identical(bbox(BITraster), bbox(points))
            )

            stop('not implemented')

            .subset_BIT_points(
              BITraster = BITraster,
              points = points0,
              variables = variables,
              years = years,
              parallel = parallel)

          })
